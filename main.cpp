#define CGAL_HAS_THREADS
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <queue>
#include <functional>
#include <set>
#include <memory>
#include <algorithm>
#include <bitset>
#include <thread>
#include "MeshKernel/Mesh.h"
#include "lib_impl.h"
#include "Eigen/Dense"
#include "Eigen/Core"
#include "Eigen/Sparse"
#include "OsqpEigen/OsqpEigen.h"
#include <vector>
#include <fstream>
#include <limits>
#include <unordered_set>
#include <sstream>
#include "DSU.h"
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Triangulation_2_projection_traits_3.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_2.h>
#include <CGAL/Triangulation_conformer_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Polygon_mesh_processing/stitch_borders.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/merge_border_vertices.h>
#include <CGAL/Polygon_mesh_processing/manifoldness.h>
#include <CGAL/Polygon_mesh_processing/intersection.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/version.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Fuzzy_sphere.h>
#include <vector>
#include <iostream>
#include <CGAL/Search_traits_3.h>
#include <atomic>
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef K::Point_3 Point;
typedef CGAL::Search_traits_3<K> STraits;
typedef CGAL::Fuzzy_sphere<STraits> Fuzzy_circle;
typedef CGAL::Kd_tree<STraits> Kd_tree;

typedef CGAL::Delaunay_triangulation_3<K> Delaunay3D;

typedef CGAL::Constrained_Delaunay_triangulation_2<K2> CDT;
//#include "BVH.h"

using namespace std;

shared_ptr <MeshKernel::SurfaceMesh> mesh;

shared_ptr<CGALPolygon>cgal_polygon;

double default_move = 0.1;
int thread_num = 12;
double max_distance_limit = 1.30;

//Thicken2 ../data/
MeshKernel::SurfaceMesh ReadObjFile(const std::string &_InputFile) {
    //std::ifstream inputfile(_InputFile, std::ios::in);
    std::vector<MeshKernel::iGameVertex> vertices;
    std::vector<std::vector<MeshKernel::iGameVertexHandle> > faces;
    std::vector<double> move_dist;
    std::vector<std::vector<double>> normals;
    std::vector<std::vector<double>> uvs;
    std::unordered_map<int, int> V2N;// vertex to normal
    std::unordered_map<int, int> V2T;// vertex to uv
    std::cout << "Reading " << _InputFile << " File" << std::endl;
    FILE *inputfile = fopen(_InputFile.c_str(), "r");
    char inputs[100];
    while (fscanf(inputfile, "%[^\n]\n", inputs) != EOF) {
        string line(inputs);
        if (line[0] == '#') {
            continue;
        }
        std::stringstream linestream;
        linestream.str(line);
        std::string flag;
        linestream >> flag;
        if (flag == "v") {
            double x, y, z;
            linestream >> x >> y >> z;
            vertices.push_back(MeshKernel::iGameVertex(x, y, z));
        } else if (flag == "f") {
            std::vector<std::string> vex;
            std::string tmp;
            while (linestream >> tmp) vex.push_back(tmp);
            std::vector<MeshKernel::iGameVertexHandle> face(3);
            for (size_t i = 0; i < 3; i++) {
                size_t idx = 0;
                while (idx < vex[i].length() && std::isdigit(vex[i][idx])) idx++;
                int vh = std::stoi(vex[i].substr(0, idx)) - 1;//  obj start v from 1
                face[i] = (MeshKernel::iGameVertexHandle)(vh);
                if (vh >= vertices.size()) {
                    std::cerr << vh << " " << vertices.size() << std::endl;
                }
            }
            if (vex.size() >= 4)
                move_dist.push_back(std::stod(vex[3]));
            else
                move_dist.push_back(default_move);
            faces.push_back(face);
        }
    }
    if (!normals.empty()) {
        int ncnt = normals.size();
        for (int i = 0; i < vertices.size(); ++i) {
            int nidx = V2N[i];
            //if (nidx < 0 || nidx >= ncnt) printf("error: nidx = %d\n", nidx);// debug 用
            assert(nidx >= 0 && nidx < ncnt);
            vertices[i].setNormal(normals[nidx]);
        }
    }
    auto mesh = MeshKernel::SurfaceMesh(vertices, faces, move_dist);
    fclose(inputfile);
    return mesh;
}

double stx, sty, stz;
double grid_len = 1e10;

struct grid {
    int x;
    int y;
    int z;
    grid(){
        x=-1;
        y=-1;
        z=-1;
    }
    grid(int x,int y,int z){
        this->x=x;
        this->y=y;
        this->z=z;
    }
    friend bool operator<(const grid &a, const grid &b) {
        if (a.x != b.x)
            return a.x < b.x;
        if (a.y != b.y)
            return a.y < b.y;
        return a.z < b.z;
    }

    friend bool operator==(const grid &a, const grid &b) {
        return a.x == b.x && a.y == b.y && a.z == b.z;
    }

    bool valid(){
        return x>=0 && y>=0 && z>=0;
    }
};

std::hash<int>int_hash;
struct grid_hash{
    size_t operator () (grid x) const {
        return int_hash(x.x) ^ int_hash(x.y<<6) ^ int_hash(x.z<<12);
    }};
struct grid_equal
{
    bool operator() (grid a,  grid b) const {
        return a.x == b.x  &&  a.y == b.y &&  a.z == b.z;
    }
};



vector <vector<int>> GridVertexDir = {{0, 0, 0}, //0
                                      {1, 0, 0}, //1
                                      {0, 1, 0}, //2
                                      {0, 0, 1}, //3
                                      {1, 1, 0}, //4
                                      {1, 0, 1}, //5
                                      {0, 1, 1}, //6
                                      {1, 1, 1}}; //7


vector <vector<int> > container_grid_face = {{0,3,5,1},
                                             {0,2,6,3},
                                             {0,1,4,2},
                                             {2,4,7,6},
                                             {3,6,7,5},
                                             {1,5,7,4}
                                             };


MeshKernel::iGameVertex getGridVertex(grid g, int k) {
    double gsx = stx + g.x * grid_len;
    double gsy = sty + g.y * grid_len;
    double gsz = stz + g.z * grid_len;
    return MeshKernel::iGameVertex(gsx + GridVertexDir[k][0] * grid_len,
                                   gsy + GridVertexDir[k][1] * grid_len,
                                   gsz + GridVertexDir[k][2] * grid_len);
}

vector <vector<int>> DirectedGridEdge = {{1, 2, 3},
                                 { 4, 5},
                                 { 4, 6},
                                 {5, 6},
                                 { 7},
                                 { 7},
                                 { 7}};

MeshKernel::iGameVertex getGridiGameVertex(const MeshKernel::iGameVertex& small, const MeshKernel::iGameVertex& big, int k) {

    double x = small.x();
    double y = small.y();
    double z = small.z();
    if(GridVertexDir[k][0] > 0 )
        x = big.x();
    if(GridVertexDir[k][1] > 0 )
        y = big.y();
    if(GridVertexDir[k][2] > 0 )
        z = big.z();
    return MeshKernel::iGameVertex(x,y,z);
}

K2::Point_3 getGridK2Vertex(const K2::Point_3& small, const K2::Point_3& big, int k) {

    CGAL::Epeck::FT x = small.x();
    CGAL::Epeck::FT y = small.y();
    CGAL::Epeck::FT z = small.z();
    if(GridVertexDir[k][0] > 0 )
        x = big.x();
    if(GridVertexDir[k][1] > 0 )
        y = big.y();
    if(GridVertexDir[k][2] > 0 )
        z = big.z();
    return K2::Point_3(x,y,z);
}

grid vertex_to_grid(MeshKernel::iGameVertex v) {
    int x = int((v.x() - stx) / grid_len + myeps); // 先不管
    int y = int((v.y() - sty) / grid_len + myeps); // 先不管
    int z = int((v.z() - stz) / grid_len + myeps); // 先不管
    return grid{x, y, z};
}

struct GlobalFace{
    int idx0,idx1,idx2;
    int field_id;
    int useful;
    K2::Point_3 center;
    set<int>special_field_id;
    // 特例： 交点是自己
    // 交点是二合一点
    //

    std::unordered_map<int,vector<pair<K2::Point_3,int> > > ray_detect_map;//field 交点 int
};


struct GridVertex {
    int grid_type;
    vector <MeshKernel::iGameFaceHandle> field_list;
    vector<MeshKernel::iGameVertex>build_v;
    vector<int>build_v_global;
    vector<vector<size_t> >build_face;
    vector<K2::Triangle_3> generate_face_list;
    map<int,vector<int> >  field_face_though_list;
    map<size_t,pair<int,int> > face_hash_id_map;
    vector<int> global_face_list;
    std::list<K2::Triangle_3 > build_aabb_tree_triangle_list;
    GridVertex() {
        grid_type = -1;
    }
};



vector<MeshKernel::iGameVertex> field_move_vertex;
vector<vector<MeshKernel::iGameVertex> > field_move_vertices;
vector<vector<MeshKernel::iGameVertex> > field_move_face;
vector<K2::Triangle_3> field_move_K2_triangle;

//vector<vector<int> >approximate_field_face_table = {{0,1,2},{3,4,5},{0,2,4},{4,3,0},{1,5,4},{4,2,1},{0,3,5},{5,1,0}};

inline K2::Point_3 iGameVertex_to_Point_K2(const MeshKernel::iGameVertex& v){
    return K2::Point_3(v.x(),v.y(),v.z());
}


inline  MeshKernel::iGameVertex Point_K2_to_iGameVertex(const K2::Point_3& v){
    return MeshKernel::iGameVertex(CGAL::to_double(v.x()),CGAL::to_double(v.y()),CGAL::to_double(v.z()));
}


bool segment_in_line(K2::Segment_3 a,K2::Segment_3  b){
    CGAL::Epeck::FT d0 = CGAL::squared_distance(a.supporting_line(),b.vertex(0));
    CGAL::Epeck::FT d1 = CGAL::squared_distance(a.supporting_line(),b.vertex(1));
    if(d0 <= CGAL::Epeck::FT(0) &&
       d1 <= CGAL::Epeck::FT(0))
        return true;
    return false;
}

bool segment_in_line(K2::Segment_2 a,K2::Segment_2  b){
    K2::FT d0 = CGAL::squared_distance(a.supporting_line(),b.vertex(0));
    K2::FT d1 = CGAL::squared_distance(a.supporting_line(),b.vertex(1));
    if(d0 <= CGAL::Epeck::FT(0) &&
       d1 <= CGAL::Epeck::FT(0))
        return true;
    return false;
}
inline K::Point_3 PointK2_Point(K2::Point_3 p){
    return K::Point_3(CGAL::to_double(p.x()),CGAL::to_double(p.y()),CGAL::to_double(p.z()));
}
inline bool is_same_triangle(K::Point_3 p0,K::Point_3 p1,K::Point_3 p2,K::Point_3 q0,K::Point_3 q1,K::Point_3 q2,double eps ){
    int check = 0;
    check+=( sqrt(CGAL::to_double(squared_distance(p0,q0))) < eps || sqrt(CGAL::to_double(squared_distance(p0,q1))) <eps || sqrt(CGAL::to_double(squared_distance(p0,q2)))<eps );
    check+=( sqrt(CGAL::to_double(squared_distance(p1,q0))) < eps || sqrt(CGAL::to_double(squared_distance(p1,q1))) <eps || sqrt(CGAL::to_double(squared_distance(p1,q2)))<eps );
    check+=( sqrt(CGAL::to_double(squared_distance(p2,q0))) < eps || sqrt(CGAL::to_double(squared_distance(p2,q1))) <eps || sqrt(CGAL::to_double(squared_distance(p2,q2)))<eps );

    check+=( sqrt(CGAL::to_double(squared_distance(q0,p0))) < eps || sqrt(CGAL::to_double(squared_distance(q0,p1))) <eps || sqrt(CGAL::to_double(squared_distance(q0,p2)))<eps );
    check+=( sqrt(CGAL::to_double(squared_distance(q1,p0))) < eps || sqrt(CGAL::to_double(squared_distance(q1,p1))) <eps || sqrt(CGAL::to_double(squared_distance(q1,p2)))<eps );
    check+=( sqrt(CGAL::to_double(squared_distance(q2,p0))) < eps || sqrt(CGAL::to_double(squared_distance(q2,p1))) <eps || sqrt(CGAL::to_double(squared_distance(q2,p2)))<eps );
    if(check)
    cout << "check "<< check <<endl;
    return check == 6;

}


CGAL::Epeck::FT triangle_squared_aspect_ratio(vector<K2::Point_2>v){
    vector<CGAL::Epeck::FT> length;
    for(int i=0;i<3;i++){
        length.push_back((v[i] - v[(i+1)%3]).squared_length());
    }
    return *max_element(length.begin(),length.end()) / *min_element(length.begin(),length.end());
}

double triangle_squared_aspect_ratio(vector<MeshKernel::iGameVertex > v){
    vector<double> length;
    for(int i=0;i<3;i++){
        length.push_back((v[i] - v[(i+1)%3]).norm2());
    }
    return *max_element(length.begin(),length.end()) / *min_element(length.begin(),length.end());
}

//


vector<vector<K2::Point_3> > CGAL_CDT(vector<K2::Point_3> sorted_bound_vertex, vector<K2::Segment_3> cs,K2::Triangle_3 origin_face) {

    //CGAL::make_conforming_Delaunay_2();
   // cout << sorted_bound_vertex.size() <<" "<< cs.size() << endl;
    K2::Point_3 base_point_3d = origin_face.vertex(0);


    K2::Vector_3 origin_face_v0 = origin_face.vertex(1) - origin_face.vertex(0);
    K2::Vector_3 origin_face_v1 = origin_face.vertex(2) - origin_face.vertex(0);
    K2::Vector_3 new_direct = CGAL::cross_product(origin_face_v0,origin_face_v1);
    new_direct = new_direct /  CGAL::Epeck::FT(Point_K2_to_iGameVertex(K2::Point_3(0,0,0) +  new_direct).norm());
    K2::Vector_3 X_axis = origin_face_v0 / CGAL::Epeck::FT(Point_K2_to_iGameVertex(K2::Point_3(0,0,0) + origin_face_v0).norm());
    K2::Vector_3 Y_axis = CGAL::cross_product(new_direct,X_axis);


    vector<CDT::Vertex_handle> sorted_bound_vertex_in_cdt(sorted_bound_vertex.size());


    CDT cdt;
    map<CDT::Vertex_handle,K2::Point_3> mp;
    for(int i=0;i<sorted_bound_vertex.size();i++){
        CGAL::Epeck::FT x = (sorted_bound_vertex[i]-base_point_3d) * X_axis;
        CGAL::Epeck::FT y = (sorted_bound_vertex[i]-base_point_3d) * Y_axis;
        auto t = Point_K2_to_iGameVertex(sorted_bound_vertex[i]);
       // cout <<"X Y : " <<CGAL::to_double(x) <<" "<<CGAL::to_double(y) << endl;
        //cout <<"t: "<< t.x() <<" "<<t.y()<<" "<<t.z() << endl;
        auto tt = base_point_3d + x*X_axis + y*Y_axis;
      //  cout <<"tt: "<< tt.x() <<" "<<tt.y()<<" "<<tt.z() << endl;

        sorted_bound_vertex_in_cdt[i] = cdt.insert(K2::Point_2(x,y));
        mp[ sorted_bound_vertex_in_cdt[i] ] = sorted_bound_vertex[i];
    }


    CGAL::Polygon_2<K2>poly;

    for(int i=0;i<sorted_bound_vertex.size();i++){
        poly.push_back(cdt.point(sorted_bound_vertex_in_cdt[i]));
    }

    CGAL::make_conforming_Gabriel_2(cdt);

    //TODO 约束裁剪 ，，新增德劳内点怎么计算 ！！！！！！这里还有问题
    vector<vector<K2::Point_2> > faces;
    for(auto fit = cdt.finite_faces_begin();fit != cdt.finite_faces_end();fit++){
        faces.push_back({cdt.point(fit->vertex(0)),cdt.point(fit->vertex(1)),cdt.point(fit->vertex(2))});
    }
    for(K2::Segment_3 seg: cs){ // 这里是人工的带约束cdt，自己通过计算求交，感觉这里可以加速
        vector<vector<K2::Point_2> > faces_new;
        CGAL::Epeck::FT x0 = (seg.vertex(0) - base_point_3d) * X_axis;
        CGAL::Epeck::FT y0 = (seg.vertex(0) - base_point_3d) * Y_axis;
        CGAL::Epeck::FT x1 = (seg.vertex(1) - base_point_3d) * X_axis;
        CGAL::Epeck::FT y1 = (seg.vertex(1) - base_point_3d) * Y_axis;

        K2::Segment_2 seg2(K2::Point_2(x0,y0),K2::Point_2(x1,y1));

        for(auto face : faces){
            K2::Triangle_2 tri(face[0],face[1],face[2]);
            CGAL::cpp11::result_of<K2::Intersect_2(K2::Segment_2 , K2::Triangle_2)>::type
                    res_st = intersection(seg2,tri);

            if (res_st) {
                if (const K2::Segment_2 *s = boost::get<K2::Segment_2>(&*res_st)) {
                    bool is_same_edge = false;
                    for(int i=0;i<3;i++){
                        K2::Segment_2 edge(face[i],face[(i+1)%3]);
                        if(segment_in_line(edge,*s))
                            is_same_edge = true;
                    }
                    if(!is_same_edge){
                        CGAL::cpp11::result_of<K2::Intersect_2(K2::Line_2 , K2::Triangle_2)>::type
                                res_lt = intersection(seg2.supporting_line(),tri);
                        if (const K2::Segment_2 *ss = boost::get<K2::Segment_2>(&*res_lt)) {
                            K2::Point_2 v0 = ss->vertex(0);
                            K2::Point_2 v1 = ss->vertex(1);

                            vector <K2::Point_2> positive_side;
                            vector <K2::Point_2> negative_side;
                            for (int i = 0; i < 3; i++) {
                                if (seg2.supporting_line().has_on_positive_side(face[i])) {
                                    positive_side.push_back(face[i]);
                                } else if (seg2.supporting_line().has_on_negative_side(face[i])) {
                                    negative_side.push_back(face[i]);
                                }
                            }
                            for (vector <K2::Point_2> vs: {positive_side, negative_side}) {
                                if (vs.size() == 1) {
                                    faces_new.push_back({v0,v1,vs[0]});
                                } else if (vs.size() == 2) {
                                    vector<K2::Point_2> ans0v0;
                                    vector<K2::Point_2> ans0v1;
                                    vector<K2::Point_2> ans1v0;
                                    vector<K2::Point_2> ans1v1;
                                    K2::Segment_2 v0s0(v0,vs[0]);
                                    K2::Segment_2 v1s1(v1,vs[1]);
                                    CGAL::cpp11::result_of<K2::Intersect_2(K2::Segment_2 , K2::Segment_2)>::type
                                            res_ss2 = intersection(v0s0,v1s1);
                                    if(!res_ss2){
                                        ans0v0 = {v0,vs[0],v1};
                                        ans1v0 = {v0,v1,vs[1]};
                                    }
                                    else{
                                        ans0v0 = {v0,vs[1],v1};
                                        ans1v0 = {v0,v1,vs[0]};
                                    }
                                    ans0v1 = {vs[0],vs[1],v1};
                                    ans1v1 = {vs[0],vs[1],v0};
                                    if(triangle_squared_aspect_ratio(ans0v0) + triangle_squared_aspect_ratio(ans0v1) <
                                            triangle_squared_aspect_ratio(ans1v0) + triangle_squared_aspect_ratio(ans1v1)
                                    ){
                                        faces_new.push_back(ans0v0);
                                        faces_new.push_back(ans0v1);
                                    }
                                    else{
                                        faces_new.push_back(ans1v0);
                                        faces_new.push_back(ans1v1);
                                    }

                                }

                            }
                            continue;
                        }
                    }
                }
            }
            faces_new.push_back(face);
        }
        swap(faces_new,faces);
    }

    vector<vector<K2::Point_3> > ret;
    K2::Vector_3 vec = cross_product(origin_face.vertex(1)-origin_face.vertex(0),
                                    origin_face.vertex(2)-origin_face.vertex(0)
                                    );
    for(auto i : faces){
        vector<K2::Point_3>tmp(3);
        for(int j=0;j<3;j++){
            tmp[j] = base_point_3d + i[j].x()*X_axis + i[j].y()*Y_axis;
        }
        K2::Vector_3 vec1 = cross_product(tmp[1] - tmp[0],tmp[2] - tmp[0]);
        if(vec * vec1 < CGAL::Epeck::FT(0))
            swap(tmp[1],tmp[2]);

        ret.push_back(tmp);
    }

    return ret;
}


vector<K2::Point_3> global_vertex_list;

struct CoverageField {
    vector<K2::Point_3 > bound_face_vertex_exact;
    vector<K::Point_3 > bound_face_vertex_inexact;
    vector<vector<int> > bound_face_id;
    vector<vector<grid> >  bound_face_cross_field_list;
    vector<vector<K2::Segment_3> > bound_face_cutting_segment;
    vector<vector<K2::Point_3> > bound_face_cutting_point;
    vector<bool> bound_face_useful;

    K2::Point_3 center;
    std::vector<std::vector<K2::Point_3> > cdt_result;
    vector<int>cdt_result_cross_field_list_id;
    vector<bool>cdt_result_cross_field_list_useful;

    vector<K2::Point_3 > renumber_bound_face_vertex;
    vector<vector<int> > renumber_bound_face_id;
    vector<vector<grid> > renumber_bound_face_cross_field_list;
    int field_id;
    void do_cdt(){
        //            vector<K2::Point_3> sorted_bound_vertex;
//            vector<K2::Segment_3> cs;
//            set<pair<int,int> >cs_set;
        for(int i=0;i<bound_face_id.size();i++){
            vector<K2::Point_3> sorted_bound_vertex{bound_face_vertex_exact[bound_face_id[i][0]],
                                                    bound_face_vertex_exact[bound_face_id[i][1]],
                                                    bound_face_vertex_exact[bound_face_id[i][2]]
            };
            if(!bound_face_useful[i]) {
                cdt_result.push_back(sorted_bound_vertex);
                cdt_result_cross_field_list_id.push_back(i);
                cdt_result_cross_field_list_useful.push_back(false);
                continue;
            }

            vector<K2::Segment_3> cs;
            for(auto j : bound_face_cutting_segment[i]){
                int cnt = 0;
                cnt += (j.vertex(0) == bound_face_vertex_exact[bound_face_id[i][0]] );
                cnt += (j.vertex(0) == bound_face_vertex_exact[bound_face_id[i][1]] );
                cnt += (j.vertex(0) == bound_face_vertex_exact[bound_face_id[i][2]] );
                cnt += (j.vertex(1) == bound_face_vertex_exact[bound_face_id[i][0]] );
                cnt += (j.vertex(1) == bound_face_vertex_exact[bound_face_id[i][1]] );
                cnt += (j.vertex(1) == bound_face_vertex_exact[bound_face_id[i][2]] );
                if(cnt <2 ){
                    cs.push_back(j);
                }
            }
            sort(cs.begin(),cs.end(),[&](const K2::Segment_3& a,const K2::Segment_3& b){
                if(a.vertex(0).x() != b.vertex(0).x()){
                    return a.vertex(0).x() < b.vertex(0).x();
                }
                else if(a.vertex(0).y() != b.vertex(0).y()){
                    return a.vertex(0).y() < b.vertex(0).y();
                }
                else if(a.vertex(0).z() != b.vertex(0).z()){
                    return a.vertex(0).z() < b.vertex(0).z();
                }
                else if(a.vertex(1).x() != b.vertex(1).x()){
                    return a.vertex(1).x() < b.vertex(1).x();
                }
                else if(a.vertex(1).y() != b.vertex(1).y()){
                    return a.vertex(1).y() < b.vertex(1).y();
                } else
                    return a.vertex(1).z() < b.vertex(1).z();
            });


            cs.resize(std::unique(cs.begin(),cs.end())-cs.begin());

            K2::Triangle_3 tri(bound_face_vertex_exact[bound_face_id[i][0]],
                               bound_face_vertex_exact[bound_face_id[i][1]],
                               bound_face_vertex_exact[bound_face_id[i][2]]);

            for(auto j : bound_face_cutting_point[i]){
                if(j != bound_face_vertex_exact[bound_face_id[i][0]] &&
                        j != bound_face_vertex_exact[bound_face_id[i][1]] &&
                        j != bound_face_vertex_exact[bound_face_id[i][2]]
                        )
                sorted_bound_vertex.push_back(j);
            }
            sort(sorted_bound_vertex.begin(),sorted_bound_vertex.end(),[&](const K2::Point_3& a, const K2::Point_3& b){
                if(a.x() != b.x()){
                    return a.x() < b.x();
                }
                else if(a.y() != b.y()){
                    return a.y() < b.y();
                }
                return a.z() < b.z();
            });
            sorted_bound_vertex.resize(std::unique(sorted_bound_vertex.begin(),sorted_bound_vertex.end())-sorted_bound_vertex.begin());

            vector<vector<K2::Point_3> > res = CGAL_CDT(sorted_bound_vertex,cs,tri);

            for(int j=0;j<res.size();j++){
                cdt_result.push_back(res[j]);
                cdt_result_cross_field_list_id.push_back(i);
                cdt_result_cross_field_list_useful.push_back(true);
            }

        }

    }

    Kd_tree * tree;
    vector<K::Vector_3>encode_vec;
    vector<int>encode_num;
    vector<int>renumber_bound_face_vertex_global_id;
    vector<int>renumber_bound_face_global_id;
    vector<bool>renumber_bound_face_useful;
    std::unordered_map<unsigned long long,int> encode_map;
    int get_kdtree_id(K::Point_3 p){
        std::vector<Point> result;
        Fuzzy_circle fs(p, myeps);
        tree->search(std::back_inserter(result), fs);
        //  cout <<"result.size() : " << result.size() << endl;
        return encode_map[CGAL::hash_value(*result.begin())];
    };
    void renumber(){
        std::vector<K::Point_3> kd_tree_points;
        DSU dsu;
        for(int i=0;i<cdt_result.size();i++){
            for(int j=0;j<3;j++){
               kd_tree_points.emplace_back(CGAL::to_double(cdt_result[i][j].x()),
                                           CGAL::to_double(cdt_result[i][j].y()),
                                           CGAL::to_double(cdt_result[i][j].z())
                                           );
            }
        }

        tree = new Kd_tree(kd_tree_points.begin(),kd_tree_points.end());
        for(int i=0;i<kd_tree_points.size();i++){
            std::vector<Point> result;
            Fuzzy_circle fs(kd_tree_points[i], myeps);
            tree->search(std::back_inserter(result), fs);
            for (const Point& p : result) {
                dsu.join(CGAL::hash_value(kd_tree_points[i]),hash_value(p));
            }
        }
        int encode_cnt = 0;
        encode_map = dsu.encode(encode_cnt);
        renumber_bound_face_vertex.resize(encode_cnt);
        encode_vec.resize(encode_cnt);
        encode_num.resize(encode_cnt);
        renumber_bound_face_vertex_global_id.resize(encode_cnt);
        for(int i=0;i<encode_cnt;i++){
            encode_vec[i] = K::Vector_3(0,0,0);
            encode_num[i] = 0;
        }

        for(int i=0;i<kd_tree_points.size();i++){
          //  cout <<"getbelong:" <<i<<" "<< encode_map[CGAL::hash_value(kd_tree_points[i])] <<endl;
            encode_vec[encode_map[CGAL::hash_value(kd_tree_points[i])]] += kd_tree_points[i] - K::Point_3 (0,0,0);
            encode_num[encode_map[CGAL::hash_value(kd_tree_points[i])]]++;
        }


        for(int i=0;i<encode_cnt;i++){
            K::Point_3 avg = K::Point_3 (0,0,0) + encode_vec[i] / encode_num[i];
            renumber_bound_face_vertex[i] = K2::Point_3(avg.x(),avg.y(),avg.z());
        }
        //重构搜索要用原来的点
        for(int i=0;i<cdt_result.size();i++){

            K::Point_3 v0(CGAL::to_double(cdt_result[i][0].x()),CGAL::to_double(cdt_result[i][0].y()),CGAL::to_double(cdt_result[i][0].z()));
            K::Point_3 v1(CGAL::to_double(cdt_result[i][1].x()),CGAL::to_double(cdt_result[i][1].y()),CGAL::to_double(cdt_result[i][1].z()));
            K::Point_3 v2(CGAL::to_double(cdt_result[i][2].x()),CGAL::to_double(cdt_result[i][2].y()),CGAL::to_double(cdt_result[i][2].z()));
            int id0 = get_kdtree_id(v0);
            int id1 = get_kdtree_id(v1);
            int id2 = get_kdtree_id(v2);
            if(set<int>{id0,id1,id2}.size() != 3)continue;
            renumber_bound_face_id.push_back({id0,id1,id2});
            renumber_bound_face_cross_field_list.push_back(bound_face_cross_field_list[cdt_result_cross_field_list_id[i]]);
            renumber_bound_face_useful.push_back(cdt_result_cross_field_list_useful[i]);

        }
        std::list<K2::Triangle_3>tri_list;
        for(int i=0;i<renumber_bound_face_id.size();i++){
            tri_list.emplace_back(renumber_bound_face_vertex[renumber_bound_face_id[i][0]],
                                  renumber_bound_face_vertex[renumber_bound_face_id[i][1]],
                                  renumber_bound_face_vertex[renumber_bound_face_id[i][2]]);
        }
        Tree aabb_tree(tri_list.begin(),tri_list.end());
        auto iter = tri_list.begin();
        for(int i=0;i<renumber_bound_face_id.size();i++,iter++){
            K2::Point_3 tri_center = CGAL::centroid(*iter);
            K2::Ray_3 ray(tri_center,iter->supporting_plane().orthogonal_vector());
            std::list< Tree::Intersection_and_primitive_id<K2::Ray_3>::Type> intersections;
            aabb_tree.all_intersections(ray,std::back_inserter(intersections));
            vector<K2::Point_3> intersection_v;
            for(auto item : intersections) {
                if(const K2::Point_3* p = boost::get<K2::Point_3>(&(item.first))){
                    if(*p != tri_center){
                        intersection_v.push_back(*p);
                    }
                }
            }
            sort(intersection_v.begin(),intersection_v.end(),[&](const K2::Point_3 &a,const K2::Point_3 &b){
                if(a.x() != b.x()){
                    return a.x() < b.x();
                }
                else if(a.y() != b.y()){
                    return a.y() < b.y();
                }
                return a.z() < b.z();
            });
            intersection_v.resize(unique(intersection_v.begin(),intersection_v.end())-intersection_v.begin());
            if(intersection_v.size() % 2 == 0) { // 调整方向向外
                swap(renumber_bound_face_id[i][1],renumber_bound_face_id[i][2]);
            }
        }
        renumber_bound_face_global_id.resize(renumber_bound_face_id.size());
    }

    CoverageField(MeshKernel::iGameFaceHandle fh) {
        bound_face_vertex_inexact.emplace_back(mesh->fast_iGameVertex[mesh->fast_iGameFace[fh].vh(0)].x(),
                                       mesh->fast_iGameVertex[mesh->fast_iGameFace[fh].vh(0)].y(),
                                       mesh->fast_iGameVertex[mesh->fast_iGameFace[fh].vh(0)].z()
                            );
        bound_face_vertex_inexact.emplace_back(mesh->fast_iGameVertex[mesh->fast_iGameFace[fh].vh(1)].x(),
                                       mesh->fast_iGameVertex[mesh->fast_iGameFace[fh].vh(1)].y(),
                                       mesh->fast_iGameVertex[mesh->fast_iGameFace[fh].vh(1)].z()
                            );
        bound_face_vertex_inexact.emplace_back(mesh->fast_iGameVertex[mesh->fast_iGameFace[fh].vh(2)].x(),
                                       mesh->fast_iGameVertex[mesh->fast_iGameFace[fh].vh(2)].y(),
                                       mesh->fast_iGameVertex[mesh->fast_iGameFace[fh].vh(2)].z()
                            );
        for(auto v: field_move_vertices[mesh->fast_iGameFace[fh].vh(0)])
            bound_face_vertex_inexact.emplace_back(v.x(),v.y(),v.z());
        for(auto v: field_move_vertices[mesh->fast_iGameFace[fh].vh(1)])
            bound_face_vertex_inexact.emplace_back(v.x(),v.y(),v.z());
        for(auto v: field_move_vertices[mesh->fast_iGameFace[fh].vh(2)])
            bound_face_vertex_inexact.emplace_back(v.x(),v.y(),v.z());


        map<size_t,int> mp;
        for(int i=0;i<bound_face_vertex_inexact.size();i++){

            mp[CGAL::hash_value(bound_face_vertex_inexact[i])] = i;
            bound_face_vertex_exact.emplace_back(bound_face_vertex_inexact[i].x(),
                                                 bound_face_vertex_inexact[i].y(),
                                                 bound_face_vertex_inexact[i].z());
        }
        Delaunay3D dt;
        dt.insert(bound_face_vertex_inexact.begin(), bound_face_vertex_inexact.end());
        std::vector<K::Triangle_3> surface_triangles;
        for (auto fit = dt.finite_cells_begin(); fit != dt.finite_cells_end(); ++fit) {
            for (int i = 0; i < 4; ++i) {
                if (dt.is_infinite(fit->neighbor(i))) {
                    surface_triangles.push_back(dt.triangle(fit, i));
                }
            }
        }


        K2::Vector_3 center_vec = {0,0,0};
        for (const auto& triangle : surface_triangles) {
            int v0_id = mp[CGAL::hash_value(triangle.vertex(0))];
            int v1_id = mp[CGAL::hash_value(triangle.vertex(1))];
            int v2_id = mp[CGAL::hash_value(triangle.vertex(2))];

            bound_face_id.push_back({v0_id, v1_id, v2_id});
            center_vec += (centroid(K2::Triangle_3(bound_face_vertex_exact[v0_id],
                                                   bound_face_vertex_exact[v1_id],
                                                   bound_face_vertex_exact[v2_id])) - K2::Point_3(0,0,0)) ;
        }
        center =  K2::Point_3(0,0,0) + (center_vec / surface_triangles.size());

        std::vector<std::vector<std::size_t> > faces_list;

        for (auto i: bound_face_id) {
            faces_list.push_back({std::size_t(i[0]), std::size_t(i[1]), std::size_t(i[2])});
            bound_face_useful.push_back(true);
        }
        bound_face_cross_field_list.resize(bound_face_id.size());
        bound_face_cutting_segment.resize(bound_face_id.size());
        bound_face_cutting_point.resize(bound_face_id.size());
        poly = new CGAL::Polyhedron_3<K2>();

        PMP::polygon_soup_to_polygon_mesh(bound_face_vertex_exact, faces_list, *poly, CGAL::parameters::all_default());
        inside_ptr = new CGAL::Side_of_triangle_mesh<CGAL::Polyhedron_3<K2>, K2>(*poly);

    }
public:
    bool in_field(K2::Point_3 v) {
        //CGAL::Side_of_triangle_mesh<CGAL::Polyhedron_3<K2>, K2> inside(*poly);
        if ((*inside_ptr)(v) == CGAL::ON_BOUNDED_SIDE)
            return true;
        return false;
    }

    bool in_or_on_field(K2::Point_3 v) {
        //CGAL::Side_of_triangle_mesh<CGAL::Polyhedron_3<K2>, K2> inside(*poly);
        auto side = (*inside_ptr)(v);
        if (side== CGAL::ON_BOUNDED_SIDE || side == CGAL::ON_BOUNDARY)
            return true;
        return false;
    }
    CGAL::Polyhedron_3<K2> * poly;
    CGAL::Side_of_triangle_mesh<CGAL::Polyhedron_3<K2>, K2> * inside_ptr;
};

vector<CoverageField> coverage_field_list;

MeshKernel::iGameVertex do_quadratic_error_metric_check(MeshKernel::iGameVertexHandle vh,vector<MeshKernel::iGameFaceHandle> neighbor_face_list,bool &is_succ,double& dist){
    MeshKernel::iGameVertex v = mesh->fast_iGameVertex[vh];
    int m = neighbor_face_list.size();
    Eigen::SparseMatrix<double> hessian(3, 3);     //P: n*n正定矩阵,必须为稀疏矩阵SparseMatrix
    hessian.setZero();
    Eigen::VectorXd gradient(3);                  //Q: n*1向量
    gradient.setZero();
    Eigen::SparseMatrix<double> linearMatrix(m, 3); //A: m*n矩阵,必须为稀疏矩阵SparseMatrix
    linearMatrix.setZero();
    Eigen::VectorXd lowerBound(m);                  //L: m*1下限向量
    lowerBound.setZero();
    Eigen::VectorXd upperBound(m);                  //U: m*1上限向量
    upperBound.setZero();
    int cnt = 0;
    dist = 0;

    double avg_move_dist=0;
    double max_move_dist = 0;
    MeshKernel::iGameVertex avg_move_vertex(0,0,0);
    for (auto f : neighbor_face_list){
        avg_move_dist += mesh->fast_iGameFace[f].move_dist;
        max_move_dist += mesh->fast_iGameFace[f].move_dist;
    }
    avg_move_dist /= (1.0*neighbor_face_list.size());

    for (auto f : neighbor_face_list) {
        MeshKernel::iGameVertex normal
                = ((mesh->fast_iGameVertex[mesh->fast_iGameFace[f].vh(1)]
                    - mesh->fast_iGameVertex[mesh->fast_iGameFace[f].vh(0)]) %
                   (mesh->fast_iGameVertex[mesh->fast_iGameFace[f].vh(2)]
                    - mesh->fast_iGameVertex[mesh->fast_iGameFace[f].vh(0)])).normalize();

        normal = normal * -1;
//        MeshKernel::iGameVertex new_v = v + normal * mesh->fast_iGameFace[f].move_dist;
//
//        MeshKernel::iGameVertex move_max_v = v + normal * mesh->fast_iGameFace[f].move_dist*1.35;
//        MeshKernel::iGameVertex move_min_v = v + normal * mesh->fast_iGameFace[f].move_dist*0.8;
        MeshKernel::iGameVertex new_v = v + normal * avg_move_dist;
        avg_move_vertex += normal * avg_move_dist;

        // MeshKernel::iGameVertex move_max_v = v + normal * avg_move_dist * (1.25+0.1*depth);
        MeshKernel::iGameVertex move_max_v = v + normal * mesh->fast_iGameFace[f].move_dist * max_distance_limit;
        MeshKernel::iGameVertex move_min_v = v + normal * mesh->fast_iGameFace[f].move_dist;

        double d = -(normal.x() * new_v.x() + normal.y() * new_v.y() +  normal.z() * new_v.z());

        Eigen::Vector3d p(normal.x(),  normal.y(), normal.z());
        Eigen::Matrix3d A = p * p.transpose();
        Eigen::Vector3d D2(2*d*normal.x(),2*d*normal.y(),2*d*normal.z());
        Eigen::Vector3d LimitV(normal.x(),normal.y(),normal.z());

        double lower = normal.x()*move_min_v.x() + normal.y()*move_min_v.y() + normal.z()*move_min_v.z();
        double upper = normal.x()*move_max_v.x() + normal.y()*move_max_v.y() + normal.z()*move_max_v.z();



//        for(int i=0;i<3;i++)
//            for(int j=0;j<3;j++)
//                hessian.coeffRef(i,j)+=A.coeff(i,j)*2;
//        for(int i=0;i<3;i++)
//            gradient.coeffRef(i) +=  D2.coeffRef(i);


        for(int i=0;i<3;i++)
            linearMatrix.coeffRef(cnt,i) = LimitV[i];
        lowerBound.coeffRef(cnt) = lower;
        upperBound.coeffRef(cnt) = upper;

        cnt++;
    }
    avg_move_vertex /= neighbor_face_list.size();//mesh->FastNeighborFhOfVertex_[vh].size();

    min_move_g[vh] =  v + avg_move_vertex;

    avg_move_vertex = v + avg_move_vertex;;
    hessian.coeffRef(0,0) += (2.0);
    hessian.coeffRef(1,1) += (2.0);
    hessian.coeffRef(2,2) += (2.0);

    gradient.coeffRef(0) -= (2.0) * v.x();
    gradient.coeffRef(1) -= (2.0) * v.y();
    gradient.coeffRef(2) -= (2.0) * v.z();

    OsqpEigen::Solver solver;
    solver.settings()->setVerbosity(false);
    solver.settings()->setWarmStart(true);
    solver.data()->setNumberOfVariables(3);   //变量数n
    solver.data()->setNumberOfConstraints(m); //约束数m
    solver.data()->setHessianMatrix(hessian);
    solver.data()->setGradient(gradient);
    solver.data()->setLinearConstraintsMatrix(linearMatrix);
    solver.data()->setLowerBound(lowerBound);
    solver.data()->setUpperBound(upperBound);
    solver.initSolver();


    if(solver.solve()) {
        // if 长度过短
        is_succ = true;
        Eigen::VectorXd QPSolution = solver.getSolution();
        //dist
        MeshKernel::iGameVertex ret {QPSolution.coeffRef(0), QPSolution.coeffRef(1), QPSolution.coeffRef(2)};
        dist = 0;
        for(auto f : neighbor_face_list){
            MeshKernel::iGameVertex normal
                    = ((mesh->fast_iGameVertex[mesh->fast_iGameFace[f].vh(1)]
                        - mesh->fast_iGameVertex[mesh->fast_iGameFace[f].vh(0)]) %
                       (mesh->fast_iGameVertex[mesh->fast_iGameFace[f].vh(2)]
                        - mesh->fast_iGameVertex[mesh->fast_iGameFace[f].vh(0)])).normalize();
            dist += abs((ret - v) * normal);
        }
        return ret;
    }
    else{
        is_succ = false;
        return {0,0,0};
    }
}


vector<int> get_sub_state(int state,int face_list_size){
    vector<int> ret;
    function<void(int,int,int)> dfs = [&](int state,int pos,int change){
        if(pos == -1 && change) {
            ret.push_back(state);
            return;
        }
        for(int i=pos;i>=0;i--){
            if(state & (1<<i)){
                dfs(state,i-1,change);
                dfs(state ^ (1<<i),i-1,1);
                break;
            }
            else if(i == 0){
                dfs(state,i-1,change);
                break;
            }
        }
    };
    dfs(state,face_list_size,0);
    return ret;
}

vector<int> solve_by_dp(MeshKernel::iGameVertexHandle vh,vector<MeshKernel::iGameFaceHandle> neighbor_face_list){
    if(neighbor_face_list.size()>=25)exit(0);
    vector<double>dp;
    vector<int>dp_source;
    dp.resize(1<<neighbor_face_list.size());
    dp_source.resize(1<<neighbor_face_list.size());
    fill(dp.begin(),dp.end(),-1);
    fill(dp_source.begin(),dp_source.end(),-1);
    function<double(int)> dfs = [&](int state){
        if(dp[state] != -1) return dp[state];
        vector<MeshKernel::iGameFaceHandle> local_neighbor_face_list;
        for(int i=0;i<neighbor_face_list.size();i++){
            if(state & (1<<i))
                local_neighbor_face_list.push_back(neighbor_face_list[i]);
        }
        bool succ = false;
        double exceed_dist = 0;
        MeshKernel::iGameVertex v_new = do_quadratic_error_metric_check(vh,local_neighbor_face_list,succ,exceed_dist);
        if(succ){
           dp[state] = (v_new - mesh->fast_iGameVertex[vh]).norm();
           dp_source[state] = 0;
           return dp[state];
        }
        vector<int>sub_state = std::move(get_sub_state(state,neighbor_face_list.size()));
        double minx = 1e100;
        int from = -1;
        for(auto next: sub_state){
            if(next == 0 || state-next == 0)continue;
            if(next > state-next)continue;
            double sub_ans = dfs(next) + dfs(state-next);
            if(sub_ans < minx){
                from = next;
                minx = sub_ans;
            }
        }
        dp[state] = minx;
        dp_source[state] = from;
        return minx;
    };
    dfs((1<<neighbor_face_list.size())-1);
    queue<int>q;
    vector<int>ret;
    q.push((1<<neighbor_face_list.size())-1);
    while(!q.empty()){
        int now = q.front();
        q.pop();
        if(dp_source[now] == 0){
            ret.push_back(now);
        }
        else{
            q.push(dp_source[now]);
            q.push(now - dp_source[now]);
        }
    }
    return ret;
}

double mix_factor = 0.5;
void mix(int T){
    for (int times = 0; times < T; times++) {
        std::vector<double> fix_move_dist;
        fix_move_dist.resize(mesh->FaceSize() + 1);
        for (auto i: mesh->allfaces()) {
            double dist = 0;
            for (auto j: mesh->NeighborFh(i.first)) {
                dist += mesh->faces(j).move_dist;
            }
            dist /= mesh->NeighborFh(i.first).size();
            fix_move_dist[i.first] = (1.0 - mix_factor) * i.second.move_dist + mix_factor * dist;
        }
        for (auto i: mesh->allfaces()) {
            mesh->faces(i.first).move_dist = fix_move_dist[i.first];
        }
    }
}

vector<vector<int> > container_grid_dir{{-1,-1,-1},{-1,-1,0},{-1,-1,1},
                                        {-1,0,-1},{-1,0,0},{-1,0,1},
                                        {-1,1,-1},{-1,1,0},{-1,1,1},
                                        {0,-1,-1},{0,-1,0},{0,-1,1},
                                        {0,0,-1},{0,0,1},
                                        {0,1,-1},{0,1,0},{0,1,1},
                                        {1,-1,-1},{1,-1,0},{1,-1,1},
                                        {1,0,-1},{1,0,0},{1,0,1},
                                        {1,1,-1},{1,1,0},{1,1,1},};


void sort_by_polar_order(vector<K2::Point_3>& v,K2::Vector_3 orthogonal_direction){
    K2::Vector_3 center_v(0,0,0);
    for(auto j: v){
        center_v += j - K2::Point_3 (0,0,0);
    }
    center_v /= v.size();
    if(v.size() <=1)
        return ;
    K2::Point_3 center =  K2::Point_3 (0,0,0) + center_v;

    function<int(CGAL::Epeck::FT,CGAL::Epeck::FT)> quadrant = [](CGAL::Epeck::FT x,CGAL::Epeck::FT y){
        auto zero = CGAL::Epeck::FT(0);
        if(x>  zero && y > zero)return 1;
        else if(x<= zero && y >  zero)return 2;
        else if(x<= zero && y <= zero)return 3;
        return 4;
    };
    Plane_3 p;
    K2::Vector_3 x_axis = (v[0] - center) / sqrt(CGAL::to_double(CGAL::squared_distance(v[0],center)));
    K2::Vector_3 y_axis = CGAL::cross_product(orthogonal_direction,x_axis);
    y_axis = y_axis / sqrt(CGAL::to_double(y_axis.squared_length()));

    sort(v.begin(),v.end(),[&](K2::Point_3 a, K2::Point_3 b){

        CGAL::Epeck::FT x1 = (a - center) * x_axis;
        CGAL::Epeck::FT y1 = (a - center) * y_axis;
        CGAL::Epeck::FT x2 = (b - center) * x_axis;
        CGAL::Epeck::FT y2 = (b - center) * y_axis;
        int q1 = quadrant(x1,y1);
        int q2 = quadrant(x2,y2);
        if(q1!=q2)return q1<q2;
        else
            return x1*y2 - x2*y1 > 0;

    });
    return ;
}

vector<vector<size_t> >face_type_012{{0,1,2}};
std::vector<std::vector<std::size_t> > each_grid_face_list;


bool check_vertex_in_grid(K2::Point_3 small, K2::Point_3 big, K2::Point_3 v){
    return small.x() <= v.x() && v.x() <= big.x() &&
           small.y() <= v.y() && v.y() <= big.y() &&
           small.z() <= v.z() && v.z() <= big.z();
};
bool check_triangle_through_grid(K2::Point_3 small, K2::Point_3 big,CGAL::Polyhedron_3<K2>& frame_poly , K2::Triangle_3 tri) {
    if(check_vertex_in_grid(small,big,centroid(tri))){
        return true;
    }
    CGAL::Polyhedron_3<K2> this_face;
    vector<K2::Point_3> vs{tri.vertex(0),tri.vertex(1),tri.vertex(2)};
    PMP::polygon_soup_to_polygon_mesh(vs, face_type_012, this_face, CGAL::parameters::all_default());
    return CGAL::Polygon_mesh_processing::do_intersect(frame_poly,this_face);
};

int main(int argc, char* argv[]) {

    cout <<"CGAL_RELEASE_DATE:" << CGAL_RELEASE_DATE << endl;
    string input_filename(argv[1]);
    FILE *file9 = fopen( (input_filename + "_9.obj").c_str(), "w");
    FILE *file10 = fopen( (input_filename + "_10.obj").c_str(), "w");
    FILE *file11 = fopen( (input_filename + "_11.obj").c_str(), "w");
//
//
     FILE *file4 = fopen( (input_filename + "_4.obj").c_str(), "w");
    FILE *file5 = fopen( (input_filename + "_5.obj").c_str(), "w");
    FILE *file55 = fopen( (input_filename + "_55.obj").c_str(), "w");
    FILE *file3 = fopen( (input_filename + "_3.obj").c_str(), "w");
    FILE *file12 = fopen( (input_filename + "_12.obj").c_str(), "w");
    FILE *file13 = fopen( (input_filename + "_13.obj").c_str(), "w");
    FILE *file14 = fopen( (input_filename + "_14.obj").c_str(), "w");
    FILE *file2 = fopen( (input_filename + "_2.obj").c_str(), "w");
    FILE *file1 = fopen( (input_filename + "_1.obj").c_str(), "w");
//    //FILE *file5_3 = fopen( (input_filename + "_5_3.obj").c_str(), "w");
    FILE *file6 = fopen( (input_filename + "_6.obj").c_str(), "w");
    FILE *file7 = fopen( (input_filename + "_7.obj").c_str(), "w");
    FILE *file8 = fopen( (input_filename + "_8.obj").c_str(), "w");
    // freopen("../debugoutput.txt","w",stdout);
    default_move = 0.01;
    grid_len = 2.5;

    mix_factor = 0.5;
  //3.59
    mesh = make_shared<MeshKernel::SurfaceMesh>(ReadObjFile(input_filename)); grid_len = 0.1;
    mesh->initBBox();
    mesh->build_fast();
    cout <<"mesh->build_fast() succ" << endl;
    double default_move_dist = 0.05;
    if(argc > 2 )
        grid_len = stod(string(argv[2]));
    else
    {
        double sum = 0;
        for(int i=0;i<mesh->FaceSize();i++){
            double minx = *set<double>{(mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(0)] - mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(1)]).norm(),
                                (mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(1)] - mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(2)]).norm(),
                                (mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(2)] - mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(0)]).norm()
                                }.begin();
            double maxx = *set<double>{(mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(0)] - mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(1)]).norm(),
                                       (mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(1)] - mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(2)]).norm(),
                                       (mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(2)] - mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(0)]).norm()
            }.rbegin();

            sum += sqrt(minx*maxx);

        }
        default_move_dist = sum/mesh->FaceSize()*1.5/4;
        double x_len = (mesh->BBoxMax - mesh->BBoxMin).x();
        double y_len = (mesh->BBoxMax - mesh->BBoxMin).y();
        double z_len = (mesh->BBoxMax - mesh->BBoxMin).z();
        double min_len = min(min(x_len,y_len),z_len);
        double rate = x_len/min_len*y_len/min_len*z_len/min_len;

        cout <<"??"<<(4*thread_num) <<" "<<rate << endl; // minlen rate 格了 4* thread_num 格子
        double need_div = max(cbrt((4*thread_num)/rate),1.0);
        grid_len = min_len / need_div;//(mesh->BBoxMax - mesh->BBoxMin).norm()/ thread_num * 2;
        cout << "grid_len "<< grid_len<<endl;

    }
    //mesh = make_shared<MeshKernel::SurfaceMesh>(ReadObjFile("../data/test_orgv2.obj2")); grid_len = 12.5; double default_move_dist = 0.8;
    if(argc > 3 ) {
        cout <<"default_move_dist : "<< stod(string(argv[3])) << endl;
        for (int i = 0; i < mesh->FaceSize(); i++) {
            mesh->fast_iGameFace[MeshKernel::iGameFaceHandle(i)].move_dist =  stod(string(argv[3]));
        }
    }
    else if(*input_filename.rbegin() != '2') {
        auto tmp = mesh->BBoxMax -  mesh->BBoxMin;
        //default_move_dist = grid_len/4;//abs(*set<double>{tmp.x(),tmp.y(),tmp.z()}.begin());
        double minx = mesh->BBoxMin.z();
        double maxx = mesh->BBoxMax.z();
        cout << "default_move_dist: "<< default_move_dist << endl;
        for(int i=0;i<mesh->FaceSize();i++){
            double this_z = ((mesh->fast_iGameVertex[mesh->fast_iGameFace[MeshKernel::iGameFaceHandle(i)].vh(0)]+
                    mesh->fast_iGameVertex[mesh->fast_iGameFace[MeshKernel::iGameFaceHandle(i)].vh(1)]+
                    mesh->fast_iGameVertex[mesh->fast_iGameFace[MeshKernel::iGameFaceHandle(i)].vh(2)])/3).z();

            mesh->fast_iGameFace[MeshKernel::iGameFaceHandle(i)].move_dist = max(default_move_dist*1.5*((maxx-this_z)/(maxx-minx)),
                                                                                 default_move_dist/2);
        }
    }


    field_move_vertex.resize(mesh->VertexSize());
    field_move_vertices.resize(mesh->VertexSize());


    field_move_face.resize(mesh->FaceSize());
    field_move_K2_triangle.resize(mesh->FaceSize());

    cout <<"st do_quadratic_error_metric" << endl;

    for(int i=0;i<mesh->VertexSize();i++){
        vector<MeshKernel::iGameFaceHandle> neighbor_list;
        for(auto j : mesh->FastNeighborFhOfVertex_[i]){
            neighbor_list.push_back(j);
        }
        vector<int> dp_result = solve_by_dp(MeshKernel::iGameVertexHandle(i),neighbor_list);
       // cout <<"dp_result.size(): " <<dp_result.size() << endl;
        for(auto tmp: dp_result){
            vector<MeshKernel::iGameFaceHandle> dp_neighbor_list;
            for(int j=0;j<neighbor_list.size();j++){
                if(tmp & (1<< j))
                    dp_neighbor_list.push_back(neighbor_list[j]);
            }
            bool is_succ = false;
            double exceed_dist;
            MeshKernel::iGameVertex t = do_quadratic_error_metric_check(MeshKernel::iGameVertexHandle(i),dp_neighbor_list,is_succ,exceed_dist);
            field_move_vertices[i].push_back(t);
        }
    }



    for(int i=0;i<mesh->FaceSize();i++){
        coverage_field_list.push_back(CoverageField(MeshKernel::iGameFaceHandle(i)));
    }

    std::vector <std::shared_ptr<std::thread> > one_ring_select_thread_pool(thread_num);
    for(int i=0;i<thread_num;i++) {
        one_ring_select_thread_pool[i] = make_shared<std::thread>([&](int now_id) {
            for(int i=0;i<mesh->FaceSize();i++) {
                if (i % thread_num != now_id)continue;
                if(i%500==0)
                    cout << "one_ring_select_thread_pool "<< i << endl;

                set<int>neighbor_field;
                for(int j=0;j<3;j++) {
                    for (auto neighbor_id: mesh->FastNeighborFhOfEdge_[mesh->fast_iGameFace[i].eh(j)]) {
                        if(neighbor_id == i)continue;
                        neighbor_field.insert(neighbor_id);
                    }
                }

                vector<K2::Triangle_3>neighbor_face;
                for(auto neighbor_id: neighbor_field){
                    for (int k = 0; k < coverage_field_list[neighbor_id].bound_face_id.size(); k++) {
                        K2::Triangle_3 tri_this(coverage_field_list[neighbor_id].bound_face_vertex_exact[coverage_field_list[neighbor_id].bound_face_id[k][0]],
                                                coverage_field_list[neighbor_id].bound_face_vertex_exact[coverage_field_list[neighbor_id].bound_face_id[k][1]],
                                                coverage_field_list[neighbor_id].bound_face_vertex_exact[coverage_field_list[neighbor_id].bound_face_id[k][2]]
                        );

                        neighbor_face.push_back(tri_this);
                    }
                }

                for(int j=0;j<coverage_field_list[i].bound_face_id.size();j++) {
                    if ((coverage_field_list[i].bound_face_id[j][0] >= 3 ||
                            coverage_field_list[i].bound_face_id[j][1] >= 3 ||
                            coverage_field_list[i].bound_face_id[j][2] >= 3)){
                        bool flag = false;

                        K2::Triangle_3 tri_this(coverage_field_list[i].bound_face_vertex_exact[coverage_field_list[i].bound_face_id[j][0]],
                                                coverage_field_list[i].bound_face_vertex_exact[coverage_field_list[i].bound_face_id[j][1]],
                                                coverage_field_list[i].bound_face_vertex_exact[coverage_field_list[i].bound_face_id[j][2]]);
                        K2::Segment_3 e0(tri_this.vertex(0), tri_this.vertex(1));
                        K2::Segment_3 e1(tri_this.vertex(1), tri_this.vertex(2));
                        K2::Segment_3 e2(tri_this.vertex(2), tri_this.vertex(0));

                        vector<K2::Segment_3> vs_tmp;
                        for (auto other: neighbor_face) { //Point_3, or Segment_3, or Triangle_3, or std::vector < Point_3 >
                            CGAL::cpp11::result_of<K2::Intersect_3(K2::Triangle_3, K2::Triangle_3)>::type
                                    res_tt = intersection(tri_this, other);
                            if (res_tt) {
                                if (const K2::Segment_3 *s = boost::get<K2::Segment_3>(&*res_tt)) {
                                    vs_tmp.push_back(*s);
                                } else if (const K2::Triangle_3 *t = boost::get<K2::Triangle_3>(&*res_tt)) {
                                    vs_tmp.emplace_back(t->vertex(0), t->vertex(1));
                                    vs_tmp.emplace_back(t->vertex(1), t->vertex(2));
                                    vs_tmp.emplace_back(t->vertex(2), t->vertex(0));
                                } else if (std::vector<K2::Point_3> *vs = boost::get<std::vector<K2::Point_3 >>(&*res_tt)) {
                                    sort_by_polar_order(*vs, tri_this.supporting_plane().orthogonal_vector());
                                    for (int k = 0; k < vs->size(); k++) {
                                        vs_tmp.emplace_back(vs->operator[](k), vs->operator[]((k + 1) % vs->size()));
                                        // cerr << "run iiiiiiiiiiiiiiit" << endl;
                                    }
                                }
                            }
                        }
                        vector<K2::Segment_3> vs;
                        for (auto se: vs_tmp) {
                            if (!segment_in_line(se, e0) && !segment_in_line(se, e1) && !segment_in_line(se, e2)) {
                                vs.push_back(se);
                            }
                        }
                        vector<vector<K2::Point_3> > res = CGAL_CDT({coverage_field_list[i].bound_face_vertex_exact[coverage_field_list[i].bound_face_id[j][0]],
                                                                     coverage_field_list[i].bound_face_vertex_exact[coverage_field_list[i].bound_face_id[j][1]],
                                                                     coverage_field_list[i].bound_face_vertex_exact[coverage_field_list[i].bound_face_id[j][2]]}, vs, tri_this);
                        //cout << res.size() <<endl;
                        for (auto each_tri: res) {
                            bool patch_flag = false;
                            K2::Point_3 center = CGAL::centroid(K2::Triangle_3(each_tri[0], each_tri[1], each_tri[2]));
                            for (auto j: neighbor_field) {
                                if (tri_this.supporting_plane().oriented_side(coverage_field_list[i].center) !=
                                    tri_this.supporting_plane().oriented_side(coverage_field_list[j].center) &&
                                    tri_this.supporting_plane().oriented_side(coverage_field_list[i].center)+
                                    tri_this.supporting_plane().oriented_side(coverage_field_list[j].center) ==0 &&
                                        coverage_field_list[j].in_or_on_field(center)) {
                                    patch_flag = true;
                                    break;
                                }
                            }
                            if(!patch_flag){
                                flag = true;
                                break;
                            }
                        }
                        coverage_field_list[i].bound_face_useful[j]=flag;
                        //cout <<i<<" "<<j <<" "<< (flag?"yes":"no") << endl;
                    }
                    else{
                        coverage_field_list[i].bound_face_useful[j]=false;
                    }

                }
            }
        },i);
    }

    for(int i=0;i<thread_num;i++)
        one_ring_select_thread_pool[i]->join();

    for(auto  each_container_face : container_grid_face){
        each_grid_face_list.push_back({(size_t)each_container_face[0],(size_t)each_container_face[1],(size_t)each_container_face[2]});
        each_grid_face_list.push_back({(size_t)each_container_face[2],(size_t)each_container_face[3],(size_t)each_container_face[0]});
    }

    cout <<"build end "<< endl;

    cgal_polygon = make_shared<CGALPolygon>(mesh);

    std::list < K2::Triangle_3> triangles;

    mesh->initBBox();

    double max_move = 0;
    double min_move = 1e10;
    for (auto i: mesh->allfaces()) {
        max_move = max(max_move, abs(i.second.move_dist));
        min_move = min(min_move, abs(i.second.move_dist));
    }

    stx = mesh->BBoxMin.x() - max_move;
    sty = mesh->BBoxMin.y() - max_move;
    stz = mesh->BBoxMin.z() - max_move;

    printf("GL %lf\n", grid_len);

    unordered_map <grid, GridVertex, grid_hash, grid_equal> frame_grid_mp;

    vector <vector<int>> bfs_dir = {{0,  0,  1},
                                    {0,  1,  0},
                                    {1,  0,  0},
                                    {0,  0,  -1},
                                    {0,  -1, 0},
                                    {-1, 0,  0}};
    std::function < vector<grid>(grid) > get_neighbor = [&](grid g) {
        vector <grid> ret;
        for (auto i: bfs_dir) {
            int xx, yy, zz;
            xx = g.x + i[0];
            yy = g.y + i[1];
            zz = g.z + i[2];
            if (xx >= 0 && yy >= 0 && zz >= 0)
                ret.push_back({xx, yy, zz});
        }
        return ret;
    };
    int fsize = mesh->FaceSize();

    cout <<"bfs start \n" << endl;
    std::mutex bfs_mutex;
    std::vector <std::shared_ptr<std::thread> > bfs_thread_pool(thread_num);
    for(int i=0;i<thread_num;i++) {
        bfs_thread_pool[i] = make_shared<std::thread>([&](int now_id) {
            for (int face_id = 0; face_id < fsize; face_id++) {
                if(face_id % thread_num !=  now_id)continue;
                if (face_id % 1000 == 0)
                    printf("%d/%d\n", face_id, fsize);
                auto fh = make_pair(MeshKernel::iGameFaceHandle(face_id),
                                    mesh->faces(MeshKernel::iGameFaceHandle(face_id)));
                MeshKernel::iGameVertex center = (mesh->fast_iGameVertex[fh.second.vh(0)] +
                                                  mesh->fast_iGameVertex[fh.second.vh(1)] +
                                                  mesh->fast_iGameVertex[fh.second.vh(2)]) / 3;

                double move_limit = 0;
                for(int j=0;j<3;j++){
                    for(int k=0;k<field_move_vertices[j].size();k++){
                        move_limit = max(move_limit,(field_move_vertices[fh.second.vh(j)][k] - mesh->fast_iGameVertex[fh.second.vh(j)]).norm());
                    }
                }

                grid now = vertex_to_grid(center);
                vector <MeshKernel::iGameFaceHandle> face_and_neighbor;
                face_and_neighbor.push_back(fh.first);
                for (auto i: mesh->NeighborFh(fh.first)) {
                    face_and_neighbor.push_back(i);
                }

                queue <grid> q;
                unordered_set <grid,grid_hash,grid_equal> is_visit;
                vector <grid> center_neighbor = get_neighbor(now);
                for (auto bfs_start_node: center_neighbor) {
                    is_visit.insert(bfs_start_node);
                    q.push(bfs_start_node);
                }
                while (!q.empty()) {
                    now = q.front();
                    q.pop();
                    //cout << now.x <<" "<< now.y <<""
                    std::unique_lock<std::mutex>lock1(bfs_mutex,std::defer_lock);
                    lock1.lock();
                    auto iter = frame_grid_mp.find(now);
                    if (iter == frame_grid_mp.end()) {
                        iter = frame_grid_mp.insert(make_pair(now, GridVertex())).first;
                    }
                    iter->second.field_list.push_back(fh.first);
                    lock1.unlock();
                    vector <grid> neighbor = get_neighbor(now);
                    for (auto j: neighbor) {
                        if (!is_visit.count(j)) {
                            double dist = cgal_vertex_triangle_dist(fh.second, getGridVertex(j, 0), mesh);

                            if (dist <  (grid_len*1.74/2+move_limit) ) { //TODO : zheli youhua cheng pianyi juli de shiji jisuan
                                q.push(j);
                                is_visit.insert(j);
                            }
                        }
                    }
                }
            }
        },i);
    }

    for(int i=0;i<thread_num;i++)
        bfs_thread_pool[i]->join();

    cout <<"bfs end \n" << endl;

    auto frame_grid_mp_bk =  frame_grid_mp;
    for(auto i : frame_grid_mp_bk) {
        for(auto j : container_grid_dir){
            int nx = i.first.x+j[0];
            int ny = i.first.y+j[1];
            int nz = i.first.z+j[2];
            //if(nx>=0 && ny>=0 && nz>=0){
                for(auto z : i.second.field_list)
                    frame_grid_mp[{nx,ny,nz}].field_list.push_back(z);
           // }
        }
    }


    std::vector <std::shared_ptr<std::thread> > each_grid_thread_pool(thread_num);
    for(int i=0;i<thread_num;i++) { //todo 存在多线程异常
        each_grid_thread_pool[i] = make_shared<std::thread>([&](int now_id) {
            int each_grid_cnt =-1;
            for (auto each_grid = frame_grid_mp.begin(); each_grid != frame_grid_mp.end(); each_grid++) { // todo 逻辑修改
                each_grid_cnt++;
                if(each_grid_cnt % thread_num != now_id)continue;
                if(each_grid_cnt % (thread_num*200) == now_id)
                    printf("each_grid_cnt %d/%d\n",each_grid_cnt,(int)frame_grid_mp.size());

                MeshKernel::iGameVertex grid_vertex = getGridVertex(each_grid->first, 0);

                sort(each_grid->second.field_list.begin(), each_grid->second.field_list.end());
                each_grid->second.field_list.resize(
                        unique(each_grid->second.field_list.begin(), each_grid->second.field_list.end()) -
                        each_grid->second.field_list.begin());


                K2::Point_3 small = iGameVertex_to_Point_K2(getGridVertex(each_grid->first, 0));
                K2::Point_3 big = iGameVertex_to_Point_K2(getGridVertex(each_grid->first, 7));

                std::vector<K2::Point_3> ps;
                for(int i=0;i<8;i++){
                    ps.push_back(getGridK2Vertex(small,big,i));
                }
                std::list<K2::Triangle_3 >  frame_faces_list;

                for(auto i : each_grid_face_list){
                    frame_faces_list.emplace_back(ps[i[0]],ps[i[1]],ps[i[2]]);
                }
                Tree frame_aabb_tree(frame_faces_list.begin(),frame_faces_list.end());
                CGAL::Polyhedron_3<K2> frame_poly;
                PMP::polygon_soup_to_polygon_mesh(ps, each_grid_face_list, frame_poly, CGAL::parameters::all_default());

                for(int j = 0; j <  each_grid->second.field_list.size() ; j++){
                    int field_id = each_grid->second.field_list[j];
                    for(int k=0;k<coverage_field_list[field_id].bound_face_id.size();k++){

                        int v0_id = coverage_field_list[field_id].bound_face_id[k][0];
                        int v1_id = coverage_field_list[field_id].bound_face_id[k][1];
                        int v2_id = coverage_field_list[field_id].bound_face_id[k][2];

                        K2::Triangle_3 this_tri(coverage_field_list[field_id].bound_face_vertex_exact[v0_id],
                                                coverage_field_list[field_id].bound_face_vertex_exact[v1_id],
                                                coverage_field_list[field_id].bound_face_vertex_exact[v2_id]
                                                );
                        if(!coverage_field_list[field_id].bound_face_useful[k])continue;
                        if(check_triangle_through_grid(small,big,frame_poly,this_tri)){

                            each_grid->second.field_face_though_list[field_id].push_back(k);
                            each_grid->second.face_hash_id_map[this_tri.id()] = {field_id,k};
                            each_grid->second.build_aabb_tree_triangle_list.push_back(this_tri);

                        }

                    }
                }
                Tree this_grid_aabb_tree( each_grid->second.build_aabb_tree_triangle_list.begin(),
                                          each_grid->second.build_aabb_tree_triangle_list.end());

                //X: 这里写查询 因为aabb树指针的bug，所以反着做，用每个格子去查询每个面

                for(auto item : each_grid->second.field_face_though_list){
                    int field_id = item.first;
                   for(int bound_face_id : item.second){
                       int v0_id = coverage_field_list[field_id].bound_face_id[bound_face_id][0];
                       int v1_id = coverage_field_list[field_id].bound_face_id[bound_face_id][1];
                       int v2_id = coverage_field_list[field_id].bound_face_id[bound_face_id][2];
                        K2::Triangle_3 this_tri(coverage_field_list[field_id].bound_face_vertex_exact[v0_id],
                                                coverage_field_list[field_id].bound_face_vertex_exact[v1_id],
                                                coverage_field_list[field_id].bound_face_vertex_exact[v2_id]
                        );

                       std::list< Tree::Intersection_and_primitive_id<K2::Triangle_3>::Type> intersections;

                       this_grid_aabb_tree.all_intersections(this_tri,std::back_inserter(intersections));

                       unique_lock<mutex> lock(mutex);
                       for(auto item : intersections) { //todo 这里改成批量插入
                           if(const K2::Segment_3 * s = boost::get<K2::Segment_3>(&(item.first))){

                               map<size_t,pair<int,int> >::iterator iter = each_grid->second.face_hash_id_map.find(item.second->id());
                               if(iter->second.first != field_id ) {
                                   coverage_field_list[field_id].bound_face_cutting_segment[bound_face_id].push_back(
                                           *s);
                               }
                           }
                           else if(const K2::Point_3 * p = boost::get<K2::Point_3>(&(item.first))){
                               map<size_t,pair<int,int> >::iterator iter = each_grid->second.face_hash_id_map.find(item.second->id());
                               if(iter->second.first != field_id ) {
                                   coverage_field_list[field_id].bound_face_cutting_point[bound_face_id].push_back(*p);
                               }
                           }
                           else if(const std::vector<K2::Point_3> * v = boost::get<std::vector<K2::Point_3> >(&(item.first)) ){
                               map<size_t,pair<int,int> >::iterator iter = each_grid->second.face_hash_id_map.find(item.second->id());
                               if(iter->second.first != field_id ){
                                   for(int j=0;j<v->size();j++){
                                       coverage_field_list[field_id].bound_face_cutting_point[bound_face_id].push_back(v->at(j));
                                   }
                               }
                           }
                       }
                       coverage_field_list[field_id].bound_face_cross_field_list[bound_face_id].push_back(each_grid->first);
                   }
                }


              }
        }, i);
    }

    for(int i=0;i<thread_num;i++)
        each_grid_thread_pool[i]->join();



    // 这里修正每个面

    // 应该先做点的序号合并


    for (auto each_grid= frame_grid_mp.begin(); each_grid != frame_grid_mp.end(); each_grid++) {
        //if(!(each_grid->first.x == 26 && each_grid->first.y == 28 && each_grid->first.z == 9  ))continue;
        //if(!(each_grid->first.x == 19 && each_grid->first.y == 19 && each_grid->first.z == 1 ))continue;
        //if(!(each_grid->first.x == 2 && each_grid->first.y == 3 && each_grid->first.z == 0 ))continue;
        auto small  = getGridVertex(each_grid->first,0);
        auto big  = getGridVertex(each_grid->first,7);
        static int f3_id = 1;
        for (int ii = 0; ii < 7; ii++) {
            for (int jj = 0; jj < DirectedGridEdge[ii].size(); jj++) {
                int from = ii;
                int to = DirectedGridEdge[ii][jj];
                MeshKernel::iGameVertex fv = getGridiGameVertex(small, big, from);
                MeshKernel::iGameVertex tv = getGridiGameVertex(small, big, to);
                fprintf(file11, "v %lf %lf %lf\n", fv.x(), fv.y(), fv.z());
                fprintf(file11, "v %lf %lf %lf\n", tv.x(), tv.y(), tv.z());
                fprintf(file11, "l %d %d\n", f3_id, f3_id + 1);
                f3_id += 2;
            }
        }
    }

    std::vector <std::shared_ptr<std::thread> > field_vertex_numbering_thread_pool(thread_num);
    std::atomic<int>global_vertex_id_sum;
    for(int i=0;i<thread_num;i++) {
        field_vertex_numbering_thread_pool[i] = make_shared<std::thread>([&](int now_id) {
            for (int field_id = 0; field_id < fsize; field_id++) {
                if (field_id % thread_num != now_id)continue;
                coverage_field_list[field_id].field_id = field_id;
                coverage_field_list[field_id].do_cdt();
                coverage_field_list[field_id].renumber();
                global_vertex_id_sum+= coverage_field_list[field_id].renumber_bound_face_vertex.size();
            }
        },i);
    }
    for(int i=0;i<thread_num;i++)
        field_vertex_numbering_thread_pool[i]->join();


    global_vertex_list.resize(global_vertex_id_sum);

    int global_cnt = 0;
    for (int field_id = 0; field_id < fsize; field_id++) {
        for (int i = 0; i < coverage_field_list[field_id].renumber_bound_face_vertex.size(); i++) {

            coverage_field_list[field_id].renumber_bound_face_vertex_global_id[i] = global_cnt + i;
            global_vertex_list[global_cnt + i] =  coverage_field_list[field_id].renumber_bound_face_vertex[i];
        }
        global_cnt += coverage_field_list[field_id].renumber_bound_face_vertex.size();
    }

    std::vector<K::Point_3> kd_tree_points;
    map<unsigned long long,int> global_kd_tree_mp;
    for(int i=0;i<global_vertex_list.size();i++){
        K::Point_3 p = PointK2_Point(global_vertex_list[i]);
        unsigned long long id = CGAL::hash_value(p);
        global_kd_tree_mp[id] = i;
        kd_tree_points.push_back(p);
    }
    DSUMultiThread dsu_multi_thread((int)kd_tree_points.size());

    std::vector <std::shared_ptr<std::thread> > global_dsu_thread_pool(thread_num);
    dsu_multi_thread.run();
    for(int i=0;i<thread_num;i++) {
        global_dsu_thread_pool[i] = make_shared<std::thread>([&](int now_id) {
            Kd_tree global_kd_tree(kd_tree_points.begin(),kd_tree_points.end());
            for(int i=0;i<kd_tree_points.size();i++){
                if (i % thread_num != now_id)continue;
                std::vector<K::Point_3> result;
                Fuzzy_circle fs(kd_tree_points[i], myeps);
                global_kd_tree.search(std::back_inserter(result), fs);
                for (const Point& p : result) {
                    dsu_multi_thread.join(i,global_kd_tree_mp[CGAL::hash_value(p)]);
                    //dsu.join(CGAL::hash_value(kd_tree_points[i]),hash_value(p));
                }
            }
        },i);
    }
    for(int i=0;i<thread_num;i++)
        global_dsu_thread_pool[i]->join();

    cout <<"wait dsu_multi_thread stop"<< endl;
    dsu_multi_thread.stop();
    cout <<"wait dsu_multi_thread.join_thread->join()"<< endl;
    dsu_multi_thread.join_thread->join();

   // exit(0);
    cout <<"global_vertex_id_sum:" << global_vertex_id_sum <<endl;
    for (int field_id = 0; field_id < fsize; field_id++) {
        for (int i = 0; i < coverage_field_list[field_id].renumber_bound_face_vertex_global_id.size(); i++) {

            coverage_field_list[field_id].renumber_bound_face_vertex_global_id[i] = dsu_multi_thread.find_root(coverage_field_list[field_id].renumber_bound_face_vertex_global_id[i]);
        }
    }
    int global_face_cnt = 0;
    for (int field_id = 0; field_id < fsize; field_id++) {
        for (int i = 0; i < coverage_field_list[field_id].renumber_bound_face_id.size(); i++) {
            int tv0 =  coverage_field_list[field_id].renumber_bound_face_vertex_global_id[coverage_field_list[field_id].renumber_bound_face_id[i][0]];
            int tv1 =  coverage_field_list[field_id].renumber_bound_face_vertex_global_id[coverage_field_list[field_id].renumber_bound_face_id[i][1]];
            int tv2 =  coverage_field_list[field_id].renumber_bound_face_vertex_global_id[coverage_field_list[field_id].renumber_bound_face_id[i][2]];
            if(set<int>{tv0,tv1,tv2}.size() < 3)continue;
            global_face_cnt++;
        }
    }
    //exit(0);
    vector<GlobalFace> global_face_list(global_face_cnt);
    global_face_cnt = 0;
    for (int field_id = 0; field_id < fsize; field_id++) {
        for (int i = 0; i < coverage_field_list[field_id].renumber_bound_face_id.size(); i++) {
            int tv0 =  coverage_field_list[field_id].renumber_bound_face_vertex_global_id[coverage_field_list[field_id].renumber_bound_face_id[i][0]];
            int tv1 =  coverage_field_list[field_id].renumber_bound_face_vertex_global_id[coverage_field_list[field_id].renumber_bound_face_id[i][1]];
            int tv2 =  coverage_field_list[field_id].renumber_bound_face_vertex_global_id[coverage_field_list[field_id].renumber_bound_face_id[i][2]];

            int useful = coverage_field_list[field_id].renumber_bound_face_useful[i];
            if(set<int>{tv0,tv1,tv2}.size() < 3)continue;
            for(int j =0 ; j < coverage_field_list[field_id].renumber_bound_face_cross_field_list[i].size();j++){
                grid cross_grid = coverage_field_list[field_id].renumber_bound_face_cross_field_list[i][j];
                frame_grid_mp[cross_grid].global_face_list.push_back(global_face_cnt);
            }
            coverage_field_list[field_id].renumber_bound_face_global_id[i] = global_face_cnt;
            global_face_list[global_face_cnt].field_id = field_id;
            global_face_list[global_face_cnt].idx0 = tv0;
            global_face_list[global_face_cnt].idx1 = tv1;
            global_face_list[global_face_cnt].idx2 = tv2;
            global_face_list[global_face_cnt].useful = useful;
            global_face_list[global_face_cnt].center = CGAL::centroid(K2::Triangle_3(global_vertex_list[tv0],global_vertex_list[tv1],global_vertex_list[tv2]));
            global_face_cnt++;

        }
    }


    //exit(0);
    // 接下去就是求交

    cout << " 接下去时求交"<<endl;
    std::vector <std::shared_ptr<std::thread> > face_generate_ray_detect_thread_pool(thread_num);
    for(int i=0;i<thread_num;i++) {
        face_generate_ray_detect_thread_pool[i] = make_shared<std::thread>([&](int now_id) {
            int each_grid_cnt = -1;
            for (auto each_grid = frame_grid_mp.begin(); each_grid != frame_grid_mp.end(); each_grid++) {
                each_grid_cnt++;
                if (each_grid_cnt % thread_num != now_id)continue; //todo 这里需要开启
                if (each_grid_cnt % (thread_num * 200) == now_id)
                    printf("face_generate_ray_detect_thread_pool %d/%d\n", each_grid_cnt, (int) frame_grid_mp.size());

                unordered_map<unsigned long long , pair<int,int> > triangle_mapping;
                std::list<K2::Triangle_3>l;
                K2::Point_3 small = iGameVertex_to_Point_K2(getGridVertex(each_grid->first, 0));
                K2::Point_3 big = iGameVertex_to_Point_K2(getGridVertex(each_grid->first, 7));
                std::vector<K2::Point_3> ps;
                for(int i=0;i<8;i++){
                    ps.push_back(getGridK2Vertex(small,big,i));
                }

                CGAL::Polyhedron_3<K2> frame_poly;
                PMP::polygon_soup_to_polygon_mesh(ps, each_grid_face_list, frame_poly, CGAL::parameters::all_default());

                std::function<bool(K2::Point_3)> vertex_in_frame = [&](K2::Point_3 v){
                    return small.x() <= v.x() && v.x() <= big.x() &&
                           small.y() <= v.y() && v.y() <= big.y() &&
                           small.z() <= v.z() && v.z() <= big.z();
                };

                for(int i=0;i<each_grid->second.field_list.size();i++){
                    int field_id = each_grid->second.field_list[i];
                    bool useful = false;
                    if(coverage_field_list[field_id].in_field(midpoint(K2::Segment_3(ps[0],ps[7])))){
                        useful = true;
                    }
                    if(!useful)
                    {
                        if(vertex_in_frame(coverage_field_list[field_id].center)){
                            useful = true;
                        }
                    }
                    if(!useful)
                        useful = CGAL::Polygon_mesh_processing::do_intersect(*coverage_field_list[field_id].poly,frame_poly);
                    if(useful){
                        for(int j=0;j<coverage_field_list[field_id].renumber_bound_face_global_id.size();j++){
                            int fid_global = coverage_field_list[field_id].renumber_bound_face_global_id[j];
                            K2::Point_3 v0 = global_vertex_list[global_face_list[fid_global].idx0];
                            K2::Point_3 v1 = global_vertex_list[global_face_list[fid_global].idx1];
                            K2::Point_3 v2 = global_vertex_list[global_face_list[fid_global].idx2];
                            K2::Triangle_3 tri(v0,v1,v2);
                            triangle_mapping[tri.id()] = {field_id,fid_global};
                            l.push_back(tri);
                        }
                    }
                }
                Tree aabb_tree(l.begin(),l.end());//todo 明天把这里打开debug一下看看？？

                for(int i=0;i<each_grid->second.global_face_list.size();i++){
                    int global_face_id = each_grid->second.global_face_list[i];
                    if(global_face_list[global_face_id].useful == false)continue;
                    K2::Point_3 v0 = global_vertex_list[global_face_list[global_face_id].idx0];
                    K2::Point_3 v1 = global_vertex_list[global_face_list[global_face_id].idx1];
                    K2::Point_3 v2 = global_vertex_list[global_face_list[global_face_id].idx2];

                    K2::Vector_3 ray_vec = K2::Triangle_3(v0,v1,v2).supporting_plane().orthogonal_vector();

                    K2::Ray_3 ray(global_face_list[global_face_id].center,-ray_vec);
                    std::list< Tree::Intersection_and_primitive_id<K2::Ray_3>::Type> intersections;
                    aabb_tree.all_intersections(ray,std::back_inserter(intersections));
                    for(auto item : intersections) {
                        pair<int,int> belong = triangle_mapping[item.second->id()];
                        int which_field = belong.first;
                        int which_id = belong.second;
                        if(which_field == global_face_list[global_face_id].field_id)continue;
                        if(const K2::Point_3* p = boost::get<K2::Point_3>(&(item.first))){
                            if(*p != global_face_list[global_face_id].center)
                                global_face_list[global_face_id].ray_detect_map[which_field].emplace_back(*p,which_id);
//                            else
//                                global_face_list[global_face_id].special_face_id.insert(which_field);
                        }
                        else{
                            global_face_list[global_face_id].special_field_id.insert(which_field);
                        }
                    }
                    // 参考3164行的写法写
                    // 注意事项：1. 建立射线的时候要所有面
                    // 注意事项：2。可能产生重复
                }

            }
        },i);
    }
    for(int i=0;i<thread_num;i++)
        face_generate_ray_detect_thread_pool[i]->join();


    cout <<"start generate"<<endl;

    std::vector <std::shared_ptr<std::thread> > global_face_final_generate_thread_pool(thread_num);
    atomic<int> flag = 0;
    for(int i=0;i<thread_num;i++) {
        global_face_final_generate_thread_pool[i] = make_shared<std::thread>([&](int now_id) {
            std::list<K2::Triangle_3>origin_face_list;
            for(int i=0;i<mesh->FaceSize();i++){
                K2::Point_3 v0(mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(0)].x(),
                               mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(0)].y(),
                               mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(0)].z());

                K2::Point_3 v1(mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(1)].x(),
                               mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(1)].y(),
                               mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(1)].z());

                K2::Point_3 v2(mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(2)].x(),
                               mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(2)].y(),
                               mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(2)].z());
                origin_face_list.emplace_back(v0,v1,v2);
            }

            Tree origin_face_tree(origin_face_list.begin(),origin_face_list.end());

            for (int i = 0; i < global_face_list.size(); i++) {
//                if(field_id %100 ==0)
//                    cout << "field_id" << field_id << endl;
                if (i % thread_num != now_id)continue;
                if(global_face_list[i].useful == false )continue;
                for(std::unordered_map<int,vector<pair<K2::Point_3,int> > >::iterator it = global_face_list[i].ray_detect_map.begin();
                it != global_face_list[i].ray_detect_map.end(); it++
                ){////field 交点 对应的面
                    unordered_set<int>se;
                    std::vector<K2::Point_3> intersection_v;
                    for(int j=0;j<it->second.size();j++){
                        if(se.count(it->second[j].second))continue;
                        se.insert(it->second[j].second);
                        intersection_v.push_back(it->second[j].first);
                    }
                    sort(intersection_v.begin(),intersection_v.end(),[&](const K2::Point_3 &a,const K2::Point_3 &b){
                        if(a.x() != b.x()){
                            return a.x() < b.x();
                        }
                        else if(a.y() != b.y()){
                            return a.y() < b.y();
                        }
                        return a.z() < b.z();
                    });
                    int old_size = intersection_v.size();
                    intersection_v.resize(unique(intersection_v.begin(),intersection_v.end())-intersection_v.begin());
                    if(old_size != intersection_v.size()){
                        global_face_list[i].special_field_id.insert(it->first);
                        continue;
                    }
                    if(intersection_v.size()%2 == 1) {

                        global_face_list[i].useful = -200;
                        break;
                    }
                    // 1 是自己  // 若全无交点，可想而知一定在外，若有交点，那不一定在外，则正方向被覆盖，然后需要探测负方向，和自己
                    // 这里加入自己射线正负方向的寻找
                    // 2 出现多次 加入
                }
                if(  global_face_list[i].useful ){
                    K2::Point_3 v0 = global_vertex_list[global_face_list[i].idx0];
                    K2::Point_3 v1 = global_vertex_list[global_face_list[i].idx1];
                    K2::Point_3 v2 = global_vertex_list[global_face_list[i].idx2];
                    if(!cgal_polygon->inMesh(centroid(K2::Triangle_3(v0,v1,v2)))){
                        global_face_list[i].useful = -300;
                    }
                    if(origin_face_tree.squared_distance(centroid(K2::Triangle_3(v0,v1,v2))) < CGAL::Epeck::FT(myeps)){
                        global_face_list[i].useful = -300;
                    }
                }
            }
        },i);
    }
    for(int i=0;i<thread_num;i++)
        global_face_final_generate_thread_pool[i]->join();

    for(int i=0;i<global_vertex_list.size();i++){
        fprintf(file6,"v %lf %lf %lf\n",CGAL::to_double(global_vertex_list[i].x()),
                CGAL::to_double(global_vertex_list[i].y()),
                CGAL::to_double(global_vertex_list[i].z()));
    }

    for (int i = 0; i < global_face_list.size(); i++) {
        if(global_face_list[i].useful>0)
        fprintf(file6,"f %d %d %d\n",global_face_list[i].idx0+1,global_face_list[i].idx1+1,global_face_list[i].idx2+1);
    }
    for(int i=0;i<global_vertex_list.size();i++){
        fprintf(file8,"v %lf %lf %lf\n",CGAL::to_double(global_vertex_list[i].x()),
                CGAL::to_double(global_vertex_list[i].y()),
                CGAL::to_double(global_vertex_list[i].z()));
    }

    for (int i = 0; i < global_face_list.size(); i++) {
        if(global_face_list[i].useful>0 && global_face_list[i].special_field_id.size())//这里还得搞的
            fprintf(file8,"f %d %d %d\n",global_face_list[i].idx0+1,global_face_list[i].idx1+1,global_face_list[i].idx2+1);
    }

    //8月22 解决漏求交的问题，先定位到具体的面，然后再修改
    return 0;
}
// 1 2 3 4 5 6
// 5 1 2 3 7 8
//




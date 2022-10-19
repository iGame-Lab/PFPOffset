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
typedef CGAL::Constrained_Delaunay_triangulation_2<K2> CDT;
//#include "BVH.h"


using namespace std;


shared_ptr <MeshKernel::SurfaceMesh> mesh;

shared_ptr<CGALPolygon>cgal_polygon;

Tree cgal_global_aabbtree;
double default_move = 0.1;
int thread_num = 14;
int file_id;
const double tolerance = 0.2;

const int out_side = -1;
const int on_side = 0;
const int in_side = 1;

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

vector <vector<int>> GridEdge = {{1, 2, 3}, //0
                                 {0, 4, 5}, //1
                                 {0, 4, 6}, //2
                                 {0, 5, 6}, //3
                                 {1, 2, 7}, //4
                                 {1, 3, 7}, //5
                                 {2, 3, 7}, //6
                                 {4, 5, 6}}; //7

vector <vector<int> > container_grid_face = {{0,3,5,1},
                                             {0,2,6,3},
                                             {0,1,4,2},
                                             {2,4,7,6},
                                             {3,6,7,5},
                                             {1,5,7,4}
                                             };
/*
f 1 4 6 2
f 1 3 7 4
f 1 2 5 3
f 3 5 8 7
f 4 7 8 6
f 2 6 8 5
 */
vector <grid> get_edge_neighbor(grid g1, grid g2) {

    grid now = min(g1,g2);
    int diff = -1;
    int nowx = now.x;
    int nowy = now.y;
    int nowz = now.z;
    if(g1.x != g2.x) {
        diff = 0;
    }
    else if(g1.y != g2.y) {
        diff = 1;
    }
    else {
        diff = 2;
    }

    if (diff == 0) {
        return {{nowx, nowy,     nowz},
                {nowx, nowy - 1, nowz},
                {nowx, nowy - 1, nowz - 1},
                {nowx, nowy,     nowz - 1}};
    } else if (diff == 1) {
        return {{nowx,     nowy, nowz},
                {nowx,     nowy, nowz - 1},
                {nowx - 1, nowy, nowz - 1},
                {nowx - 1, nowy, nowz}};
    } else {
        assert(diff == 2);
        return {{nowx,     nowy,     nowz},
                {nowx - 1, nowy,     nowz},
                {nowx - 1, nowy - 1, nowz},
                {nowx,     nowy - 1, nowz}};
    }
}




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

MeshKernel::iGameVertex getTinyGridVertex(const MeshKernel::iGameVertex& small, const MeshKernel::iGameVertex& big,int k) {

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





grid getGridFrameVertex(grid g, int k) {
    grid v = g;
    v.x += GridVertexDir[k][0];
    v.y += GridVertexDir[k][1];
    v.z += GridVertexDir[k][2];
    return v;
}

grid vertex_to_grid(MeshKernel::iGameVertex v) {
    int x = int((v.x() - stx) / grid_len + myeps);
    int y = int((v.y() - sty) / grid_len + myeps);
    int z = int((v.z() - stz) / grid_len + myeps);
    return grid{x, y, z};
}




bool vertex_in_grid(grid g, MeshKernel::iGameVertex v) {
    auto v1 = getGridVertex(g, 0);
    auto v2 = getGridVertex(g, 7);
    return v1.x() <= v.x() && v.x() <= v2.x() &&
           v1.y() <= v.y() && v.y() <= v2.y() &&
           v1.z() <= v.z() && v.z() <= v2.z();
}

struct GridVertex {
    int grid_type;
    vector <MeshKernel::iGameFaceHandle> face_list;

    GridVertex() {
        grid_type = -1;
    }
};


int vertex_state(std::vector<MeshKernel::iGameFaceHandle>face_list, MeshKernel::iGameVertex v,
                 std::shared_ptr<std::vector<MeshKernel::iGameFaceHandle> > field_face_list = nullptr){


    int  side_type=-2;
    vector<MeshKernel::iGameFace>vf;
    for(auto i : face_list){
        vf.push_back(mesh->fast_iGameFace.at(i));
    }
    cgal_aabbtree_query(vf, v, mesh, side_type, cgal_polygon);
    if(side_type == 0)
        return 0;
    bool in_field = false;
    //  cout <<"pos:"<< v.x()<<" "<< v.y()<<" "<< v.z()<<endl;
    for(MeshKernel::iGameFaceHandle fh : face_list) {
        double dist = cgal_vertex_triangle_dist(mesh->fast_iGameFace.at(fh), v, mesh);
        //   cout << "fh:"<< fh<<" dist:"<<dist <<endl;
        if(vertex_field_state(fh,v,mesh)){
            in_field = true;
            if(field_face_list){
                field_face_list->push_back(fh);
            }
        }

//        if(dist <= mesh->fast_iGameFace.at(fh).move_dist) {
//            in_field = true;
//            if(field_face_list){
//                field_face_list->push_back(fh);
//            }
//        }

    }

    if(in_field) {
        return 1;
    }
    return 2;
}


//struct ThickenCubeEdge{
//    grid from;
//    grid to;
//    int type;
//    MeshKernel::iGameVertex edge_mesh_intersect_vertex;
//    MeshKernel::iGameVertex field_far_vertex;
//};

//map<grid,vector<ThickenCubeEdge> > thicken_cube_edge_mp;
//
//double ternary_search(MeshKernel::iGameVertex v1,MeshKernel::iGameVertex v2,MeshKernel::iGameFaceHandle fh){
//    double l=0;
//    double r=1;
//    const double eps = 0.001;
//    Point_3 a(mesh->fast_iGameVertex[mesh->fast_iGameFace[fh].vh(0)].x(),
//              mesh->fast_iGameVertex[mesh->fast_iGameFace[fh].vh(0)].y(),
//              mesh->fast_iGameVertex[mesh->fast_iGameFace[fh].vh(0)].z());
//    Point_3 b(mesh->fast_iGameVertex[mesh->fast_iGameFace[fh].vh(1)].x(),
//              mesh->fast_iGameVertex[mesh->fast_iGameFace[fh].vh(1)].y(),
//              mesh->fast_iGameVertex[mesh->fast_iGameFace[fh].vh(1)].z());
//    Point_3 c(mesh->fast_iGameVertex[mesh->fast_iGameFace[fh].vh(2)].x(),
//              mesh->fast_iGameVertex[mesh->fast_iGameFace[fh].vh(2)].y(),
//              mesh->fast_iGameVertex[mesh->fast_iGameFace[fh].vh(2)].z());
//    Triangle tri (a,b,c);
//    while(r-l>eps)
//    {
//        double mid1=(2*l+r)/3;
//        double mid2=(l+2*r)/3;
//        MeshKernel::iGameVertex Vmid1 = v1*(1-mid1) + v2*mid1;
//        Point Pmid1(Vmid1.x(),Vmid1.y(),Vmid1.z());
//        MeshKernel::iGameVertex Vmid2 = v1*(1-mid2) + v2*mid2;
//        Point Pmid2(Vmid2.x(),Vmid2.y(),Vmid2.z());
//        auto dist1 =  CGAL::squared_distance(Pmid1,tri);
//        auto dist2 =  CGAL::squared_distance(Pmid2,tri);
//       if(dist1>dist2) l=mid1;
//       else r=mid2;
//    }
//    return l;
//}
//
//void solve_grid(grid g,std::vector<MeshKernel::iGameFaceHandle>face_handle_list){
//    sort(face_handle_list.begin(), face_handle_list.end());
//    face_handle_list.resize(unique(face_handle_list.begin(), face_handle_list.end()) - face_handle_list.begin());
//    vector<MeshKernel::iGameFace>face_list;
//    for(auto i : face_handle_list){
//        face_list.push_back(mesh->fast_iGameFace[i]);
//    }
//    std::vector<BVH::BVH_Face>bvh_faces;
//    for(MeshKernel::iGameFaceHandle i: face_handle_list){
//        BVH::BVH_Face tmp;
//        for(int id = 0 ; id < 3 ;id++) {
//            tmp.vertices.emplace_back(mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(id)].x(),
//                                      mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(id)].y(),
//                                      mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(id)].z());
//        }
//        tmp.index = i;
//        bvh_faces.push_back(tmp);
//    }
//    BVH::BVH_Tree bvhTree;
//    bvhTree.buildBVH_Tree(bvh_faces);
//
//    for (int ev1 = 0; ev1 < GridEdge.size(); ev1++) {
//        for (auto ev2: GridEdge[ev1]) {
//            if (ev1 < ev2) {
//                for(MeshKernel::iGameFaceHandle i: face_handle_list) {
//                    grid gv1 = getGridFrameVertex(g, ev1);
//                    grid gv2 = getGridFrameVertex(g, ev2);
//                    MeshKernel::iGameVertex iv1 = getGridVertex(g, ev1);
//                    MeshKernel::iGameVertex iv2 = getGridVertex(g, ev2);
//                    double near_pos_div = ternary_search(iv1, iv2, i);
//                    MeshKernel::iGameVertex near_pos = iv1 * (1 - near_pos_div) + iv2 * near_pos_div;
//                    double dist = cgal_vertex_triangle_dist(mesh->fast_iGameFace[i],near_pos,mesh);
//                    if(dist > mesh->fast_iGameFace[i].move_dist)
//                        continue;
//
//
//
//                    MeshKernel::iGameVertex left_start;
//                    MeshKernel::iGameVertex right_start;
//                    bool need_left_start;
//                    bool need_right_start;
//                    int side = -1;
//                    cgal_aabbtree_query(face_list,near_pos,mesh,side,cgal_polygon);
//                    if(side ==0){
//
//                    }
//                    else{
//
//                    }
//                }
//            }
//        }
//    }
//
///*
//MeshKernel::iGameVertex cgal_aabbtree_query(std::vector<MeshKernel::iGameFace>& f_list,
//                                            MeshKernel::iGameVertex v,
//                                            std::shared_ptr<MeshKernel::SurfaceMesh>mesh,
//                                            int & side,std::shared_ptr<CGALPolygon>cgalinmesh);
// */
//
//
////    for(MeshKernel::iGameFaceHandle i: face_list){
////
////
////    }
//
//
//}
//
//vector<pair<MeshKernel::iGameVertex,MeshKernel::iGameVertex> >
//        solve_edge_field(MeshKernel::iGameVertex v1,MeshKernel::iGameVertex v2,std::vector<MeshKernel::iGameFaceHandle>face_list
//                ,BVH::BVH_Tree& local_bvh_tree){
//    vector<pair<MeshKernel::iGameVertex,MeshKernel::iGameVertex>  >segments_in_field;
//    for(MeshKernel::iGameFaceHandle fh : face_list){
//        // 先计算最近的点555
//
//
//
//
//        // 在计算平面与交点
//
//
//    }
//
//    return segments_in_field;
//};


vector<MeshKernel::iGameVertex> field_move_vertex;
vector<vector<MeshKernel::iGameVertex> > field_move_face;
vector<K2::Triangle_3> field_move_K2_triangle;

//vector<vector<int> >approximate_field_face_table = {{0,1,2},{3,4,5},{0,2,4},{4,3,0},{1,5,4},{4,2,1},{0,3,5},{5,1,0}};

inline K::Point_3 iGameVertex_to_Point(const MeshKernel::iGameVertex& v){
    return K::Point_3(v.x(),v.y(),v.z());
}

inline K2::Point_3 iGameVertex_to_Point_K2(const MeshKernel::iGameVertex& v){
    return K2::Point_3(v.x(),v.y(),v.z());
}

inline  MeshKernel::iGameVertex Point_K_to_iGameVertex(const K::Point_3& v){
    return MeshKernel::iGameVertex(CGAL::to_double(v.x()),CGAL::to_double(v.y()),CGAL::to_double(v.z()));
}

inline  MeshKernel::iGameVertex Point_K2_to_iGameVertex(const K2::Point_3& v){
    return MeshKernel::iGameVertex(CGAL::to_double(v.x()),CGAL::to_double(v.y()),CGAL::to_double(v.z()));
}

bool in_triangle_positive_side(MeshKernel::iGameFaceHandle fh,const MeshKernel::iGameVertex& v){
    MeshKernel::iGameVertex v0 = mesh->fast_iGameVertex[mesh->fast_iGameFace[fh].vh(0)];
    MeshKernel::iGameVertex v1 = mesh->fast_iGameVertex[mesh->fast_iGameFace[fh].vh(1)];
    MeshKernel::iGameVertex v2 = mesh->fast_iGameVertex[mesh->fast_iGameFace[fh].vh(2)];
    MeshKernel::iGameVertex normal = ((v1 - v0) % (v2 - v0)).normalize();
    MeshKernel::iGameVertex center = (v0 + v1 + v2)/3;
    if(((v - center) * normal)>0){
        return true;
    }
    return false;
}



//vector<K2::Triangle_3> intersection_divide(vector<K2::Triangle_3> v){
//
//
//    vector<K2::Triangle_3> ret;
//    for(int i=0;i<v.size();i++){
//        vector<K2::Triangle_3>tri_list;
//        tri_list.push_back(v[i]);
//        for(int j=0;i<v.size();j++){
//            if(j == i)continue;
//            for(int k=0;k<tri_list.size();k++){
//                CGAL::cpp11::result_of<K2::Intersect_3(K2::Triangle_3 , K2::Triangle_3)>::type
//                                res_tt = intersection(tri_list[k],v[j]);
//                    if (res_tt) {
//
//                        if (const K2::Segment_3 *s = boost::get<K2::Segment_3>(&*res_tt)) {
//
//                        }
//                        // TODO : 后面再补交成面
//                    }
//            }
//        }
//
//
//        for(int j=0;j<tri_list.size();j++){
//            if(){
//
//            }
//            else{
//
//            }
//
//
//        }
//
//
//        for(int j=0;j<v.size();j++){
//            if(j!= i){
//
//            }
//        }
//
//
//    }
//
//
//
//}

bool segment_in_line(K2::Segment_3 a,K2::Segment_3  b){
    CGAL::Epeck::FT d0 = CGAL::squared_distance(a.supporting_line(),b.vertex(0));
    CGAL::Epeck::FT d1 = CGAL::squared_distance(a.supporting_line(),b.vertex(1));
    if(d0 < CGAL::Epeck::FT(myeps) &&
       d1 < CGAL::Epeck::FT(myeps) &&
       abs(CGAL::to_double(d1-d0)) < myeps)
        return true;
    return false;
}

bool segment_in_line(K2::Segment_2 a,K2::Segment_2  b){
    K2::FT d0 = CGAL::squared_distance(a.supporting_line(),b.vertex(0));
    K2::FT d1 = CGAL::squared_distance(a.supporting_line(),b.vertex(1));
    if(d0 < CGAL::Epeck::FT(myeps) &&
       d1 < CGAL::Epeck::FT(myeps) &&
       abs(CGAL::to_double(d1-d0)) < myeps)
        return true;
    return false;
}




// 0 不相交 1 线交   2 面交
pair<int,K2::Segment_3> face_cross_face(const K2::Triangle_3 &tri1 ,const K2::Triangle_3 &tri2){
    K2::Segment_3 e0(tri1.vertex(0),
                     tri1.vertex(1));
    K2::Segment_3 e1(tri1.vertex(1),
                     tri1.vertex(2));
    K2::Segment_3 e2(tri1.vertex(2),
                     tri1.vertex(0));
    CGAL::cpp11::result_of<K2::Intersect_3(K2::Triangle_3, K2::Triangle_3)>::type
            res = intersection(tri1,tri2);
    if (res) {
        if (const K2::Segment_3 *s = boost::get<K2::Segment_3>(&*res)) {
            bool same_edge = false;
            for(K2::Segment_3 edge : {e0,e1,e2}){
                if(segment_in_line(edge,*s))
                    same_edge = true;
            }
            if(!same_edge){
                return std::make_pair(1,*s);
            }
            //判断交线是不是边6
        }// TODO : 面交再说！；；；；

    }
    return std::make_pair(0,K2::Segment_3());
}


CGAL::Epeck::FT triangle_squared_aspect_ratio(vector<K2::Point_2>v){
    vector<CGAL::Epeck::FT> length;
    for(int i=0;i<3;i++){
        length.push_back((v[i] - v[(i+1)%3]).squared_length());
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
    for(K2::Segment_3 seg: cs){
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

    for(auto i : faces){
        vector<K2::Point_3>tmp(3);
        for(int j=0;j<3;j++){
            tmp[j] = base_point_3d + i[j].x()*X_axis + i[j].y()*Y_axis;
        }
        ret.push_back(tmp);
    }


    return ret;
    for(auto fit = cdt.finite_faces_begin();fit != cdt.finite_faces_end();fit++){
        bool flag = true;
        for(int vid = 0 ;vid < 3 ;vid++) {
            auto v1 = cdt.point(fit->vertex(vid));
            if(CGAL::oriented_side(v1,poly) ==  CGAL::NEGATIVE ){
                flag = false;
            }
        }
        if(flag){
            vector<K2::Point_3>face;
            for(int vid = 0 ;vid < 3 ;vid++){

                if(mp.count(fit->vertex(vid))) {
                    cout <<"preis" <<CGAL::to_double(mp[fit->vertex(vid)].x()) <<" "<< CGAL::to_double(mp[fit->vertex(vid)].y())
                    <<" "<< CGAL::to_double(mp[fit->vertex(vid)].z()) << endl;
                    auto now = base_point_3d + cdt.point(fit->vertex(vid)).x()*X_axis + cdt.point(fit->vertex(vid)).y()*Y_axis;
                    cout <<"nowis" <<CGAL::to_double(now.x()) <<" "<< CGAL::to_double(now.y())
                         <<" "<< CGAL::to_double(now.z()) << endl;

                }
                face.push_back(base_point_3d + cdt.point(fit->vertex(vid)).x()*X_axis + cdt.point(fit->vertex(vid)).y()*Y_axis);
                cout <<"$$$$$" << (CGAL::squared_distance(
                                                          base_point_3d + cdt.point(fit->vertex(vid)).x()*X_axis + cdt.point(fit->vertex(vid)).y()*Y_axis,
                                                  origin_face
                        ) /*== CGAL::Epick::FT(0)*/ ) << std::endl;
            }
            ret.push_back(face);
        }
    }



    return ret;















    vector<K2::Segment_3>cutting_segment;
    for(int i=0; i<cs.size();i++){
        vector<K2::Segment_3> subdivide;
        subdivide.push_back(cs[i]);
        for(int j=0;j<sorted_bound_vertex.size();j++){
            vector<K2::Segment_3> subdivide_new;
            K2::Segment_3 seg (sorted_bound_vertex[j], sorted_bound_vertex[(j+1) % sorted_bound_vertex.size()]);
            for(auto k:subdivide){
                bool flag = true;
                CGAL::cpp11::result_of<K2::Intersect_3(K2::Segment_3 , K2::Segment_3)>::type
                        res_ss = intersection(k,seg);
                if (res_ss) {
                    if (const K2::Point_3 *p = boost::get<K2::Point_3>(&*res_ss)) {
                        if(*p != k.vertex(0) && *p != k.vertex(1)){
                            subdivide_new.push_back(K2::Segment_3(k.vertex(0),*p));
                            subdivide_new.push_back(K2::Segment_3(*p,k.vertex(1)));
                            flag = false;
                        }
                    }
                }
                if(flag)
                    subdivide_new.push_back(k);
            }
            swap(subdivide_new,subdivide);
        }

        for(int j=0;j<cs.size();j++){
            if(i==j)continue;
            vector<K2::Segment_3> subdivide_new;
            for(auto k:subdivide){
                bool flag = true;
                CGAL::cpp11::result_of<K2::Intersect_3(K2::Segment_3 , K2::Segment_3)>::type
                        res_ss = intersection(k,cs[j]);
                if (res_ss) {
                    if (const K2::Point_3 *p = boost::get<K2::Point_3>(&*res_ss)) {
                        if(*p != k.vertex(0) && *p != k.vertex(1)){
                            subdivide_new.push_back(K2::Segment_3(k.vertex(0),*p));
                            subdivide_new.push_back(K2::Segment_3(*p,k.vertex(1)));
                            flag = false;
                        }
                    }
                }
                if(flag)
                    subdivide_new.push_back(k);
            }
            swap(subdivide_new,subdivide);
        }
        for(auto j: subdivide)
            cutting_segment.push_back(j);
    }

    for(int i=0;i<cutting_segment.size();i++){
        for(int j=i+1;j<cutting_segment.size();j++){
            CGAL::cpp11::result_of<K2::Intersect_3(K2::Segment_3 , K2::Segment_3)>::type
                    res_ss = intersection(cutting_segment[i],cutting_segment[j]);
            if (res_ss) {
                if (const K2::Point_3 *p = boost::get<K2::Point_3>(&*res_ss)) {
                    if (*p != cutting_segment[i].vertex(0) && *p != cutting_segment[i].vertex(1)) {
                        cout <<"ERROR " << endl;
                    }
                }
            }
        }
    }

    for(auto i = cdt.all_vertices_begin();i!= cdt.all_vertices_end();i++){
        printf("v %lf %lf %lf\n",CGAL::to_double(i->point().x()),CGAL::to_double(i->point().y()),0);
    }




    vector<pair<CDT::Vertex_handle,CDT::Vertex_handle> > cutting_segment_in_cdt(cutting_segment.size());









    for(int i=0;i<cutting_segment.size();i++){
        CGAL::Epeck::FT x = (cutting_segment[i].vertex(0)-base_point_3d)*X_axis;
        CGAL::Epeck::FT y = (cutting_segment[i].vertex(0)-base_point_3d)*Y_axis;
        cutting_segment_in_cdt[i].first = cdt.insert(K2::Point_2(x,y));
        mp[cutting_segment_in_cdt[i].first] = cutting_segment[i].vertex(0);
        x = (cutting_segment[i].vertex(1)-base_point_3d)*X_axis;
        y = (cutting_segment[i].vertex(1)-base_point_3d)*Y_axis;
        cutting_segment_in_cdt[i].second = cdt.insert(K2::Point_2(x,y));
    }



    for(int i=0;i<sorted_bound_vertex.size();i++){
        cdt.insert_constraint(sorted_bound_vertex_in_cdt[i],sorted_bound_vertex_in_cdt[(i+1)%sorted_bound_vertex.size()]);
    }


    for(int i=0;i<cutting_segment.size();i++){
        cdt.insert_constraint(cutting_segment_in_cdt[i].first,cutting_segment_in_cdt[i].second);
    }




    return ret;

//    p.push_back(K2::Point_2(0,0));
//    p.push_back(K2::Point_2(1,0));
//    p.push_back(K2::Point_2(2,1));
//    p.push_back(K2::Point_2(1,2));
//    p.push_back(K2::Point_2(0,2));
//    cout << CGAL::ON_ORIENTED_BOUNDARY << endl;
//    cout << CGAL::POSITIVE << endl;
//    cout << CGAL::NEGATIVE << endl;
//    cout << "***************" << endl;
//    cout << CGAL::oriented_side(K2::Point_2(2, 0),p) << endl;
//    cout << CGAL::oriented_side(K2::Point_2(2, 1),p) << endl;
//    cout << CGAL::oriented_side(K2::Point_2(1, 1),p) << endl;
//    cout << CGAL::oriented_side(K2::Point_2(-1, 1),p) << endl;
//    cout << CGAL::oriented_side(K2::Point_2(0.5, 0.5),p) << endl;



    //for(int i=0;i<)




    std::cerr << "succ " << std::endl;
}



struct ApproximateField {
    vector<MeshKernel::iGameVertex> origin_vertices;
    vector<MeshKernel::iGameVertex> extend_vertices;
    vector<K2::Tetrahedron_3>tet_list;
    vector<vector<MeshKernel::iGameVertex> > outer_face;
    vector<vector<MeshKernel::iGameVertex> > inner_face;
    vector<vector<MeshKernel::iGameVertex> > side_face;
    vector<MeshKernel::iGameVertex > bound_face_vertex;
    vector<vector<int> > bound_face_id;
    K2::Point_3 center;

   // vector<vector<MeshKernel::iGameVertex> > bounded_face;
   // vector<K2::Triangle_3 > bounded_face_k2;
    MeshKernel::iGameFaceHandle fh;
    ApproximateField(){}
    ApproximateField(MeshKernel::iGameFaceHandle fh) {
        this->fh=fh;
        origin_vertices.push_back(mesh->fast_iGameVertex[mesh->fast_iGameFace[fh].vh(0)]);
        origin_vertices.push_back(mesh->fast_iGameVertex[mesh->fast_iGameFace[fh].vh(1)]);
        origin_vertices.push_back(mesh->fast_iGameVertex[mesh->fast_iGameFace[fh].vh(2)]);
        //TODO 注意分解
        MeshKernel::iGameVertex normal = ((origin_vertices[1] - origin_vertices[0]) % (origin_vertices[2] - origin_vertices[0])).normalize();
        for(int i=0;i<3;i++){
            //extend_vertices.push_back(field_move_vertex[mesh->fast_iGameFace[fh].vh(i)]);
            //if(!in_triangle_positive_side(this->fh,field_move_vertex[mesh->fast_iGameFace[fh].vh(i)])){
                extend_vertices.push_back(field_move_vertex[mesh->fast_iGameFace[fh].vh(i)]);
            //}
           /* else{
                extend_vertices.push_back(origin_vertices[i] + normal * mesh->fast_iGameFace[fh].move_dist);
            }*/
        }
        MeshKernel::iGameVertex new_normal = ((extend_vertices[1] - extend_vertices[0]) % (extend_vertices[2] - extend_vertices[0])).normalize();
        if(new_normal * normal <0)
            swap(extend_vertices[0],extend_vertices[1]);

        for(int i=0;i<3;i++) {
            tet_list.emplace_back(iGameVertex_to_Point_K2(origin_vertices[0]),
                                  iGameVertex_to_Point_K2(origin_vertices[1]),
                                  iGameVertex_to_Point_K2(origin_vertices[2]),
                                  iGameVertex_to_Point_K2(extend_vertices[i]));

        }
        for(int i=0;i<3;i++) {
            tet_list.emplace_back(iGameVertex_to_Point_K2(extend_vertices[0]),
                                  iGameVertex_to_Point_K2(extend_vertices[1]),
                                  iGameVertex_to_Point_K2(extend_vertices[2]),
                                  iGameVertex_to_Point_K2(origin_vertices[i]));
        }
        std::vector<K2::Triangle_3 >face_list;
//        for(int i=0;i<3;i++){
//            auto iv0 = origin_vertices[i];
//            auto iv1 = origin_vertices[(i+1)%3];
//            for(int j=0;j<3;j++) {
//                auto ov = extend_vertices[j];
//                K2::Triangle_3 this_tri(iGameVertex_to_Point_K2(iv0),iGameVertex_to_Point_K2(ov),iGameVertex_to_Point_K2(iv1));
//                face_list.push_back(this_tri);
//            }
//            iv0 = extend_vertices[i];
//            iv1 = extend_vertices[(i+2)%3];
//            for(int j=0;j<3;j++) {
//                auto ov = extend_vertices[j];
//                K2::Triangle_3 this_tri(iGameVertex_to_Point_K2(iv0),iGameVertex_to_Point_K2(ov),iGameVertex_to_Point_K2(iv1));
//                face_list.push_back(this_tri);
//            }
//        }
        vector<vector<int> >side_face_generate;
        side_face_generate.push_back({0,1,2});
        side_face_generate.push_back({0+3,2+3,1+3});

        for(int i=0;i<3;i++)
            bound_face_vertex.push_back(origin_vertices[i]);
        for(int i=0;i<3;i++)
            bound_face_vertex.push_back(extend_vertices[i]);
        MeshKernel::iGameVertex c(0,0,0);
        for(int i=0;i<3;i++){
            c+= origin_vertices[i];
            c+= extend_vertices[i];
        }
        c/=6;
        center = iGameVertex_to_Point_K2(c);


        int cnt = 0;
        for(int i=0;i<3;i++) {
            auto iv0 = origin_vertices[(i+1)%3];
            auto iv1 = origin_vertices[i];
            for(int j=0;j<3;j++){
                auto ov = extend_vertices[j];
                vector<MeshKernel::iGameVertex> new_face{iv0,ov,iv1};
                K2::Triangle_3 this_tri(iGameVertex_to_Point_K2(iv0),iGameVertex_to_Point_K2(ov),
                                        iGameVertex_to_Point_K2(iv1));

                K::Plane_3 this_plane_K1(iGameVertex_to_Point(iv0),iGameVertex_to_Point(ov),
                                          iGameVertex_to_Point(iv1));
                bool positive_side = false;
                bool negative_side = false;
                for(int k=0;k<tet_list.size();k++){
                    positive_side |= this_tri.supporting_plane().has_on_positive_side(centroid(tet_list[k]));
                    negative_side |= this_tri.supporting_plane().has_on_negative_side(centroid(tet_list[k]));
                }
                if( positive_side ^ negative_side){
                    bool flag = true;
                    for(auto k : side_face_generate)
                    {

                        K::Point_3 other_center = iGameVertex_to_Point((bound_face_vertex[k[0]] + bound_face_vertex[k[1]] + bound_face_vertex[k[2]]) / 3);
                        if(CGAL::squared_distance(other_center,this_plane_K1) < myeps){
                            K2::Point_3 other_v0 = iGameVertex_to_Point_K2(bound_face_vertex[k[0]]);
                            K2::Point_3 other_v1 = iGameVertex_to_Point_K2(bound_face_vertex[k[1]]);
                            K2::Point_3 other_v2 = iGameVertex_to_Point_K2(bound_face_vertex[k[2]]);
                            K2::Triangle_3 other_tri(other_v0,other_v1,other_v2);

                            CGAL::cpp11::result_of<K2::Intersect_3(K2::Triangle_3 , K2::Triangle_3)>::type
                                    res_tt = intersection(this_tri,other_tri);
                            if (res_tt) {
                                if (const K2::Point_3 *p = boost::get<K2::Point_3>(&*res_tt)) {
                                    continue;
                                }
                                if (const K2::Segment_3 *p = boost::get<K2::Segment_3>(&*res_tt)) {
                                    continue;
                                }
                                flag = false;
                                break;
                            }
                        }
                    }
                    if(flag){
                        side_face_generate.push_back({(i+1)%3,i,j+3});
                    }
                    //side_face_generate.push_back(new_face);
                }
            }
        }



        for(int i=0;i<3;i++){
            auto iv0 = extend_vertices[i];
            auto iv1 = extend_vertices[(i+1)%3];
            for(int j=0;j<3;j++){
                auto ov = origin_vertices[j];
                vector<MeshKernel::iGameVertex> new_face{iv0,ov,iv1};
                K2::Triangle_3 this_tri(iGameVertex_to_Point_K2(iv0),iGameVertex_to_Point_K2(ov),
                                        iGameVertex_to_Point_K2(iv1));
                K::Plane_3 this_plane_K1(iGameVertex_to_Point(iv0),iGameVertex_to_Point(ov),
                                          iGameVertex_to_Point(iv1));
                bool positive_side = false;
                bool negative_side = false;
                for(int k=0;k<tet_list.size();k++){
                    positive_side |= this_tri.supporting_plane().has_on_positive_side(centroid(tet_list[k]));
                    negative_side |= this_tri.supporting_plane().has_on_negative_side(centroid(tet_list[k]));
                }
                if( positive_side ^ negative_side){
                    bool flag = true;
                    for(auto k : side_face_generate)
                    {
                        K::Point_3 other_center = iGameVertex_to_Point((bound_face_vertex[k[0]] + bound_face_vertex[k[1]] + bound_face_vertex[k[2]]) / 3);
                        if(CGAL::squared_distance(other_center,this_plane_K1) < myeps){
                            K2::Point_3 other_v0 = iGameVertex_to_Point_K2(bound_face_vertex[k[0]]);
                            K2::Point_3 other_v1 = iGameVertex_to_Point_K2(bound_face_vertex[k[1]]);
                            K2::Point_3 other_v2 = iGameVertex_to_Point_K2(bound_face_vertex[k[2]]);
                            K2::Triangle_3 other_tri(other_v0,other_v1,other_v2);

                            CGAL::cpp11::result_of<K2::Intersect_3(K2::Triangle_3 , K2::Triangle_3)>::type
                                    res_tt = intersection(this_tri,other_tri);
                            if (res_tt) {
                                if (const K2::Point_3 *p = boost::get<K2::Point_3>(&*res_tt)) {
                                    continue;
                                }
                                if (const K2::Segment_3 *p = boost::get<K2::Segment_3>(&*res_tt)) {
                                    continue;
                                }
                                flag = false;
                                break;
                            }
                        }
                    }
                    if(flag){
                        side_face_generate.push_back({i+3,(i+1)%3+3,j});
                    }
                    //side_face_generate.push_back(new_face);
                }
            }
        }



        //int xx=1;

        for(int i=0;i<side_face_generate.size();i++) {
//            if(i<2){
//                bound_face_id.push_back({side_face_generate[i][0],side_face_generate[i][1],side_face_generate[i][2]});
//                continue;
//            }


            K2::Point_3 this_v0 = iGameVertex_to_Point_K2(bound_face_vertex[side_face_generate[i][0]]);
            K2::Point_3 this_v1 = iGameVertex_to_Point_K2(bound_face_vertex[side_face_generate[i][1]]);
            K2::Point_3 this_v2 = iGameVertex_to_Point_K2(bound_face_vertex[side_face_generate[i][2]]);
            K2::Triangle_3 this_tri(this_v0,this_v1,this_v2);
            K2::Ray_3 ray(CGAL::centroid(this_tri),CGAL::centroid(this_tri) + this_tri.supporting_plane().orthogonal_vector());
            K2::Ray_3 ray2(CGAL::centroid(this_tri),CGAL::centroid(this_tri) - this_tri.supporting_plane().orthogonal_vector());
//            auto v0 = Point_K2_to_iGameVertex(CGAL::centroid(this_tri));
//            auto v1 = Point_K2_to_iGameVertex(CGAL::centroid(this_tri) + this_tri.supporting_plane().orthogonal_vector());
//            auto v2 = Point_K2_to_iGameVertex(CGAL::centroid(this_tri) - this_tri.supporting_plane().orthogonal_vector());
//            cout << "v "<< v0.x() <<" "<< v0.y()<<" "<< v0.z() << endl;
//            cout << "v "<< v1.x() <<" "<< v1.y()<<" "<< v1.z() << endl;
//            cout << "v "<< v2.x() <<" "<< v2.y()<<" "<< v2.z() << endl;
//            cout << "l "<< xx <<" "<< xx+1 << endl;
//            cout << "l "<< xx <<" "<< xx+2 << endl;

            bool flag1 = true;
            bool flag2 = true;

            int yy= 0;
            for(int j=0;j<side_face_generate.size();j++) {
                if(i==j)continue;
                K2::Point_3 other_v0 = iGameVertex_to_Point_K2(bound_face_vertex[side_face_generate[j][0]]);
                K2::Point_3 other_v1 = iGameVertex_to_Point_K2(bound_face_vertex[side_face_generate[j][1]]);
                K2::Point_3 other_v2 = iGameVertex_to_Point_K2(bound_face_vertex[side_face_generate[j][2]]);
                K2::Triangle_3 other_tri(other_v0,other_v1,other_v2);
                CGAL::cpp11::result_of<K2::Intersect_3(K2::Triangle_3 , K2::Ray_3)>::type
                        res_tt = intersection(other_tri,ray);
                if(res_tt) {
                    flag1 = false;
                }
                CGAL::cpp11::result_of<K2::Intersect_3(K2::Triangle_3 , K2::Ray_3)>::type
                        res_tt2 = intersection(other_tri,ray2);
                if(res_tt2) {
                    flag2 = false;
                }
//                if(xx == 10){
//                    auto vvv0 = bound_face_vertex[side_face_generate[j][0]];
//                    auto vvv1 = bound_face_vertex[side_face_generate[j][1]];
//                    auto vvv2 = bound_face_vertex[side_face_generate[j][2]];
//                    printf("v %lf %lf %lf\n",vvv0.x(),vvv0.y(),vvv0.z());
//                    printf("v %lf %lf %lf\n",vvv1.x(),vvv1.y(),vvv1.z());
//                    printf("v %lf %lf %lf\n",vvv2.x(),vvv2.y(),vvv2.z());
//                    printf("f %d %d %d\n",yy+1,yy+2,yy+3);
//                    yy+=3;
//                }
            }
            if(flag1 || flag2 ){
                bound_face_id.push_back({side_face_generate[i][0],side_face_generate[i][1],side_face_generate[i][2]});
            }
            //cout << xx <<" "<< flag1 <<" "<< flag2 << endl; //10 坏了
            //xx+=3;
        }
        set<int >se;
        queue<int>q;
        q.push(0);
        while(!q.empty()){
            int id = q.front();
            q.pop();
            if(se.count(id))continue;
            se.insert(id);
            map<int,int> mp_id_vertex;
            for(int i=0;i<3;i++)
                mp_id_vertex[bound_face_id[id][i]] = i;

            for(int j=0;j<bound_face_id.size();j++){ // 找其他所有
                if(!se.count(j) && mp_id_vertex.count(bound_face_id[j][0])
                                   + mp_id_vertex.count(bound_face_id[j][1]) + mp_id_vertex.count(bound_face_id[j][2]) >= 2  ){
                    for(int k=0;k<3;k++){
                        if(mp_id_vertex.count(bound_face_id[j][k]) && mp_id_vertex.count(bound_face_id[j][(k + 1) % 3])) {
                           int in_id_k = mp_id_vertex[bound_face_id[j][k]];
                           int in_id_kp1 = mp_id_vertex[bound_face_id[j][(k+1)%3]];
                           if( (in_id_kp1 + 1) % 3 != in_id_k){
                               swap(bound_face_id[j][k],bound_face_id[j][(k+1)%3]);
                               break;
                           }
                        }
                    }
                    q.push(j);
                }
            }
        }






//
//
//        for(int i=0;i<3;i++){
//            auto iv0 = extend_vertices[i];
//            auto iv1 = extend_vertices[(i+1)%3];
//            for(int j=0;j<3;j++){
//                auto ov = origin_vertices[j];
//                vector<MeshKernel::iGameVertex> new_face{iv0,iv1,ov};
//                K2::Triangle_3 this_tri(iGameVertex_to_Point_K2(iv0),iGameVertex_to_Point_K2(iv1),
//                                        iGameVertex_to_Point_K2(ov));
//                bool positive_side = false;
//                bool negative_side = false;
//                for(int k=0;k<tet_list.size();k++){
//                    positive_side |= this_tri.supporting_plane().has_on_positive_side(centroid(tet_list[k]));
//                    negative_side |= this_tri.supporting_plane().has_on_negative_side(centroid(tet_list[k]));
//                }
//                if( positive_side ^ negative_side){
//                    bool redundancy = false;
//                    for(auto k : bounded_face){
//                        K2::Triangle_3 other_tri(iGameVertex_to_Point_K2(k[0]) , iGameVertex_to_Point_K2(k[1]), iGameVertex_to_Point_K2(k[2]));
//
//
//                        CGAL::cpp11::result_of<K2::Intersect_3(K2::Triangle_3 , K2::Triangle_3)>::type
//                                res_tt = intersection(this_tri,other_tri);
//                        vector<K2::Point_3> intersect_res;
//                        if (res_tt) {
//                            bool flag = true;
//                            if (const K2::Point_3 *p = boost::get<K2::Point_3>(&*res_tt)) {
//                                flag = false;
//                            }
//                            else if (const K2::Segment_3 *s = boost::get<K2::Segment_3>(&*res_tt)) {
//                                flag = false;
//                            }
//                            if(flag){
//                               // cout << "redundancy" << endl;
//                                redundancy = true;
//                            }
//
//                        }
//                    }
//                    if(redundancy) {
//                        cout<<"RE GG???GGG!!!?G" << std::endl;
//                    }
//
//                    redundancy = false;
//                    if(!redundancy){
//                        bounded_face.push_back(new_face);
//                    }
//                }
//            }
//        }


//        for(int i=0;i<3;i++){
//            auto iv0 = origin_vertices[i];
//            auto iv1 = origin_vertices[(i+1)%3];
//            for(int j=0;j<3;j++){
//                auto ov = extend_vertices[j];
//                vector<MeshKernel::iGameVertex> new_face{iv0,iv1,ov};
//                bounded_face.push_back(new_face);
//            }
//        }
//
//
//        for(int i=0;i<3;i++){
//            auto iv0 = extend_vertices[i];
//            auto iv1 = extend_vertices[(i+1)%3];
//            for(int j=0;j<3;j++){
//                auto ov = origin_vertices[j];
//                vector<MeshKernel::iGameVertex> new_face{iv0,iv1,ov};
//                bounded_face.push_back(new_face);
//            }
//        }



       // bounded_face.push_back({extend_vertices[0],extend_vertices[2],extend_vertices[1]});

        outer_face.push_back({extend_vertices[0],extend_vertices[2],extend_vertices[1]});
        inner_face.push_back({origin_vertices[0],origin_vertices[2],origin_vertices[1]});
       // bounded_face.push_back(outer_face[0]);

//        K2::Triangle_3 inner_face_K2(iGameVertex_to_Point_K2(extend_vertices[0]),
//                                     iGameVertex_to_Point_K2(extend_vertices[0]),
//                                     iGameVertex_to_Point_K2(extend_vertices[0])
//                                     );
//
//        for(vector<MeshKernel::iGameVertex> i : side_face){
//            K2::Triangle_3 tri_i(iGameVertex_to_Point_K2(i[0]),
//                                 iGameVertex_to_Point_K2(i[1]),
//                                 iGameVertex_to_Point_K2(i[2]));
//            if(face_cross_face(tri_i,inner_face_K2).first == 1){
//                bounded_face.push_back(i);
//            }
//        }
//
//
//
//        for(vector<MeshKernel::iGameVertex> i : bounded_face) {
//            bounded_face_k2.emplace_back(iGameVertex_to_Point_K2(i[0]),
//                                      iGameVertex_to_Point_K2(i[1]),
//                                      iGameVertex_to_Point_K2(i[2]));
//        }
       // bounded_face.push_back({origin_vertices[0],origin_vertices[1],origin_vertices[2]});


    // 6个四面体法，防止相交处理麻烦 并且用tet 来判断内外这样每一个面的偏移就是6个tet，；
    }

    bool in_field(K2::Point_3 v){
        CGAL::Polyhedron_3<K2>  poly;
        vector<K2::Point_3> vertex_list;
        std::vector<std::vector<std::size_t> > face_list;




        for(int i=0;i<6;i++){
            vertex_list.push_back(iGameVertex_to_Point_K2(bound_face_vertex[i]));
        }
        for(auto i : bound_face_id){
            face_list.push_back({std::size_t(i[0]),std::size_t(i[1]),std::size_t(i[2])});
        }


//        for(auto i : bound_face_id) {
//            K2::Triangle_3 tri(vertex_list[i[0]],vertex_list[i[1]],vertex_list[i[2]]);
//            if(CGAL::squared_distance(v,tri) < CGAL::Epick::FT(myeps)){
//                return false;
//            }
//        }

//        for(int k=0;k<bound_face_vertex.size();k++){
//            printf("v %lf %lf %lf\n",bound_face_vertex[k].x(),bound_face_vertex[k].y(),bound_face_vertex[k].z());
//        }
//        for(int k=0;k<bound_face_id.size();k++){
//            printf("f %d %d %d\n",bound_face_id[k][0]+1,
//                    bound_face_id[k][1]+1,
//                    bound_face_id[k][2]+1);
//
//        }


        PMP::polygon_soup_to_polygon_mesh(vertex_list, face_list, poly, CGAL::parameters::all_default());
        CGAL::Side_of_triangle_mesh<CGAL::Polyhedron_3<K2>, K2> inside(poly);
        CGAL::Bounded_side res = (inside)(v);
        if (res == CGAL::ON_BOUNDED_SIDE)
            return true;
        return false;
    }




    vector<K2::Point_3> get_near_v(K2::Point_3 v){
        vector<K2::Point_3> ret;
        for(auto i : bound_face_id){

            K2::Triangle_3 tri(iGameVertex_to_Point_K2(bound_face_vertex[i[0]]),
                               iGameVertex_to_Point_K2(bound_face_vertex[i[1]]),
                               iGameVertex_to_Point_K2(bound_face_vertex[i[2]]));
            if(CGAL::squared_distance(tri,v) < CGAL::Epeck::FT(myeps/1000)) {
                K2::Vector_3 vec = (v - center) * CGAL::Epeck::FT(1+myeps);
                K2::Point_3 new_p = center + vec;
                ret.push_back(new_p);
            }

        }
        return ret;
    }
    ApproximateField(MeshKernel::iGameFaceHandle fh,vector<MeshKernel::iGameVertex>moved_vertex){
        this->fh=fh;
        origin_vertices.push_back(mesh->fast_iGameVertex[mesh->fast_iGameFace[fh].vh(0)]);
        origin_vertices.push_back(mesh->fast_iGameVertex[mesh->fast_iGameFace[fh].vh(1)]);
        origin_vertices.push_back(mesh->fast_iGameVertex[mesh->fast_iGameFace[fh].vh(2)]);
        //TODO 注意分解
        MeshKernel::iGameVertex normal = ((origin_vertices[1] - origin_vertices[0]) % (origin_vertices[2] - origin_vertices[0])).normalize();
        for(int i=0;i<3;i++){
            extend_vertices.push_back(moved_vertex[i]);
        }
        MeshKernel::iGameVertex new_normal = ((extend_vertices[1] - extend_vertices[0]) % (extend_vertices[2] - extend_vertices[0])).normalize();
        if(new_normal * normal <0)
            swap(extend_vertices[0],extend_vertices[1]);

        for(int i=0;i<3;i++) {
            tet_list.emplace_back(iGameVertex_to_Point_K2(origin_vertices[0]),
                                  iGameVertex_to_Point_K2(origin_vertices[1]),
                                  iGameVertex_to_Point_K2(origin_vertices[2]),
                                  iGameVertex_to_Point_K2(extend_vertices[i]));

        }
        for(int i=0;i<3;i++) {
            tet_list.emplace_back(iGameVertex_to_Point_K2(extend_vertices[0]),
                                  iGameVertex_to_Point_K2(extend_vertices[1]),
                                  iGameVertex_to_Point_K2(extend_vertices[2]),
                                  iGameVertex_to_Point_K2(origin_vertices[i]));
        }

        // 6个四面体法，防止相交处理麻烦 并且用tet 来判断内外这样每一个面的偏移就是6个tet，；
    }
};
vector<ApproximateField>faces_approximate_field;



bool check_in_approximate_field_list_debug(const vector<K2::Tetrahedron_3>& tet_list ,K2::Point_3 v){
    //todo : 这里可以改成aabbtree 版本 6 tet 升级为一颗aabbtree ！！！！！！！！！！1  重大优化



    vector<K2::Triangle_3 >near_face;
    vector<K2::Point_3 > tet_center_list;

    MeshKernel::iGameVertex vv = Point_K2_to_iGameVertex(v);
    cout <<"K2::POINT3V  v " <<vv.x() <<" "<< vv.y() <<" "<< vv.z() << endl;
    int cnt=0;
    for (auto j: tet_list) {
        FILE *f = fopen(("../data/output" + to_string(file_id) + "_tet" + to_string(cnt)+".obj").c_str(), "w");
        int x = 0;
        for (vector<int> id: vector<vector<int> >{{0, 2, 1},{1, 3, 0},{1, 2, 3},{0, 3, 2}}) {
            K2::Triangle_3 tri(j.vertex(id[0]),j.vertex(id[1]),j.vertex(id[2]));
            if(CGAL::squared_distance(tri,v) < CGAL::Epeck::FT(myeps)) {
                K2::Point_3 center = centroid(j);
                near_face.push_back(tri);
                tet_center_list.push_back(center);
                cout <<"ne id :"<< cnt << endl;
            }
            MeshKernel::iGameVertex v0 = Point_K2_to_iGameVertex(j.vertex(id[0]));
            MeshKernel::iGameVertex v1 = Point_K2_to_iGameVertex(j.vertex(id[1]));
            MeshKernel::iGameVertex v2 = Point_K2_to_iGameVertex(j.vertex(id[2]));
            fprintf(f, "v %lf %lf %lf\n", v0.x(), v0.y(),v0.z());
            fprintf(f, "v %lf %lf %lf\n", v1.x(), v1.y(),v1.z());
            fprintf(f, "v %lf %lf %lf\n", v2.x(), v2.y(),v2.z());
            fprintf(f, "f %d %d %d\n",x+1,x+2,x+3);
            x+=3;
        }
        cnt++;
    }

    cout <<"near face :!!!" <<near_face.size() << endl;
//
//
//
//    exit(1);
//    cout <<"tet_size : " <<" "<< nearby_field_id_list.size() << endl;
//    cout << near_face.size()<< " center size "<< tet_center_list.size() << endl;
//    int xxx =0;
//    for(auto i : near_face){
//        auto v0 = Point_K2_to_iGameVertex(i.vertex(0));
//        auto v1 = Point_K2_to_iGameVertex(i.vertex(1));
//        auto v2 = Point_K2_to_iGameVertex(i.vertex(2));
//        cout << "v "<< v0.x() <<" "<< v0.y() <<" "<< v0.z() << endl;
//        cout << "v "<< v1.x() <<" "<< v1.y() <<" "<< v1.z() << endl;
//        cout << "v "<< v2.x() <<" "<< v2.y() <<" "<< v2.z() << endl;
//        cout <<"f "<< xxx+1 <<" "<< xxx+2 <<" "<< xxx+3 << endl;
//        xxx+=3;
//    }
//    cout <<"*******" << endl;
//    for(auto i : tet_center_list){
//        auto v0 = Point_K2_to_iGameVertex(i);
//        cout << "v "<< v0.x() <<" "<< v0.y() <<" "<< v0.z() << endl;
//    }

    if(!near_face.empty()) {// 这里只能ret true ，false还得下面继续判断后才能判断
        for(auto i : near_face) {
            bool negative = false;
            bool positive = false;
            for(auto j : tet_center_list) {
                if(i.supporting_plane().has_on_negative_side(j))
                    negative = true;
                else if(i.supporting_plane().has_on_positive_side(j))
                    positive = true;
            }
            cout << negative <<" *  "<< positive << endl;
            if(negative && positive)return true;
        }
    }

    for (auto j: tet_list) {
        if(j.has_on_bounded_side(v)){
            return true;
        }
    }
    cout << "false" << endl;
    return false;
}


bool check_in_approximate_tet_list(const vector<K2::Tetrahedron_3>& tet_list, K2::Point_3 v,vector<MeshKernel::iGameFaceHandle> debug_tet_list_belong_face) {
    FILE *file1 = fopen(("../data/in1tet" + to_string(file_id) + ".obj").c_str(), "w");
    FILE *file2 = fopen(("../data/in2tet" + to_string(file_id) + ".obj").c_str(), "w");
    FILE *file3 = fopen(("../data/in3tet" + to_string(file_id) + ".obj").c_str(), "w");
    FILE *file4 = fopen(("../data/in4tet" + to_string(file_id) + ".obj").c_str(), "w");

    FILE *filebound1 = fopen(("../data/inbound1tet" + to_string(file_id) + ".obj").c_str(), "w");
    FILE *filebound2 = fopen(("../data/inbound2tet" + to_string(file_id) + ".obj").c_str(), "w");
    FILE *filebound3 = fopen(("../data/inbound3tet" + to_string(file_id) + ".obj").c_str(), "w");

    int fis1id=0;
    for(auto i : tet_list){

        for (vector<int> id: vector<vector<int> >{{0, 2, 1},{1, 3, 0},{1, 2, 3},{0, 3, 2}}){
            MeshKernel::iGameVertex v0 = Point_K2_to_iGameVertex(i.vertex(id[0]));
            MeshKernel::iGameVertex v1 = Point_K2_to_iGameVertex(i.vertex(id[1]));
            MeshKernel::iGameVertex v2 = Point_K2_to_iGameVertex(i.vertex(id[2]));

            fprintf(file1, "v %lf %lf %lf\n", v0.x(), v0.y(),v0.z());
            fprintf(file1, "v %lf %lf %lf\n", v1.x(), v1.y(),v1.z());
            fprintf(file1, "v %lf %lf %lf\n", v2.x(), v2.y(),v2.z());
            fprintf(file1, "f %d %d %d\n",fis1id+1,fis1id+2,fis1id+3);
            fis1id+=3;
        }
    }


    vector<K2::Point_3> check_point_list;
    for(const K2::Tetrahedron_3& i : tet_list){
        for(vector<int> id : vector<vector<int> >{{0,2,1},{1,3,0},{1,2,3},{0,3,2}}){
            K2::Triangle_3 tri(i.vertex(id[0]),i.vertex(id[1]),i.vertex(id[2]));
            if(CGAL::squared_distance(tri,v) < CGAL::Epeck::FT(myeps/1000)) {
                K2::Point_3 center = centroid(i);
                K2::Vector_3 vec = (v - center) * CGAL::Epeck::FT(1+myeps);
                K2::Point_3 new_p = center + vec;
                check_point_list.push_back(new_p);
            }
        }
    }
    for(const K2::Tetrahedron_3& i : tet_list) {
        for(const K2::Point_3& j : check_point_list) {
            if(i.has_on_bounded_side(j)) {
                return true;
            }
        }
    }
    auto v0 = Point_K2_to_iGameVertex(v);
    fprintf(file2, "v %lf %lf %lf\n", v0.x(), v0.y(),v0.z());
    cout <<"tet_list.size() : " <<tet_list.size() << endl;
    cout << "method2 res :"<< check_in_approximate_field_list_debug(tet_list,v) <<" "<<endl;
    int x = 0;
    for(auto i : faces_approximate_field[debug_tet_list_belong_face[34]].tet_list){
        cout <<"34 bside" <<i.bounded_side(v) <<" "<< i.oriented_side(v)<< endl;
        for (vector<int> id: vector<vector<int> >{{0, 2, 1},{1, 3, 0},{1, 2, 3},{0, 3, 2}}) {
            MeshKernel::iGameVertex v0 = Point_K2_to_iGameVertex(i.vertex(id[0]));
            MeshKernel::iGameVertex v1 = Point_K2_to_iGameVertex(i.vertex(id[1]));
            MeshKernel::iGameVertex v2 = Point_K2_to_iGameVertex(i.vertex(id[2]));
            fprintf(file3, "v %lf %lf %lf\n", v0.x(), v0.y(),v0.z());
            fprintf(file3, "v %lf %lf %lf\n", v1.x(), v1.y(),v1.z());
            fprintf(file3, "v %lf %lf %lf\n", v2.x(), v2.y(),v2.z());
            fprintf(file3, "f %d %d %d\n",x+1,x+2,x+3);
            x+=3;
        }
    }
    x = 0;
    for(auto i : faces_approximate_field[debug_tet_list_belong_face[5]].tet_list){
        cout <<"5 bside" <<i.bounded_side(v) <<" "<<i.oriented_side(v)<< endl;
        for (vector<int> id: vector<vector<int> >{{0, 2, 1},{1, 3, 0},{1, 2, 3},{0, 3, 2}}) {
            MeshKernel::iGameVertex v0 = Point_K2_to_iGameVertex(i.vertex(id[0]));
            MeshKernel::iGameVertex v1 = Point_K2_to_iGameVertex(i.vertex(id[1]));
            MeshKernel::iGameVertex v2 = Point_K2_to_iGameVertex(i.vertex(id[2]));
            fprintf(file4, "v %lf %lf %lf\n", v0.x(), v0.y(),v0.z());
            fprintf(file4, "v %lf %lf %lf\n", v1.x(), v1.y(),v1.z());
            fprintf(file4, "v %lf %lf %lf\n", v2.x(), v2.y(),v2.z());
            fprintf(file4, "f %d %d %d\n",x+1,x+2,x+3);
            x+=3;
        }
    }
    x = 0;
    for(auto i : faces_approximate_field[debug_tet_list_belong_face[34]].inner_face){
        MeshKernel::iGameVertex v0 = i[0];
        MeshKernel::iGameVertex v1 = i[1];
        MeshKernel::iGameVertex v2 = i[2];
        fprintf(filebound1, "v %lf %lf %lf\n", v0.x(), v0.y(),v0.z());
        fprintf(filebound1, "v %lf %lf %lf\n", v1.x(), v1.y(),v1.z());
        fprintf(filebound1, "v %lf %lf %lf\n", v2.x(), v2.y(),v2.z());
        fprintf(filebound1, "f %d %d %d\n",x+1,x+2,x+3);
        x+=3;
    } x = 0;
    for(auto i : faces_approximate_field[debug_tet_list_belong_face[34]].side_face){
        MeshKernel::iGameVertex v0 = i[0];
        MeshKernel::iGameVertex v1 = i[1];
        MeshKernel::iGameVertex v2 = i[2];
        fprintf(filebound2, "v %lf %lf %lf\n", v0.x(), v0.y(),v0.z());
        fprintf(filebound2, "v %lf %lf %lf\n", v1.x(), v1.y(),v1.z());
        fprintf(filebound2, "v %lf %lf %lf\n", v2.x(), v2.y(),v2.z());
        fprintf(filebound2, "f %d %d %d\n",x+1,x+2,x+3);
        x+=3;
    } x = 0;
    for(auto i : faces_approximate_field[debug_tet_list_belong_face[34]].outer_face){
        MeshKernel::iGameVertex v0 = i[0];
        MeshKernel::iGameVertex v1 = i[1];
        MeshKernel::iGameVertex v2 = i[2];
        fprintf(filebound3, "v %lf %lf %lf\n", v0.x(), v0.y(),v0.z());
        fprintf(filebound3, "v %lf %lf %lf\n", v1.x(), v1.y(),v1.z());
        fprintf(filebound3, "v %lf %lf %lf\n", v2.x(), v2.y(),v2.z());
        fprintf(filebound3, "f %d %d %d\n",x+1,x+2,x+3);
        x+=3;
    }


    exit(111);

    return false;
}

//todo : 这里可以改成aabbtree 版本 6 tet 升级为一颗aabbtree ！！！！！！！！！！1  重大优化
bool check_in_approximate_field_list(const set<int> &nearby_field_id_list ,K2::Point_3 v) {
    // 检查接近
    //vector<ApproximateField>nearest_face;
    vector<K2::Point_3> check_point_list;
    check_point_list.push_back(v);
    for(auto i : nearby_field_id_list) {
        vector<K2::Point_3> points = faces_approximate_field[i].get_near_v(v);
        for(auto j : points)
            check_point_list.push_back(j);
    }
    for(auto i : nearby_field_id_list) {
        for(auto j : check_point_list){
            if(faces_approximate_field[i].in_field(j)) {
                {
                    return true;
                }
            }
        }
    }

    return false;



//        for(auto i : bounded_face_k2) {
//            if(CGAL::squared_distance(p,i) < CGAL::Epeck::FT(myeps)){
//                return true;
//            }
//        }
}



bool check_in_approximate_field_list(const set<int> &nearby_field_id_list ,K2::Triangle_3 tri){
    // 检查接近
    //vector<ApproximateField>nearest_face;
    K2::Point_3 v = CGAL::centroid(tri);
    vector<K2::Point_3> check_point_list;
    check_point_list.push_back(v);
    for(auto i : nearby_field_id_list) {
        vector<K2::Point_3> points = faces_approximate_field[i].get_near_v(v);
        for(auto j : points)
            check_point_list.push_back(j);
    }
    for(auto i : nearby_field_id_list) {
        for(auto j : check_point_list){
            if(faces_approximate_field[i].in_field(j)) {
                {
                    printf("ddd\n");
                    static int xx = 0;
                    FILE *filexx = fopen(("../data/ddd/output" + to_string(file_id) + "_dianyun"+ to_string(xx) + ".obj").c_str(), "w");
                    xx++;
                    int id = 1;
                    fprintf(filexx,"v %lf %lf %lf\n",CGAL::to_double(tri.vertex(0).x()),CGAL::to_double(tri.vertex(0).y()),CGAL::to_double(tri.vertex(0).z()));
                    fprintf(filexx,"v %lf %lf %lf\n",CGAL::to_double(tri.vertex(1).x()),CGAL::to_double(tri.vertex(1).y()),CGAL::to_double(tri.vertex(1).z()));
                    fprintf(filexx,"v %lf %lf %lf\n",CGAL::to_double(tri.vertex(2).x()),CGAL::to_double(tri.vertex(2).y()),CGAL::to_double(tri.vertex(2).z()));
                    fprintf(filexx,"f %d %d %d\n",id,id+1,id+2);
                    id =4;
                    for(int k=0;k<faces_approximate_field[i].bound_face_vertex.size();k++){
                        fprintf(filexx,"v %lf %lf %lf\n",faces_approximate_field[i].bound_face_vertex[k].x(),faces_approximate_field[i].bound_face_vertex[k].y(),faces_approximate_field[i].bound_face_vertex[k].z());
                    }
                    for(int k=0;k<faces_approximate_field[i].bound_face_id.size();k++){
                        fprintf(filexx,"f %d %d %d\n",faces_approximate_field[i].bound_face_id[k][0]+4,
                                faces_approximate_field[i].bound_face_id[k][1]+4,
                                faces_approximate_field[i].bound_face_id[k][2]+4);
                    }
                    return true;
                }
            }
        }
    }

    return false;



//        for(auto i : bounded_face_k2) {
//            if(CGAL::squared_distance(p,i) < CGAL::Epeck::FT(myeps)){
//                return true;
//            }
//        }
}











vector<MeshKernel::iGameVertex> min_move_g;
vector<MeshKernel::iGameVertex> max_move_g;
MeshKernel::iGameVertex do_quadratic_error_metric(MeshKernel::iGameVertexHandle vh){
    MeshKernel::iGameVertex v = mesh->fast_iGameVertex[vh];
    int m = mesh->FastNeighborFhOfVertex_[vh].size();
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

    double avg_move_dist=0;
    MeshKernel::iGameVertex avg_move_vertex(0,0,0);
    for (auto f : mesh->FastNeighborFhOfVertex_[vh]){
        avg_move_dist += mesh->fast_iGameFace[f].move_dist;

    }
    avg_move_dist /= (1.0*mesh->FastNeighborFhOfVertex_[vh].size());

    for (auto f : mesh->FastNeighborFhOfVertex_[vh]) {
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

        MeshKernel::iGameVertex move_max_v = v + normal * avg_move_dist * 1.15;
        MeshKernel::iGameVertex move_min_v = v + normal * avg_move_dist * 0.8;

        double d = -(normal.x() * new_v.x() + normal.y() * new_v.y() +  normal.z() * new_v.z());

        Eigen::Vector3d p(normal.x(),  normal.y(), normal.z());
        Eigen::Matrix3d A = p * p.transpose();
        Eigen::Vector3d D2(2*d*normal.x(),2*d*normal.y(),2*d*normal.z());
        Eigen::Vector3d LimitV(normal.x(),normal.y(),normal.z());

        double lower = normal.x()*move_min_v.x() + normal.y()*move_min_v.y() + normal.z()*move_min_v.z();
        double upper = normal.x()*move_max_v.x() + normal.y()*move_max_v.y() + normal.z()*move_max_v.z();



        for(int i=0;i<3;i++)
            for(int j=0;j<3;j++)
                hessian.coeffRef(i,j)+=A.coeff(i,j)*2;
        for(int i=0;i<3;i++)
            gradient.coeffRef(i) +=  D2.coeffRef(i);


        for(int i=0;i<3;i++)
            linearMatrix.coeffRef(cnt,i) = LimitV[i];
        lowerBound.coeffRef(cnt) = lower;
        upperBound.coeffRef(cnt) = upper;

        cnt++;
    }
    avg_move_vertex/= mesh->FastNeighborFhOfVertex_[vh].size();

    min_move_g[vh] =  v + avg_move_vertex;

    avg_move_vertex = v + avg_move_vertex;;
    hessian.coeffRef(0,0) += (2.0)/1000;
    hessian.coeffRef(1,1) += (2.0)/1000;
    hessian.coeffRef(2,2) += (2.0)/1000;

    gradient.coeffRef(0) -= (2.0/1000) * v.x();
    gradient.coeffRef(1) -= (2.0/1000) * v.y();
    gradient.coeffRef(2) -= (2.0/1000) * v.z();

//    MeshKernel::iGameVertex assist_vertex = avg_vertex/sum_move;

//    ErrorMat(0,0) += assist_vertex.x()/10;
//    ErrorMat(1,1) += assist_vertex.y()/10;
//    ErrorMat(2,2) += assist_vertex.z()/10;

//    Eigen::Matrix4d Q  =  ErrorMat;
//   // cout << Q << endl;
//
//    Eigen::Matrix4d A;
//    A << Q(0, 0), Q(0, 1), Q(0, 2), Q(0, 3),
//            Q(0, 1), Q(1, 1), Q(1, 2), Q(1, 3),
//            Q(0, 2), Q(1, 2), Q(2, 2), Q(2, 3),
//            0.f, 0.f, 0.f, 1.f;
//
//
//
//    if (abs(A.determinant())>myeps) {
//        Eigen::Vector4d b(0.f, 0.f, 0.f, 1.f);
//
//        //Eigen::Vector4d X = A.colPivHouseholderQr().solve(b);
//        Eigen::Vector4d X = A.inverse() * b;
//        //cout << "cost:" <<X.transpose() * Q * X << endl;
//        //cout << v.x() <<" "<< v.y()<<" "<< v.z() <<" "<<X[0]<<" "<< X[1]<<" "<< X[2] << endl;
//        return {X[0] , X[1] , X[2]};
//    }
//    else{
//        cout <<"err\n";
//        //return avg_vertex/sum_move;
//        return {1000,1000,1000};
//    }
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
        Eigen::VectorXd QPSolution = solver.getSolution();
        return {QPSolution.coeffRef(0), QPSolution.coeffRef(1), QPSolution.coeffRef(2)};
    }
    else{
        puts("use avg_move_vertex instead");
        return avg_move_vertex;
    }
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




bool vertex_in_tiny_grid(const MeshKernel::iGameVertex& small,const MeshKernel::iGameVertex& big,
                         const MeshKernel::iGameVertex& v) {
    return small.x() <= v.x()+myeps  && v.x()-myeps <= big.x() &&
    small.y() <= v.y()+myeps  && v.y()-myeps <= big.y() &&
    small.z() <= v.z()+myeps  && v.z()-myeps <= big.z();
}

K2::Point_3 point_k_to_k2(K::Point_3 p){
    return K2::Point_3 (p.x(),p.y(),p.z());
}


bool face_through_grid(const MeshKernel::iGameVertex& small, const MeshKernel::iGameVertex& big,
                       const vector<MeshKernel::iGameVertex>& v){
    for(int i=0;i<v.size();i++){
        if(vertex_in_tiny_grid(small,big,v[i]))
            return true;
    }
    K::Triangle_3 tri(K::Point_3(v[0].x(),v[0].y(),v[0].z()),
                       K::Point_3(v[1].x(),v[1].y(),v[1].z()),
                       K::Point_3(v[2].x(),v[2].y(),v[2].z()));
    vector<K::Point_3>grid_vertex(8);
    for(int i=0;i<8;i++)
        grid_vertex[i] = iGameVertex_to_Point(getTinyGridVertex(small,big,i));

    for(auto  each_container_face : container_grid_face){
        K::Triangle_3 tri1(grid_vertex[each_container_face[0]],
                           grid_vertex[each_container_face[1]],
                           grid_vertex[each_container_face[2]]);
        K::Triangle_3 tri2(grid_vertex[each_container_face[2]],
                           grid_vertex[each_container_face[3]],
                           grid_vertex[each_container_face[0]]);
        if(sqrt(squared_distance(tri1,tri)) < myeps){
            return true;
        }
        if(sqrt(squared_distance(tri2,tri)) < myeps){
            return true;
        }
    }


//    for(auto  each_container_face : container_grid_face){
//        K2::Triangle_3 tri(K2::Point_3(v[0].x(),v[0].y(),v[0].z()),
//                          K2::Point_3(v[1].x(),v[1].y(),v[1].z()),
//                          K2::Point_3(v[2].x(),v[2].y(),v[2].z()));
//        K2::Triangle_3 tri1(point_k_to_k2(grid_vertex[each_container_face[0]]),
//                            point_k_to_k2(grid_vertex[each_container_face[1]]),
//                            point_k_to_k2(grid_vertex[each_container_face[2]]));
//        K2::Triangle_3 tri2(point_k_to_k2(grid_vertex[each_container_face[1]]),
//                            point_k_to_k2(grid_vertex[each_container_face[2]]),
//                            point_k_to_k2(grid_vertex[each_container_face[3]]));
//        CGAL::cpp11::result_of<K2::Intersect_3(K2::Triangle_3, K2::Triangle_3)>::type result = intersection(tri1, tri);
//        if(result) {
//            cout <<"GGGGGG" << endl;
//        }
//        CGAL::cpp11::result_of<K2::Intersect_3(K2::Triangle_3, K2::Triangle_3)>::type result2 = intersection(tri2, tri);
//        if(result2) {
//            cout <<"GGGGGG" << endl;
//        }
//
//    }




    return false;

//    for(int i=0;i<6;i++){ // 用K1核就放宽成本，相差一点点的也算在里面，避免后期BFS误差，然后具体的用K2核心生成时修复;;;;;;
//        MeshKernel::iGameVertex v1 = getTinyGridVertex(small,big,i);
//        for(int j=0;j<DirectedGridEdge[i].size();j++){
//            MeshKernel::iGameVertex v2 = getTinyGridVertex(small,big,DirectedGridEdge[i][j]);
//            K::Segment_3 se(K::Point_3(v1.x(),v1.y(),v1.z()),
//                             K::Point_3(v2.x(),v2.y(),v2.z()));
//            if(sqrt(squared_distance(K::Point_3(v1.x(),v1.y(),v1.z()),tri)  ) < myeps ){
//                return true;
//            }
//            if(sqrt(squared_distance(K::Point_3(v2.x(),v2.y(),v2.z()),tri)  ) < myeps ){
//                return true;
//            }
//
//            CGAL::cpp11::result_of<K::Intersect_3(K::Segment_3, K::Triangle_3)>::type
//                    result = intersection(se, tri);
//            if(result)
//                return true;
//        }
//    }
    return false;
}




bool tet_through_grid(const MeshKernel::iGameVertex& small, const MeshKernel::iGameVertex& big,
                       const K2::Tetrahedron_3& tet){

    vector<MeshKernel::iGameVertex> v;
    for(int i=0;i<4;i++){
        v.push_back(Point_K2_to_iGameVertex(tet.vertex(i)));
    }
    for(vector<int> i : vector<vector<int> > {{0,1,2},{0,2,3},{1,2,3},{0,1,3}}){
        if(face_through_grid(small,big,{v[i[0]],v[i[1]],v[i[2]]}))
            return true;
    }
    return false;
}





bool face_through_grid222debug(const MeshKernel::iGameVertex& small, const MeshKernel::iGameVertex& big,
                       const vector<MeshKernel::iGameVertex>& v) {
    for(int i=0;i<v.size();i++){
        if(vertex_in_tiny_grid(small,big,v[i]))
            return true;
    }
    K::Triangle_3 tri(K::Point_3(v[0].x(),v[0].y(),v[0].z()),
                      K::Point_3(v[1].x(),v[1].y(),v[1].z()),
                      K::Point_3(v[2].x(),v[2].y(),v[2].z()));
    vector<K::Point_3>grid_vertex(8);
    for(int i=0;i<8;i++)
        grid_vertex[i] = iGameVertex_to_Point(getTinyGridVertex(small,big,i));

    for(auto  each_container_face : container_grid_face){
        K::Triangle_3 tri1(grid_vertex[each_container_face[0]],
                           grid_vertex[each_container_face[1]],
                           grid_vertex[each_container_face[2]]);
        K::Triangle_3 tri2(grid_vertex[each_container_face[1]],
                           grid_vertex[each_container_face[2]],
                           grid_vertex[each_container_face[3]]);
        printf("v %lf %lf %lf\n",tri1.vertex(0).x(),tri1.vertex(0).y(),tri1.vertex(0).z());
        printf("v %lf %lf %lf\n",tri1.vertex(1).x(),tri1.vertex(1).y(),tri1.vertex(1).z());
        printf("v %lf %lf %lf\n",tri1.vertex(2).x(),tri1.vertex(2).y(),tri1.vertex(2).z());
        cout<<"dddd1: "<< sqrt(squared_distance(tri1,tri)) << endl;
        printf("v %lf %lf %lf\n",tri2.vertex(0).x(),tri2.vertex(0).y(),tri2.vertex(0).z());
        printf("v %lf %lf %lf\n",tri2.vertex(1).x(),tri2.vertex(1).y(),tri2.vertex(1).z());
        printf("v %lf %lf %lf\n",tri2.vertex(2).x(),tri2.vertex(2).y(),tri2.vertex(2).z());
        cout<<"dddd2: "<< sqrt(squared_distance(tri2,tri)) << endl;
        if(sqrt(squared_distance(tri1,tri)) < myeps){
            return true;
        }
        if(sqrt(squared_distance(tri2,tri)) < myeps){
            return true;
        }
    }


//    for(auto  each_container_face : container_grid_face){
//        K2::Triangle_3 tri(K2::Point_3(v[0].x(),v[0].y(),v[0].z()),
//                          K2::Point_3(v[1].x(),v[1].y(),v[1].z()),
//                          K2::Point_3(v[2].x(),v[2].y(),v[2].z()));
//        K2::Triangle_3 tri1(point_k_to_k2(grid_vertex[each_container_face[0]]),
//                            point_k_to_k2(grid_vertex[each_container_face[1]]),
//                            point_k_to_k2(grid_vertex[each_container_face[2]]));
//        K2::Triangle_3 tri2(point_k_to_k2(grid_vertex[each_container_face[1]]),
//                            point_k_to_k2(grid_vertex[each_container_face[2]]),
//                            point_k_to_k2(grid_vertex[each_container_face[3]]));
//        CGAL::cpp11::result_of<K2::Intersect_3(K2::Triangle_3, K2::Triangle_3)>::type result = intersection(tri1, tri);
//        if(result) {
//            cout <<"GGGGGG" << endl;
//        }
//        CGAL::cpp11::result_of<K2::Intersect_3(K2::Triangle_3, K2::Triangle_3)>::type result2 = intersection(tri2, tri);
//        if(result2) {
//            cout <<"GGGGGG" << endl;
//        }
//
//    }




    return false;

//    for(int i=0;i<6;i++){ // 用K1核就放宽成本，相差一点点的也算在里面，避免后期BFS误差，然后具体的用K2核心生成时修复;;;;;;
//        MeshKernel::iGameVertex v1 = getTinyGridVertex(small,big,i);
//        for(int j=0;j<DirectedGridEdge[i].size();j++){
//            MeshKernel::iGameVertex v2 = getTinyGridVertex(small,big,DirectedGridEdge[i][j]);
//            K::Segment_3 se(K::Point_3(v1.x(),v1.y(),v1.z()),
//                             K::Point_3(v2.x(),v2.y(),v2.z()));
//            if(sqrt(squared_distance(K::Point_3(v1.x(),v1.y(),v1.z()),tri)  ) < myeps ){
//                return true;
//            }
//            if(sqrt(squared_distance(K::Point_3(v2.x(),v2.y(),v2.z()),tri)  ) < myeps ){
//                return true;
//            }
//
//            CGAL::cpp11::result_of<K::Intersect_3(K::Segment_3, K::Triangle_3)>::type
//                    result = intersection(se, tri);
//            if(result)
//                return true;
//        }
//    }
    return false;
}

/*
 *  N*N*N 来确定点
 *   然后保证后，通过格点表面采样法，然后连接边还原
 *    正负不管了
 */

vector<vector<int> > container_grid_dir{{-1,-1,-1},{-1,-1,0},{-1,-1,1},
                                        {-1,0,-1},{-1,0,0},{-1,0,1},
                                        {-1,1,-1},{-1,1,0},{-1,1,1},
                                        {0,-1,-1},{0,-1,0},{0,-1,1},
                                        {0,0,-1},{0,0,1},
                                        {0,1,-1},{0,1,0},{0,1,1},
                                        {1,-1,-1},{1,-1,0},{1,-1,1},
                                        {1,0,-1},{1,0,0},{1,0,1},
                                        {1,1,-1},{1,1,0},{1,1,1},};


struct SharpPoint{
    K2::Point_3 p;
    vector<int> source_face_local_id;
};

int furthest_vertex_in_vector(int x , const vector<MeshKernel::iGameVertex>& v){
    int ans = x;
    for(int i=0;i<v.size();i++){
        if( (v[i] - v[x]).norm() > (v[ans] - v[x]).norm()){
            ans = i;
        }
    }

    return ans;
}

int furthest_K2_point_in_vector(int x , const vector<K2::Point_3>& v){
    int ans = x;
    for(int i=0;i<v.size();i++){
        if( (v[i] - v[x]).squared_length() > (v[ans] - v[x]).squared_length()){
            ans = i;
        }
    }
    return ans;
}


vector<K2::Point_3> delete_same_K2_point_in_vector( const vector<K2::Point_3>& v){
    vector<K2::Point_3> ret;
    for(int i=0;i<v.size();i++){
        bool flag = false;
        for(int j=0;j<ret.size();j++){
            if( (v[i] - ret[j]).squared_length() < CGAL::Epeck::FT(myeps)){
                flag = true;
            }
        }
        if(!flag){
            ret.push_back(v[i]);
        }
    }
    return ret;
}



void delete_same_vertex(vector<MeshKernel::iGameVertex>& v){
    vector<MeshKernel::iGameVertex> tmp;
    for(auto i : v){
        bool flag = false;
        for(auto j : tmp){
            if( (i-j).norm() < myeps)
                flag = true;
        }
        if(!flag)
            tmp.push_back(i);
    }
    swap(tmp,v);
}



MeshKernel::iGameVertex sort_by_polar_order(vector<MeshKernel::iGameVertex>& v,MeshKernel::iGameVertex orthogonal_direction){
    MeshKernel::iGameVertex center(0,0,0);
    for(auto i : v){
        center = center + i;
    }
    center /= v.size();
    if(v.size() <=1)
        return center;

    function<int(double,double)> quadrant = [](double x,double y){
        if(x>0 && y > 0)return 1;
        else if(x<0 && y > 0)return 2;
        else if(x<0 && y < 0)return 3;
        return 4;
    };
    Plane_3 p;
    MeshKernel::iGameVertex x_axis = (v[0] - center).normalize();
    MeshKernel::iGameVertex y_axis = (orthogonal_direction % x_axis).normalize();
//    FILE *file11 = fopen(("../data/output" + to_string(file_id) + "_dianyun11.obj").c_str(), "w");
//    FILE *file12 = fopen(("../data/output" + to_string(file_id) + "_dianyun12.obj").c_str(), "w");
//    FILE *file13 = fopen(("../data/output" + to_string(file_id) + "_dianyun13.obj").c_str(), "w");
//    FILE *file14 = fopen(("../data/output" + to_string(file_id) + "_dianyun14.obj").c_str(), "w");
//    FILE *file15 = fopen(("../data/output" + to_string(file_id) + "_dianyun15.obj").c_str(), "w");
//    fprintf(file15,"v %lf %lf %lf\n",center.x(),center.y(),center.z());
//    fprintf(file15,"v %lf %lf %lf\n",(center + x_axis).x(),(center + x_axis).y(),(center + x_axis).z());
//    fprintf(file15,"v %lf %lf %lf\n",(center + y_axis).x(),(center + y_axis).y(),(center + y_axis).z());
//    fprintf(file15,"l 1 2\n");
//    fprintf(file15,"l 1 3\n");
//    for(int i=0;i<v.size();i++){
//        double x1 = (v[i] - center) * x_axis;
//        double y1 = (v[i] - center) * y_axis;
//        int q1 = quadrant(x1,y1);
//        if(q1 == 1){
//            fprintf(file11,"v %lf %lf %lf\n",v[i].x(),v[i].y(),v[i].z());
//        }
//        if(q1 == 2){
//            fprintf(file12,"v %lf %lf %lf\n",v[i].x(),v[i].y(),v[i].z());
//        }
//        if(q1 == 3){
//            fprintf(file13,"v %lf %lf %lf\n",v[i].x(),v[i].y(),v[i].z());
//        }
//        if(q1 == 4){
//            fprintf(file14,"v %lf %lf %lf\n",v[i].x(),v[i].y(),v[i].z());
//        }
//
//    }

    sort(v.begin(),v.end(),[&](MeshKernel::iGameVertex v1 , MeshKernel::iGameVertex v2){
        double x1 = (v1 - center) * x_axis;
        double y1 = (v1 - center) * y_axis;
        double x2 = (v2 - center) * x_axis;
        double y2 = (v2 - center) * y_axis;
        int q1 = quadrant(x1,y1);
        int q2 = quadrant(x2,y2);
        if(q1!=q2)return q1<q2;
        else
            return x1*y2 - x2*y1 > 0;

    });
//    for(auto i : v){
//        double x1 = (i - center) * x_axis;
//        double y1 = (i - center) * y_axis;
//        cout << quadrant(x1,y1) << "??\n";
//    }
    return center;
}


void sort_by_polar_order(vector<K2::Point_3>& v,MeshKernel::iGameVertex orthogonal_direction){
    MeshKernel::iGameVertex center(0,0,0);
    for(auto i : v){
        center = center + Point_K2_to_iGameVertex(i);
    }
    center /= v.size();
    if(v.size() <=1)
        return ;

    function<int(double,double)> quadrant = [](double x,double y){
        if(x>0 && y > 0)return 1;
        else if(x<0 && y > 0)return 2;
        else if(x<0 && y < 0)return 3;
        return 4;
    };
    Plane_3 p;
    MeshKernel::iGameVertex x_axis = (Point_K2_to_iGameVertex(v[0]) - center).normalize();
    MeshKernel::iGameVertex y_axis = (orthogonal_direction % x_axis).normalize();
//    FILE *file11 = fopen(("../data/output" + to_string(file_id) + "_dianyun11.obj").c_str(), "w");
//    FILE *file12 = fopen(("../data/output" + to_string(file_id) + "_dianyun12.obj").c_str(), "w");
//    FILE *file13 = fopen(("../data/output" + to_string(file_id) + "_dianyun13.obj").c_str(), "w");
//    FILE *file14 = fopen(("../data/output" + to_string(file_id) + "_dianyun14.obj").c_str(), "w");
//    FILE *file15 = fopen(("../data/output" + to_string(file_id) + "_dianyun15.obj").c_str(), "w");
//    fprintf(file15,"v %lf %lf %lf\n",center.x(),center.y(),center.z());
//    fprintf(file15,"v %lf %lf %lf\n",(center + x_axis).x(),(center + x_axis).y(),(center + x_axis).z());
//    fprintf(file15,"v %lf %lf %lf\n",(center + y_axis).x(),(center + y_axis).y(),(center + y_axis).z());
//    fprintf(file15,"l 1 2\n");
//    fprintf(file15,"l 1 3\n");
//    for(int i=0;i<v.size();i++){
//        double x1 = (v[i] - center) * x_axis;
//        double y1 = (v[i] - center) * y_axis;
//        int q1 = quadrant(x1,y1);
//        if(q1 == 1){
//            fprintf(file11,"v %lf %lf %lf\n",v[i].x(),v[i].y(),v[i].z());
//        }
//        if(q1 == 2){
//            fprintf(file12,"v %lf %lf %lf\n",v[i].x(),v[i].y(),v[i].z());
//        }
//        if(q1 == 3){
//            fprintf(file13,"v %lf %lf %lf\n",v[i].x(),v[i].y(),v[i].z());
//        }
//        if(q1 == 4){
//            fprintf(file14,"v %lf %lf %lf\n",v[i].x(),v[i].y(),v[i].z());
//        }
//
//    }

    sort(v.begin(),v.end(),[&](K2::Point_3 a, K2::Point_3 b){


        MeshKernel::iGameVertex v1 = Point_K2_to_iGameVertex(a);
        MeshKernel::iGameVertex v2 = Point_K2_to_iGameVertex(b);
        double x1 = (v1 - center) * x_axis;
        double y1 = (v1 - center) * y_axis;
        double x2 = (v2 - center) * x_axis;
        double y2 = (v2 - center) * y_axis;
        int q1 = quadrant(x1,y1);
        int q2 = quadrant(x2,y2);
        if(q1!=q2)return q1<q2;
        else
            return x1*y2 - x2*y1 > 0;

    });
//    for(auto i : v){
//        double x1 = (i - center) * x_axis;
//        double y1 = (i - center) * y_axis;
//        cout << quadrant(x1,y1) << "??\n";
//    }
    return ;
}





void de1014bug(){
    CGAL::Polygon_2<K2>p;
    p.push_back(K2::Point_2(0,0));
    p.push_back(K2::Point_2(1,0));
    p.push_back(K2::Point_2(2,1));
    p.push_back(K2::Point_2(1,2));
    p.push_back(K2::Point_2(0,2));
    cout << CGAL::ON_ORIENTED_BOUNDARY << endl;
    cout << CGAL::POSITIVE << endl;
    cout << CGAL::NEGATIVE << endl;
    cout << "***************" << endl;
    cout << CGAL::oriented_side(K2::Point_2(2, 0),p) << endl;
    cout << CGAL::oriented_side(K2::Point_2(2, 1),p) << endl;
    cout << CGAL::oriented_side(K2::Point_2(1, 1),p) << endl;
    cout << CGAL::oriented_side(K2::Point_2(-1, 1),p) << endl;
    cout << CGAL::oriented_side(K2::Point_2(0.5, 0.5),p) << endl;

}

int main() {

   // d1013bug();

//    CGAL::Epeck::FT c1 = 2;
//    CGAL::Epeck::FT c2 = 5;
//    auto c3 = c1*c2;
//    cout << CGAL::to_double(c3) << endl;
//    return 0;
    file_id = 4015;//54.489426 73.433586 -17.990524
    FILE *file = fopen(("../data/output" + to_string(file_id) + ".obj").c_str(), "w");

//    vector<FILE*> part_debug_file_list(10);
//
//    vector<int > part_id;
//    for(int i=0;i<10;i++){
//        part_debug_file_list[i] =  fopen(("../data/output" + to_string(file_id) + "_part"+ to_string(i)+".obj").c_str(), "w");
//        part_id.push_back(0 );
//    }

    FILE *file0 = fopen(("../data/output" + to_string(file_id) + "_dianyun0.obj").c_str(), "w");
    FILE *file1 = fopen(("../data/output" + to_string(file_id) + "_dianyun1.obj").c_str(), "w");
    FILE *file2 = fopen(("../data/output" + to_string(file_id) + "_dianyun2.obj").c_str(), "w");
    FILE *file3 = fopen(("../data/output" + to_string(file_id) + "_dianyun3.obj").c_str(), "w");
    FILE *file4 = fopen(("../data/output" + to_string(file_id) + "_dianyun4.obj").c_str(), "w");
    FILE *file4_5 = fopen(("../data/output" + to_string(file_id) + "_dianyun4.5.obj").c_str(), "w");
    FILE *file5 = fopen(("../data/output" + to_string(file_id) + "_dianyun5.obj").c_str(), "w");
    FILE *file6 = fopen(("../data/output" + to_string(file_id) + "_dianyun6.obj").c_str(), "w");
    FILE *file7 = fopen(("../data/output" + to_string(file_id) + "_dianyun7.obj").c_str(), "w");
    FILE *file8 = fopen(("../data/output" + to_string(file_id) + "_dianyun8.obj").c_str(), "w");
    FILE *file9 = fopen(("../data/output" + to_string(file_id) + "_dianyun9.obj").c_str(), "w");
    FILE *file10 = fopen(("../data/output" + to_string(file_id) + "_dianyun10.obj").c_str(), "w");
    FILE *file11 = fopen(("../data/output" + to_string(file_id) + "_dianyun11.obj").c_str(), "w");
    FILE *file12 = fopen(("../data/output" + to_string(file_id) + "_dianyun12.obj").c_str(), "w");
    FILE *file13 = fopen(("../data/output" + to_string(file_id) + "_dianyun13.obj").c_str(), "w");
    FILE *file14 = fopen(("../data/output" + to_string(file_id) + "_dianyun14.obj").c_str(), "w");
    FILE *file15 = fopen(("../data/output" + to_string(file_id) + "_dianyun15.obj").c_str(), "w");
    FILE *file16 = fopen(("../data/output" + to_string(file_id) + "_dianyun16.obj").c_str(), "w");

    FILE *file21 = fopen(("../data/output" + to_string(file_id) + "_dianyun21.obj").c_str(), "w");
    FILE *file22 = fopen(("../data/output" + to_string(file_id) + "_dianyun22.obj").c_str(), "w");
    FILE *file23 = fopen(("../data/output" + to_string(file_id) + "_dianyun23.obj").c_str(), "w");


    FILE *fileis1 = fopen(("../data/output" + to_string(file_id) + "_dianyunis1.obj").c_str(), "w");
    FILE *fileis2 = fopen(("../data/output" + to_string(file_id) + "_dianyunis2.obj").c_str(), "w");
    FILE *fileis3 = fopen(("../data/output" + to_string(file_id) + "_dianyunis3.obj").c_str(), "w");

    FILE *filedebug2_1 = fopen(("../data/output" + to_string(file_id) + "_dianyundebug2_1.obj").c_str(), "w");
    FILE *filedebug2_2 = fopen(("../data/output" + to_string(file_id) + "_dianyundebug2_2.obj").c_str(), "w");

    // freopen("../debugoutput.txt","w",stdout);
    default_move = 0.01;
    grid_len = 2.5;
    // 2.5 552531870234
    // 1.7 550170141470
    cout << grid_len <<endl;
    mix_factor = 0.5;
    //mesh = make_shared<MeshKernel::SurfaceMesh>(ReadObjFile("../data/debug4.obj2")); grid_len = 0.05;double default_move_dist = 0.1;

//    MeshKernel::iGameVertex v1(100,-1,-1);
//    MeshKernel::iGameVertex v2(1,-1,-1);
//    double t = ternary_search(v1,v2,MeshKernel::iGameFaceHandle(0)) ;
//    cout << t << endl;
//
//    return 0;


 //    mesh = make_shared<MeshKernel::SurfaceMesh>(ReadObjFile("../data/Armadillo.obj")); grid_len = 4.5; double default_move_dist =2;


    mesh = make_shared<MeshKernel::SurfaceMesh>(ReadObjFile("../data/test_orgv2.obj2")); grid_len = 10.5; double default_move_dist = 0.8;

//    for(int i=0;i<mesh->FaceSize();i++){
//        mesh->faces(MeshKernel::iGameFaceHandle(i)).move_dist = 0.05;
//    }

    for(int i=0;i<mesh->FaceSize();i++){
        mesh->faces(MeshKernel::iGameFaceHandle(i)).move_dist = default_move_dist;
    }

    mesh->build_fast();

    //mix(30);


//    BVH::BVH_Tree bvhTree;
//    bvhTree.buildBVH_Tree(*mesh);
//    auto res = bvhTree.getIntersection(BVH::Ray(Vector3d(1,1,1),Vector3d(-1,-1,-1)));
//    cout << res.pos[0]<<" "<< res.pos[1]<<" "<< res.pos[2] << endl;


    //只动xy
    faces_approximate_field.resize(mesh->FaceSize());
    field_move_vertex.resize(mesh->VertexSize());
    min_move_g.resize(mesh->VertexSize());
    max_move_g.resize(mesh->VertexSize());

    field_move_face.resize(mesh->FaceSize());
    field_move_K2_triangle.resize(mesh->FaceSize());

    for(int i=0;i<mesh->VertexSize();i++){
        field_move_vertex[i] = do_quadratic_error_metric(MeshKernel::iGameVertexHandle(i));
    }
    cout <<"build st "<< endl;
    std::vector <std::shared_ptr<std::thread> > build_thread_pool(thread_num);
    for(int i=0;i<thread_num;i++) {
        build_thread_pool[i] = make_shared<std::thread>([&](int now_id) {
            for(int i=0;i<mesh->FaceSize();i++) {
                if (i % thread_num != now_id)continue;
                if(i%500==0)
                    cout << "build "<< i << endl;
                faces_approximate_field[i] = ApproximateField(MeshKernel::iGameFaceHandle(i));

                //continue;

                MeshKernel::iGameVertex v0 = field_move_vertex[mesh->fast_iGameFace[i].vh(0)];
                MeshKernel::iGameVertex v1 = field_move_vertex[mesh->fast_iGameFace[i].vh(1)];
                MeshKernel::iGameVertex v2 = field_move_vertex[mesh->fast_iGameFace[i].vh(2)];

                MeshKernel::iGameVertex ov0 = mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(0)];
                MeshKernel::iGameVertex ov1 = mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(1)];
                MeshKernel::iGameVertex ov2 = mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(2)];

                MeshKernel::iGameVertex normal = (v1 - v0) % (v2 - v0);
                MeshKernel::iGameVertex normal_o = (ov1 - ov0) % (ov2 - ov0);
                if(normal * normal_o <0){
                    field_move_face[i]=vector<MeshKernel::iGameVertex>{v0,v2,v1};
                    //field_move_face[i]=vector<MeshKernel::iGameVertex>{v0,v1,v2};
                    auto center = (v0 + v1 + v2)/3;
                }
                else{
                    field_move_face[i]=vector<MeshKernel::iGameVertex>{v0,v1,v2};
                }
                field_move_K2_triangle[i] = K2::Triangle_3(iGameVertex_to_Point_K2(field_move_face[i][0]),
                                                           iGameVertex_to_Point_K2(field_move_face[i][1]),
                                                           iGameVertex_to_Point_K2(field_move_face[i][2]));


            }
        },i);
    }

    for(int i=0;i<thread_num;i++)
        build_thread_pool[i]->join();


    cout <<"build end "<< endl;
   //return 0;
//    for(int i=0;i<6;i++)
//    {
//        fprintf(file23,"v %lf %lf %lf\n",faces_approximate_field[1419].bound_face_vertex[i].x(),
//                faces_approximate_field[1419].bound_face_vertex[i].y(),
//                faces_approximate_field[1419].bound_face_vertex[i].z());
//    }
//    for(auto i : faces_approximate_field[1419].bound_face_id)
//        fprintf(file23,"f %d %d %d\n",i[0]+1,i[1]+1,i[2]+1);
//
//    exit(0);



//    for(int i=0;i<mesh->FaceSize();i++){
//        for(int j=0;j<faces_approximate_field[i].outer_face.size();j++) {
//            MeshKernel::iGameVertex v0 = faces_approximate_field[i].outer_face[j][0];
//            MeshKernel::iGameVertex v1 = faces_approximate_field[i].outer_face[j][1];
//            MeshKernel::iGameVertex v2 = faces_approximate_field[i].outer_face[j][2];
//            static int f11id=0;
//            fprintf(file11, "v %lf %lf %lf\n", v0.x(), v0.y(),v0.z());
//            fprintf(file11, "v %lf %lf %lf\n", v1.x(), v1.y(),v1.z());
//            fprintf(file11, "v %lf %lf %lf\n", v2.x(), v2.y(),v2.z());
//            fprintf(file11, "f %d %d %d\n",f11id+1,f11id+2,f11id+3);
//            f11id+=3;
//        }
//    }



    //return 0;
//    for(int i=0;i<mesh->FaceSize();i++){
//        for(int j=0;j<6;j++) {
//            for(int k=0;k<4;k++) {
//                K2::Point_3 p = faces_approximate_field[i].tet_list[j].vertex(k);
//                MeshKernel::iGameVertex v = Point_K2_to_iGameVertex(p);
//                fprintf(file8, "v %lf %lf %lf\n", v.x() , v.y(), v.z());
//            }
//            static int cid =1;
//            for(vector<int> k : vector<vector<int> >{{0,1,3},{1,2,3},{0,3,2},{0,1,2}}) {
//                //auto Point_K_to_iGameVertex(faces_approximate_field[i].tet_list[j].vertex(k[0]));
//                fprintf(file8,"f %d %d %d\n",cid+k[0],cid+k[1],cid+k[2]);
//            }
//            cid +=4;
//        }
//
//    }



//


//    int cnt=1;
//
//    for(int i=0;i<mesh->FaceSize();i++){
////        auto v0 = faces_approximate_field[i].extend_vertices[0];
////        auto v1 = faces_approximate_field[i].extend_vertices[1];
////        auto v2 = faces_approximate_field[i].extend_vertices[2];
//        //field_move_vertex[i]
//        auto v0 = field_move_face[i][0];
//        auto v1 = field_move_face[i][1];
//        auto v2 = field_move_face[i][2];
//        fprintf(file0, "v %lf %lf %lf\n", v0.x(), v0.y(),v0.z());
//        fprintf(file0, "v %lf %lf %lf\n", v1.x(), v1.y(),v1.z());
//        fprintf(file0, "v %lf %lf %lf\n", v2.x(), v2.y(),v2.z());
//        //fprintf(file0, "f %d %d %d\n", cnt, cnt+1,cnt+2);
//    }
//    for(int i=0;i<mesh->FaceSize();i++){
//        fprintf(file0, "f %d %d %d\n", cnt, cnt+1,cnt+2);
//        cnt+=3;
//
//    }
//    cnt=0;
//    for(int  i =0 ; i< mesh->VertexSize() ; i++){
//    //    if(flag) {
////            cout << "origin " <<mesh->fast_iGameVertex[i].x()<<" "<< mesh->fast_iGameVertex[i].y()<<" "<<
////                    mesh->fast_iGameVertex[i].z()<<endl;
////            cout << "updated  " <<field_move_vertex[i].x()<<" "<< field_move_vertex[i].y()<<" "<<
////                    field_move_vertex[i].z()<<endl;
////            cout <<"***************"<< endl;
//
////            cout <<"***************"<< endl;
//            fprintf(file1, "v %lf %lf %lf\n", mesh->fast_iGameVertex[i].x(), mesh->fast_iGameVertex[i].y(),
//                   mesh->fast_iGameVertex[i].z());
//            fprintf(file1, "v %lf %lf %lf\n", field_move_vertex[i].x(), field_move_vertex[i].y(),
//                    field_move_vertex[i].z());
//
////            fprintf(file0, "v %lf %lf %lf\n", min_move_g[i].x(), min_move_g[i].y(),
////                    min_move_g[i].z());
//
//            fprintf(file1,"l %d %d\n",cnt+1,cnt+2);
//          //  fprintf(file0,"l %d %d\n",cnt+1,cnt+3);
//            cnt+=2;
//
//       // }
//    }
//



    cgal_polygon = make_shared<CGALPolygon>(mesh);

    std::list <Triangle> triangles;
    unordered_map <Triangle, MeshKernel::iGameFaceHandle, triangle_hash, triangle_equal> triangle_map;
    for (int i = 0; i < mesh->FaceSize(); i++) {
        Point_3 p0(mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(0)].x(),
                   mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(0)].y(),
                   mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(0)].z());

        Point_3 p1(mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(1)].x(),
                   mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(1)].y(),
                   mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(1)].z());

        Point_3 p2(mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(2)].x(),
                   mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(2)].y(),
                   mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(2)].z());
        Triangle t(p0, p1, p2);
        triangle_map[t] = MeshKernel::iGameFaceHandle(i);
        triangles.push_back(t);
    }
    cgal_global_aabbtree = Tree(triangles.begin(), triangles.end());

    //Armadillo.obj
    // freopen("../log.txt","w",stdout);

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



    // mesh->init_plane_point();
    //  grid_len = (mesh->BBoxMax.x() - mesh->BBoxMin.x())/50;
    // grid_len = (mesh->BBoxMax.x() - mesh->BBoxMin.x())/50;
    // grid_len = min_move / 2.5;

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



//    double avg_area=0;
//    vector<bool>need_tiny_face;
//    need_tiny_face.resize(mesh->FaceSize());
//    for(auto i : mesh->allfaces()){
//        MeshKernel::iGameVertex v0 = mesh->fast_iGameVertex[i.second.vh(0)];
//        MeshKernel::iGameVertex v1 = mesh->fast_iGameVertex[i.second.vh(1)];
//        MeshKernel::iGameVertex v2 = mesh->fast_iGameVertex[i.second.vh(2)];
//        avg_area+=((v1 - v0)%(v2 - v0)).norm()/2;
//    }
//    avg_area/=mesh->FaceSize();
//    for(auto i : mesh->allfaces()){
//        MeshKernel::iGameVertex v0 = mesh->fast_iGameVertex[i.second.vh(0)];
//        MeshKernel::iGameVertex v1 = mesh->fast_iGameVertex[i.second.vh(1)];
//        MeshKernel::iGameVertex v2 = mesh->fast_iGameVertex[i.second.vh(2)];
//        if(((v1 - v0)%(v2 - v0)).norm()/2 < avg_area/10){
//            need_tiny_face[i.first] = true;
//        }
//        else
//            need_tiny_face[i.first] = false;
//    }



/*
    std::vector<std::shared_ptr<std::thread> > bfs_thread_pool(thread_num);
    for(int i=0;i<thread_num;i++){
        bfs_thread_pool[i] = make_shared<std::thread>([&](int now_id) {
            for(auto iter : build_mesh_vertex){
                if(iter.second%thread_num == now_id){
                    qef_result[iter.second] = qef_grid_pos(iter.first);
                    if(iter.second%100==0) // cut the output time usage
                        printf("succ %d/%d \n",iter.second, build_mesh_vertex.size());
                }
            }
        },i);
    }
    for(int i=0;i<thread_num;i++)
        bfs_thread_pool[i]->join();
 */

    //1.463933 19.793690 -4.492545
    MeshKernel::iGameVertex debug_v(54.489426,73.433586,-17.990524);
    grid debug_g =  vertex_to_grid(debug_v);

    cout <<"v to g :" <<debug_g.x <<" "<< debug_g.y <<" "<<debug_g.z << endl;

    auto debugcenter = (getGridVertex(debug_g,0) + getGridVertex(debug_g,7))/2;
    fprintf(file7,"v %lf %lf %lf\n",debugcenter.x(),debugcenter.y(),debugcenter.z());

    cout <<"bfs end \n" << endl;
    std::mutex bfs_mutex;
    std::vector <std::shared_ptr<std::thread> > bfs_thread_pool(thread_num);
    for(int i=0;i<thread_num;i++) {
        bfs_thread_pool[i] = make_shared<std::thread>([&](int now_id) {
            for (int face_id = 0; face_id < fsize; face_id++) {
                if(face_id%thread_num !=  now_id)continue;
                if (face_id % 1000 == 0)
                    printf("%d/%d\n", face_id, fsize);
                auto fh = make_pair(MeshKernel::iGameFaceHandle(face_id),
                                    mesh->faces(MeshKernel::iGameFaceHandle(face_id)));
                MeshKernel::iGameVertex center = (mesh->fast_iGameVertex[fh.second.vh(0)] +
                                                  mesh->fast_iGameVertex[fh.second.vh(1)] +
                                                  mesh->fast_iGameVertex[fh.second.vh(2)]) / 3;
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
                    std::unique_lock<std::mutex>lock1(bfs_mutex,std::defer_lock);
                    lock1.lock();
                    auto iter = frame_grid_mp.find(now);
                    if (iter == frame_grid_mp.end()) {
                        iter = frame_grid_mp.insert(make_pair(now, GridVertex())).first;
                    }
                    iter->second.face_list.push_back(fh.first);
                    lock1.unlock();
                    vector <grid> neighbor = get_neighbor(now);
                    for (auto j: neighbor) {
                        if (!is_visit.count(j)) {
                            double dist = cgal_vertex_triangle_dist(fh.second, getGridVertex(j, 0), mesh);
                            double move_limit = *set<double>
                                    {(field_move_vertex[fh.second.vh(0)]-mesh->fast_iGameVertex[fh.second.vh(0)]).norm(),
                                     (field_move_vertex[fh.second.vh(1)]-mesh->fast_iGameVertex[fh.second.vh(1)]).norm(),
                                     (field_move_vertex[fh.second.vh(2)]-mesh->fast_iGameVertex[fh.second.vh(2)]).norm()
                                     }.rbegin();
                            if (dist <  max(grid_len,fh.second.move_dist)*1.01 ) { //TODO : zheli youhua cheng pianyi juli de shiji jisuan
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
            if(nx>=0 && ny>=0 && nz>=0){
                for(auto z : i.second.face_list)
                    frame_grid_mp[{nx,ny,nz}].face_list.push_back(z);
            }
        }
    }




    std::vector <std::shared_ptr<std::thread> > each_grid_thread_pool(thread_num);
    for(int i=0;i<thread_num;i++) {
        each_grid_thread_pool[i] = make_shared<std::thread>([&](int now_id) {
            int each_grid_cnt =-1;
            for (auto each_grid = frame_grid_mp.begin(); each_grid != frame_grid_mp.end(); each_grid++) {
                each_grid_cnt++;
                if(each_grid_cnt % thread_num != now_id)continue;
                if(each_grid_cnt % (thread_num*1000) == now_id)
                    printf("each_grid_cnt %d/%d\n",each_grid_cnt,(int)frame_grid_mp.size());


                MeshKernel::iGameVertex grid_vertex = getGridVertex(each_grid->first, 0);
                MeshKernel::iGameFaceHandle near_face_handle =
                        triangle_map[*cgal_global_aabbtree.closest_point_and_primitive(
                                Point_3(grid_vertex.x(), grid_vertex.y(), grid_vertex.z())).second];


                each_grid->second.face_list.push_back(near_face_handle);
                for (auto i: mesh->fast_iGameFace[near_face_handle].getSortedVertexHandle()) {
                    for (auto j: mesh->FastNeighborFhOfVertex_[i])
                        each_grid->second.face_list.push_back(j);
                }
                sort(each_grid->second.face_list.begin(), each_grid->second.face_list.end());
                each_grid->second.face_list.resize(
                        unique(each_grid->second.face_list.begin(), each_grid->second.face_list.end()) -
                        each_grid->second.face_list.begin());
//                each_grid->second.grid_type = vertex_state(each_grid->second.face_list,
//                                                           getGridVertex(each_grid->first, 0));
            }
        }, i);
    }


    for(int i=0;i<thread_num;i++)
        each_grid_thread_pool[i]->join();

    cout << "each_grid_cnt succ " <<endl;

    std::function<bool(Plane_3 plane,MeshKernel::iGameVertex,double)> vertex_in_plane
    = [&](Plane_3 pla,MeshKernel::iGameVertex v,double eps){
                if(sqrt(squared_distance(pla,Point_3(v.x(),v.y(),v.z())))<eps){
                    return true;
                }
                return false;
    };// 20 10 10 5

    long long  sum_grid = 0;
    for (auto each_grid= frame_grid_mp.begin(); each_grid != frame_grid_mp.end(); each_grid++) {
       // if(!(each_grid->first.x == 26 && each_grid->first.y == 28 && each_grid->first.z == 9  ))continue;
       // if(!(each_grid->first.x == 14 && each_grid->first.y == 11 && each_grid->first.z == 24  ))continue;
        auto small  = getGridVertex(each_grid->first,0);
        auto big  = getGridVertex(each_grid->first,7);
        static int f3_id = 1;
        for (int ii = 0; ii < 7; ii++) {
            for (int jj = 0; jj < DirectedGridEdge[ii].size(); jj++) {
                int from = ii;
                int to = DirectedGridEdge[ii][jj];
                MeshKernel::iGameVertex fv = getTinyGridVertex(small, big, from);
                MeshKernel::iGameVertex tv = getTinyGridVertex(small, big, to);
                fprintf(file9, "v %lf %lf %lf\n", fv.x(), fv.y(), fv.z());
                fprintf(file9, "v %lf %lf %lf\n", tv.x(), tv.y(), tv.z());
                fprintf(file9, "l %d %d\n", f3_id, f3_id + 1);
                f3_id += 2;
            }
        }
    }

    // 上述代码完成距离场建格子的过程 8 ;
    atomic<int>sum_face_size(0);
    atomic<int>maxx_face_size(0);
    std::mutex mu;
    std::vector<std::shared_ptr<std::thread> > each_frame_thread(thread_num);
    for(int i=0;i<thread_num;i++){
        each_frame_thread[i] = make_shared<thread>([&](int id){
            int tt=0;
            for (auto each_grid = frame_grid_mp.begin(); each_grid != frame_grid_mp.end(); each_grid++){
                tt++;
                if (tt % 200 == 0) { //2520
                    cout <<id <<" "<< tt << " // " <<" "<< frame_grid_mp.size() << endl;
                }
                if(tt % thread_num != id)continue;

               // if(tt % 50 != 0 )continue;




                function<int(MeshKernel::iGameVertex, MeshKernel::iGameFaceHandle,
                             const vector<MeshKernel::iGameVertex>)> vertex_in_tet_field_state
                        = [&](MeshKernel::iGameVertex v, MeshKernel::iGameFaceHandle fh,
                              const vector<MeshKernel::iGameVertex> moved_v) {
                            K2::Triangle_3 tri(iGameVertex_to_Point_K2(moved_v[0]),
                                               iGameVertex_to_Point_K2(moved_v[1]),
                                               iGameVertex_to_Point_K2(moved_v[2]));
                            K2::Point_3 v_k2 = iGameVertex_to_Point_K2(v);
                            if (sqrt(CGAL::to_double(CGAL::squared_distance(tri, v_k2))) < myeps*10000) {
                                return on_side;
                            }

                            ApproximateField approximateField(fh, moved_v);
                            for (int each_tet_id = 0; each_tet_id < 6; each_tet_id++) {
                                if (approximateField.tet_list[each_tet_id].has_on_bounded_side(v_k2)) {
                                    return in_side;
                                }
                            }
                            return out_side;
                        };





                set < MeshKernel::iGameFaceHandle > face_set;
                vector<MeshKernel::iGameFaceHandle > face_list;
               // vector<MeshKernel::iGameFaceHandle > possible_face_list;
                vector<MeshKernel::iGameFace> vf;
                for (int j = 0; j < 8; j++) {
                    grid g = getGridFrameVertex(each_grid->first, j);
                    unordered_map<grid, GridVertex, grid_hash, grid_equal>::iterator it = frame_grid_mp.find(g);
                    if (it == frame_grid_mp.end())continue;
                    for (auto k: it->second.face_list) {
                        if (!face_set.count(k)) {
                            face_set.insert(k);
                            face_list.push_back(k);
                        }
                    }
                }

               // if(!(each_grid->first.x == 17 && each_grid->first.y == 27 && each_grid->first.z == 10  ))continue;
               // if(!(each_grid->first.x == 26 && each_grid->first.y == 28 && each_grid->first.z == 9  ))continue; // 15 22 21 // 26 31 7




//                for (MeshKernel::iGameFaceHandle i: face_list) {
//
//                    static int f5id=0;
//                    auto v0 = field_move_face[i][0];
//                    auto v1 = field_move_face[i][1];
//                    auto v2 = field_move_face[i][2];
//                    fprintf(file5, "v %lf %lf %lf\n", v0.x(), v0.y(),v0.z());
//                    fprintf(file5, "v %lf %lf %lf\n", v1.x(), v1.y(),v1.z());
//                    fprintf(file5, "v %lf %lf %lf\n", v2.x(), v2.y(),v2.z());
//                    fprintf(file5, "f %d %d %d\n",f5id+1,f5id+2,f5id+3);
//                    f5id+=3;
//                }


                function<int(MeshKernel::iGameVertex)> get_side_type = [&](MeshKernel::iGameVertex igame_p2) {
                    int side_type = 0;
                    //cout <<"******\n" << vf.size() << " "<<face_list.size() << endl;
                    //cgal_aabbtree_query(vf, igame_p2, mesh, side_type, cgal_polygon);

                    if(cgal_polygon->inMesh(igame_p2))
                        side_type = 1;
                    else
                        side_type = 0;

                    return side_type;
                };
                MeshKernel::iGameVertex small = getGridVertex(each_grid->first, 0);
                MeshKernel::iGameVertex big = getGridVertex(each_grid->first, 7);

//                if(each_grid->first.x == 14 && each_grid->first.y == 16 && each_grid->first.z == 12 ){
//                    auto center = (small + big)/2;
//                    fprintf(file6,"v %lf %lf %lf\n",center.x(),center.y(),center.z());
//                }
//                else
//                    continue;


//                int all_in = 0;
//                for (MeshKernel::iGameFaceHandle i: face_list){
//
//                    for(int j=0;j<6;j++) {
//                        faces_approximate_field[i].tet_list[j].
//                    }
//                    vector<MeshKernel::iGameVertex> v;
//
//                }


                bool all_in = false;
                for (MeshKernel::iGameFaceHandle i: face_list) {
                    vector<MeshKernel::iGameVertex> v;
                    for (int j = 0; j < 3; j++)
                        v.push_back(field_move_vertex[mesh->fast_iGameFace[i].vh(j)]);

                        int cnt_in = 0;
                        for (int j = 0; j < 8; j++) {

                            // if(getGridFrameVertex(each_grid->first,j))
                            auto vertex = getGridVertex(each_grid->first, j);
                            int state = vertex_in_tet_field_state(vertex, i, {field_move_face[i][0],
                                                                              field_move_face[i][1],
                                                                              field_move_face[i][2]});
                            if (state == in_side)
                                cnt_in++;

                        }
                        if (cnt_in == 8) {
                            all_in = true;
                        }
                }
                if(all_in)continue;





                for (auto i: face_list) {
                    vf.push_back(mesh->fast_iGameFace.at(i));
                }
              //  cout << possible_face_list.size() <<" "<<all_in << " "<< face_list.size() << endl;


                vector<vector<MeshKernel::iGameVertex> > maybe_used_face;
                vector<vector<MeshKernel::iGameVertex> > maybe_used_side_face;
                //set<int>maybe_used_face_field;
                for (MeshKernel::iGameFaceHandle i: face_list) {
                    for(auto j : faces_approximate_field[i].outer_face) {
                        if (face_through_grid(small, big, j)) {
                            maybe_used_face.push_back(j);
                            //maybe_used_face_field.insert(i);
//                            static int f21id=0;
//                            auto v0 = j[0];
//                            auto v1 = j[1];
//                            auto v2 = j[2];
//                            fprintf(file21, "v %lf %lf %lf\n", v0.x(), v0.y(),v0.z());
//                            fprintf(file21, "v %lf %lf %lf\n", v1.x(), v1.y(),v1.z());
//                            fprintf(file21, "v %lf %lf %lf\n", v2.x(), v2.y(),v2.z());
//                            fprintf(file21, "f %d %d %d\n",f21id+1,f21id+2,f21id+3);
//                            f21id+=3;
                        }
                        else{
//                            static int f22id=0;
//                            auto v0 = j[0];
//                            auto v1 = j[1];
//                            auto v2 = j[2];
//                            fprintf(file22, "v %lf %lf %lf\n", v0.x(), v0.y(),v0.z());
//                            fprintf(file22, "v %lf %lf %lf\n", v1.x(), v1.y(),v1.z());
//                            fprintf(file22, "v %lf %lf %lf\n", v2.x(), v2.y(),v2.z());
//                            fprintf(file22, "f %d %d %d\n",f22id+1,f22id+2,f22id+3);
//                            f22id+=3;
                        }
                    }
                }
                set<int> debug_face_id;

               // std::vector<K2::Tetrahedron_3 >field_tet_list;
                std::vector<MeshKernel::iGameFaceHandle> debug_tet_list_belong_face;
                set<int> field_through_set;
                for (MeshKernel::iGameFaceHandle i: face_list) {
                    for(auto j : faces_approximate_field[i].tet_list){
                        if(tet_through_grid(small,big,j)){
                            //field_tet_list.push_back(j);
                            field_through_set.insert(i);
                            break;
                            //debug_tet_list_belong_face.push_back(i);
                        }
                    }
                }
//                for(auto i : field_tet_list) {
//                    static int fis1id=0;
//                    for (vector<int> id: vector<vector<int> >{{0, 2, 1},{1, 3, 0},{1, 2, 3},{0, 3, 2}}){
//                        MeshKernel::iGameVertex v0 = Point_K2_to_iGameVertex(i.vertex(id[0]));
//                        MeshKernel::iGameVertex v1 = Point_K2_to_iGameVertex(i.vertex(id[1]));
//                        MeshKernel::iGameVertex v2 = Point_K2_to_iGameVertex(i.vertex(id[2]));
//
//                        fprintf(fileis1, "v %lf %lf %lf\n", v0.x(), v0.y(),v0.z());
//                        fprintf(fileis1, "v %lf %lf %lf\n", v1.x(), v1.y(),v1.z());
//                        fprintf(fileis1, "v %lf %lf %lf\n", v2.x(), v2.y(),v2.z());
//                        fprintf(fileis1, "f %d %d %d\n",fis1id+1,fis1id+2,fis1id+3);
//                        fis1id+=3;
//                    }
//                }

                // 旧的maybe used face
//                for (MeshKernel::iGameFaceHandle i: face_list) {
//                    //maybe_used_face_field.insert(i);
//                    //if( i == 1419)continue;
//                    for(vector<MeshKernel::iGameVertex> j : faces_approximate_field[i].side_face){
//                        K2::Triangle_3 tri1(iGameVertex_to_Point_K2(j[0]),
//                                        iGameVertex_to_Point_K2(j[1]),
//                                        iGameVertex_to_Point_K2(j[2])
//                                        );
//                        for(vector<MeshKernel::iGameVertex> k: maybe_used_face){
//                            K2::Triangle_3 tri2(iGameVertex_to_Point_K2(k[0]),
//                                                iGameVertex_to_Point_K2(k[1]),
//                                                iGameVertex_to_Point_K2(k[2])
//                            );
//                            if(face_cross_face(tri1,tri2).first == 1){
//                                if(face_through_grid(small, big, j)) {
//                                    maybe_used_side_face.push_back(j);
//
////                                    if( i == 1419) { //TODO THIS CODE IS DEBUG CODE
////                                        debug_face_id.insert( maybe_used_face.size());
////                                        cout << "debug_face_id??" << maybe_used_face.size() << endl;
////                                    }
//                                    maybe_used_face.push_back(j);
//                                }
//                                break;
//                            }
//                        }
//                    }
//                }
                //新的maybe used face
                for (MeshKernel::iGameFaceHandle i: face_list){
                    for(auto j : faces_approximate_field[i].bound_face_id){
                        vector<MeshKernel::iGameVertex> tmp{faces_approximate_field[i].bound_face_vertex[j[0]],
                                                        faces_approximate_field[i].bound_face_vertex[j[1]],
                                                        faces_approximate_field[i].bound_face_vertex[j[2]]};
                        if(j[0] >= 3 || j[1] >= 3 || j[2] >= 3) {
                            if(face_through_grid(small, big, tmp))
                                maybe_used_face.push_back(tmp);
                        }
                    }
                }





                // 9.22 删除 疑似和 上面的代码可替代
//                for(auto j : maybe_used_side_face){
//                    if( i == 1419) {
//                        debug_face_id = maybe_used_face.size();
//                        cout << "debug_face_id??" << debug_face_id << endl;
//                    }
//                    maybe_used_face.push_back(j);
//                }



                // todo: 只需要判断内层面和外层面相交是不是就能做




//                sort(maybe_used_face_field.begin(),maybe_used_face_field.end());
//                maybe_used_face_field.resize(unique(maybe_used_face_field.begin(),maybe_used_face_field.end())
//                - maybe_used_face_field.begin());

//                for(auto j : maybe_used_face) {
//                    static int f6id=0;
//                    auto v0 = j[0];
//                    auto v1 = j[1];
//                    auto v2 = j[2];
//                    fprintf(file6, "v %lf %lf %lf\n", v0.x(), v0.y(),v0.z());
//                    fprintf(file6, "v %lf %lf %lf\n", v1.x(), v1.y(),v1.z());
//                    fprintf(file6, "v %lf %lf %lf\n", v2.x(), v2.y(),v2.z());
//                    fprintf(file6, "f %d %d %d\n",f6id+1,f6id+2,f6id+3);
//                    f6id+=3;
//                }


//                for (int i: maybe_used_face_field) {
//                    for(int j =0; j<faces_approximate_field[i].side_face.size(); j++) {
//                        static int f5id=0;
//                        auto v0 = faces_approximate_field[i].side_face[j][0];
//                        auto v1 = faces_approximate_field[i].side_face[j][1];
//                        auto v2 = faces_approximate_field[i].side_face[j][2];
//                        MeshKernel::iGameVertex debug_v(-0.853005,-3.295282,48.552765);
//
//                        K::Triangle_3 tri(iGameVertex_to_Point(v0),iGameVertex_to_Point(v1),iGameVertex_to_Point(v2));
//
//                        if(CGAL::squared_distance(iGameVertex_to_Point(debug_v),tri) < myeps){
//                            cout <<"*****" << i<<" "<<j<< endl;
//                        }
//
//                        fprintf(file5, "v %lf %lf %lf\n", v0.x(), v0.y(),v0.z());
//                        fprintf(file5, "v %lf %lf %lf\n", v1.x(), v1.y(),v1.z());
//                        fprintf(file5, "v %lf %lf %lf\n", v2.x(), v2.y(),v2.z());
//                        fprintf(file5, "f %d %d %d\n",f5id+1,f5id+2,f5id+3);
//                        f5id+=3;
//                    }
//                    for(auto j : faces_approximate_field[i].outer_face) {
//                        static int f4id=0;
//                        auto v0 = j[0];
//                        auto v1 = j[1];
//                        auto v2 = j[2];
//                        fprintf(file4, "v %lf %lf %lf\n", v0.x(), v0.y(),v0.z());
//                        fprintf(file4, "v %lf %lf %lf\n", v1.x(), v1.y(),v1.z());
//                        fprintf(file4, "v %lf %lf %lf\n", v2.x(), v2.y(),v2.z());
//                        fprintf(file4, "f %d %d %d\n",f4id+1,f4id+2,f4id+3);
//                        f4id+=3;
//                    }
//                    for(auto j : faces_approximate_field[i].inner_face) {
//                        static int f3id=0;
//                        auto v0 = j[0];
//                        auto v1 = j[1];
//                        auto v2 = j[2];
//                        fprintf(file3, "v %lf %lf %lf\n", v0.x(), v0.y(),v0.z());
//                        fprintf(file3, "v %lf %lf %lf\n", v1.x(), v1.y(),v1.z());
//                        fprintf(file3, "v %lf %lf %lf\n", v2.x(), v2.y(),v2.z());
//                        fprintf(file3, "f %d %d %d\n",f3id+1,f3id+2,f3id+3);
//                        f3id+=3;
//                    }
//                }



              //  cout <<"********444****" << endl;
                std::function<vector<K2::Point_3 >(
                        MeshKernel::iGameVertex, MeshKernel::iGameVertex, MeshKernel::iGameVertex,
                        MeshKernel::iGameVertex, MeshKernel::iGameVertex)> get_grid_intersect_triangle
                        = [&](MeshKernel::iGameVertex v0, MeshKernel::iGameVertex v1, MeshKernel::iGameVertex v2,
                              MeshKernel::iGameVertex tinygrid_st, MeshKernel::iGameVertex tinygrid_end) {
                            vector<K2::Point_3 > ret;
                            vector<MeshKernel::iGameVertex> vs;

                            for (int j = 0; j < 8; j++) {
                                vs.push_back(getTinyGridVertex(tinygrid_st, tinygrid_end, j));
                            }
                            // 这里采用大三角形法： ;
                            K2::Triangle_3 pla(iGameVertex_to_Point_K2(v0),
                                               iGameVertex_to_Point_K2(v1),
                                               iGameVertex_to_Point_K2(v2));
                            for (auto i: container_grid_face) {
                                vector<K2::Point_3> cross_vertex;
                                vector<K2::Point_3> type3_cross_vertex;
                                for (auto j: vector<vector<int> >{{0, 1, 2},
                                                                  {2, 3, 0}}) {
                                    K2::Triangle_3 grid_tri(iGameVertex_to_Point_K2(vs[i[j[0]]]),
                                                            iGameVertex_to_Point_K2(vs[i[j[1]]]),
                                                            iGameVertex_to_Point_K2(vs[i[j[2]]]));

                                    const auto result = intersection(pla, grid_tri);
                                    if (result) {
                                        if (const K2::Segment_3 *s = boost::get<K2::Segment_3>(&*result)) {

//                                    MeshKernel::iGameVertex cross_product1 = (v1_f - v0) % (v2_f - v0);
//                                    MeshKernel::iGameVertex cross_product2 = (v1 - v0) % (v2 - v0) ;
//                                    if(cross_product1 * cross_product2 < 0)
//                                        swap(v1_f,v2_f);

                                            cross_vertex.push_back(s->vertex(0));
                                            cross_vertex.push_back(s->vertex(1));

                                        }//todo : fix with point cut 这里还有bug !!!!!!!!!!!1; 如果点面交怎么办;;;
                                        else if(const K2::Triangle_3 *t = boost::get<K2::Triangle_3>(&*result)){
                                            type3_cross_vertex.push_back(t->vertex(0));
                                            type3_cross_vertex.push_back(t->vertex(1));
                                            type3_cross_vertex.push_back(t->vertex(2));
                                        }
                                        else if(const std::vector<K2::Point_3> *v = boost::get<std::vector<K2::Point_3>>(&*result)){
                                            for(auto ii : *v)
                                                type3_cross_vertex.push_back(ii);
                                        }
                                    }
                                }
                                if(!type3_cross_vertex.empty()){
                                    return type3_cross_vertex;
                                }
                                if (cross_vertex.size() >= 2) {
                                    int id1 = furthest_K2_point_in_vector(0, cross_vertex);
                                    int id2 = furthest_K2_point_in_vector(id1, cross_vertex);
                                    K2::Point_3 v1_f = cross_vertex[id1];
                                    K2::Point_3 v2_f = cross_vertex[id2];
//                            MeshKernel::iGameVertex cross_product1 = (v1_f - v0) % (v2_f - v0);此处不需要;
//                            MeshKernel::iGameVertex cross_product2 = (v1 - v0) % (v2 - v0) ;
//                            if(cross_product1 * cross_product2 < 0)
//                                swap(v1_f,v2_f);
                                    ret.push_back(v1_f);
                                    ret.push_back(v2_f);
                                } else if (cross_vertex.size()) {
                                    ret.push_back(*cross_vertex.begin());
                                }


                            }


                            //getTinyGridVertex
                            // container_grid_face


                            return ret;
                        };

                vector<K2::Triangle_3 >possible_face_inner_part;
               // cout <<"********1****" << endl;
                static int vid = 1;
                int xx=0;
                vector<vector<K2::Triangle_3> > face_inner_grid_part(maybe_used_face.size());
                vector<vector<K2::Point_3> > face_inner_grid_polygon(maybe_used_face.size());

                // 下面这段是和边框切割的代码
                for(int maybe_used_face_id=0; maybe_used_face_id < maybe_used_face.size(); maybe_used_face_id++) {

                    vector<K2::Point_3 > res = get_grid_intersect_triangle(maybe_used_face[maybe_used_face_id][0],maybe_used_face[maybe_used_face_id][1],maybe_used_face[maybe_used_face_id][2],
                                                                           small,big);
                    for(int j=0;j<3;j++){
                        MeshKernel::iGameVertex v = maybe_used_face[maybe_used_face_id][j];
                        if(vertex_in_tiny_grid(small,big,v) ){
                            res.push_back(iGameVertex_to_Point_K2(v));
                        }
                    }
                   // cout <<"************" << endl;
                    res = delete_same_K2_point_in_vector(res);

//                    for(auto j: res) {
//                        auto tt = Point_K2_to_iGameVertex(j);
//                        fprintf(file7, "v %lf %lf %lf\n",tt.x(),tt.y(),tt.z());
//                    }
//                    continue;

                    MeshKernel::iGameVertex center(0,0,0);
                    vector<K2::Point_3 > triangle_vertex_list;
                    for(auto j: res){
                        center += Point_K2_to_iGameVertex(j);
                        triangle_vertex_list.push_back(j);
                    }
                    center /= res.size();
                    auto origin_direct = (maybe_used_face[maybe_used_face_id][1] - maybe_used_face[maybe_used_face_id][0]) %
                                         (maybe_used_face[maybe_used_face_id][2] - maybe_used_face[maybe_used_face_id][0]);

                    sort_by_polar_order(
                            triangle_vertex_list, origin_direct);
                    //*****************************************************
                    // 修改triangle list
                    face_inner_grid_polygon[maybe_used_face_id] = triangle_vertex_list;



                    if(triangle_vertex_list.size() <2 ) {
                        continue;
                    }
                    MeshKernel::iGameVertex v0 = Point_K2_to_iGameVertex(triangle_vertex_list[0]);
                    for (int i = 0; i < triangle_vertex_list.size()-2; i++) {
                        MeshKernel::iGameVertex v1 = Point_K2_to_iGameVertex(triangle_vertex_list[i + 1]);
                        MeshKernel::iGameVertex v2 = Point_K2_to_iGameVertex(triangle_vertex_list[i + 2]);
                        auto triangle_vertex_list_i1 = triangle_vertex_list[i + 1];
                        auto triangle_vertex_list_i2 = triangle_vertex_list[i + 2];
                        MeshKernel::iGameVertex now_direct = (v1 - v0) % (v2 - v0);
                        if (now_direct * origin_direct < 0) {
                            swap(v1, v2);
                            swap(triangle_vertex_list_i1, triangle_vertex_list_i2);
                        }//10 11 15
                        if(!(each_grid->first.x == 17 && each_grid->first.y == 28 && each_grid->first.z == 10 ))continue;

                        face_inner_grid_part[maybe_used_face_id].push_back(K2::Triangle_3(triangle_vertex_list[0],
                                                                                          triangle_vertex_list_i1,
                                                                                          triangle_vertex_list_i2));

//                        std::unique_lock<std::mutex> lock(mu);
//                        static int xxx = 0;
//                        fprintf(file6, "v %lf %lf %lf\n", v0.x(), v0.y(), v0.z());
//                        fprintf(file6, "v %lf %lf %lf\n", v1.x(), v1.y(), v1.z());
//                        fprintf(file6, "v %lf %lf %lf\n", v2.x(), v2.y(), v2.z());
//                        fprintf(file6, "f %d %d %d\n", xxx + 1, xxx + 2, xxx + 3);
//                        xxx += 3;

//                        if(i<10){
//                            std::unique_lock<std::mutex>lock(mu);
//                            fprintf(part_debug_file_list[i],"v %lf %lf %lf\n",v0.x(),v0.y(),v0.z());
//                            fprintf(part_debug_file_list[i],"v %lf %lf %lf\n",v1.x(),v1.y(),v1.z());
//                            fprintf(part_debug_file_list[i],"v %lf %lf %lf\n",v2.x(),v2.y(),v2.z());
//                            fprintf(part_debug_file_list[i],"f %d %d %d\n",part_id[i]+1,part_id[i]+2,part_id[i]+3);
//                            part_id[i]+=3;
//
//                        }
                    }

                    if(!(each_grid->first.x == 17 && each_grid->first.y == 28 && each_grid->first.z == 10  ))continue;


                    //todo : debug code 2

//                    static int file6id = 0;
//                    fprintf(file6,"v %lf %lf %lf\n", field_move_face[fh][0].x(),
//                            field_move_face[fh][0].y(),
//                            field_move_face[fh][0].z());
//                    fprintf(file6,"v %lf %lf %lf\n", field_move_face[fh][1].x(),
//                            field_move_face[fh][1].y(),
//                            field_move_face[fh][1].z());
//                    fprintf(file6,"v %lf %lf %lf\n", field_move_face[fh][2].x(),
//                            field_move_face[fh][2].y(),
//                            field_move_face[fh][2].z());
//                    fprintf(file6,"f %d %d %d\n",file6id+1,file6id+2,file6id+3);
//                    file6id+=3;



                    //todo : debug code and recover before



//                    if(triangle_vertex_list.size() <3 )continue;
//                    MeshKernel::iGameVertex cc(0,0,0);
//                    for(int iii=0;iii<triangle_vertex_list.size();iii++)
//                        cc+=Point_K2_to_iGameVertex(triangle_vertex_list[iii]);
//                    cc /= triangle_vertex_list.size();
//
//                    for (int i = 0; i < triangle_vertex_list.size(); i++) {
//                        MeshKernel::iGameVertex v1 = Point_K2_to_iGameVertex(triangle_vertex_list[i]);
//                        MeshKernel::iGameVertex v2 = Point_K2_to_iGameVertex(triangle_vertex_list[(i+1) % triangle_vertex_list.size()]);
//                        MeshKernel::iGameVertex now_direct = (cc - v1) % (cc - v2);
//                        if (now_direct * origin_direct < 0)
//                            swap(v1, v2);
//                        static int file5id = 0;
//                        std::unique_lock<std::mutex>lock(mu);
//                        fprintf(file4,"v %lf %lf %lf\n",cc.x(),cc.y(),cc.z());
//                        fprintf(file4,"v %lf %lf %lf\n",v1.x(),v1.y(),v1.z());
//                        fprintf(file4,"v %lf %lf %lf\n",v2.x(),v2.y(),v2.z());
//                        fprintf(file4,"f %d %d %d\n",file5id+1,file5id+2,file5id+3);
//                        file5id+=3;
//                    }







                    //*****************************************************

//                    for(int tt1 = 0;tt1<  face_inner_grid_part[maybe_used_face_id].size(); tt1++){
//                        for(int tt2 = tt1+1;tt2<  face_inner_grid_part[maybe_used_face_id].size(); tt2++){
//                            auto x1 = face_inner_grid_part[maybe_used_face_id][tt1];
//                            int cnt=0;
//                            for(int kk=0;kk<3;kk++) {
//                                auto x2 = face_inner_grid_part[maybe_used_face_id][tt2].vertex(kk);
//                                if(CGAL::squared_distance(x1,x2) <CGAL::Epeck::FT(myeps)){
//                                    cnt++;
//                                }
//                            }
//                            if(cnt==3){
//                                std::unique_lock<std::mutex>lock(mu);
//                                static int vid5=1;
//                                auto ver10 = Point_K2_to_iGameVertex(face_inner_grid_part[maybe_used_face_id][tt1].vertex(0));
//                                auto ver11 = Point_K2_to_iGameVertex(face_inner_grid_part[maybe_used_face_id][tt1].vertex(1));
//                                auto ver12 = Point_K2_to_iGameVertex(face_inner_grid_part[maybe_used_face_id][tt1].vertex(2));
//                                auto ver20 = Point_K2_to_iGameVertex(face_inner_grid_part[maybe_used_face_id][tt2].vertex(0));
//                                auto ver21 = Point_K2_to_iGameVertex(face_inner_grid_part[maybe_used_face_id][tt2].vertex(1));
//                                auto ver22 = Point_K2_to_iGameVertex(face_inner_grid_part[maybe_used_face_id][tt2].vertex(2));
//                                fprintf(file5,"v %lf %lf %lf\n" ,ver10.x(),ver10.y(),ver10.z());
//                                fprintf(file5,"v %lf %lf %lf\n" ,ver11.x(),ver11.y(),ver11.z());
//                                fprintf(file5,"v %lf %lf %lf\n" ,ver12.x(),ver12.y(),ver12.z());
//                                fprintf(file5,"v %lf %lf %lf\n" ,ver20.x(),ver20.y(),ver20.z());
//                                fprintf(file5,"v %lf %lf %lf\n" ,ver21.x(),ver21.y(),ver21.z());
//                                fprintf(file5,"v %lf %lf %lf\n" ,ver22.x(),ver22.y(),ver22.z());
//                                fprintf(file5,"f %d %d %d\n",vid5,vid5+1,vid5+2);
//                                fprintf(file5,"f %d %d %d\n",vid5+3,vid5+4,vid5+5);
//                                vid5+=6;
//                                cout <<"GGGG EBBEBEBEBE" << endl;
//                            }
//
//                        }
//                    }








//                    for(int i=0;i<first_divided_face.size();i++){
//                        for(int j=0;j<first_divided_face.size();j++)
//                    }
/*
2.5 8
963992
 */

//                    static int vid = 1;
//
//                    std::unique_lock<std::mutex>lock(mu);
//                    fprintf(file10, "v %lf %lf %lf\n", center.x(), center.y(),
//                            center.z());
//                    int cid = vid;
//                    vid++;
//                    for (int i = 0; i < triangle_vertex_list.size(); i++) {
//
//                        MeshKernel::iGameVertex v1 = triangle_vertex_list[i];
//                        MeshKernel::iGameVertex v2 = triangle_vertex_list[(i + 1) % triangle_vertex_list.size()];
//                        MeshKernel::iGameVertex now_direct = (v1 - center) % (v2 - center);
//                        if (now_direct * origin_direct < 0)
//                            swap(v1, v2);
//                        fprintf(file10, "v %lf %lf %lf \n", v1.x(), v1.y(), v1.z());
//                        fprintf(file10, "v %lf %lf %lf \n", v2.x(), v2.y(), v2.z());
//                        fprintf(file10, "f %d %d %d \n", cid, vid, vid + 1);
//
//                        vid += 2;
//                    }
//                    lock.unlock();
//                    continue;

                }
                // 理论上面的都是对的

               // vector<vector<K2::Triangle_3> >second_div_face(maybe_used_face.size());



//                if(maybe_used_face.size() > 10){
//                    cout << maybe_used_face.size() << endl;
//                    std::unique_lock<std::mutex>lock(mu);
//                    static int vid = 1;
//                    for(auto i : maybe_used_face){
//                        auto v0 = i[0];
//                        auto v1 = i[1];
//                        auto v2 = i[2];
//
//                        fprintf(file10, "v %lf %lf %lf \n", v0.x(), v0.y(), v0.z());
//                        fprintf(file10, "v %lf %lf %lf \n", v1.x(), v1.y(), v1.z());
//                        fprintf(file10, "v %lf %lf %lf \n", v2.x(), v2.y(), v2.z());
//                        fprintf(file10, "f %d %d %d \n", vid, vid +1 , vid + 2);
//                        vid+=3;
//                    }
//                }
//                continue;

                //todo ： 这里加个continue 看看效率
                for(int i=0;i<maybe_used_face.size();i++) { // 处理单片面的切割


                    list<K2::Triangle_3 >now_tri_list;
                    K2::Segment_3 e0(iGameVertex_to_Point_K2(maybe_used_face[i][0]),
                                     iGameVertex_to_Point_K2(maybe_used_face[i][1]));
                    K2::Segment_3 e1(iGameVertex_to_Point_K2(maybe_used_face[i][1]),
                                     iGameVertex_to_Point_K2(maybe_used_face[i][2]));
                    K2::Segment_3 e2(iGameVertex_to_Point_K2(maybe_used_face[i][2]),
                                     iGameVertex_to_Point_K2(maybe_used_face[i][0]));
                    vector<pair<K2::Segment_3,K2::Triangle_3> >cutting_segment;
                    K2::Triangle_3  tri_i(iGameVertex_to_Point_K2(maybe_used_face[i][0]),
                                          iGameVertex_to_Point_K2(maybe_used_face[i][1]),
                                          iGameVertex_to_Point_K2(maybe_used_face[i][2]));
                    //下面这个for 计算出所有的交线  直接用原来的面去切，不要用inner part 切 ，加快效率
                    for(int j=0;j<maybe_used_face.size();j++) {
                        if(i==j)continue;
                        K2::Triangle_3  tri_j(iGameVertex_to_Point_K2(maybe_used_face[j][0]),
                                              iGameVertex_to_Point_K2(maybe_used_face[j][1]),
                                              iGameVertex_to_Point_K2(maybe_used_face[j][2]));
                        CGAL::cpp11::result_of<K2::Intersect_3(K2::Triangle_3, K2::Triangle_3)>::type
                                res = intersection(tri_i,tri_j);
                        if (res) {
                            if (const K2::Segment_3 *s = boost::get<K2::Segment_3>(&*res)) {
                                bool same_edge = false;
                                for(K2::Segment_3 edge : {e0,e1,e2}){
                                    if(segment_in_line(edge,*s))
                                        same_edge = true;
                                }
                                if(!same_edge){
                                    cutting_segment.emplace_back(*s,tri_j);
//                                    std::unique_lock<std::mutex>lock(mu);
//
//                                    static int xxx=1;
//                                    auto vv0 = Point_K2_to_iGameVertex(s->vertex(0));
//                                    auto vv1 = Point_K2_to_iGameVertex(s->vertex(1));
//
//                                    fprintf(file9,"v %lf %lf %lf\n",vv0.x(),vv0.y(),vv0.z());
//                                    fprintf(file9,"v %lf %lf %lf\n",vv1.x(),vv1.y(),vv1.z());
//                                    fprintf(file9,"l %d %d\n",xxx,xxx+1);
//                                    {
//                                        CGAL::cpp11::result_of<K2::Intersect_3(K2::Triangle_3, K2::Line_3)>::type
//                                                res22 = intersection(field_move_K2_triangle[possible_face_list[i]],
//                                                                     s->supporting_line());
//                                        if(res22){
//                                            if(const K2::Segment_3 *ss2 = boost::get<K2::Segment_3>(&*res22)){
//                                                auto ss20 = Point_K2_to_iGameVertex(ss2->vertex(0));
//                                                auto ss21 = Point_K2_to_iGameVertex(ss2->vertex(1));
//                                                fprintf(file7,"v %lf %lf %lf\n",ss20.x(),ss20.y(),ss20.z());
//                                                fprintf(file7,"v %lf %lf %lf\n",ss21.x(),ss21.y(),ss21.z());
//                                                fprintf(file7,"l %d %d\n",xxx,xxx+1);
//                                            }
//                                            else
//                                                cout <<"GGG" << endl;
//                                        }
//                                        else
//                                            cout <<"GGG" << endl;
//                                    }
//                                    xxx+=2;
                                }
                                //判断交线是不是边6

                            }//面交再说！；；；；
                        }

                    }
                    if(face_inner_grid_polygon[i].size() <3)continue;

                    for(int j=0;j<face_inner_grid_part[i].size();j++) {
                        /******debug******/
                        static int f12id = 0;
                        MeshKernel::iGameVertex v0 = Point_K2_to_iGameVertex(face_inner_grid_part[i][j].vertex(0));
                        MeshKernel::iGameVertex v1 = Point_K2_to_iGameVertex(face_inner_grid_part[i][j].vertex(1));
                        MeshKernel::iGameVertex v2 = Point_K2_to_iGameVertex(face_inner_grid_part[i][j].vertex(2));

                        fprintf(file12,"v %lf %lf %lf\n",v0.x(),v0.y(),v0.z());
                        fprintf(file12,"v %lf %lf %lf\n",v1.x(),v1.y(),v1.z());
                        fprintf(file12,"v %lf %lf %lf\n",v2.x(),v2.y(),v2.z());
                        fprintf(file12,"f %d %d %d\n",f12id+1,f12id+2,f12id+3);
                        f12id+=3;
                        now_tri_list.push_back(face_inner_grid_part[i][j]);
                        /*************/
                    }

                   // now_tri_list.push_back(field_move_K2_triangle[possible_face_list[i]]);
                    MeshKernel::iGameVertex direct = (maybe_used_face[i][1] - maybe_used_face[i][0]) %
                            (maybe_used_face[i][2] - maybe_used_face[i][0]);

                    vector<K2::Segment_3>cs;
                    for(auto j : cutting_segment){
                        cs.push_back(j.first);
                    }

                    //下面这块用德劳内代替

                  //  cout<<face_inner_grid_part[i].size() <<" *** "<< face_inner_grid_polygon[i].size() << endl;

                    K2::Vector_3 origin_normal = CGAL::cross_product(tri_i.vertex(1) -tri_i.vertex(0),
                                                                     tri_i.vertex(2) -tri_i.vertex(0));

                    vector<vector<K2::Point_3> > cdt_res = CGAL_CDT(face_inner_grid_polygon[i],cs,tri_i);

                    static int file15id = 0;

                    for(auto j : cdt_res){
                        auto v0 = Point_K2_to_iGameVertex(j[0]);
                        auto v1 = Point_K2_to_iGameVertex(j[1]);
                        auto v2 = Point_K2_to_iGameVertex(j[2]);
                        //cout << "F15ISFACE"<<CGAL::to_double(CGAL::squared_distance(j[0],tri_i.supporting_plane())) << std::endl;
                        //cout << "F15ISFACE"<< CGAL::to_double(CGAL::squared_distance(j[1],tri_i.supporting_plane())) << std::endl;
                       // cout << "F15ISFACE"<< CGAL::to_double(CGAL::squared_distance(j[2],tri_i.supporting_plane())) << std::endl;
                        fprintf(file15, "v %lf %lf %lf\n", v0.x(), v0.y(), v0.z());
                        fprintf(file15, "v %lf %lf %lf\n", v1.x(), v1.y(), v1.z());
                        fprintf(file15, "v %lf %lf %lf\n", v2.x(), v2.y(), v2.z());
                        fprintf(file15, "f %d %d %d\n", file15id + 1, file15id + 2, file15id + 3);
                        file15id += 3;
                    }
                    //cout  <<"** CDT: "<< cdt_res.size() << endl;




// 在这里启动老式切割:::

//                    for(int j=0;j<cutting_segment.size();j++) {
//
//                        vector<K2::Triangle_3>new_triangles;
//
//                        for (list<K2::Triangle_3>::iterator it = now_tri_list.begin();it != now_tri_list.end(); ) {
//                            //break;
//                            vector<K2::Point_3>on_side_vertices;
//                            vector<K2::Point_3>positive_side_vertices;
//                            vector<K2::Point_3>negative_side_vertices;
//
//
//                            for(int k=0;k<3;k++){
//                                if(cutting_segment[j].second.supporting_plane().has_on(it->vertex(k))){
//                                    on_side_vertices.push_back(it->vertex(k));
//                                }
//                                else if(cutting_segment[j].second.supporting_plane().has_on_positive_side(it->vertex(k))){
//                                    positive_side_vertices.push_back(it->vertex(k));
//                                }
//                                else {
//                                    negative_side_vertices.push_back(it->vertex(k));
//                                }
//                            }
//
//                            if(!positive_side_vertices.empty() && !negative_side_vertices.empty()){
//                                K2::Line_3 line(cutting_segment[j].first.vertex(0),cutting_segment[j].first.vertex(1));
//
//                                CGAL::cpp11::result_of<K2::Intersect_3(K2::Triangle_3, K2::Line_3)>::type
//                                        res = intersection(*it,line);
//                                if (res) {
//                                    if (const K2::Segment_3 *s = boost::get<K2::Segment_3>(&*res)) {
//
//
//                                        K2::Point_3 intersect_v0 = s->vertex(0);
//                                        K2::Point_3 intersect_v1 = s->vertex(1);
//                                        for(const vector<K2::Point_3>& v_side: {positive_side_vertices,negative_side_vertices}){
//                                       // for(const vector<K2::Point_3>& v_side: {positive_side_vertices}){
//                                            if(v_side.size() == 1){
//                                                //if(CGAL::cross_product(it->supporting_plane().orthogonal_vector()),)
//
//
//                                                new_triangles.emplace_back(intersect_v0,intersect_v1,v_side[0]);
//                                                auto it2 = new_triangles.rbegin();
//                                                if(it->supporting_plane().orthogonal_vector()*it2->supporting_plane().orthogonal_vector() < CGAL::Epeck::FT(0)){
//                                                    *it2 = K2::Triangle_3 (it2->vertex(0),it2->vertex(2),it2->vertex(1));
//                                                }
//                                            }
//                                            else if(v_side.size() == 2){
//                                                vector<K2::Point_3>vec{intersect_v0,intersect_v1,v_side[0],v_side[1]};
//                                                sort_by_polar_order(vec,direct);
//                                                new_triangles.emplace_back(vec[0],vec[1],vec[2]);
//                                                auto it2 = new_triangles.rbegin();
//                                                if(it->supporting_plane().orthogonal_vector()*it2->supporting_plane().orthogonal_vector() < CGAL::Epeck::FT(0)){
//                                                    *it2 = K2::Triangle_3 (it2->vertex(0),it2->vertex(2),it2->vertex(1));
//                                                }
//                                                new_triangles.emplace_back(vec[0],vec[2],vec[3]);
//                                                it2 = new_triangles.rbegin();
//                                                if(it->supporting_plane().orthogonal_vector()*it2->supporting_plane().orthogonal_vector() < CGAL::Epeck::FT(0)){
//                                                    *it2 = K2::Triangle_3 (it2->vertex(0),it2->vertex(2),it2->vertex(1));
//                                                }
//                                            }
//                                        }
//                                        it = now_tri_list.erase(it);
//
////                                        std::unique_lock<std::mutex>lock(mu);
////
////
////
////                                        static int xxx=1;
////                                        auto vv0 = Point_K2_to_iGameVertex(s->vertex(0));
////                                        auto vv1 = Point_K2_to_iGameVertex(s->vertex(1));
////                                        fprintf(file8,"v %lf %lf %lf\n",vv0.x(),vv0.y(),vv0.z());
////                                        fprintf(file8,"v %lf %lf %lf\n",vv1.x(),vv1.y(),vv1.z());
////                                        fprintf(file8,"l %d %d\n",xxx,xxx+1);
////                                        xxx+=2;
//
//                                        continue;
//                                    }
//                                }
//                            }
//                            it++;
//                        }
//                        for(auto tmp:new_triangles) {
//
//
//                            now_tri_list.push_back(tmp);
//                        }
//                    }
//                    cout <<"** now_tri_list: "<< now_tri_list.size() << endl;


                    static int vid = 1;
                   // cout << i <<" "<< now_tri_list.size() << endl;
//                   if(i!=15)continue;
//                   else{
//                       for(auto tri: now_tri_list){
//                           static int f12id = 0;
//                           auto v0 = Point_K2_to_iGameVertex(tri.vertex(0));
//                           auto v1 = Point_K2_to_iGameVertex(tri.vertex(1));
//                           auto v2 = Point_K2_to_iGameVertex(tri.vertex(2));
//
//                           fprintf(file12, "v %lf %lf %lf\n", v0.x(), v0.y(),v0.z());
//                           fprintf(file12, "v %lf %lf %lf\n", v1.x(), v1.y(),v1.z());
//                           fprintf(file12, "v %lf %lf %lf\n", v2.x(), v2.y(),v2.z());
//                           fprintf(file12, "f %d %d %d\n",  f12id+1,f12id+2,f12id+3);
//                           f12id+=3;
//                       }
//                   }


//FIXME: 在这里启用cdt
                    //*********************

                    now_tri_list.clear();

                    for(int j=0;j<cdt_res.size();j++){
                        K2::Vector_3 this_normal = CGAL::cross_product(cdt_res[j][1] - cdt_res[j][0],
                                                                       cdt_res[j][2] - cdt_res[j][0]);
                        if(this_normal * origin_normal < CGAL::Epeck::FT(0)){
                            swap(cdt_res[j][1],cdt_res[j][2]);
                        }

                        now_tri_list.emplace_back(cdt_res[j][0],cdt_res[j][1],cdt_res[j][2]);
                    }

                    //*********************




                   for(auto tri: now_tri_list) {
                        // todo : open this code this code is check weather is need generate

                        //todo :   注意这里没有和外部面判断相交进行裁切 说不定有问题!!!!!!!!!；
/*****************************/
                       // cout <<"start check "<< endl;

                       // bool flag = check_in_approximate_field_list(maybe_used_face_field ,CGAL::centroid(tri));;
                        bool flag = check_in_approximate_field_list(field_through_set ,CGAL::centroid(tri));

                            // 这个地方是不是可以去掉底相交呢？？？？？
                            if(flag == 0) {
                                static int f13id = 0;
                                auto v0 = Point_K2_to_iGameVertex(tri.vertex(0));
                                auto v1 = Point_K2_to_iGameVertex(tri.vertex(1));
                                auto v2 = Point_K2_to_iGameVertex(tri.vertex(2));

                                fprintf(file13, "v %lf %lf %lf\n", v0.x(), v0.y(), v0.z());
                                fprintf(file13, "v %lf %lf %lf\n", v1.x(), v1.y(), v1.z());
                                fprintf(file13, "v %lf %lf %lf\n", v2.x(), v2.y(), v2.z());
                                fprintf(file13, "f %d %d %d\n", f13id + 1, f13id + 2, f13id + 3);
                                f13id += 3;
                            }
                            else {
                                static int f14id = 0;
                                auto v0 = Point_K2_to_iGameVertex(tri.vertex(0));
                                auto v1 = Point_K2_to_iGameVertex(tri.vertex(1));
                                auto v2 = Point_K2_to_iGameVertex(tri.vertex(2));

                                fprintf(file14, "v %lf %lf %lf\n", v0.x(), v0.y(), v0.z());
                                fprintf(file14, "v %lf %lf %lf\n", v1.x(), v1.y(), v1.z());
                                fprintf(file14, "v %lf %lf %lf\n", v2.x(), v2.y(), v2.z());
                                fprintf(file14, "f %d %d %d\n", f14id + 1, f14id + 2, f14id + 3);
                                f14id += 3;
                            }

//                        if(!flag){
//
//
//                            bool fff = check_in_approximate_field_list_debug(maybe_used_face_field ,CGAL::centroid(tri));;
//                            if(!fff) {
//                                static int f13id = 0;
//                                auto v0 = Point_K2_to_iGameVertex(tri.vertex(0));
//                                auto v1 = Point_K2_to_iGameVertex(tri.vertex(1));
//                                auto v2 = Point_K2_to_iGameVertex(tri.vertex(2));
//
//                                fprintf(file13, "v %lf %lf %lf\n", v0.x(), v0.y(), v0.z());
//                                fprintf(file13, "v %lf %lf %lf\n", v1.x(), v1.y(), v1.z());
//                                fprintf(file13, "v %lf %lf %lf\n", v2.x(), v2.y(), v2.z());
//                                fprintf(file13, "f %d %d %d\n", f13id + 1, f13id + 2, f13id + 3);
//                                f13id += 3;
//                            }
//                        }


                        /*if(debug_face_id.count(i)) {
                            cout << i <<" ???? "<< flag << endl;
                        };*/

                        /*for(int j=0;j<possible_face_list.size();j++){
                            for (int each_tet_id = 0; each_tet_id < 6 && !flag; each_tet_id++) {

                                K2::Tetrahedron_3 tet = faces_approximate_field[possible_face_list[j]].tet_list[each_tet_id];
                                CGAL::Epeck::FT dist = CGAL::squared_distance(CGAL::centroid(tri),field_move_K2_triangle[possible_face_list[j]]);
                                std::unique_lock<std::mutex>lock(mu);
                                static int f11id = 0;
                                fprintf(file11,"v %lf %lf %lf\n", field_move_face[possible_face_list[j]][0].x(),field_move_face[possible_face_list[j]][0].y(),field_move_face[possible_face_list[j]][0].z());
                                fprintf(file11,"v %lf %lf %lf\n", field_move_face[possible_face_list[j]][1].x(),field_move_face[possible_face_list[j]][1].y(),field_move_face[possible_face_list[j]][1].z());
                                fprintf(file11,"v %lf %lf %lf\n", field_move_face[possible_face_list[j]][2].x(),field_move_face[possible_face_list[j]][2].y(),field_move_face[possible_face_list[j]][2].z());
                                fprintf(file11, "f %d %d %d\n", f11id +1 , f11id + 2, f11id +3);
                                f11id+=3;


                                for(int k=0;k<4;k++) {
                                    K2::Point_3 p = tet.vertex(k);
                                    MeshKernel::iGameVertex v = Point_K2_to_iGameVertex(p);
                                    fprintf(file2, "v %lf %lf %lf\n", v.x() , v.y(), v.z());
                                }
                                static int cid =1;
                                for(vector<int> k : vector<vector<int> >{{0,1,3},{1,2,3},{0,3,2},{0,1,2}}) {
                                    //auto Point_K_to_iGameVertex(faces_approximate_field[i].tet_list[j].vertex(k[0]));
                                    fprintf(file2,"f %d %d %d\n",cid+k[0],cid+k[1],cid+k[2]);
                                }
                                cid +=4;

                                if(tet.has_on_bounded_side(CGAL::centroid(tri)) && dist > CGAL::Epeck::FT(myeps)){
                                    flag = true;

                                    // 这里话四面体 找bug 0914 ！！！！！！！！！！！！！

//                                    for(int k=0;k<4;k++) {
//                                            K2::Point_3 p = tet.vertex(k);
//                                            MeshKernel::iGameVertex v = Point_K2_to_iGameVertex(p);
//                                            fprintf(file2, "v %lf %lf %lf\n", v.x() , v.y(), v.z());
//                                        }
//                                        static int cid =1;
//                                        for(vector<int> k : vector<vector<int> >{{0,1,3},{1,2,3},{0,3,2},{0,1,2}}) {
//                                            //auto Point_K_to_iGameVertex(faces_approximate_field[i].tet_list[j].vertex(k[0]));
//                                            fprintf(file2,"f %d %d %d\n",cid+k[0],cid+k[1],cid+k[2]);
//                                        }
//                                        cid +=4;

                                }
//                                else if(tet.has_on_bounded_side(tri.vertex(0))){
//                                    flag = true;
//                                }
//                                else if(tet.has_on_bounded_side(tri.vertex(1))){
//                                    flag = true;
//                                }
//                                else if(tet.has_on_bounded_side(tri.vertex(2))){
//                                    flag = true;
//                                }
//                                CGAL::cpp11::result_of<K2::Intersect_3(K2::Tetrahedron_3 , K2::Triangle_3)>::type
//                                        res_tet_tri = intersection(tet,
//                                                           tri);
//                                vector<K2::Point_3> intersect_res;
//                                if (res_tet_tri) {
//                                    if (const K2::Segment_3 *s = boost::get<K2::Segment_3>(&*res_tet_tri)) {
//                                        intersect_res.push_back(s->vertex(0));
//                                        intersect_res.push_back(s->vertex(1));
//
//                                    }
//                                    if (const K2::Triangle_3 *t = boost::get<K2::Triangle_3>(&*res_tet_tri)) {
//                                        intersect_res.push_back(t->vertex(0));
//                                        intersect_res.push_back(t->vertex(1));
//                                        intersect_res.push_back(t->vertex(2));
//                                    }
//                                    if (const std::vector <K2::Point_3> *v = boost::get<std::vector <K2::Point_3> >(&*res_tet_tri)) {
//                                        for(const auto& v_i: *v)
//                                            intersect_res.push_back(v_i);
//                                    }
//                                }
//                                K2::Triangle_3 moved_face = field_move_K2_triangle[possible_face_list[j]];
//
//                                for(auto item : intersect_res) {
//                                    if (CGAL::squared_distance(item, moved_face) > CGAL::Epeck::FT(myeps)) {
//                                        flag = true;
//                                    }
//                                }

                            }
                        }*/
                        if(!flag) {
                            int side_type = //get_side_type(Point_K2_to_iGameVertex(CGAL::centroid(tri)));
                                    cgal_polygon->inMesh(Point_K2_to_iGameVertex(CGAL::centroid(tri)));
                            if(side_type !=1) {
                                flag = true;
                                cout <<"GGST" << endl;
                            }
                        }



                        if(flag) {

 //                           static int dvid = 1;


//                            std::unique_lock<std::mutex>lock(mu);
//
//
//
//
//ø
//
//
//                            fprintf(file4_5,"v %lf %lf %lf\n",vvv0.x(),vvv0.y(),vvv0.z());
//                            fprintf(file4_5,"v %lf %lf %lf\n",vvv1.x(),vvv1.y(),vvv1.z());
//                            fprintf(file4_5,"v %lf %lf %lf\n",vvv2.x(),vvv2.y(),vvv2.z());
//                            fprintf(file4_5,"f %d %d %d\n",dvid,dvid+1,dvid+2);
//                            dvid+=3;
                          //  printf("GGG4455\n");
                            continue;
                        }



/*****************************/





                        auto vvv0 = Point_K2_to_iGameVertex(tri.vertex(0));
                        auto vvv1 = Point_K2_to_iGameVertex(tri.vertex(1));
                        auto vvv2 = Point_K2_to_iGameVertex(tri.vertex(2));
                        auto ddd = (vvv1 - vvv0) % (vvv2- vvv0);
                        //if(ddd*direct <0)
                           // swap(vvv1,vvv2);
                        std::unique_lock<std::mutex>lock(mu);
                        fprintf(file10,"v %lf %lf %lf\n",vvv0.x(),vvv0.y(),vvv0.z());
                        fprintf(file10,"v %lf %lf %lf\n",vvv2.x(),vvv2.y(),vvv2.z());
                        fprintf(file10,"v %lf %lf %lf\n",vvv1.x(),vvv1.y(),vvv1.z());
                        fprintf(file10,"f %d %d %d\n",vid,vid+1,vid+2);
                        vid+=3;
                    }





//                    for(auto it1 = now_tri_list.begin();it1!=now_tri_list.end();it1++){
//
//                        for(auto it2 = it1;it2!=now_tri_list.end();it2++){
//                            bool flag = false;
//                            if(it1 != it2){
//                                int cnt = 0;
//                                cnt += CGAL::squared_distance(it1->vertex(0),*it2) < CGAL::Epeck::FT(myeps);
//                                cnt += CGAL::squared_distance(it1->vertex(1),*it2) < CGAL::Epeck::FT(myeps);
//                                cnt += CGAL::squared_distance(it1->vertex(2),*it2) < CGAL::Epeck::FT(myeps);
//                                if(cnt ==3)
//                                    flag = true;
//                            }
//                            if(flag){
//
//                                static int vvvvid = 1;
//                                std::unique_lock<std::mutex>lock(mu);
//                                auto vvvvv0 = Point_K2_to_iGameVertex(it1->vertex(0));
//                                auto vvvvv1 = Point_K2_to_iGameVertex(it1->vertex(1));
//                                auto vvvvv2 = Point_K2_to_iGameVertex(it1->vertex(2));
//                                auto vvvvv3 = Point_K2_to_iGameVertex(it2->vertex(0));
//                                auto vvvvv4 = Point_K2_to_iGameVertex(it2->vertex(1));
//                                auto vvvvv5 = Point_K2_to_iGameVertex(it2->vertex(2));
//
//                                fprintf(file2,"v %lf %lf %lf\n",vvvvv0.x(),vvvvv0.y(),vvvvv0.z());
//                                fprintf(file2,"v %lf %lf %lf\n",vvvvv1.x(),vvvvv1.y(),vvvvv1.z());
//                                fprintf(file2,"v %lf %lf %lf\n",vvvvv2.x(),vvvvv2.y(),vvvvv2.z());
//                                fprintf(file2,"v %lf %lf %lf\n",vvvvv3.x(),vvvvv3.y(),vvvvv3.z());
//                                fprintf(file2,"v %lf %lf %lf\n",vvvvv4.x(),vvvvv4.y(),vvvvv4.z());
//                                fprintf(file2,"v %lf %lf %lf\n",vvvvv5.x(),vvvvv5.y(),vvvvv5.z());
//                                fprintf(file2,"f %d %d %d\n",vvvvid,vvvvid+1,vvvvid+2);
//                                fprintf(file2,"f %d %d %d\n",vvvvid+3,vvvvid+4,vvvvid+5);
//                                vvvvid+=6;
//                            }
//                        }
//                    }
//
//
//                    for(auto it1 = now_tri_list.begin();it1!=now_tri_list.end();it1++){
//                        K2::Segment_3 e0(it1->vertex(0),
//                                         it1->vertex(1));
//                        K2::Segment_3 e1(it1->vertex(1),
//                                         it1->vertex(2));
//                        K2::Segment_3 e2(it1->vertex(2),
//                                         it1->vertex(0));
//                        for(auto it2 = it1;it2!=now_tri_list.end();it2++){
//                            bool flag = false;
//                            MeshKernel::iGameVertex tv1,tv2;
//                            if(it1 != it2){
//                                CGAL::cpp11::result_of<K2::Intersect_3(K2::Triangle_3 , K2::Triangle_3)>::type
//                                        res_tt = intersection(*it1,*it2);
//                                vector<K2::Point_3> intersect_res;
//
//                                if (res_tt) {
//
//                                    if (const K2::Segment_3 *s = boost::get<K2::Segment_3>(&*res_tt)) {
//                                        bool same_edge = false;
//                                        for(auto item  : {e0,e1,e2})
//                                            if(segment_in_line(*s,item))
//                                                same_edge = true;
//                                        tv1 = Point_K2_to_iGameVertex(s->vertex(0));
//                                        tv2 = Point_K2_to_iGameVertex(s->vertex(1));
//                                        if(!same_edge)
//                                            flag = true;
//                                    }
//
//                                }
//                            }
//                            if(flag){
//                                cout << "GGF" << endl;
//                                static int vvvvid = 1;
//                                std::unique_lock<std::mutex>lock(mu);
//                                auto vvvvv0 = Point_K2_to_iGameVertex(it1->vertex(0));
//                                auto vvvvv1 = Point_K2_to_iGameVertex(it1->vertex(1));
//                                auto vvvvv2 = Point_K2_to_iGameVertex(it1->vertex(2));
//                                auto vvvvv3 = Point_K2_to_iGameVertex(it2->vertex(0));
//                                auto vvvvv4 = Point_K2_to_iGameVertex(it2->vertex(1));
//                                auto vvvvv5 = Point_K2_to_iGameVertex(it2->vertex(2));
//
//                                fprintf(file3,"v %lf %lf %lf\n",vvvvv0.x(),vvvvv0.y(),vvvvv0.z());
//                                fprintf(file3,"v %lf %lf %lf\n",vvvvv1.x(),vvvvv1.y(),vvvvv1.z());
//                                fprintf(file3,"v %lf %lf %lf\n",vvvvv2.x(),vvvvv2.y(),vvvvv2.z());
//                                fprintf(file3,"v %lf %lf %lf\n",vvvvv3.x(),vvvvv3.y(),vvvvv3.z());
//                                fprintf(file3,"v %lf %lf %lf\n",vvvvv4.x(),vvvvv4.y(),vvvvv4.z());
//                                fprintf(file3,"v %lf %lf %lf\n",vvvvv5.x(),vvvvv5.y(),vvvvv5.z());
//                                fprintf(file3,"v %lf %lf %lf\n",tv1.x(),tv1.y(),tv1.z());
//                                fprintf(file3,"v %lf %lf %lf\n",tv2.x(),tv2.y(),tv2.z());
//                                fprintf(file3,"f %d %d %d\n",vvvvid,vvvvid+1,vvvvid+2);
//                                fprintf(file3,"f %d %d %d\n",vvvvid+3,vvvvid+4,vvvvid+5);
//                                vvvvid+=8;
//                            }
//                        }
//                    }
                }
                //for(int other_face_id = 0 ; other_face_id < possible_face_list.size())

            }

        },i);
    }
    for(int i=0;i<thread_num;i++)
        each_frame_thread[i]->join();
    cout << maxx_face_size << endl;
    cout << sum_face_size << endl;

    // 二类面判断交点在什么方向;
    return 0;
}
// 1 2 3 4 5 6
// 5 1 2 3 7 8
//
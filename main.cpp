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
#include "DSU.h"
//#include "BVH.h"


using namespace std;


shared_ptr <MeshKernel::SurfaceMesh> mesh;

shared_ptr<CGALPolygon>cgal_polygon;

Tree cgal_global_aabbtree;
double default_move = 0.1;
int thread_num = 16;

const double tolerance = 0.2;


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



vector <vector<int>> GridVertexDir = {{0, 0, 0},
                                      {1, 0, 0},
                                      {0, 1, 0},
                                      {0, 0, 1},
                                      {1, 1, 0},
                                      {1, 0, 1},
                                      {0, 1, 1},
                                      {1, 1, 1}};

vector <vector<int>> GridEdge = {{1, 2, 3},
                                 {0, 4, 5},
                                 {0, 4, 6},
                                 {0, 5, 6},
                                 {1, 2, 7},
                                 {1, 3, 7},
                                 {2, 3, 7},
                                 {4, 5, 6}};

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

    return MeshKernel::iGameVertex(small.x() + GridVertexDir[k][0] * (big.x()-small.x()),
                                   small.y() + GridVertexDir[k][1] * (big.y()-small.y()),
                                   small.z() + GridVertexDir[k][2] * (big.z()-small.z()));
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


vector<MeshKernel::iGameVertex>field_move_vertex;

//vector<vector<int> >approximate_field_face_table = {{0,1,2},{3,4,5},{0,2,4},{4,3,0},{1,5,4},{4,2,1},{0,3,5},{5,1,0}};

inline Point iGameVertex_to_Point(const MeshKernel::iGameVertex& v){
    return Point(v.x(),v.y(),v.z());
}

inline K2::Point_3 iGameVertex_to_Point_K2(const MeshKernel::iGameVertex& v){
    return K2::Point_3(v.x(),v.y(),v.z());
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

struct ApproximateField{
    vector<MeshKernel::iGameVertex> origin_vertices;
    vector<MeshKernel::iGameVertex> extend_vertices;
    vector<K2::Tetrahedron_3>tet_list;
    MeshKernel::iGameFaceHandle fh;
    ApproximateField(){}
    ApproximateField(MeshKernel::iGameFaceHandle fh){
        this->fh=fh;
        origin_vertices.push_back(mesh->fast_iGameVertex[mesh->fast_iGameFace[fh].vh(0)]);
        origin_vertices.push_back(mesh->fast_iGameVertex[mesh->fast_iGameFace[fh].vh(1)]);
        origin_vertices.push_back(mesh->fast_iGameVertex[mesh->fast_iGameFace[fh].vh(2)]);
        //TODO 注意分解
        MeshKernel::iGameVertex normal = ((origin_vertices[1] - origin_vertices[0]) % (origin_vertices[2] - origin_vertices[0])).normalize();
        for(int i=0;i<3;i++){
            //extend_vertices.push_back(field_move_vertex[mesh->fast_iGameFace[fh].vh(i)]);
            if(!in_triangle_positive_side(this->fh,field_move_vertex[mesh->fast_iGameFace[fh].vh(i)])){
                extend_vertices.push_back(field_move_vertex[mesh->fast_iGameFace[fh].vh(i)]);
            }
            else{
                extend_vertices.push_back(origin_vertices[i] + normal * mesh->fast_iGameFace[fh].move_dist);
            }
        }

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
struct ApproximateFieldTet{
    MeshKernel::iGameFaceHandle fh;
    int tet_id;
    friend bool operator < (const ApproximateFieldTet &a,const ApproximateFieldTet &b){
        if(a.fh != b.fh)
            return a.fh < b.fh;
        else
            return a.tet_id < b.tet_id;
    }
};


vector<ApproximateField>faces_approximate_field;
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

        MeshKernel::iGameVertex move_max_v = v + normal * avg_move_dist * 1.35;
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


struct Fraction{
    int numerator;
    int denominator;
    friend bool operator == (const Fraction &a,const Fraction &b){
        return a.numerator == b.numerator && a.denominator == b.denominator;
    }

};

bool vertex_in_tiny_grid(const MeshKernel::iGameVertex& small,const MeshKernel::iGameVertex& big,
                         const MeshKernel::iGameVertex& v){
    return small.x() <= v.x()  && v.x() <= big.x() &&
    small.y() <= v.y()  && v.y() <= big.y() &&
    small.z() <= v.z()  && v.z() <= big.z();
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

    for(int i=0;i<7;i++){
        MeshKernel::iGameVertex v1 = getTinyGridVertex(small,big,i);
        for(int j=0;j<DirectedGridEdge[i].size();j++){
            MeshKernel::iGameVertex v2 = getTinyGridVertex(small,big,DirectedGridEdge[i][j]);
            K::Segment_3 se(K::Point_3(v1.x(),v1.y(),v1.z()),
                             K::Point_3(v2.x(),v2.y(),v2.z()));

            CGAL::cpp11::result_of<K::Intersect_3(K::Segment_3, K::Triangle_3)>::type
                    result = intersection(se, tri);
            if(result)
                return true;
        }
    }
    return false;
}

/*
 *  N*N*N 来确定点
 *   然后保证后，通过格点表面采样法，然后连接边还原
 *    正负不管了
 */


struct SharpPoint{
    K2::Point_3 p;
    vector<int> source_face_local_id;
};

int main() {


    // freopen("../debugoutput.txt","w",stdout);
    default_move = 1;
    grid_len = 2.5;
    cout << grid_len <<endl;
    mix_factor = 0.5;
    //mesh = make_shared<MeshKernel::SurfaceMesh>(ReadObjFile("../data/debug5.obj2"));

//    MeshKernel::iGameVertex v1(100,-1,-1);
//    MeshKernel::iGameVertex v2(1,-1,-1);
//    double t = ternary_search(v1,v2,MeshKernel::iGameFaceHandle(0)) ;
//    cout << t << endl;
//
//    return 0;


    // mesh = make_shared<MeshKernel::SurfaceMesh>(ReadObjFile("../data/Armadillo.obj"));


    mesh = make_shared<MeshKernel::SurfaceMesh>(ReadObjFile("../data/test_orgv2.obj2"));

    for(int i=0;i<mesh->FaceSize();i++){
        mesh->faces(MeshKernel::iGameFaceHandle(i)).move_dist = 0.8;
    }

//    for(int i=0;i<mesh->FaceSize();i++){
//        mesh->faces(MeshKernel::iGameFaceHandle(i)).move_dist = 0.8;
//    }

    mesh->build_fast();


  //  mix(30);

//    BVH::BVH_Tree bvhTree;
//    bvhTree.buildBVH_Tree(*mesh);
//    auto res = bvhTree.getIntersection(BVH::Ray(Vector3d(1,1,1),Vector3d(-1,-1,-1)));
//    cout << res.pos[0]<<" "<< res.pos[1]<<" "<< res.pos[2] << endl;


    //只动xy
    faces_approximate_field.resize(mesh->FaceSize());
    field_move_vertex.resize(mesh->VertexSize());
    min_move_g.resize(mesh->VertexSize());
    max_move_g.resize(mesh->VertexSize());



    for(int i=0;i<mesh->VertexSize();i++){
        field_move_vertex[i] = do_quadratic_error_metric(MeshKernel::iGameVertexHandle(i));
    }

    for(int i=0;i<mesh->FaceSize();i++){
        faces_approximate_field[i] = ApproximateField(MeshKernel::iGameFaceHandle(i));
    }
//

    int file_id = 3013;
    FILE *file = fopen(("../data/output" + to_string(file_id) + ".obj").c_str(), "w");
    FILE *file0 = fopen(("../data/output" + to_string(file_id) + "_dianyun0.obj").c_str(), "w");
    FILE *file1 = fopen(("../data/output" + to_string(file_id) + "_dianyun1.obj").c_str(), "w");
    FILE *file2 = fopen(("../data/output" + to_string(file_id) + "_dianyun2.obj").c_str(), "w");
    FILE *file10 = fopen(("../data/output" + to_string(file_id) + "_dianyun10.obj").c_str(), "w");

    int cnt=1;
    for(int i=0;i<mesh->FaceSize();i++){
//        auto v0 = faces_approximate_field[i].extend_vertices[0];
//        auto v1 = faces_approximate_field[i].extend_vertices[1];
//        auto v2 = faces_approximate_field[i].extend_vertices[2];
        //field_move_vertex[i]
        auto v0 = field_move_vertex[mesh->fast_iGameFace[i].vh(0)];
        auto v1 = field_move_vertex[mesh->fast_iGameFace[i].vh(1)];
        auto v2 = field_move_vertex[mesh->fast_iGameFace[i].vh(2)];
        fprintf(file0, "v %lf %lf %lf\n", v0.x(), v0.y(),v0.z());
        fprintf(file0, "v %lf %lf %lf\n", v1.x(), v1.y(),v1.z());
        fprintf(file0, "v %lf %lf %lf\n", v2.x(), v2.y(),v2.z());
        //fprintf(file0, "f %d %d %d\n", cnt, cnt+1,cnt+2);
    }
    for(int i=0;i<mesh->FaceSize();i++){
        fprintf(file0, "f %d %d %d\n", cnt, cnt+1,cnt+2);
        cnt+=3;

    }
    cnt=0;
    for(int  i =0 ; i< mesh->VertexSize() ; i++){
    //    if(flag) {
//            cout << "origin " <<mesh->fast_iGameVertex[i].x()<<" "<< mesh->fast_iGameVertex[i].y()<<" "<<
//                    mesh->fast_iGameVertex[i].z()<<endl;
//            cout << "updated  " <<field_move_vertex[i].x()<<" "<< field_move_vertex[i].y()<<" "<<
//                    field_move_vertex[i].z()<<endl;
//            cout <<"***************"<< endl;

//            cout <<"***************"<< endl;
            fprintf(file1, "v %lf %lf %lf\n", mesh->fast_iGameVertex[i].x(), mesh->fast_iGameVertex[i].y(),
                   mesh->fast_iGameVertex[i].z());
            fprintf(file1, "v %lf %lf %lf\n", field_move_vertex[i].x(), field_move_vertex[i].y(),
                    field_move_vertex[i].z());

//            fprintf(file0, "v %lf %lf %lf\n", min_move_g[i].x(), min_move_g[i].y(),
//                    min_move_g[i].z());

            fprintf(file1,"l %d %d\n",cnt+1,cnt+2);
          //  fprintf(file0,"l %d %d\n",cnt+1,cnt+3);
            cnt+=2;

       // }
    }

    for(int  i =0 ; i< mesh->VertexSize() ; i++){

       if((mesh->fast_iGameVertex[i] - field_move_vertex[i]).norm()>100){
           fprintf(file1, "v %lf %lf %lf\n", mesh->fast_iGameVertex[i].x(), mesh->fast_iGameVertex[i].y(),mesh->fast_iGameVertex[i].z());
       }

    }


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

    vector <vector<int>> dir = {{0,  0,  1},
                                {0,  1,  0},
                                {1,  0,  0},
                                {0,  0,  -1},
                                {0,  -1, 0},
                                {-1, 0,  0}};
    std::function < vector<grid>(grid) > get_neighbor = [&](grid g) {
        vector <grid> ret;
        for (auto i: dir) {
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

    std::mutex bfs_mutex;
    std::vector <std::shared_ptr<std::thread> > bfs_thread_pool(thread_num);
    for(int i=0;i<thread_num;i++){
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
                set <grid> is_visit;
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
                            if (dist <  max(grid_len ,fh.second.move_dist)*1.1 ) { //TODO : zheli youhua cheng pianyi juli de shiji jisuan
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

    std::function<bool(Plane_3 plane,MeshKernel::iGameVertex,double)> vertex_in_plane
    = [&](Plane_3 pla,MeshKernel::iGameVertex v,double eps){
                if(sqrt(squared_distance(pla,Point_3(v.x(),v.y(),v.z())))<eps){
                    return true;
                }
                return false;
    };

    // 上述代码完成距离场建格子的过程 8 ;
    int tt=0;
    int sum_grid = 0;
    for (auto each_grid = frame_grid_mp.begin(); each_grid != frame_grid_mp.end(); each_grid++) {
        if(tt++%5000==0){
            cout << tt <<" // "<<  frame_grid_mp.size() << endl;
        }
        set<MeshKernel::iGameFaceHandle> face_set;
        vector<MeshKernel::iGameFaceHandle> face_list;
        vector<MeshKernel::iGameFaceHandle> possible_face_list;
        vector<MeshKernel::iGameFace>vf;


        for (int j = 0; j < 8; j++) {
            grid g = getGridFrameVertex(each_grid->first, j);
            unordered_map<grid, GridVertex, grid_hash, grid_equal>::iterator it = frame_grid_mp.find(g);
            if (it == frame_grid_mp.end())continue;
            for (auto k: it->second.face_list) {
                if(!face_set.count(k)){
                    face_set.insert(k);
                    face_list.push_back(k);
                }
            }
        }
        MeshKernel::iGameVertex small = getGridVertex(each_grid->first,0);
        MeshKernel::iGameVertex big = getGridVertex(each_grid->first,7);

        for(MeshKernel::iGameFaceHandle i : face_list){
            vector<MeshKernel::iGameVertex> v;
            for(int j=0;j<3;j++)
                v.push_back(field_move_vertex[mesh->fast_iGameFace[i].vh(j)]);
            if(face_through_grid(small,big,v))
                possible_face_list.push_back(i);
        }

        for(auto i : face_list){
            vf.push_back(mesh->fast_iGameFace.at(i));
        }
       // cout << possible_face_list.size() <<"/" << face_list.size() << endl;

        //face_through_grid


        DSU dsu(possible_face_list.size());
        vector<vector<int> >plane_cross;
        plane_cross.resize(possible_face_list.size());
        for(int i=0;i<possible_face_list.size();i++){
            for(int j=i+1;j<possible_face_list.size();j++){
                int root_i = dsu.find_root(i);
                MeshKernel::iGameVertex v0i = field_move_vertex[mesh->fast_iGameFace[possible_face_list[i]].vh(0)];
                MeshKernel::iGameVertex v1i = field_move_vertex[mesh->fast_iGameFace[possible_face_list[i]].vh(1)];
                MeshKernel::iGameVertex v2i = field_move_vertex[mesh->fast_iGameFace[possible_face_list[i]].vh(2)];
                Plane_3 pla_i (Point_3(v0i.x(),v0i.y(),v0i.z()),
                             Point_3(v1i.x(),v1i.y(),v1i.z()),
                             Point_3(v2i.x(),v2i.y(),v2i.z()));
                MeshKernel::iGameVertex v0j = field_move_vertex[mesh->fast_iGameFace[possible_face_list[j]].vh(0)];
                MeshKernel::iGameVertex v1j = field_move_vertex[mesh->fast_iGameFace[possible_face_list[j]].vh(1)];
                MeshKernel::iGameVertex v2j = field_move_vertex[mesh->fast_iGameFace[possible_face_list[j]].vh(2)];

                K::Triangle_3 tri_face1(Point_3(v0i.x(),v0i.y(),v0i.z()),
                                        Point_3(v1i.x(),v1i.y(),v1i.z()),
                                        Point_3(v2i.x(),v2i.y(),v2i.z()));
                K::Triangle_3 tri_face2(Point_3(v0j.x(),v0j.y(),v0j.z()),
                                        Point_3(v1j.x(),v1j.y(),v1j.z()),
                                        Point_3(v2j.x(),v2j.y(),v2j.z()));



                double edge_eps = sqrt(max({(v0i-v1i).norm(),(v2i-v1i).norm(),(v0i-v2i).norm()}))/15;
                if(vertex_in_plane(pla_i,v0j,edge_eps) &&
                        vertex_in_plane(pla_i,v1j,edge_eps) &&
                        vertex_in_plane(pla_i,v2j,edge_eps) &&
                        sqrt(squared_distance(tri_face1,tri_face2)) < myeps) { //并且直接相连 ;;;;
                    dsu.join(i,j);
                }
            }
        }

//        for(int i=0;i<face_list.size();i++){
//            cout <<"dsu.find_root(i):" << dsu.find_root(i) << endl;
//        }
        vector<SharpPoint> sharp_point_list;
        for(int i=0;i<possible_face_list.size();i++) {
            MeshKernel::iGameVertex v0i = field_move_vertex[mesh->fast_iGameFace[possible_face_list[i]].vh(0)];
            MeshKernel::iGameVertex v1i = field_move_vertex[mesh->fast_iGameFace[possible_face_list[i]].vh(1)];
            MeshKernel::iGameVertex v2i = field_move_vertex[mesh->fast_iGameFace[possible_face_list[i]].vh(2)];
            K2::Triangle_3 tri_i (K2::Point_3(v0i.x(),v0i.y(),v0i.z()),
                                  K2::Point_3(v1i.x(),v1i.y(),v1i.z()),
                                  K2::Point_3(v2i.x(),v2i.y(),v2i.z()));

            for (int j = i + 1; j < possible_face_list.size(); j++) {
                if(dsu.find_root(i) != dsu.find_root(j)){

                    MeshKernel::iGameVertex v0j = field_move_vertex[mesh->fast_iGameFace[possible_face_list[j]].vh(0)];
                    MeshKernel::iGameVertex v1j = field_move_vertex[mesh->fast_iGameFace[possible_face_list[j]].vh(1)];
                    MeshKernel::iGameVertex v2j = field_move_vertex[mesh->fast_iGameFace[possible_face_list[j]].vh(2)];
                    K2::Triangle_3 tri_j (K2::Point_3(v0j.x(),v0j.y(),v0j.z()),
                                          K2::Point_3(v1j.x(),v1j.y(),v1j.z()),
                                          K2::Point_3(v2j.x(),v2j.y(),v2j.z()));

                    CGAL::cpp11::result_of<K2::Intersect_3(K2::Triangle_3, K2::Triangle_3)>::type
                            result = intersection(tri_i, tri_j);
                    if (result) {

                        if (const K2::Point_3 * p = boost::get<K2::Point_3>(&*result)) {
                            // std::cout << (*p) << std::endl;
                            // cout <<"poi: "<< *p << endl;
                            for (int k = j + 1; k < possible_face_list.size(); k++){
                                //cout << "PPPPPPPP " << endl;
                                if(dsu.find_root(i) == dsu.find_root(k) ||
                                    dsu.find_root(j) == dsu.find_root(k))continue;

                                MeshKernel::iGameVertex v0k = field_move_vertex[mesh->fast_iGameFace[possible_face_list[k]].vh(0)];
                                MeshKernel::iGameVertex v1k = field_move_vertex[mesh->fast_iGameFace[possible_face_list[k]].vh(1)];
                                MeshKernel::iGameVertex v2k = field_move_vertex[mesh->fast_iGameFace[possible_face_list[k]].vh(2)];
                                K2::Triangle_3 tri_k (K2::Point_3(v0k.x(),v0k.y(),v0k.z()),
                                                    K2::Point_3(v1k.x(),v1k.y(),v1k.z()),
                                                    K2::Point_3(v2k.x(),v2k.y(),v2k.z()));
                                CGAL::cpp11::result_of<K2::Intersect_3(K2::Point_3, K2::Triangle_3)>::type
                                        result2 = intersection(*p, tri_k);
                                if (result2) {
                                    //  cout << "GGGGGGG2 " << endl;
                                    if (const K2::Point_3 * p2 = boost::get<K2::Point_3>(&*result2)) {
                                        MeshKernel::iGameVertex igame_p2(CGAL::to_double(p2->x()),
                                         CGAL::to_double(p2->y()),
                                         CGAL::to_double(p2->z()));
                                        if(!vertex_in_grid(each_grid->first,igame_p2))
                                            continue;

                                        bool in_tet_field = false;
                                        for(int field_face_id=0;field_face_id<face_list.size();field_face_id++){
                                            for(int each_tet_id =0; each_tet_id<6;each_tet_id++ ){
                                                if(faces_approximate_field[face_list[field_face_id]].tet_list[each_tet_id].has_on_bounded_side(*p2)){
                                                    in_tet_field = true;
                                                    break;
                                                }
                                            }
                                            if(in_tet_field) break;
                                        }
                                        int side_type = 0;
                                        //cout <<"******\n" << vf.size() << " "<<face_list.size() << endl;
                                        cgal_aabbtree_query(vf, igame_p2, mesh, side_type, cgal_polygon);

                                        if ( (!in_tet_field) && side_type ==1 ) {
                                            SharpPoint sp;
                                            sp.p = *p2;
                                            sp.source_face_local_id.push_back(i);
                                            sp.source_face_local_id.push_back(j);
                                            sp.source_face_local_id.push_back(k);
                                            sharp_point_list.push_back(sp);
                                        }
                                    }
                                }
                            }
                        }
                        if (const K2::Segment_3 * s = boost::get<K2::Segment_3 >(&*result)) {
                            plane_cross[dsu.find_root(i)].push_back(dsu.find_root(j));
                            plane_cross[dsu.find_root(j)].push_back(dsu.find_root(i));


                            // std::cout << (*p) << std::endl;
//                            cout <<"tri:" <<tri_i.has_on(s->point(0))<<" "<< tri_i.has_on(s->point(1))<<" "
//                            <<  tri_j.has_on(s->point(0))<<" "<< tri_j.has_on(s->point(1)) << endl;
//                             cout << "tri i "<< tri_i << endl;
//                            cout << "tri j "<< tri_j << endl;
                            for (int k = j + 1; k < possible_face_list.size(); k++){
                                if(dsu.find_root(i) == dsu.find_root(k) ||
                                   dsu.find_root(j) == dsu.find_root(k))continue;

                                //cout << "GGGGGGG " << endl;

                                MeshKernel::iGameVertex v0k = field_move_vertex[mesh->fast_iGameFace[possible_face_list[k]].vh(0)];
                                MeshKernel::iGameVertex v1k = field_move_vertex[mesh->fast_iGameFace[possible_face_list[k]].vh(1)];
                                MeshKernel::iGameVertex v2k = field_move_vertex[mesh->fast_iGameFace[possible_face_list[k]].vh(2)];
                                K2::Triangle_3 tri_k (K2::Point_3(v0k.x(),v0k.y(),v0k.z()),
                                                      K2::Point_3(v1k.x(),v1k.y(),v1k.z()),
                                                      K2::Point_3(v2k.x(),v2k.y(),v2k.z()));

                               // cout << "seg s " << *s << endl;
                               // cout << "tri_k  " << tri_k << endl;
                                CGAL::cpp11::result_of<K2::Intersect_3(K2::Segment_3, K2::Triangle_3)>::type
                                        result2 = intersection(*s, tri_k);
                                if (result2) {
                                  //  cout << "GGGGGGG2 " << endl;
                                    if (const K2::Point_3 * p2 = boost::get<K2::Point_3>(&*result2)) {

                                        MeshKernel::iGameVertex igame_p2(CGAL::to_double(p2->x()),
                                                                         CGAL::to_double(p2->y()),
                                                                         CGAL::to_double(p2->z()));

                                        if(!vertex_in_grid(each_grid->first,igame_p2))
                                            continue;

                                        bool in_tet_field = false;
                                        for(int field_face_id=0;field_face_id<face_list.size();field_face_id++){
                                            for(int each_tet_id =0; each_tet_id<6;each_tet_id++ ){
                                                if(faces_approximate_field[face_list[field_face_id]].tet_list[each_tet_id].has_on_bounded_side(*p2)){
                                                    in_tet_field = true;
                                                    break;
                                                }
                                            }
                                            if(in_tet_field) break;
                                        }
                                        int side_type = 0;
                                        //cout <<"******\n" << vf.size() << " "<<face_list.size() << endl;
                                        cgal_aabbtree_query(vf, igame_p2, mesh, side_type, cgal_polygon);

                                        if ( (!in_tet_field) && side_type ==1 ) {
                                            SharpPoint sp;
                                            sp.p = *p2;
                                            sp.source_face_local_id.push_back(i);
                                            sp.source_face_local_id.push_back(j);
                                            sp.source_face_local_id.push_back(k);
                                            sharp_point_list.push_back(sp);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        //cout << "sharp_point_list.size(): "<< sharp_point_list.size() << endl;

        vector<SharpPoint> sharp_point_list_merge;
        for(int i=0;i<sharp_point_list.size();i++){
            bool flag = false;
            for(int j=0;j<sharp_point_list_merge.size() && !flag; j++){
                if(sqrt(CGAL::to_double(squared_distance(sharp_point_list_merge[j].p,sharp_point_list[i].p)))<myeps){
                    flag = true;
                    for(auto k : sharp_point_list[i].source_face_local_id)
                        sharp_point_list_merge[j].source_face_local_id.push_back(k);
                }
            }
            if(!flag){
                sharp_point_list_merge.push_back(sharp_point_list[i]);
            }
        }
        for(int i=0;i<sharp_point_list_merge.size();i++){
            sort(sharp_point_list_merge[i].source_face_local_id.begin(),sharp_point_list_merge[i].source_face_local_id.end());
            sharp_point_list_merge[i].source_face_local_id.resize(
                    unique(sharp_point_list_merge[i].source_face_local_id.begin(),
                           sharp_point_list_merge[i].source_face_local_id.end())-sharp_point_list_merge[i].source_face_local_id.begin());
        }

        swap(sharp_point_list_merge,sharp_point_list);

        for(auto it : sharp_point_list){
            fprintf(file2,"v %lf %lf %lf\n",CGAL::to_double(it.p.x()),CGAL::to_double(it.p.y()),CGAL::to_double(it.p.z()));
        }


        // 注意合并处理 ;



        double maxlen = grid_len;
        for(int i=0;i<sharp_point_list.size();i++){
            for(int j=i+1;j<sharp_point_list.size();j++){
                double dist = sqrt(CGAL::to_double(squared_distance(sharp_point_list[i].p,sharp_point_list[j].p))/3.0);
                maxlen = min(maxlen , dist);
            }
            for (int j = 0; j < possible_face_list.size(); j++) {
                bool flag = true;
                for(int k=0;k<sharp_point_list[i].source_face_local_id.size();k++){
                    if( dsu.find_root(sharp_point_list[i].source_face_local_id[k]) == dsu.find_root(j))
                        flag = false;
                }
                if (flag) {
                    MeshKernel::iGameVertex v0j = field_move_vertex[mesh->fast_iGameFace[possible_face_list[j]].vh(0)];
                    MeshKernel::iGameVertex v1j = field_move_vertex[mesh->fast_iGameFace[possible_face_list[j]].vh(1)];
                    MeshKernel::iGameVertex v2j = field_move_vertex[mesh->fast_iGameFace[possible_face_list[j]].vh(2)];
                    K2::Triangle_3 tri(K2::Point_3(v0j.x(),v0j.y(),v0j.z()),
                                   K2::Point_3(v1j.x(),v1j.y(),v1j.z()),
                                   K2::Point_3(v2j.x(),v2j.y(),v2j.z()));
                    double dist = sqrt(CGAL::to_double(squared_distance(sharp_point_list[i].p,tri))/3.0);
                    maxlen = min(maxlen , dist);
                }
            }
        }
        int num = 1 ;
        if(sharp_point_list.size() > 1 ){
            //cout << sharp_point_list.size() <<" "<< maxlen <<  endl;
            int div = int((grid_len / maxlen)+1);
            int d = 0;
            int tmp = div -1 ;
            while(tmp){
                tmp>>=1;d++;
            }
            num  = 1 << d;
           // cout << num << endl;
            sum_grid+= num;
        }
        else
            sum_grid+= 1;
        double tiny_grid_len = grid_len/num;



        /*
         * 三类格子
         *  1356 行的地方先确定相交确定为两个
         *  两条然后相交线延长出去
         *  做格子外壳与三角面片求交
         *  然后连线
         *  正负方向通过在tet内外搞
         *  能写！
         *  二类格子
         *
         */







//        for(int i=0;i<sharp_point_list.size();i++){
//            for(int k1 = 0;k1<sharp_point_list[i].source_face_local_id.size();k1++)
//            {
//                for(int k2 = 0;k2<sharp_point_list[i].source_face_local_id.size();k2++){
//                    if(k1!=k2){
//
//                    }
//                }
//            }
//        }

//        for(int x=0;x<num;x++)
//            for(int y=0;y<num;y++)
//                for(int z=0;z<num;z++){
//                    MeshKernel::iGameVertex start_pos = small + MeshKernel::iGameVertex(x*tiny_grid_len,
//                                                                 y*tiny_grid_len,
//                                                                 z*tiny_grid_len);
//                    MeshKernel::iGameVertex end_pos = small + MeshKernel::iGameVertex(x * (tiny_grid_len+1),
//                                                                                     y * (tiny_grid_len+1),
//                                                                                     z * (tiny_grid_len+1));
//
//                    int tiny_grid_type = 0;
//                    for(int i = 0 ;i < sharp_point_list.size();i++){
//                        MeshKernel::iGameVertex sharp_vertex(CGAL::to_double(sharp_point_list[i].p.x()),
//                                                        CGAL::to_double(sharp_point_list[i].p.y()),
//                                                        CGAL::to_double(sharp_point_list[i].p.z()));
//                        if(vertex_in_tiny_grid(start_pos,end_pos,sharp_vertex)){
//                            tiny_grid_type = 3;
//                            for(int j=0;j<sharp_point_list[i].source_face_local_id.size();j++) {
//
//                            }
//                            break;
//                        }
//                    }
//                    if(tiny_grid_type == 3)
//                        continue;
//
//                    // check is type 3 !!!!!!
//
//
//
//
//                    /*
//                     *         for(MeshKernel::iGameFaceHandle i : face_list){
//            vector<MeshKernel::iGameVertex> v;
//            for(int j=0;j<3;j++)
//                v.push_back(field_move_vertex[mesh->fast_iGameFace[i].vh(j)]);
//            if(face_through_grid(small,big,v))
//                possible_face_list.push_back(i);
//        }
//                     */
//
//                }


//        if(num >60){
//           // cout <<"sharp_point_list.size()" <<sharp_point_list.size() << "\n";
//            for(auto it : sharp_point_list){
//                fprintf(file10,"v %lf %lf %lf\n",CGAL::to_double(it.p.x()),CGAL::to_double(it.p.y()),CGAL::to_double(it.p.z()));
//            }
//        }

        //vertex_in_grid


//        set<int>sse;
//        for(int i=0;i<face_list.size();i++){
//            sse.insert(dsu.find_root(i));
//        }

        //cout << "sse.size(): " << sse.size() << endl;
    }
    cout << sum_grid << endl;
    // 二类面判断交点在什么方向;














    return 0;


    vector<vector<int> >tet_face_order{{0,1,2},{0,1,3},{0,2,3},{1,2,3}};
    // 用cgal



//    for(auto each_dir : dir){
//        ThickenCubeEdge edge_info;
//        edge_info.from = each_grid->first;
//        edge_info.to = grid(edge_info.from.x+each_dir[0],edge_info.from.y+each_dir[1],edge_info.from.z+each_dir[2]);
//        if(!frame_grid_mp.count(edge_info.to))continue;
//        edge_info.type = each_grid->second.grid_type;
//        vector<MeshKernel::iGameFaceHandle> face_list;
//        for(auto i : each_grid->second.face_list) {
//            if ()// zaili mian 不对 这个地方 123
//            {
//            }
//        }
//        if(edge_info.type == 0){
//
//        }
//        else if(edge_info.type ==1){
//
//        }
//
//    }



    //thicken_cube_edge_list.resize()

    //DSU
    //TODO 可改多线程 555




























    return 0;



    function<MeshKernel::iGameVertex(grid)>qef_grid_pos_all_consider = [&](grid g) {
        vector<double>b;
        vector<double>sx,sy,sz,nx,ny,nz;
        MeshKernel::iGameVertex center_of_weight(0,0,0);
        //cout<< "************************************************** "<< endl;
        //cout<< "grid g : "<<g.x<<" "<<g.y<<" "<< g.z << endl;
        // todo: 改这里，改成所有面能力的叠加 20220623 搞这里是时刻
        int center_cnt = 0;
        for (int ev1 = 0; ev1 < GridEdge.size(); ev1++) {
            for (auto ev2: GridEdge[ev1]) {
                if (ev1 < ev2) {
                    grid v1 = getGridFrameVertex(g,ev1);
                    grid v2 = getGridFrameVertex(g,ev2);
                    auto itv1 = frame_grid_mp.find(v1);
                    auto itv2 = frame_grid_mp.find(v2);
                    if(itv1 == frame_grid_mp.end() || itv2 == frame_grid_mp.end()){
                        continue;
                    }
                    if(itv1->second.grid_type > itv2->second.grid_type) {
                        swap(itv1, itv2);
                        swap(v1,v2);
                    }
                    // 这里确实有问题！ 4454545 先不管了 小丑 55555

                    if( (itv1->second.grid_type ==1||itv1->second.grid_type ==0) && itv2->second.grid_type ==2) { //判断1 2 的face list
                        //  cout << "itv1 : "<<itv1->first.x<<" "<<  itv1->first.y<<" "<<  itv1->first.z<<" "<<itv1->second.grid_type<< endl;
                        // cout << "itv2 : "<<itv2->first.x<<" "<<  itv2->first.y<<" "<<  itv2->first.z<<" "<< itv2->second.grid_type<<endl;

                        //cout <<"itv1fl size:"<< itv1->second.face_list.size() <<endl;


                        MeshKernel::iGameVertex v0_vertex = getGridVertex(v1, 0);
                        MeshKernel::iGameVertex v1_vertex = getGridVertex(v2, 0);
                        MeshKernel::iGameVertex pi = v0_vertex;
                        if(itv1->second.grid_type ==0)
                            pi = (v0_vertex + v1_vertex)/2;
                        MeshKernel::iGameVertex ni;
                        int tmp;
                        for (MeshKernel::iGameFaceHandle fht: itv1->second.face_list) {

//                            ni += (((mesh->fast_iGameVertex[mesh->fast_iGameFace.at(fht).vh(1)]
//                            - mesh->fast_iGameVertex[mesh->fast_iGameFace.at(fht).vh(0)]) %
//                                    (mesh->fast_iGameVertex[mesh->fast_iGameFace.at(fht).vh(2)]
//                                     - mesh->fast_iGameVertex[mesh->fast_iGameFace.at(fht).vh(0)]))*-1).normalize();

                            ni += (pi - cgal_aabbtree_query({mesh->fast_iGameFace.at(fht)}, pi, mesh, tmp,
                                                            nullptr)).normalize();
                        }
                        ni /= itv1->second.face_list.size();

                        double l = 0, r = 1;
                        double limit = 0.01;
                        //    cout << "st erfen"<< endl;
                        while (l + limit < r) {
                            double mid = (l + r) / 2;
                            MeshKernel::iGameVertex cut_point = (v0_vertex + (v1_vertex - v0_vertex) * mid);
                            bool flag = false;
                            std::shared_ptr<vector<MeshKernel::iGameFaceHandle> > field_face_list =
                                    std::make_shared<vector<MeshKernel::iGameFaceHandle> >();
                            int state = vertex_state(itv1->second.face_list,cut_point,field_face_list);
                            if (state == 2)
                                r = mid;
                            else
                                l = mid;

                            if (state == 1) {
                                pi = cut_point;
                                ni = (pi - cgal_aabbtree_query({mesh->fast_iGameFace.at(*field_face_list->begin())}, pi, mesh,tmp,
                                                               nullptr)).normalize();

//                                ni = (((mesh->fast_iGameVertex[mesh->fast_iGameFace.at(*field_face_list->begin()).vh(1)]
//                                         - mesh->fast_iGameVertex[mesh->fast_iGameFace.at(*field_face_list->begin()).vh(0)]) %
//                                        (mesh->fast_iGameVertex[mesh->fast_iGameFace.at(*field_face_list->begin()).vh(2)]
//                                         - mesh->fast_iGameVertex[mesh->fast_iGameFace.at(*field_face_list->begin()).vh(0)]))*-1).normalize();
                                flag = true;
                            }

//                            cout <<CGALPlaneDist(mesh->faces(fh),v0_vertex)<<" "<< CGALPlaneDist(mesh->faces(fh),v1_vertex)<<
//                            " "<<  CGALPlaneDist(mesh->faces(fh),pi)<<endl;
                            //   cout <<l<< "  en erfen"<< endl;
                            //   cout << "pi.:" << pi.x() << " " << pi.y() << " " << pi.z() << endl;
                            sx.push_back(pi.x());
                            sy.push_back(pi.y());
                            sz.push_back(pi.z());
                            //   cout << "ni.:" << ni.x() << " " << ni.y() << " " << ni.z() << endl;

                            nx.push_back(ni.x());
                            ny.push_back(ni.y());
                            nz.push_back(ni.z());
                            center_of_weight = center_of_weight + pi;

                            // cout <<"pi "<<pi.x() <<" "<< pi.y()<<" "<<pi.z() << endl;
//                            fprintf(file1, "v %lf %lf %lf\n",pi.x(),pi.y(),pi.z());

                            center_cnt++;
                            //    cout <<"sxsize:" <<sx.size() << endl;
                        }
                    }

                }
            }
        }
        // cout <<"center_cnt:" <<center_cnt << endl;

        vector<double>res(3);
        QEF_calculate({getGridVertex(g,0).x(),getGridVertex(g,0).y(),getGridVertex(g,0).z()},
                      {getGridVertex(g,7).x(),getGridVertex(g,7).y(),getGridVertex(g,7).z()},
                      sx.size(),sx,sy,sz,nx,ny,nz,res);
        //  cout << "res: "<<res[0]<<" "<<res[1]<<" "<<res[2]<<endl;
        // cout << "qefinner? "<< vertex_in_grid(g,MeshKernel::iGameVertex(res[0],res[1],res[2])) << endl;
        // cout << center_cnt << endl;
        //cout <<"jiuji bug :"<< vertex_in_grid(g,center_of_weight/center_cnt)<<endl;
        //  cout << (center_of_weight/center_cnt).x() <<" "<<  (center_of_weight/center_cnt).y()<<" "<<
        //                                                                                      (center_of_weight/center_cnt).z()<<endl;
        //  cout << getGirdVertex(g,0).x()<<" "<<getGirdVertex(g,0).y()<<" "<<getGirdVertex(g,0).z()<<endl;
        //   cout << getGirdVertex(g,7).x()<<" "<<getGirdVertex(g,7).y()<<" "<<getGirdVertex(g,7).z()<<endl;
        // return center_of_weight/center_cnt;


        //TODO : debug code need del

        if(vertex_in_grid(g,MeshKernel::iGameVertex(res[0],res[1],res[2]))){
            return MeshKernel::iGameVertex(res[0],res[1],res[2]);
        }
        else{
            // cout<< "out"<<endl;
            return center_of_weight/center_cnt;
        }
    };




    function<MeshKernel::iGameVertex(grid)>qef_grid_pos = [&](grid g) {

        vector<double>b;
        vector<double>sx,sy,sz,nx,ny,nz;
        MeshKernel::iGameVertex center_of_weight(0,0,0);
        //cout<< "************************************************** "<< endl;
        //cout<< "grid g : "<<g.x<<" "<<g.y<<" "<< g.z << endl;
        // todo: 改这里，改成所有面能力的叠加 20220623 搞这里是时刻
        int center_cnt = 0;
        for (int ev1 = 0; ev1 < GridEdge.size(); ev1++) {
            for (auto ev2: GridEdge[ev1]) {
                if (ev1 < ev2) {
                    grid v1 = getGridFrameVertex(g,ev1);
                    grid v2 = getGridFrameVertex(g,ev2);
                    auto itv1 = frame_grid_mp.find(v1);
                    auto itv2 = frame_grid_mp.find(v2);
                    if(itv1 == frame_grid_mp.end() || itv2 == frame_grid_mp.end()){
                        continue;
                    }
                    if(itv1->second.grid_type > itv2->second.grid_type) {
                        swap(itv1, itv2);
                        swap(v1,v2);
                    }


                    if( (itv1->second.grid_type ==1||itv1->second.grid_type ==0) && itv2->second.grid_type ==2) { //判断1 2 的face list
                        MeshKernel::iGameVertex v0_vertex = getGridVertex(v1, 0);
                        MeshKernel::iGameVertex v1_vertex = getGridVertex(v2, 0);
                        MeshKernel::iGameVertex edge_pi ;
                        MeshKernel::iGameVertex edge_ni ;
                        int edge_cnt = 0;
                        for (MeshKernel::iGameFaceHandle fht: itv1->second.face_list) {
                            int state_v0 = vertex_state({fht}, v0_vertex);
                            int state_v1 = vertex_state({fht}, v1_vertex);
                            if (state_v0 != state_v1 && state_v1 == 2) {
                                double l = 0, r = 1;
                                double limit = 0.01;
                                bool flag = false;
                                MeshKernel::iGameVertex pi,ni;
                                while (l + limit < r) {
                                    double mid = (l + r) / 2;
                                    MeshKernel::iGameVertex cut_point = (v0_vertex + (v1_vertex - v0_vertex) * mid);
                                    int state_mid = vertex_state({fht}, cut_point);
                                    if (state_mid == 2)
                                        r = mid;
                                    else
                                        l = mid;
                                    if (state_mid == 1){
                                        pi = (v0_vertex + (v1_vertex - v0_vertex) * l);
                                        ni = (((mesh->fast_iGameVertex[mesh->fast_iGameFace.at(fht).vh(1)]
                                                  - mesh->fast_iGameVertex[mesh->fast_iGameFace.at(fht).vh(0)]) %
                                                 (mesh->fast_iGameVertex[mesh->fast_iGameFace.at(fht).vh(2)]
                                                  - mesh->fast_iGameVertex[mesh->fast_iGameFace.at(fht).vh(0)])) *
                                                -1).normalize();
                                        flag = true;
                                    }
                                }
                                if(flag) {
                                    edge_pi+=pi;
                                    edge_ni+=ni;
                                    edge_cnt++;
                                }
                                else{
                                    edge_ni+= (((mesh->fast_iGameVertex[mesh->fast_iGameFace.at(fht).vh(1)]
                                            - mesh->fast_iGameVertex[mesh->fast_iGameFace.at(fht).vh(0)]) %
                                           (mesh->fast_iGameVertex[mesh->fast_iGameFace.at(fht).vh(2)]
                                            - mesh->fast_iGameVertex[mesh->fast_iGameFace.at(fht).vh(0)])) *
                                          -1).normalize();
                                    edge_pi+=v1_vertex;
                                    edge_cnt++;
                                }
                            }
                        }
                        if(edge_cnt>0){
                            edge_pi/=edge_cnt;
                            edge_ni/=edge_cnt;
                            sx.push_back(edge_pi.x());
                            sy.push_back(edge_pi.y());
                            sz.push_back(edge_pi.z());
                            nx.push_back(edge_ni.x());
                            ny.push_back(edge_ni.y());
                            nz.push_back(edge_ni.z());
                            center_of_weight+=edge_pi;
                            center_cnt++;
                        }

                    }
                }//end  if (ev1 < ev2)
            }
        }
//        if(sx.size()==0 || center_cnt ==0){
//            return (MeshKernel::iGameVertex{getGridVertex(g,0).x(),getGridVertex(g,0).y(),getGridVertex(g,0).z()}+
//                    MeshKernel::iGameVertex{getGridVertex(g,7).x(),getGridVertex(g,7).y(),getGridVertex(g,7).z()})/2;
//
//        }
        // cout <<"center_cnt:" <<center_cnt << endl;

        vector<double>res(3);
        QEF_calculate({getGridVertex(g,0).x(),getGridVertex(g,0).y(),getGridVertex(g,0).z()},
                      {getGridVertex(g,7).x(),getGridVertex(g,7).y(),getGridVertex(g,7).z()},
                      sx.size(),sx,sy,sz,nx,ny,nz,res);



        if(vertex_in_grid(g,MeshKernel::iGameVertex(res[0],res[1],res[2]))){
            return MeshKernel::iGameVertex(res[0],res[1],res[2]);
        }
        else if(center_cnt>0){
            // cout<< "out"<<endl;
            return center_of_weight/center_cnt;
        }
        else{
            return qef_grid_pos_all_consider(g);
            /*cout <<"error"<<endl;
            return MeshKernel::iGameVertex{getGridVertex(g,0).x(),getGridVertex(g,0).y(),getGridVertex(g,0).z()}+
                   MeshKernel::iGameVertex{getGridVertex(g,7).x(),getGridVertex(g,7).y(),getGridVertex(g,7).z()};*/
        }
    };


    for(auto i : frame_grid_mp){
        if(i.second.grid_type==0){
            fprintf(file0, "v %lf %lf %lf\n", getGridVertex(i.first, 0).x(), getGridVertex(i.first, 0).y(),
                    getGridVertex(i.first, 0).z());
        }
    }

    for(auto i : frame_grid_mp){
        if(i.second.grid_type==1){
            fprintf(file1, "v %lf %lf %lf\n", getGridVertex(i.first, 0).x(), getGridVertex(i.first, 0).y(),
                    getGridVertex(i.first, 0).z());
        }
    }

    for(auto i : frame_grid_mp){
        if(i.second.grid_type==2){
            fprintf(file2, "v %lf %lf %lf\n", getGridVertex(i.first, 0).x(), getGridVertex(i.first, 0).y(),
                    getGridVertex(i.first, 0).z());
        }
    }

    unordered_map<grid,int ,grid_hash,grid_equal> build_mesh_vertex;
    vector<MeshKernel::iGameVertex>qef_result;
    std::vector<vector<int> > quad_face;

    string ss;
    int frame_grid_mp_cnt = 0;
    for(auto i : frame_grid_mp) {
        frame_grid_mp_cnt++;
        if(frame_grid_mp_cnt%1000==0)
            printf("%d/%d\n",frame_grid_mp_cnt,frame_grid_mp.size());
        auto v = getGridVertex(i.first, 0);
        if(i.second.grid_type==1 || i.second.grid_type==0) {
            auto tmp = get_neighbor(i.first);
            for(auto j : tmp){
                auto it = frame_grid_mp.find(j);
                if(it!= frame_grid_mp.end() && it->second.grid_type == 2){

                    auto v_gird = get_edge_neighbor(i.first,j);
                    for(auto vh : v_gird){
                        if(!build_mesh_vertex.count(vh)){
                            int id = (int)build_mesh_vertex.size()+1;
                            build_mesh_vertex.insert({vh ,id});
                        }
                    }
                    if (i.first < j) {
                        std::swap(v_gird[1], v_gird[3]);
                    }
                    quad_face.push_back({build_mesh_vertex[v_gird[0]],build_mesh_vertex[v_gird[1]] ,build_mesh_vertex[v_gird[2]],build_mesh_vertex[v_gird[3]]});
                }
            }
        }
    }
    qef_result.resize(build_mesh_vertex.size()+1);

    std::cout <<"start multi thread "<< endl;


    std::vector<std::shared_ptr<std::thread> > thread_pool(thread_num);
    for(int i=0;i<thread_num;i++){
        thread_pool[i] = make_shared<std::thread>([&](int now_id) {
            for(auto iter : build_mesh_vertex){
                if(iter.second%thread_num == now_id){
                    qef_result[iter.second] = qef_grid_pos(iter.first);
                    if(iter.second%500==0) // cut the output time usage
                        printf("succ %d/%d \n",iter.second, build_mesh_vertex.size());
                }
            }
        },i);
    }
    for(int i=0;i<thread_num;i++)
        thread_pool[i]->join();

    for(int i=1 ;i<=build_mesh_vertex.size();i++){
        fprintf(file,"v %lf %lf %lf\n",qef_result[i].x(),qef_result[i].y(),qef_result[i].z());
    }



    for(auto i : quad_face){
        if((qef_result[i[0]] -
            qef_result[i[2]]).norm() <
           (qef_result[i[1]] -
            qef_result[i[3]]).norm()) {
            ss += "f";
            for (int k: {0, 1, 2})
                ss += (" " + to_string(i[k]));
            ss += "\n";
            ss += "f";
            for (int k: {2, 3, 0})
                ss += (" " + to_string(i[k]));
            ss += "\n";
        }
        else{
            ss += "f";
            for (int k: {1, 2, 3})
                ss += (" " + to_string(i[k]));
            ss += "\n";
            ss += "f";
            for (int k: {3, 0, 1})
                ss += (" " + to_string(i[k]));
            ss += "\n";
        }

    }


    fprintf(file,"%s",ss.c_str());



    fclose(file0);
    fclose(file1);
    fclose(file2);

    unordered_map<grid,int,grid_hash,grid_equal> build_debug_grid;
    for(auto i : frame_grid_mp){
        fprintf(file10, "v %lf %lf %lf\n", getGridVertex(i.first, 0).x(), getGridVertex(i.first, 0).y(),
                getGridVertex(i.first, 0).z());
        int id = build_debug_grid.size()+1;
        build_debug_grid[i.first] = id;
    }
    for(auto i : frame_grid_mp){
        auto nei = get_neighbor(i.first);
        for(auto j : nei){
            if(i.first < j && frame_grid_mp.count(j)) {
                if(build_debug_grid.count(j)) {
                    fprintf(file10, "l %d %d\n", build_debug_grid[i.first], build_debug_grid[j]);
                }
            }
        }
    }
    fclose(file10);

    return 0;
}
// 1 2 3 4 5 6
// 5 1 2 3 7 8
//
//
// Created by te1t0ch1phead on 2022/6/6.
//

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
#include "MeshKernel/Mesh.h"


#include "GLKLib/GLKMatrixLib.h"


using namespace std;

#include "lib_impl.h"

MeshKernel::iGameVertex cgal_aabbtree_query(const std::vector<MeshKernel::iGameFace>& f_list,
                                            MeshKernel::iGameVertex v,
                                            std::shared_ptr<MeshKernel::SurfaceMesh>mesh,
                                            int & side,std::shared_ptr<CGALPolygon>cgalinmesh) {

    std::list<Triangle> triangles;

    for(auto f : f_list) {
        Point a(mesh->fast_iGameVertex[f.vh(0)].x(), mesh->fast_iGameVertex[f.vh(0)].y(), mesh->fast_iGameVertex[f.vh(0)].z());
        Point b(mesh->fast_iGameVertex[f.vh(1)].x(), mesh->fast_iGameVertex[f.vh(1)].y(), mesh->fast_iGameVertex[f.vh(1)].z());
        Point c(mesh->fast_iGameVertex[f.vh(2)].x(), mesh->fast_iGameVertex[f.vh(2)].y(), mesh->fast_iGameVertex[f.vh(2)].z());
        triangles.push_back(Triangle(a, b, c));
    }
    side = 1;
    Tree tree(triangles.begin(), triangles.end());
    Point_3 point_query(v.x(), v.y(), v.z());
//    cout <<"********************************\n";
//    cout <<"quv : "<< v.x()<<" "<< v.y()<<" "<< v.z()<<endl;
    auto pp = tree.closest_point_and_primitive(point_query);

    // TODO  太近了就搞一下 5555
//    if( (MeshKernel::iGameVertex(pp.first.x(), pp.first.y(), pp.first.z()) - v).norm() < myeps){
//        side = 1;
//        return MeshKernel::iGameVertex(pp.first.x(), pp.first.y(), pp.first.z());
//    }
    // TODO 引入判断主面
    vector<MeshKernel::iGameFace> near_face;
//    vector<double>dist_assist;
//    double dist0 = (MeshKernel::iGameVertex(pp.first.x(), pp.first.y(), pp.first.z()) - v).norm();
    bool is0 = false;
    bool is1 = false;
    for(int i=0;i<f_list.size();i++){
        if(cgal_vertex_triangle_dist(f_list[i], {pp.first.x(), pp.first.y(), pp.first.z()}, mesh) < myeps){
//            Point_3 fvh0 (mesh->fast_iGameVertex[f_list[i].vh(0)].x(),
//                       mesh->fast_iGameVertex[f_list[i].vh(0)].y(),
//                       mesh->fast_iGameVertex[f_list[i].vh(0)].z());
//            Point_3 fvh1 (mesh->fast_iGameVertex[f_list[i].vh(1)].x(),
//                          mesh->fast_iGameVertex[f_list[i].vh(1)].y(),
//                          mesh->fast_iGameVertex[f_list[i].vh(1)].z());
//            Point_3 fvh2 (mesh->fast_iGameVertex[f_list[i].vh(2)].x(),
//                          mesh->fast_iGameVertex[f_list[i].vh(2)].y(),
//                          mesh->fast_iGameVertex[f_list[i].vh(2)].z());
//
//            K::FT sqd = squared_distance(Plane_3(fvh0,fvh1,fvh2),
//                                         Point_3(v.x(),v.y(),v.z()));
//            double res = sqrt(sqd);
//            dist_assist.push_back(res);
//            near_face.push_back(f_list[i]);
//            MeshKernel::iGameVertex center = (mesh->fast_iGameVertex[f_list[i].vh(0)]
//                    + mesh->fast_iGameVertex[f_list[i].vh(1)] +  mesh->fast_iGameVertex[f_list[i].vh(2)])/3;
//            MeshKernel::iGameVertex direct = (center - MeshKernel::iGameVertex(pp.first.x(), pp.first.y(), pp.first.z())).normalize();
//            direct *= dist0;
//            MeshKernel::iGameVertex assist_v = MeshKernel::iGameVertex(pp.first.x(), pp.first.y(), pp.first.z()) + direct;

             MeshKernel::iGameVertex normal = ((mesh->fast_iGameVertex[f_list[i].vh(1)] - mesh->fast_iGameVertex[f_list[i].vh(0)])
                       % (mesh->fast_iGameVertex[f_list[i].vh(2)] - mesh->fast_iGameVertex[f_list[i].vh(0)])).normalize();
            if(normal * (v - MeshKernel::iGameVertex(pp.first.x(),pp.first.y(),pp.first.z()) ).normalize() >= 0 )
                is0 = true;
            else
                is1 = true;
//            auto nown = ((mesh->vertices(f_list[i].vh(1)) - mesh->vertices(f_list[i].vh(0)))
//                         % (mesh->vertices(f_list[i].vh(2)) - mesh->vertices(f_list[i].vh(0)))).normalize();            cout << mesh->vertices(f_list[i].vh(0)).x()<<" "<< mesh->vertices(f_list[i].vh(0)).y()<<" "
//            <<mesh->vertices(f_list[i].vh(0)).z()<<endl;
//            cout << mesh->vertices(f_list[i].vh(1)).x()<<" "<< mesh->vertices(f_list[i].vh(1)).y()<<" "
//                 <<mesh->vertices(f_list[i].vh(1)).z()<<endl;
//            cout << mesh->vertices(f_list[i].vh(2)).x()<<" "<< mesh->vertices(f_list[i].vh(2)).y()<<" "
//                 <<mesh->vertices(f_list[i].vh(2)).z()<<endl;
//            cout <<"nown: " <<nown.x()<<" "<< nown.y()<<" "<< nown.z() << endl;
            //cnt++;
//            near_face.push_back(f_list[i]);
//            dist_assist.push_back((assist_v-v).norm());
        }
    }
    if(is0 && is1){
        if(cgalinmesh->inMesh(v))
            side = 1;
        else
            side = 0;
    }
    else if(is0){
        side = 0;
    }
    else{
        side = 1;
    }



    // TODO 内外debug


//    MeshKernel::iGameFace nearest_face = near_face[min_element(dist_assist.begin(),dist_assist.end()) - dist_assist.begin()];
//
//    MeshKernel::iGameVertex normal = ((mesh->fast_iGameVertex[nearest_face.vh(1)] - mesh->fast_iGameVertex[nearest_face.vh(0)])
//            % (mesh->fast_iGameVertex[nearest_face.vh(2)] - mesh->fast_iGameVertex[nearest_face.vh(0)])).normalize();
//
//    if(normal * (v - MeshKernel::iGameVertex(pp.first.x(),pp.first.y(),pp.first.z()) ).normalize() >= 0 ){
//        side = 0;
//        if(cgalinmesh){
//            if(cgalinmesh->inMesh(v)){
//                cout <<"BUG v : v "<< v.x()<<" "<< v.y()<<" "<<v.z()<<endl;
//                cout <<"BUG near : v "<< pp.first.x()<<" "<< pp.first.y()<<" "<<pp.first.z()<<endl;
//            }
//        }
//    }
//    else {
//        side = 1;
//    }
    //    int near_face_id = 0;
//    for(int i=1;i<near_face.size();i++){
//            [i] =
//    }


    //cout <<" cnt " <<cnt << " "<< pp.first.x()<<" "<< pp.first.y()<<" "<< pp.first.z() << endl;
//    if(cnt>=3){
//        side = 0;
//    }
//
//    normal = normal /cnt;
//
//    // cout <<"normal:" <<normal.x()<<" "<< normal.y() <<" "<< normal.z()<<endl;
//    if(normal * (v - MeshKernel::iGameVertex(pp.first.x(),pp.first.y(),pp.first.z()) ).normalize() > 0 ){
//        //       cout <<"side "<< 0 << endl;
//        side = 0;
//    }
//    else {
////        if(v.z()>=1.0) {
////          //  cout << "side " << 1 << endl;
////            cout << v.x() << " " << v.y() << " " << v.z() << endl;
////        }
//    }
    //side = pp.second->supporting_plane().oriented_side(point_query);

    return MeshKernel::iGameVertex(pp.first.x(), pp.first.y(), pp.first.z());
}


double cgal_vertex_triangle_dist(MeshKernel::iGameFace f, MeshKernel::iGameVertex v, std::shared_ptr<MeshKernel::SurfaceMesh>mesh) {
    Point a(mesh->fast_iGameVertex.at(f.vh(0)).x(), mesh->fast_iGameVertex.at(f.vh(0)).y(), mesh->fast_iGameVertex.at(f.vh(0)).z());
    Point b(mesh->fast_iGameVertex.at(f.vh(1)).x(), mesh->fast_iGameVertex.at(f.vh(1)).y(), mesh->fast_iGameVertex.at(f.vh(1)).z());
    Point c(mesh->fast_iGameVertex.at(f.vh(2)).x(), mesh->fast_iGameVertex.at(f.vh(2)).y(), mesh->fast_iGameVertex.at(f.vh(2)).z());
    Point point_query(v.x(), v.y(), v.z());
    K::FT sqd = squared_distance(Triangle(a, b, c),point_query);
    double res = sqrt(sqd);
    return res;

}


double cgal_vertex_plane_dist(MeshKernel::iGameFace f, MeshKernel::iGameVertex v, std::shared_ptr<MeshKernel::SurfaceMesh>mesh) {
    Point a(mesh->fast_iGameVertex.at(f.vh(0)).x(), mesh->fast_iGameVertex.at(f.vh(0)).y(), mesh->fast_iGameVertex.at(f.vh(0)).z());
    Point b(mesh->fast_iGameVertex.at(f.vh(1)).x(), mesh->fast_iGameVertex.at(f.vh(1)).y(), mesh->fast_iGameVertex.at(f.vh(1)).z());
    Point c(mesh->fast_iGameVertex.at(f.vh(2)).x(), mesh->fast_iGameVertex.at(f.vh(2)).y(), mesh->fast_iGameVertex.at(f.vh(2)).z());
    Point point_query(v.x(), v.y(), v.z());
    K::FT sqd = squared_distance(Plane_3(a, b, c),point_query);
    double res = sqrt(sqd);
    return res;
}




bool CGALPlaneSegIntersect(MeshKernel::iGameFace f, MeshKernel::iGameVertex v1,MeshKernel::iGameVertex v2,std::shared_ptr<MeshKernel::SurfaceMesh>mesh){

    Segment_3 seg(Point_3(v1.x(),v1.y(),v1.z()), Point_3(v2.x(),v2.y(),v2.z()));
    Point_3 p0(mesh->fast_iGameVertex.at(f.vh(0)).x(),mesh->fast_iGameVertex.at(f.vh(0)).y(),mesh->fast_iGameVertex.at(f.vh(0)).z());
    Point_3 p1(mesh->fast_iGameVertex.at(f.vh(1)).x(),mesh->fast_iGameVertex.at(f.vh(1)).y(),mesh->fast_iGameVertex.at(f.vh(1)).z());
    Point_3 p2(mesh->fast_iGameVertex.at(f.vh(2)).x(),mesh->fast_iGameVertex.at(f.vh(2)).y(),mesh->fast_iGameVertex.at(f.vh(2)).z());
    Plane_3 plane3 (p0,p1,p2);
    CGAL::cpp11::result_of<Intersect_3(Segment_3, Segment_3)>::type
            result = intersection(seg, plane3);
    if(result)
        return true;
    else
        return false;
}

void QEF_calculate(std::vector<double>minP,std::vector<double>maxP, int pntNum, vector<double> sx, vector<double> sy, vector<double> sz,
                   vector<double> nx, vector<double> ny, vector<double> nz, vector<double> &pp) {
    double proj, scale;
    double criterion = 0.05;
    int i, j, k;
    double **A, **UU, **VV, **UUT, **VVT;
    double *B, *X;


    //---------------------------------------------------------------------------
    //	Preparation

    //---------------------------------------------------------------------------
    pp[0] = pp[1] = pp[2] = 0.0;
    for (int index = 0; index < pntNum; index++) {
        pp[0] += sx[index];
        pp[1] += sy[index];
        pp[2] += sz[index];
    }
    pp[0] = pp[0] / (double) pntNum;
    pp[1] = pp[1] / (double) pntNum;
    pp[2] = pp[2] / (double) pntNum;
    //---------------------------------------------------------------------------
    GLKMatrixLib::CreateMatrix(A, 3, 3);
    B = new double[3];
    X = new double[3];
    GLKMatrixLib::CreateMatrix(UU, 3, 3);
    GLKMatrixLib::CreateMatrix(VV, 3, 3);
    GLKMatrixLib::CreateMatrix(UUT, 3, 3);
    GLKMatrixLib::CreateMatrix(VVT, 3, 3);
    //---------------------------------------------------------------------------
    B[0] = B[1] = B[2] = X[0] = X[1] = X[2] = 0.0;
    for (k = 0; k < pntNum; k++) {
        proj = (sx[k] - pp[0]) * nx[k] + (sy[k] - pp[1]) * ny[k] + (sz[k] - pp[2]) * nz[k];
        B[0] += proj * nx[k];
        B[1] += proj * ny[k];
        B[2] += proj * nz[k];

        A[0][0] += nx[k] * nx[k];
        A[0][1] += nx[k] * ny[k];
        A[0][2] += nx[k] * nz[k];
        A[1][0] += ny[k] * nx[k];
        A[1][1] += ny[k] * ny[k];
        A[1][2] += ny[k] * nz[k];
        A[2][0] += nz[k] * nx[k];
        A[2][1] += nz[k] * ny[k];
        A[2][2] += nz[k] * nz[k];
    }

    //---------------------------------------------------------------------------
    //	Singular Value Decomposition
    GLKMatrixLib::SingularValueDecomposition(A, 3, 3, UU, VVT);
    GLKMatrixLib::Transpose(UU, 3, 3, UUT);
    GLKMatrixLib::Transpose(VVT, 3, 3, VV);
    double maxFactor = (fabs(A[0][0]) > fabs(A[1][1])) ? (A[0][0]) : (A[1][1]);
    maxFactor = (fabs(maxFactor) > fabs(A[2][2])) ? (maxFactor) : (A[2][2]);
    if (fabs(maxFactor) < 1.0e-6) {
        for (i = 0; i < 3; i++) {
            for (j = 0; j < 3; j++) {
                A[i][j] = 0.0;
            }
        }
    } else {
        for (i = 0; i < 3; i++) {
            for (j = 0; j < 3; j++) {
                if (i != j) {
                    A[i][j] = 0.0;
                } else {
                    if (fabs(A[i][j] / maxFactor) < criterion)
                        A[i][j] = 0.0;
                    else
                        A[i][j] = 1.0 / A[i][j];
                }
            }
        }
    }
    GLKMatrixLib::Mul(UUT, B, 3, 3, X);
    GLKMatrixLib::Mul(A, X, 3, 3, B);
    GLKMatrixLib::Mul(VV, B, 3, 3, X);
    //-----------------------------------------------------------------
    //	truncate the update vector and update node position
    scale = 1.0;
    if (fabs(X[0]) > 1.0e-5 && (pp[0] + X[0] * scale) > maxP[0]) scale = (maxP[0] - pp[0]) / X[0];
    if (fabs(X[1]) > 1.0e-5 && (pp[1] + X[1] * scale) > maxP[1]) scale = (maxP[1] - pp[1]) / X[1];
    if (fabs(X[2]) > 1.0e-5 && (pp[2] + X[2] * scale) > maxP[2]) scale = (maxP[2] - pp[2]) / X[2];
    if (fabs(X[0]) > 1.0e-5 && (pp[0] + X[0] * scale) < minP[0]) scale = (minP[0] - pp[0]) / X[0];
    if (fabs(X[1]) > 1.0e-5 && (pp[1] + X[1] * scale) < minP[1]) scale = (minP[1] - pp[1]) / X[1];
    if (fabs(X[2]) > 1.0e-5 && (pp[2] + X[2] * scale) < minP[2]) scale = (minP[2] - pp[2]) / X[2];
    pp[0] = pp[0] + X[0] * scale;
    pp[1] = pp[1] + X[1] * scale;
    pp[2] = pp[2] + X[2] * scale;

    //---------------------------------------------------------------------------
    //	Free the memory
    GLKMatrixLib::DeleteMatrix(A, 3, 3);
    delete[]B;
    delete[]X;
    GLKMatrixLib::DeleteMatrix(UU, 3, 3);
    GLKMatrixLib::DeleteMatrix(VV, 3, 3);
    GLKMatrixLib::DeleteMatrix(UUT, 3, 3);
    GLKMatrixLib::DeleteMatrix(VVT, 3, 3);
    double blendingFactor = 0.25;
    /*pp[0] = (1.0 - blendingFactor)*pp[0] + blendingFactor*g.get_center().x();
    pp[1] = (1.0 - blendingFactor)*pp[1] + blendingFactor*g.get_center().y();
    pp[2] = (1.0 - blendingFactor)*pp[2] + blendingFactor*g.get_center().z();*/
}


bool vertex_field_state(MeshKernel::iGameFaceHandle fh , MeshKernel::iGameVertex v,std::shared_ptr<MeshKernel::SurfaceMesh>mesh){
    MeshKernel::iGameFace f = mesh->fast_iGameFace.at(fh);
    vector<Point>cgal_vertex_list;
    cgal_vertex_list.emplace_back(mesh->fast_iGameVertex.at(f.vh(0)).x(), mesh->fast_iGameVertex.at(f.vh(0)).y(), mesh->fast_iGameVertex.at(f.vh(0)).z());
    cgal_vertex_list.emplace_back(mesh->fast_iGameVertex.at(f.vh(1)).x(), mesh->fast_iGameVertex.at(f.vh(1)).y(), mesh->fast_iGameVertex.at(f.vh(1)).z());
    cgal_vertex_list.emplace_back(mesh->fast_iGameVertex.at(f.vh(2)).x(), mesh->fast_iGameVertex.at(f.vh(2)).y(), mesh->fast_iGameVertex.at(f.vh(2)).z());

    Point point_query(v.x(), v.y(), v.z());
    K::FT sqd = squared_distance(Triangle(cgal_vertex_list[0], cgal_vertex_list[1], cgal_vertex_list[2]),point_query);
    double res = sqrt(sqd);
    K::FT sqd2 = squared_distance(Triangle(cgal_vertex_list[0], cgal_vertex_list[1], cgal_vertex_list[2]).supporting_plane(),point_query);
    double res2 = sqrt(sqd2);
    return res <= mesh->faces(fh).move_dist;


    if(abs(res - res2 ) < myeps){
        return res <= mesh->faces(fh).move_dist;
    }






    for(int i=0;i<3;i++){
        double dist_vtov = (mesh->fast_iGameVertex.at(f.vh(i)) - v).norm();
        if(abs(dist_vtov -res) < myeps){ //  ru guo shi yi ge dian

            std::unordered_set<MeshKernel::iGameFaceHandle> neighbor_face = mesh->FastNeighborFhOfVertex_[f.vh(i)];
            bool in_field = true;

            for(auto nei_fh : neighbor_face){
                if(nei_fh !=fh){

                    MeshKernel::iGameVertex nei_fv0= mesh->fast_iGameVertex[mesh->fast_iGameFace[nei_fh].vh(0)];
                    MeshKernel::iGameVertex nei_fv1= mesh->fast_iGameVertex[mesh->fast_iGameFace[nei_fh].vh(1)];
                    MeshKernel::iGameVertex nei_fv2= mesh->fast_iGameVertex[mesh->fast_iGameFace[nei_fh].vh(2)];

                    Point_3 cgal_nei_fv0(nei_fv0.x(),nei_fv0.y(),nei_fv0.z());
                    Point_3 cgal_nei_fv1(nei_fv1.x(),nei_fv1.y(),nei_fv1.z());
                    Point_3 cgal_nei_fv2(nei_fv2.x(),nei_fv2.y(),nei_fv2.z());

                    Triangle tr(cgal_nei_fv0,cgal_nei_fv1,cgal_nei_fv2);

                    if(abs(sqrt(squared_distance(tr,point_query)) -  sqrt(squared_distance(tr.supporting_plane(),point_query)))<myeps)
                        return false;
                }
            }


            for(auto nei_fh : neighbor_face){
                MeshKernel::iGameVertex normal =
                        ((mesh->fast_iGameVertex[mesh->fast_iGameFace[nei_fh].vh(1)]
                        - mesh->fast_iGameVertex[mesh->fast_iGameFace[nei_fh].vh(0)])
                        % (mesh->fast_iGameVertex[mesh->fast_iGameFace[nei_fh].vh(2)]
                        - mesh->fast_iGameVertex[mesh->fast_iGameFace[nei_fh].vh(0)])).normalize();
                MeshKernel::iGameVertex center = (mesh->fast_iGameVertex[mesh->fast_iGameFace[nei_fh].vh(0)]
                        + mesh->fast_iGameVertex[mesh->fast_iGameFace[nei_fh].vh(1)]
                        + mesh->fast_iGameVertex[mesh->fast_iGameFace[nei_fh].vh(2)])/3;

                if(normal * (v - MeshKernel::iGameVertex(center.x(),center.y(),center.z()) ).normalize() < 0 ){ // yaoqiu juli < x
                    Point_3 neighbor_face_v0(mesh->fast_iGameVertex[mesh->fast_iGameFace[nei_fh].vh(0)].x(),
                                             mesh->fast_iGameVertex[mesh->fast_iGameFace[nei_fh].vh(0)].y(),
                                             mesh->fast_iGameVertex[mesh->fast_iGameFace[nei_fh].vh(0)].z());
                    Point_3 neighbor_face_v1(mesh->fast_iGameVertex[mesh->fast_iGameFace[nei_fh].vh(1)].x(),
                                             mesh->fast_iGameVertex[mesh->fast_iGameFace[nei_fh].vh(1)].y(),
                                             mesh->fast_iGameVertex[mesh->fast_iGameFace[nei_fh].vh(1)].z());
                    Point_3 neighbor_face_v2(mesh->fast_iGameVertex[mesh->fast_iGameFace[nei_fh].vh(2)].x(),
                                             mesh->fast_iGameVertex[mesh->fast_iGameFace[nei_fh].vh(2)].y(),
                                             mesh->fast_iGameVertex[mesh->fast_iGameFace[nei_fh].vh(2)].z());
                    Plane_3 neighbor_plane(neighbor_face_v0,neighbor_face_v1,neighbor_face_v2);
                    double dist_ptov = sqrt(squared_distance(neighbor_plane,point_query));
                    if(dist_ptov > mesh->faces(nei_fh).move_dist){
                        in_field = false;
                    }
                }
            }
            return in_field;
        }
    }






    for(int i=0;i<3;i++){
        MeshKernel::iGameEdge e =  mesh->fast_iGameEdge.at(f.eh(i));

        Point_3 cgalv0( mesh->fast_iGameVertex[e.vh(0)].x(),
                        mesh->fast_iGameVertex[e.vh(0)].y(),
                        mesh->fast_iGameVertex[e.vh(0)].z());
        Point_3 cgalv1(mesh->fast_iGameVertex[e.vh(1)].x(),
                       mesh->fast_iGameVertex[e.vh(1)].y(),
                       mesh->fast_iGameVertex[e.vh(1)].z());

        Segment_3 segment3(cgalv0,cgalv1);

        double dist_vtoe = squared_distance(segment3,point_query);
        if(abs(dist_vtoe -res) < myeps) { //  ru guo shi yi ge dian
            std::unordered_set<MeshKernel::iGameFaceHandle> neighbor_face = mesh->FastNeighborFhOfEdge_[f.eh(i)];

            for(auto nei_fh : neighbor_face){
                if(nei_fh !=fh){

                    MeshKernel::iGameVertex nei_fv0= mesh->fast_iGameVertex[mesh->fast_iGameFace[nei_fh].vh(0)];
                    MeshKernel::iGameVertex nei_fv1= mesh->fast_iGameVertex[mesh->fast_iGameFace[nei_fh].vh(1)];
                    MeshKernel::iGameVertex nei_fv2= mesh->fast_iGameVertex[mesh->fast_iGameFace[nei_fh].vh(2)];

                    Point_3 cgal_nei_fv0(nei_fv0.x(),nei_fv0.y(),nei_fv0.z());
                    Point_3 cgal_nei_fv1(nei_fv1.x(),nei_fv1.y(),nei_fv1.z());
                    Point_3 cgal_nei_fv2(nei_fv2.x(),nei_fv2.y(),nei_fv2.z());

                    Triangle tr(cgal_nei_fv0,cgal_nei_fv1,cgal_nei_fv2);

                    if(abs(sqrt(squared_distance(tr,point_query)) -  sqrt(squared_distance(tr.supporting_plane(),point_query)))<myeps)
                        return false;
                }
            }

            bool in_field = true;
            for(auto nei_fh : neighbor_face){
                MeshKernel::iGameVertex normal =
                        ((mesh->fast_iGameVertex[mesh->fast_iGameFace[nei_fh].vh(1)]
                          - mesh->fast_iGameVertex[mesh->fast_iGameFace[nei_fh].vh(0)])
                         % (mesh->fast_iGameVertex[mesh->fast_iGameFace[nei_fh].vh(2)]
                            - mesh->fast_iGameVertex[mesh->fast_iGameFace[nei_fh].vh(0)])).normalize();

                MeshKernel::iGameVertex center = (mesh->fast_iGameVertex[mesh->fast_iGameFace[nei_fh].vh(0)]
                                                  + mesh->fast_iGameVertex[mesh->fast_iGameFace[nei_fh].vh(1)]
                                                  + mesh->fast_iGameVertex[mesh->fast_iGameFace[nei_fh].vh(2)])/3;

                if(normal * (v -  center ).normalize() < 0 ){ // yaoqiu juli < x
                    Point_3 neighbor_face_v0(mesh->fast_iGameVertex[mesh->fast_iGameFace[nei_fh].vh(0)].x(),
                                             mesh->fast_iGameVertex[mesh->fast_iGameFace[nei_fh].vh(0)].y(),
                                             mesh->fast_iGameVertex[mesh->fast_iGameFace[nei_fh].vh(0)].z());
                    Point_3 neighbor_face_v1(mesh->fast_iGameVertex[mesh->fast_iGameFace[nei_fh].vh(1)].x(),
                                             mesh->fast_iGameVertex[mesh->fast_iGameFace[nei_fh].vh(1)].y(),
                                             mesh->fast_iGameVertex[mesh->fast_iGameFace[nei_fh].vh(1)].z());
                    Point_3 neighbor_face_v2(mesh->fast_iGameVertex[mesh->fast_iGameFace[nei_fh].vh(2)].x(),
                                             mesh->fast_iGameVertex[mesh->fast_iGameFace[nei_fh].vh(2)].y(),
                                             mesh->fast_iGameVertex[mesh->fast_iGameFace[nei_fh].vh(2)].z());
                    Plane_3 neighbor_plane(neighbor_face_v0,neighbor_face_v1,neighbor_face_v2);
                    double dist_ptov = sqrt(squared_distance(neighbor_plane,point_query));
                    if(dist_ptov >  mesh->faces(nei_fh).move_dist){
                        in_field = false;
                    }
                }
            }
            return in_field;
        }
    }

    return false;
   // return res <= mesh->faces(fh).move_dist;
}

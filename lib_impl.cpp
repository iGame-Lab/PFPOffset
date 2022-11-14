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


using namespace std;

#include "lib_impl.h"


double cgal_vertex_triangle_dist(MeshKernel::iGameFace f, MeshKernel::iGameVertex v, std::shared_ptr<MeshKernel::SurfaceMesh>mesh) {
    Point a(mesh->fast_iGameVertex.at(f.vh(0)).x(), mesh->fast_iGameVertex.at(f.vh(0)).y(), mesh->fast_iGameVertex.at(f.vh(0)).z());
    Point b(mesh->fast_iGameVertex.at(f.vh(1)).x(), mesh->fast_iGameVertex.at(f.vh(1)).y(), mesh->fast_iGameVertex.at(f.vh(1)).z());
    Point c(mesh->fast_iGameVertex.at(f.vh(2)).x(), mesh->fast_iGameVertex.at(f.vh(2)).y(), mesh->fast_iGameVertex.at(f.vh(2)).z());
    Point point_query(v.x(), v.y(), v.z());
    K::FT sqd = squared_distance(K::Triangle_3(a, b, c),point_query);
    double res = sqrt(sqd);
    return res;

}




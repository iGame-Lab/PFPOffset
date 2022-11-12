//
// Created by te1t0ch1phead on 2022/6/6.
//

#ifndef THICKEN2_LIB_IMPL_H
#define THICKEN2_LIB_IMPL_H

#define myeps 1e-6


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



#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>

#include <CGAL/point_generators_3.h>
#include <CGAL/squared_distance_3.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>
#include <CGAL/enum.h>
#include <vector>
#include <fstream>
#include <limits>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point_3;
typedef K::Segment_3 Segment_3;
typedef K::Plane_3 Plane_3;
typedef K::Intersect_3 Intersect_3;
typedef K::Point_3 Point;
typedef CGAL::Polyhedron_3<K> Polyhedron;
typedef K::Ray_3 Ray;
typedef K::Line_3 Line;
typedef K::Triangle_3 Triangle;
typedef std::list<Triangle>::iterator Iterator;
typedef CGAL::AABB_triangle_primitive<K, Iterator> Primitive;
typedef CGAL::AABB_traits<K, Primitive> AABB_triangle_traits;
typedef CGAL::AABB_tree<AABB_triangle_traits> Tree;
typedef K::Tetrahedron_3  Tetrahedron;
namespace PMP = CGAL::Polygon_mesh_processing;
typedef CGAL::Exact_predicates_exact_constructions_kernel K2;


using namespace std;

struct triangle_hash{
    size_t operator () (Triangle x) const {
        return hash_value(x.vertex(0)) ^ (hash_value(x.vertex(1))) ^ (hash_value(x.vertex(2)));
    }};
struct triangle_equal
{
    bool operator() (Triangle a,  Triangle b) const {
        return a.vertex(0) == b.vertex(0)  &&  a.vertex(1) == b.vertex(1) &&  a.vertex(2) == b.vertex(2);
    }
};

class CGALPolygon {
    CGAL::Polyhedron_3<K2> * poly;
    CGAL::Side_of_triangle_mesh<CGAL::Polyhedron_3<K2>, K2 > * inside;

public:
    CGALPolygon() {}

    CGALPolygon(const shared_ptr <MeshKernel::SurfaceMesh> &mesh) {
        std::vector<K2::Point_3> ps;
        std::vector<std::vector<std::size_t> > fs;
        for (int i = 0; i < mesh->VertexSize(); i++) {
            ps.emplace_back(mesh->vertices(MeshKernel::iGameVertexHandle(i)).x(),
                            mesh->vertices(MeshKernel::iGameVertexHandle(i)).y(),
                            mesh->vertices(MeshKernel::iGameVertexHandle(i)).z());

        }
        for (int i = 0; i < mesh->FaceSize(); i++) {
            fs.push_back({static_cast<unsigned long long>(mesh->faces(MeshKernel::iGameFaceHandle(i)).vh(0).idx()),
                          static_cast<unsigned long long>(mesh->faces(MeshKernel::iGameFaceHandle(i)).vh(1).idx()),
                          static_cast<unsigned long long>(mesh->faces(MeshKernel::iGameFaceHandle(i)).vh(2).idx())});
        }
        poly = new CGAL::Polyhedron_3<K2>();
        PMP::polygon_soup_to_polygon_mesh(ps, fs, *poly, CGAL::parameters::all_default());
        inside = new CGAL::Side_of_triangle_mesh<CGAL::Polyhedron_3<K2>, K2>(*poly);



    };

    bool inMesh(MeshKernel::iGameVertex v) {
       // CGAL::Side_of_triangle_mesh<CGAL::Polyhedron_3<K2>, K2 > inside(poly);
        CGAL::Bounded_side res = (*inside)(K2::Point_3(v.x(), v.y(), v.z()));
        if (res == CGAL::ON_BOUNDED_SIDE)
            return true;
        return false;
    }
};

MeshKernel::iGameVertex cgal_aabbtree_query(const std::vector<MeshKernel::iGameFace>& f_list,
                                            MeshKernel::iGameVertex v,
                                            std::shared_ptr<MeshKernel::SurfaceMesh>mesh,
                                            int & side,std::shared_ptr<CGALPolygon>cgalinmesh);

double cgal_vertex_triangle_dist(MeshKernel::iGameFace f, MeshKernel::iGameVertex v, std::shared_ptr<MeshKernel::SurfaceMesh>mesh) ;


double cgal_vertex_plane_dist(MeshKernel::iGameFace f, MeshKernel::iGameVertex v, std::shared_ptr<MeshKernel::SurfaceMesh>mesh);

bool CGALPlaneSegIntersect(MeshKernel::iGameFace f, MeshKernel::iGameVertex v1,MeshKernel::iGameVertex v2,std::shared_ptr<MeshKernel::SurfaceMesh>mesh);

void QEF_calculate(std::vector<double>minP,std::vector<double>maxP, int pntNum, vector<double> sx, vector<double> sy, vector<double> sz,
                   vector<double> nx, vector<double> ny, vector<double> nz, vector<double> &pp);



bool vertex_field_state(MeshKernel::iGameFaceHandle fh , MeshKernel::iGameVertex v,std::shared_ptr<MeshKernel::SurfaceMesh>mesh);


#endif //THICKEN2_LIB_IMPL_H

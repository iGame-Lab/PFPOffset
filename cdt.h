//
// Created by rainbowwing on 2023/8/25.
//

#ifndef THICKEN2_CDT_H
#define THICKEN2_CDT_H

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Polygon_2.h>
#include <vector>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Fuzzy_iso_box.h>
typedef CGAL::Arr_segment_traits_2<K2> Traits_2;
typedef CGAL::Arrangement_2<Traits_2> Arrangement_2;

typedef CGAL::Constrained_Delaunay_triangulation_2<K2> CDT;
typedef CGAL::Polygon_2<K2> Polygon_2;


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef K::Point_3 Point;
typedef CGAL::Search_traits_3<K> STraits;
typedef CGAL::Fuzzy_sphere<STraits> Fuzzy_circle;
typedef CGAL::Kd_tree<STraits> Kd_tree;

typedef CGAL::Search_traits_3<K2> STraitsK2;
typedef CGAL::Fuzzy_sphere<STraitsK2> Fuzzy_circle_K2;
typedef CGAL::Kd_tree<STraitsK2> Kd_tree_K2;
typedef CGAL::Fuzzy_iso_box<STraitsK2> FuzzyBoxK2;

typedef CGAL::Delaunay_triangulation_3<K> Delaunay3D;

typedef CGAL::Constrained_Delaunay_triangulation_2<K2> CDT;



vector<vector<K2::Point_3> > CGAL_CDT_NEW(vector<K2::Point_3> sorted_bound_vertex, vector<K2::Segment_3> cs,K2::Triangle_3 origin_face) {

    //CGAL::make_conforming_Delaunay_2();
    // cout << sorted_bound_vertex.size() <<" "<< cs.size() << endl;
    K2::Point_3 base_point_3d = origin_face.vertex(0);
    CDT cdt;

    K2::Vector_3 origin_face_v0 = origin_face.vertex(1) - origin_face.vertex(0);
    K2::Vector_3 origin_face_v1 = origin_face.vertex(2) - origin_face.vertex(0);
    K2::Vector_3 new_direct = CGAL::cross_product(origin_face_v0,origin_face_v1);
    new_direct = new_direct /  CGAL::Epeck::FT(Point_K2_to_iGameVertex(K2::Point_3(0,0,0) +  new_direct).norm());
    K2::Vector_3 X_axis = origin_face_v0 / CGAL::Epeck::FT(Point_K2_to_iGameVertex(K2::Point_3(0,0,0) + origin_face_v0).norm());
    K2::Vector_3 Y_axis = CGAL::cross_product(new_direct,X_axis);

    std::vector<K2::Point_2> polygon_vertices;
    Arrangement_2 arrangement;
    for(int i=0;i<sorted_bound_vertex.size();i++){
        CGAL::Epeck::FT x = (sorted_bound_vertex[i]-base_point_3d) * X_axis;
        CGAL::Epeck::FT y = (sorted_bound_vertex[i]-base_point_3d) * Y_axis;
        //CGAL::insert(arrangement,K2::Point_2(x,y));
        cdt.insert(K2::Point_2(x,y));
        polygon_vertices.emplace_back(x,y);
    }

    Polygon_2 polygon(polygon_vertices.begin(), polygon_vertices.end());

    std::vector<K2::Point_2> bounded_polygon_vertices;
    for(int i=0;i<3;i++){
        CGAL::Epeck::FT x = (origin_face.vertex(i)-base_point_3d) * X_axis;
        CGAL::Epeck::FT y = (origin_face.vertex(i)-base_point_3d) * Y_axis;
        bounded_polygon_vertices.emplace_back(x,y);
    }
    Polygon_2 bounded_polygon(bounded_polygon_vertices.begin(),bounded_polygon_vertices.end());


    for(int i=0;i<3;i++){
        CGAL::Epeck::FT x0 = (origin_face.vertex(i)-base_point_3d) * X_axis;
        CGAL::Epeck::FT y0 = (origin_face.vertex(i)-base_point_3d) * Y_axis;
        CGAL::Epeck::FT x1 = (origin_face.vertex((i+1)%3)-base_point_3d) * X_axis;
        CGAL::Epeck::FT y1 = (origin_face.vertex((i+1)%3)-base_point_3d) * Y_axis;
        K2::Segment_2 seg(K2::Point_2(x0,y0) ,K2::Point_2(x1,y1));
        CGAL::insert(arrangement,seg);
        //cdt.insert_constraint(K2::Point_2(x0,y0) ,K2::Point_2(x1,y1));
    }
    for(auto s: cs){
        CGAL::Epeck::FT x0 = (s.vertex(0)-base_point_3d) * X_axis;
        CGAL::Epeck::FT y0 = (s.vertex(0)-base_point_3d) * Y_axis;
        CGAL::Epeck::FT x1 = (s.vertex(1)-base_point_3d) * X_axis;
        CGAL::Epeck::FT y1 = (s.vertex(1)-base_point_3d) * Y_axis;
        K2::Segment_2 seg(K2::Point_2(x0,y0) ,K2::Point_2(x1,y1));
        CGAL::insert(arrangement,seg);
        //cdt.insert_constraint(K2::Point_2(x0,y0) ,K2::Point_2(x1,y1));
    }
    for (auto eit = arrangement.edges_begin(); eit != arrangement.edges_end(); ++eit) {
        K2::Segment_2 seg = eit->curve();
        K2::Point_2 source = seg.source();
        K2::Point_2 target = seg.target();
        cdt.insert_constraint(source,target);
//        std::cout << "Segment: ((" << source.x() << ", " << source.y() << "), ("
//                  << target.x() << ", " << target.y() << "))\n";
    }


    K2::Vector_3 vec = cross_product(origin_face.vertex(1)-origin_face.vertex(0),
                                     origin_face.vertex(2)-origin_face.vertex(0)
    );

    vector<vector<K2::Point_3> > ret;
    for (CDT::Finite_faces_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) {
        K2::Triangle_2 t = cdt.triangle(fit);
        if(bounded_polygon.has_on_bounded_side(CGAL::centroid(t))) {
            vector<K2::Point_3>tmp(3);
            for(int j=0;j<3;j++){
                tmp[j] = base_point_3d + t[j].x()*X_axis + t[j].y()*Y_axis;
            }
            K2::Vector_3 vec1 = cross_product(tmp[1] - tmp[0],tmp[2] - tmp[0]);
            if(vec * vec1 < CGAL::Epeck::FT(0))
                swap(tmp[1],tmp[2]);

            ret.push_back(tmp);
        }
    }


    return ret;
}

#endif //THICKEN2_CDT_H

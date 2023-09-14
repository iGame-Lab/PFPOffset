//
// Created by rainbowwing on 2023/8/25.
//

#ifndef THICKEN2_GEOMETRIC_DATA_H
#define THICKEN2_GEOMETRIC_DATA_H
double cgal_vertex_triangle_dist(MeshKernel::iGameFace f, MeshKernel::iGameVertex v, std::shared_ptr<MeshKernel::SurfaceMesh>mesh) {
    Point a(mesh->fast_iGameVertex.at(f.vh(0)).x(), mesh->fast_iGameVertex.at(f.vh(0)).y(), mesh->fast_iGameVertex.at(f.vh(0)).z());
    Point b(mesh->fast_iGameVertex.at(f.vh(1)).x(), mesh->fast_iGameVertex.at(f.vh(1)).y(), mesh->fast_iGameVertex.at(f.vh(1)).z());
    Point c(mesh->fast_iGameVertex.at(f.vh(2)).x(), mesh->fast_iGameVertex.at(f.vh(2)).y(), mesh->fast_iGameVertex.at(f.vh(2)).z());
    Point point_query(v.x(), v.y(), v.z());
    K::FT sqd = squared_distance(K::Triangle_3(a, b, c),point_query);
    double res = sqrt(sqd);
    return res;

}

inline K2::Point_3 iGameVertex_to_Point_K2(const MeshKernel::iGameVertex& v){
    return K2::Point_3(v.x(),v.y(),v.z());
}


inline  MeshKernel::iGameVertex Point_K2_to_iGameVertex(const K2::Point_3& v){
    return MeshKernel::iGameVertex(CGAL::to_double(v.x()),CGAL::to_double(v.y()),CGAL::to_double(v.z()));
}

inline  K::Point_3 Point_K2_to_Point_K(const K2::Point_3& v){
    return K::Point_3 (CGAL::to_double(v.x()),CGAL::to_double(v.y()),CGAL::to_double(v.z()));
}

inline  MeshKernel::iGameVertex Point_K_to_iGameVertex(const K::Point_3& v){
    return MeshKernel::iGameVertex((v.x()),(v.y()),(v.z()));
}


bool segment_in_line(K2::Segment_3 a,K2::Segment_3  b){
    CGAL::Epeck::FT d0 = CGAL::squared_distance(a.supporting_line(),b.vertex(0));
    CGAL::Epeck::FT d1 = CGAL::squared_distance(a.supporting_line(),b.vertex(1));
    if(d0 <= CGAL::Epeck::FT(0) &&
       d1 <= CGAL::Epeck::FT(0))
        return true;
    return false;
}

bool segment_coincide_triangle(K2::Segment_3 a,K2::Triangle_3 tri){
    for(int i=0;i<3;i++){
        K2::Segment_3 c(tri.vertex(i),tri.vertex((i+1)%3));
        if(segment_in_line(a,c))return true;
    }
    return false;
}

bool point_coincide_triangle(K2::Point_3 a,K2::Triangle_3 tri){
    for(int i=0;i<3;i++){
        K2::Segment_3 c(tri.vertex(i),tri.vertex((i+1)%3));
        CGAL::Epeck::FT d = CGAL::squared_distance(c.supporting_line(),a);
        if(d <= CGAL::Epeck::FT(0))return true;
    }
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





#endif //THICKEN2_GEOMETRIC_DATA_H

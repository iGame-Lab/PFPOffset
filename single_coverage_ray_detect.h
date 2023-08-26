//
// Created by rainbowwing on 2023/8/26.
//

#ifndef THICKEN2_SINGLE_COVERAGE_RAY_DETECT_H
#define THICKEN2_SINGLE_COVERAGE_RAY_DETECT_H

bool in_single_coverage_field_ray_detect(int field_id,K2::Triangle_3& face){
    K2::Vector_3 ray_vec = face.supporting_plane().orthogonal_vector();
    K2::Point_3 center = centroid(face);
    vector< K2::Triangle_3> field_face_list;
    for(int i=0;i< coverage_field_list[field_id].renumber_bound_face_global_id.size();i++){
        int global_face_id = coverage_field_list[field_id].renumber_bound_face_global_id[i];
        int tv0 = global_face_list[global_face_id].idx0;
        int tv1 = global_face_list[global_face_id].idx1;
        int tv2 = global_face_list[global_face_id].idx2;
        K2::Triangle_3 tri(global_vertex_list[tv0],global_vertex_list[tv1],global_vertex_list[tv2]);
        field_face_list.push_back(tri);
    }

    for(int i=0;i<field_face_list.size();i++){
        K2::Vector_3 ray_detect = field_face_list[i].supporting_plane().orthogonal_vector();
        if(ray_detect*ray_vec>=CGAL::Epeck::FT(0)){
            ray_detect = -ray_detect;
        }
        K2::Ray_3 ray(center,ray_detect);
        bool flag = true;
        vector<K2::Point_3> v;
        for(int j=0;j<field_face_list.size();j++){
            if(CGAL::squared_distance(field_face_list[j],center) > CGAL::Epeck::FT(0)){
                CGAL::cpp11::result_of<K2::Intersect_3(K2::Ray_3 , K2::Triangle_3)>::type
                        res_lt = intersection(ray,field_face_list[j]);
                if(res_lt) {
                    if (const K2::Point_3 *ss = boost::get<K2::Point_3>(&*res_lt)) {
                        v.push_back(*ss);
                    } else {
                        flag = false;
                    }
                }
            }
        }
        if(flag){
            sort(v.begin(),v.end(),[&](const K2::Point_3 &a,const K2::Point_3 &b){
                if(a.x() != b.x()){
                    return a.x() < b.x();
                }
                else if(a.y() != b.y()){
                    return a.y() < b.y();
                }
                return a.z() < b.z();
            });
            int vsize_pre = v.size();
            v.resize(unique(v.begin(),v.end())-v.begin());
            if(v.size() == vsize_pre){//没有重复交点
                if(v.size() % 2 == 1)return true;
                //return false;
            }
        }
    }
    return false;


}

#endif //THICKEN2_SINGLE_COVERAGE_RAY_DETECT_H

//
// Created by rainbowwing on 2023/9/10.
//

#ifndef THICKEN2_MERGE_INITIAL_H
#define THICKEN2_MERGE_INITIAL_H
//
void subdivide(vector<int>v,std::vector<K::Point_3>& kd_tree_points_new,int cnt){
    if(v.size() == 0)return;
    //float delta_x = -1;
    double move_dist_avg = 0;
    K::Vector_3 avg_vec(0,0,0);
    for(int i=0;i<v.size();i++){
        move_dist_avg += merge_limit[v[i]];
        avg_vec += kd_tree_points_new[v[i]] - K::Point_3(0,0,0);
    }
    avg_vec /= v.size();
    double max_x = kd_tree_points_new[v[0]].x();
    double min_x = kd_tree_points_new[v[0]].x();
    double max_y = kd_tree_points_new[v[0]].y();
    double min_y = kd_tree_points_new[v[0]].y();
    double max_z = kd_tree_points_new[v[0]].z();
    double min_z = kd_tree_points_new[v[0]].z();
    for(int i=0;i<v.size();i++){
        max_x = max(max_x, kd_tree_points_new[v[i]].x());
        min_x = min(min_x, kd_tree_points_new[v[i]].x());
        max_y = max(max_y, kd_tree_points_new[v[i]].y());
        min_y = min(min_y, kd_tree_points_new[v[i]].y());
        max_z = max(max_z, kd_tree_points_new[v[i]].z());
        min_z = min(min_z, kd_tree_points_new[v[i]].z());
    }
    double dist = sqrt((max_x-min_x)*(max_x-min_x) + (max_y-min_y)*(max_y-min_y) + (max_z-min_z)*(max_z-min_z));
    if(dist > move_dist_avg  && cnt < 4){
        vector<vector<int> >subdivide_next(8);
        for(int i=0;i<=1;i++){
            for(int j=0;j<=1;j++){
                for(int k=0;k<=1;k++){
                    double this_min_x = min_x + (max_x - min_x)*i;
                    double this_min_y = min_y + (max_y - min_y)*j;
                    double this_min_z = min_z + (max_z - min_z)*k;
                    double this_max_x = min_x + (max_x - min_x)*(i+1);
                    double this_max_y = min_y + (max_y - min_y)*(j+1);
                    double this_max_z = min_z + (max_z - min_z)*(k+1);
                    for(int id : v){
                        if(this_min_x<= kd_tree_points_new[id].x() && kd_tree_points_new[id].x() <= this_max_x &&
                                this_min_y<= kd_tree_points_new[id].y() && kd_tree_points_new[id].y() <= this_max_y &&
                                this_min_z<= kd_tree_points_new[id].z() && kd_tree_points_new[id].z() <= this_max_z){
                            subdivide_next[i*4+j*2+k].push_back(id);
                        }
                    }
                }
            }
        }
        for(int i=0;i<8;i++){
            subdivide(subdivide_next[i],kd_tree_points_new,cnt+1);
        }
    }
    else{
        for(int i=0;i<v.size();i++){
            kd_tree_points_new[v[i]] = K::Point_3(0,0,0) + avg_vec;
        }
        return ;
    }

}



void merge_initial(){
    std::vector<K::Point_3> kd_tree_points;
    std::vector<K::Point_3> kd_tree_points_new;
//    std::vector<K::Vector_3> kd_tree_vec;
    std::vector<vector<int> > kd_tree_list;
    std::vector<pair<int,int> > kd_tree_which_source;
    std::vector<double>kd_tree_points_move_limit;
    unordered_map<unsigned long long, int > mp;
    for(int i=0;i<mesh->VertexSize();i++){
        for(int j=0;j<field_move_vertices[i].size();j++){
            mp[unique_hash_value(field_move_vertices[i][j])] = kd_tree_which_source.size();
            kd_tree_points.push_back(field_move_vertices[i][j]);
            kd_tree_points_move_limit.push_back(merge_limit[i]);
            kd_tree_list.emplace_back();
            kd_tree_which_source.emplace_back(i,j);
            //cout <<"前置map"<< i <<" "<<j <<" 第"<< mp[unique_hash_value(field_move_vertices[i][j])] <<" "<<  unique_hash_value(field_move_vertices[i][j]) << endl;
        }
    }
    kd_tree_points_new.resize(kd_tree_points.size());
    DSU dsu;

    Kd_tree tree(kd_tree_points.begin(),kd_tree_points.end());
    for(int i=0;i<kd_tree_points.size();i++){
        std::vector<Point> result;
        Fuzzy_circle fs(kd_tree_points[i],kd_tree_points_move_limit[i]);
        tree.search(std::back_inserter(result), fs);
        for (const Point& p : result) {
            int source_id = mp[unique_hash_value(p)];
            pair<int,int>which = kd_tree_which_source[source_id];
            if(sqrt(CGAL::squared_distance(p,kd_tree_points[i])) < merge_limit[which.first]){
                auto a = unique_hash_value(kd_tree_points[i]);
                auto b = unique_hash_value(p);
                if(a!=b){
                    //cout <<"join: "<<a<<" "<<b<< endl;
                    dsu.join(a,b);
                }
            }
        }
    }
    FILE *file16_6 = fopen( ("../data/debug_16.6.obj"), "w");cout <<"open succ"<<endl;
    for(int i=0;i<kd_tree_points.size();i++){
        unsigned long long hash_value = unique_hash_value(kd_tree_points[i]);
        //cout << "differ equal: "<<(hash_value != dsu.find_root(hash_value)) << endl;
        int source_id = mp[dsu.find_root(hash_value)];
        pair<int,int>which = kd_tree_which_source[source_id];
        if(hash_value != dsu.find_root(hash_value)){
            fprintf(file16_6,"v %lf %lf %lf\n",mesh->fast_iGameVertex[which.first].x(),
                    mesh->fast_iGameVertex[which.first].y(),
                    mesh->fast_iGameVertex[which.first].z());
        }
        kd_tree_list[source_id].push_back(i);
        kd_tree_points_new[i] = kd_tree_points[i];
    }
    kd_tree_points_new.resize(kd_tree_points.size());
    for(int i=0;i<kd_tree_list.size();i++){
        subdivide(kd_tree_list[i],kd_tree_points_new,0);
    }

    for(int i=0;i<kd_tree_points.size();i++){
        // cout << i<<" "<< unique_hash_value(kd_tree_points[i]) << " "<< unique_hash_value(kd_tree_points[i]) <<" "<< unique_hash_value(kd_tree_points[i]) <<endl;
        // cout << i<<" "<<CGAL::hash_value(kd_tree_points[i]) <<" "<<CGAL::hash_value(kd_tree_points[i])<<" "<<CGAL::hash_value(kd_tree_points[i])<<endl;
        unsigned long long hash_value = unique_hash_value(kd_tree_points[i]);
        int source_id = mp[dsu.find_root(hash_value)];
        pair<int,int>which = kd_tree_which_source[i];
//        cout <<i<<"point"<< kd_tree_points[i].x() <<" "<< kd_tree_points[i].y()<<" "<< kd_tree_points[i].z()<<" "<<i<<"hashvalue"<<hash_value <<"root"<<dsu.find_root(hash_value) <<" "<<source_id <<" "<<which.first <<" "<< which.second << endl;
        field_move_vertices[which.first][which.second] =  kd_tree_points_new[source_id] ;
    }
}

#endif //THICKEN2_MERGE_INITIAL_H

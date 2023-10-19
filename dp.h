//
// Created by rainbowwing on 2023/8/25.
//

#ifndef THICKEN2_DP_H
#define THICKEN2_DP_H
#include <random>
double avg_edge_limit = 1e-10;
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
vector<MeshKernel::iGameVertex> run(MeshKernel::iGameVertexHandle vh,vector<MeshKernel::iGameFaceHandle> neighbor_face_list){
//    if(neighbor_face_list.size()>=21) {
//        cout << "neighbor_face_list.size()=" <<neighbor_face_list.size() <<">21 : this mesh please do remeshing before" << endl;
//        exit(0);
//    }
    vector<double>dp;
    vector<int>dp_source;
    vector<MeshKernel::iGameVertex>dp_osqp_answer;
    dp.resize(1<<neighbor_face_list.size());
    dp_source.resize(1<<neighbor_face_list.size());
    dp_osqp_answer.resize(1<<neighbor_face_list.size());
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
        double self_value = 1e100;
        for(int times=0;times<4;times++) { //避免osqp求解器的不稳定性,多尝试几次
            MeshKernel::iGameVertex v_new = do_quadratic_error_metric_check(vh, local_neighbor_face_list, succ,
                                                                            exceed_dist);
            if (succ) {
                if((v_new - mesh->fast_iGameVertex[vh]).norm() < self_value){
                    self_value = (v_new - mesh->fast_iGameVertex[vh]).norm();
                }
                dp[state] = (v_new - mesh->fast_iGameVertex[vh]).norm();
                dp_osqp_answer[state] = v_new;
                dp_source[state] = 0;
                self_value = dp[state];
            }
        }
        double minx = 1e100;
        if(neighbor_face_list.size() >=3) {
            vector<int> sub_state = std::move(get_sub_state(state, neighbor_face_list.size()));
            int from = -1;
            for (auto next: sub_state) {
                if (next == 0 || state - next == 0)continue;
                if (next > state - next)continue;
                double sub_ans = dfs(next) + dfs(state - next);
                if (sub_ans < minx) {
                    from = next;
                    minx = sub_ans;
                }
            }
            if(self_value > minx) {
                dp[state] = minx;
                dp_source[state] = from;
                return minx;
            }
        }
        return min(minx,self_value);
    };
    dfs((1<<neighbor_face_list.size())-1);
    queue<int>q;
    vector<MeshKernel::iGameVertex>ret;

    q.push((1<<neighbor_face_list.size())-1);
    while(!q.empty()){
        int now = q.front();
        q.pop();
        if(dp_source[now] == 0){
            //ret.push_back(now);
            ret.push_back(dp_osqp_answer[now]);
        }
        else{
            q.push(dp_source[now]);
            q.push(now - dp_source[now]);
        }
    }
    return ret;
}

std::mt19937 mt;


vector<MeshKernel::iGameVertex> solve_by_dp(MeshKernel::iGameVertexHandle vh,vector<MeshKernel::iGameFaceHandle> neighbor_face_list){
    vector<MeshKernel::iGameFaceHandle> neighbor_face_list_tmp = neighbor_face_list;
    //cout <<"determins " <<neighbor_face_list.size() <<":"<<avg_edge_limit / 1000<< endl;
    neighbor_face_list.clear();
    for(int i=0;i<neighbor_face_list_tmp.size();i++){
        //cout <<"************************************" << endl;
        MeshKernel::iGameFaceHandle f = neighbor_face_list[i];
        bool flag = true;
        vector<double>len_list;
        for(int j=0;j<3;j++) {
            len_list.push_back((mesh->fast_iGameVertex[mesh->fast_iGameFace[f].vh(j)]
                                - mesh->fast_iGameVertex[mesh->fast_iGameFace[f].vh((j+1)%3)]).norm());
        }
        sort(len_list.begin(),len_list.end());
        if(len_list[0] + len_list[1]> len_list[2] + avg_edge_limit/100){
            neighbor_face_list.push_back(f);
        }
//
//
//        if(flag) {
//
//
//
////            cout <<"v "<< mesh->fast_iGameVertex[mesh->fast_iGameFace[f].vh(0)].x()
////            <<" "<<mesh->fast_iGameVertex[mesh->fast_iGameFace[f].vh(0)].y() <<" "<<
////            mesh->fast_iGameVertex[mesh->fast_iGameFace[f].vh(0)].z() << endl;
////            cout <<"v "<< mesh->fast_iGameVertex[mesh->fast_iGameFace[f].vh(1)].x()
////                 <<" "<<mesh->fast_iGameVertex[mesh->fast_iGameFace[f].vh(1)].y() <<" "<<
////                 mesh->fast_iGameVertex[mesh->fast_iGameFace[f].vh(1)].z() << endl;
////            cout <<"v "<< mesh->fast_iGameVertex[mesh->fast_iGameFace[f].vh(2)].x()
////                 <<" "<<mesh->fast_iGameVertex[mesh->fast_iGameFace[f].vh(2)].y() <<" "<<
////                 mesh->fast_iGameVertex[mesh->fast_iGameFace[f].vh(2)].z() << endl;
////            cout <<"f "<<1 <<" "<<2<<" "<<3<< endl;
////            MeshKernel::iGameVertex normal
////                    = ((mesh->fast_iGameVertex[mesh->fast_iGameFace[f].vh(1)]
////                        - mesh->fast_iGameVertex[mesh->fast_iGameFace[f].vh(0)]) %
////                       (mesh->fast_iGameVertex[mesh->fast_iGameFace[f].vh(2)]
////                        - mesh->fast_iGameVertex[mesh->fast_iGameFace[f].vh(0)])).normalize();
////            cout <<"normal "<< normal.x() <<" "<< normal.y() <<" "<< normal.z() << endl;
////            K2::Triangle_3 tri(iGameVertex_to_Point_K2(mesh->fast_iGameVertex[mesh->fast_iGameFace[f].vh(0)]),
////                                iGameVertex_to_Point_K2(mesh->fast_iGameVertex[mesh->fast_iGameFace[f].vh(1)]),
////                                iGameVertex_to_Point_K2(mesh->fast_iGameVertex[mesh->fast_iGameFace[f].vh(2)])
////                               );
////            double xx = CGAL::to_double(tri.supporting_plane().orthogonal_vector().x());
////            double yy = CGAL::to_double(tri.supporting_plane().orthogonal_vector().y());
////            double zz = CGAL::to_double(tri.supporting_plane().orthogonal_vector().z());
////            cout << xx <<" "<< yy <<" "<< zz << endl;
//            neighbor_face_list.push_back(f);
//        }
    }
    //exit(0);
    //cout <<"determine "<< neighbor_face_list.size() << endl;


    //cout << vh << endl;
    if(neighbor_face_list.size()<=16) {
        //cout <<"?? now : " << neighbor_face_list.size() << endl;
        return run(vh,neighbor_face_list);
    }
    else{
        //cout<<"1-ring size: "<< neighbor_face_list.size() <<" meeting limit 16 do random subdivide"<< endl;

        int cnt = ((int)neighbor_face_list.size() - 1)/16 + 1;
        int each = neighbor_face_list.size() / cnt + 1;
        //if(each > 16)exit(0);
        //cout <<"each "<< each<< endl;
        double ans = 1e100;
        vector<MeshKernel::iGameVertex> ret;
        for(int times=0;times<10;times++){
            vector<MeshKernel::iGameVertex> now;
            std::shuffle(neighbor_face_list.begin(), neighbor_face_list.end(), mt);
            vector<MeshKernel::iGameFaceHandle> que;
            for(int i=0;i<neighbor_face_list.size();i++){
                //cout<<"times: " << times <<" i:"<< i << endl;
                que.push_back(neighbor_face_list[i]);
                if(que.size() == each){
                    vector<MeshKernel::iGameVertex> tmp = run(vh, que);
                    for(auto item : tmp) now.push_back(item);
                    que.clear();
                }
            }
            if(que.size()){
                //cout << "last run:"<<que.size() << endl;
                vector<MeshKernel::iGameVertex> tmp = run(vh, que);
                //cout << "last run end:"<<tmp.size() << endl;
                for(auto item : tmp) now.push_back(item);
                que.clear();
            }
            double now_ans = 0;
            for(int i=0;i<now.size();i++){
                now_ans += now[i].dist(mesh->fast_iGameVertex[vh]);
            }
            if(now_ans < ans){
                ans = now_ans;
                ret = now;
            }
        }
        return ret;
    }
}



#endif //THICKEN2_DP_H

//
// Created by rainbowwing on 2023/8/25.
//

#ifndef THICKEN2_DP_H
#define THICKEN2_DP_H
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

vector<MeshKernel::iGameVertex> solve_by_dp(MeshKernel::iGameVertexHandle vh,vector<MeshKernel::iGameFaceHandle> neighbor_face_list){
    if(neighbor_face_list.size()>=21) {
        cout << "neighbor_face_list.size()=" <<neighbor_face_list.size() <<">21 : this mesh please do remeshing before" << endl;
        exit(0);
    }
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


#endif //THICKEN2_DP_H

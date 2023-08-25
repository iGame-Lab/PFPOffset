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

vector<int> solve_by_dp(MeshKernel::iGameVertexHandle vh,vector<MeshKernel::iGameFaceHandle> neighbor_face_list){
    if(neighbor_face_list.size()>=25)exit(0);
    vector<double>dp;
    vector<int>dp_source;
    dp.resize(1<<neighbor_face_list.size());
    dp_source.resize(1<<neighbor_face_list.size());
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
        MeshKernel::iGameVertex v_new = do_quadratic_error_metric_check(vh,local_neighbor_face_list,succ,exceed_dist);
        if(succ){
            dp[state] = (v_new - mesh->fast_iGameVertex[vh]).norm();
            dp_source[state] = 0;
            return dp[state];
        }
        vector<int>sub_state = std::move(get_sub_state(state,neighbor_face_list.size()));
        double minx = 1e100;
        int from = -1;
        for(auto next: sub_state){
            if(next == 0 || state-next == 0)continue;
            if(next > state-next)continue;
            double sub_ans = dfs(next) + dfs(state-next);
            if(sub_ans < minx){
                from = next;
                minx = sub_ans;
            }
        }
        dp[state] = minx;
        dp_source[state] = from;
        return minx;
    };
    dfs((1<<neighbor_face_list.size())-1);
    queue<int>q;
    vector<int>ret;
    q.push((1<<neighbor_face_list.size())-1);
    while(!q.empty()){
        int now = q.front();
        q.pop();
        if(dp_source[now] == 0){
            ret.push_back(now);
        }
        else{
            q.push(dp_source[now]);
            q.push(now - dp_source[now]);
        }
    }
    return ret;
}


#endif //THICKEN2_DP_H

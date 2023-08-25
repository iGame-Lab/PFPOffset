//
// Created by teitoku on 2022-03-20.
//
#include <Eigen/Dense>
#include <map>
#include <thread>

#include "DSU.h"
#include <chrono>

unsigned long long DSU::find_root(unsigned long long id) {
    if(pre.find(id) == pre.end())
        pre[id]=id;


    unsigned long long r = id;
    while (pre[r] != r) {
        r = pre[r];
    }
    while (pre[id] != r) {
        unsigned long long y = pre[id];
        pre[id] = r;
        id = y;
    }
    return r;
}

void DSU::join(unsigned long long x, unsigned long long y) {
    if (find_root(x) != find_root(y)) {
        pre[find_root(x)] = find_root(y);
    }
    return;
}
std::unordered_map<unsigned long long,int> DSU::encode(int &encode_cnt){
    std::unordered_map<unsigned long long,int> ret;
    encode_cnt = 0;
    for(auto i : pre){
        if(find_root(i.first) == i.first){
            ret[i.first]=encode_cnt++;
        }
    }
    for(auto i : pre){
        if(find_root(i.first) != i.first){
            ret[i.first]=ret[find_root(i.first)];
        }
    }
    return ret;
}

void DSUMultiThread::inner_join(int x,int y) {
    if (find_root(x) != find_root(y)) {
        pre[find_root(x)] = find_root(y);
    }
    return;
}

int DSUMultiThread::find_root(int id) {

    int r = id;
    while (pre[r] != r) {
        r = pre[r];
    }
    while (pre[id] != r) {
        int y = pre[id];
        pre[id] = r;
        id = y;
    }
    return r;
}


DSUMultiThread::DSUMultiThread(int x){
    pre.resize(x+10);
    for(int i=0;i<x+10;i++){
        pre[i] = i;
    }
}
void DSUMultiThread::run(){
   join_thread = std::make_shared<std::thread>([&](){
       while(1){
           std::unique_lock<std::mutex> lck(mu2);
           cv.wait(lck, [&]() {
               return !q.empty();
           });
           std::unique_lock<std::mutex> lck2(mu1);
           std::pair<int,int> p = q.front();
           q.pop();
           if(p.first<0 && p.second < 0){
               break;
           }
           inner_join(p.first,p.second);
       }
       return ;
   });
}
void DSUMultiThread::join(int x, int y) {
    std::unique_lock<std::mutex> lck2(mu1);
    q.push({x,y});
    cv.notify_one();
}
void DSUMultiThread::stop() {
    std::unique_lock<std::mutex> lck2(mu1);
    q.push({-1,-1});
    cv.notify_one();
}
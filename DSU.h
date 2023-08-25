//
// Created by teitoku on 2022-03-20.
//

#ifndef IGAMEVIEW_DSU_H
#define IGAMEVIEW_DSU_H
#include <unordered_map>
#include <queue>
#include <vector>


class DSU {
    std::unordered_map<unsigned long long,unsigned long long>pre;
public:
    void join(unsigned long long x,unsigned long long y);
    unsigned long long find_root(unsigned long long x);
    std::unordered_map<unsigned long long,int>encode(int &encode_cnt);

};


class DSUMultiThread{
public:
    std::vector<int>pre;
    DSUMultiThread(){}
    DSUMultiThread(int x);
    int find_root(int x);
    void join(int x,int y);
    int find_root();
    void run();
    void stop();
    std::shared_ptr<std::thread> join_thread;
private:
    void inner_join(int x,int y);
    std::queue<std::pair<int,int> >q;
    std::mutex mu1,mu2;
    std::condition_variable cv;

};
#endif //IGAMEVIEW_DSU_H

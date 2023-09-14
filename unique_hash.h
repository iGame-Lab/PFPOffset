//
// Created by rainbowwing on 2023/9/8.
//

#ifndef THICKEN2_UNIQUE_HASH_H
#define THICKEN2_UNIQUE_HASH_H
struct Point_K_hash{
    size_t operator () (K::Point_3 x) const {
        return CGAL::hash_value(x);
    }
};
struct Point_K_equal{
    bool operator() (K::Point_3 a,  K::Point_3  b) const {
        return a==b;
    }
};

const int hash_factor = 5000011;
vector<unordered_map<K::Point_3,unsigned long long> >point_hash_vector(hash_factor);
unsigned long long unique_hash_value(K::Point_3 p){
    unsigned long long value = CGAL::hash_value(p)&((1LL<<32)-1);
//    value += *(unsigned long long *)(&p.x());
//    value += (*(unsigned long long *)(&p.y())) << (*(unsigned long long *)(&p.x())%10);
//    value += (*(unsigned long long *)(&p.z())) << (*(unsigned long long *)(&p.y())%10);
//    return value;

 //   cout<<"query" << value<<" "<<value%hash_factor<< endl;
    unordered_map<K::Point_3,unsigned long long>::iterator it = point_hash_vector[value%hash_factor].find(p);
    if(it == point_hash_vector[value%hash_factor].end()){
        //cout <<"meeting"<<endl;
        unsigned long long conflict_size = point_hash_vector[value%hash_factor].size()+1;
//        for(auto i:point_hash_vector[value%hash_factor]){
//            cout <<"config "<< i.first.x() <<" "<< i.first.y() <<" "<< i.first.z() << endl;
//        }
       // cout <<"meeting"<< p.x()<<" "<<p.y()<<" "<<p.z()<<" "<<value <<" "<<  conflict_size<<endl;
        unsigned long long new_value = (conflict_size<<33LL)|value;
        point_hash_vector[value%hash_factor][p]=new_value;
        return new_value;
    }
    else{
        return it->second;
    }
}



#endif //THICKEN2_UNIQUE_HASH_H

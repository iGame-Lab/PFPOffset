//
// Created by rainbowwing on 2023/8/25.
//

#ifndef THICKEN2_SORT_BY_POLAR_ORDER_H
#define THICKEN2_SORT_BY_POLAR_ORDER_H

void sort_by_polar_order(vector<K2::Point_3>& v,K2::Vector_3 orthogonal_direction){
    K2::Vector_3 center_v(0,0,0);
    for(auto j: v){
        center_v += j - K2::Point_3 (0,0,0);
    }
    center_v /= v.size();
    if(v.size() <=1)
        return ;
    K2::Point_3 center =  K2::Point_3 (0,0,0) + center_v;

    function<int(CGAL::Epeck::FT,CGAL::Epeck::FT)> quadrant = [](CGAL::Epeck::FT x,CGAL::Epeck::FT y){
        auto zero = CGAL::Epeck::FT(0);
        if(x>  zero && y > zero)return 1;
        else if(x<= zero && y >  zero)return 2;
        else if(x<= zero && y <= zero)return 3;
        return 4;
    };
    Plane_3 p;
    K2::Vector_3 x_axis = (v[0] - center) / sqrt(CGAL::to_double(CGAL::squared_distance(v[0],center)));
    K2::Vector_3 y_axis = CGAL::cross_product(orthogonal_direction,x_axis);
    y_axis = y_axis / sqrt(CGAL::to_double(y_axis.squared_length()));

    sort(v.begin(),v.end(),[&](K2::Point_3 a, K2::Point_3 b){

        CGAL::Epeck::FT x1 = (a - center) * x_axis;
        CGAL::Epeck::FT y1 = (a - center) * y_axis;
        CGAL::Epeck::FT x2 = (b - center) * x_axis;
        CGAL::Epeck::FT y2 = (b - center) * y_axis;
        int q1 = quadrant(x1,y1);
        int q2 = quadrant(x2,y2);
        if(q1!=q2)return q1<q2;
        else
            return x1*y2 - x2*y1 > 0;

    });
    return ;
}

#endif //THICKEN2_SORT_BY_POLAR_ORDER_H

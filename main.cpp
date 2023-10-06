#define CGAL_HAS_THREADS
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <queue>
#include <functional>
#include <set>
#include <memory>
#include <algorithm>
#include <bitset>
#include <thread>
#include "MeshKernel/Mesh.h"
#include "CGALPolygon.h"
#include "Eigen/Dense"
#include "Eigen/Core"
#include "Eigen/Sparse"
#include "OsqpEigen/OsqpEigen.h"
#include <vector>
#include <fstream>
#include <limits>
#include <unordered_set>
#include <sstream>
#include "DSU.h"
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Triangulation_2_projection_traits_3.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_2.h>
#include <CGAL/Triangulation_conformer_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Polygon_mesh_processing/stitch_borders.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/merge_border_vertices.h>
#include <CGAL/Polygon_mesh_processing/manifoldness.h>
#include <CGAL/Polygon_mesh_processing/intersection.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/version.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Fuzzy_sphere.h>
#include <vector>
#include <iostream>
#include <CGAL/Search_traits_3.h>
#include <atomic>
#include "obj_input.h"
#include "grid.h"
#include "global_face.h"
#include "geometric_data.h"
#include "cdt.h"
#include "unique_hash.h"
#include "CoverageField.h"
#include "osqp.h"
#include "dp.h"
#include "sort_by_polar_order.h"
#include "remeshing.h"
#include "single_coverage_ray_detect.h"
#include "flag_parser.h"
#include "merge_initial.h"



using namespace std;
int main(int argc, char* argv[]) {
    google::ParseCommandLineFlags(&argc, &argv, true);
    flag_parser();

    cout <<"CGAL_RELEASE_DATE:" << CGAL_RELEASE_DATE << endl;

    FILE *file11 = fopen( (input_filename + "_grid.obj").c_str(), "w");
    FILE *file6 = fopen( (input_filename + "_tmp.obj").c_str(), "w");

    FILE *file12 = fopen( (input_filename + "_12.obj").c_str(), "w");
    FILE *file14 = fopen( (input_filename + "_14.obj").c_str(), "w");

    default_move = 0.01;
    grid_len = 2.5;

  //3.59
    mesh = make_shared<MeshKernel::SurfaceMesh>(ReadObjFile(input_filename)); grid_len = 0.1;
    mesh->initBBox();
    mesh->build_fast();
    cout <<"mesh->build_fast() succ" << endl;
    double default_move_dist = 0.05;
    {
        double sum = 0;
        for(int i=0;i<mesh->FaceSize();i++){
            double minx = *set<double>{(mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(0)] - mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(1)]).norm(),
                                (mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(1)] - mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(2)]).norm(),
                                (mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(2)] - mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(0)]).norm()
                                }.begin();
            double maxx = *set<double>{(mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(0)] - mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(1)]).norm(),
                                       (mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(1)] - mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(2)]).norm(),
                                       (mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(2)] - mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(0)]).norm()
            }.rbegin();

            sum += sqrt(minx*maxx);

        }
        default_move_dist = sum/mesh->FaceSize()*1.5/4*10;
        double x_len = (mesh->BBoxMax - mesh->BBoxMin).x();
        double y_len = (mesh->BBoxMax - mesh->BBoxMin).y();
        double z_len = (mesh->BBoxMax - mesh->BBoxMin).z();
        double min_len = min(min(x_len,y_len),z_len);
        double rate = x_len/min_len*y_len/min_len*z_len/min_len;

        cout <<"??"<<(4*thread_num) <<" "<<rate << endl; // minlen rate 格了 4* thread_num 格子
        double need_div = max(cbrt((4*thread_num)/rate),1.0);
        grid_len = min_len / need_div;//(mesh->BBoxMax - mesh->BBoxMin).norm()/ thread_num * 2;
        cout << "grid_len "<< grid_len<<endl;

    }
    for (int i = 0; i < mesh->FaceSize(); i++) {
        mesh->fast_iGameFace[MeshKernel::iGameFaceHandle(i)].move_dist =  default_move_dist;
    }




    field_move_vertex.resize(mesh->VertexSize());
    field_move_vertices.resize(mesh->VertexSize());


    field_move_face.resize(mesh->FaceSize());
    field_move_K2_triangle.resize(mesh->FaceSize());

    cout <<"st do_quadratic_error_metric" << endl;


    merge_limit.resize(mesh->VertexSize());
    std::vector <std::shared_ptr<std::thread> > dp_thread_pool(thread_num);
    for(int i=0;i<thread_num;i++) {
        dp_thread_pool[i] = make_shared<std::thread>([&](int now_id) {
            for (int i = 0; i < mesh->VertexSize(); i++) {
                if (i % thread_num != now_id)continue;
                if (i % 20 == 0)
                    cout << "dp: " << i << "/"<<mesh->VertexSize()<< endl;
                vector<MeshKernel::iGameFaceHandle> neighbor_list;
                double min_move_dist = 1e30;
                for(auto j : mesh->FastNeighborFhOfVertex_[i]){
                    neighbor_list.push_back(j);
                }
                // cout <<i<<"//"<<mesh->VertexSize()<<" dp_result.size(): " <<"start" << endl;
                vector<MeshKernel::iGameVertex> dp_result = solve_by_dp(MeshKernel::iGameVertexHandle(i),neighbor_list);
                //cout <<i<<"//"<<mesh->VertexSize()<<" dp_result.size(): " <<dp_result.size() << endl;
                for(auto t: dp_result){
                    field_move_vertices[i].emplace_back(t.x(),t.y(),t.z());
                    min_move_dist = min(min_move_dist,mesh->fast_iGameVertex[i].dist(t));
                }

                merge_limit[i] = min_move_dist/10;
            }
        }, i);
    }
    for(int i=0;i<thread_num;i++)
        dp_thread_pool[i]->join();


    merge_initial();


    for(int i=0;i<mesh->FaceSize();i++){
        coverage_field_list.push_back(CoverageField(MeshKernel::iGameFaceHandle(i)));
    }
    int pre_cnt = 0;
    for(int i=0;i<mesh->FaceSize();i++){
//        if(coverage_field_list[i].bound_face_vertex_inexact.size() == 7 ) {
        for (int j = 0; j < coverage_field_list[i].bound_face_vertex_inexact.size(); j++) {
            fprintf(file14,"v %lf %lf %lf\n",CGAL::to_double(coverage_field_list[i].bound_face_vertex_inexact[j].x()),
                    CGAL::to_double(coverage_field_list[i].bound_face_vertex_inexact[j].y()),
                    CGAL::to_double(coverage_field_list[i].bound_face_vertex_inexact[j].z()));
        }
        for(int j=0;j<coverage_field_list[i].bound_face_id.size();j++){
            fprintf(file14,"f %d %d %d\n",coverage_field_list[i].bound_face_id[j][0]+1+pre_cnt,
                    coverage_field_list[i].bound_face_id[j][1]+1+pre_cnt,
                    coverage_field_list[i].bound_face_id[j][2]+1+pre_cnt);
        }
        pre_cnt += coverage_field_list[i].bound_face_vertex_inexact.size();
//            break;
//        }
    }
    //exit(0);



    std::vector <std::shared_ptr<std::thread> > one_ring_select_thread_pool(thread_num);
    for(int i=0;i<thread_num;i++) {
        one_ring_select_thread_pool[i] = make_shared<std::thread>([&](int now_id) {
            for(int i=0;i<mesh->FaceSize();i++) {
                if (i % thread_num != now_id)continue;
                if(i%500==0)
                    cout << "one_ring_select_thread_pool "<< i << endl;

                set<int>neighbor_field;
                for(int j=0;j<3;j++) {
                    for (auto neighbor_id: mesh->FastNeighborFhOfEdge_[mesh->fast_iGameFace[i].eh(j)]) {
                        if(neighbor_id == i)continue;
                        neighbor_field.insert(neighbor_id);
                    }
                }

                vector<K2::Triangle_3>neighbor_face;
                for(auto neighbor_id: neighbor_field){
                    for (int k = 0; k < coverage_field_list[neighbor_id].bound_face_id.size(); k++) {
                        K2::Triangle_3 tri_this(coverage_field_list[neighbor_id].bound_face_vertex_exact[coverage_field_list[neighbor_id].bound_face_id[k][0]],
                                                coverage_field_list[neighbor_id].bound_face_vertex_exact[coverage_field_list[neighbor_id].bound_face_id[k][1]],
                                                coverage_field_list[neighbor_id].bound_face_vertex_exact[coverage_field_list[neighbor_id].bound_face_id[k][2]]
                        );

                        neighbor_face.push_back(tri_this);
                    }
                }

                for(int j=0;j<coverage_field_list[i].bound_face_id.size();j++) {
                    if ((coverage_field_list[i].bound_face_id[j][0] >= 3 ||
                            coverage_field_list[i].bound_face_id[j][1] >= 3 ||
                            coverage_field_list[i].bound_face_id[j][2] >= 3)){
                        bool flag = false;

                        K2::Triangle_3 tri_this(coverage_field_list[i].bound_face_vertex_exact[coverage_field_list[i].bound_face_id[j][0]],
                                                coverage_field_list[i].bound_face_vertex_exact[coverage_field_list[i].bound_face_id[j][1]],
                                                coverage_field_list[i].bound_face_vertex_exact[coverage_field_list[i].bound_face_id[j][2]]);
                        K2::Segment_3 e0(tri_this.vertex(0), tri_this.vertex(1));
                        K2::Segment_3 e1(tri_this.vertex(1), tri_this.vertex(2));
                        K2::Segment_3 e2(tri_this.vertex(2), tri_this.vertex(0));

                        vector<K2::Segment_3> vs_tmp;
                        for (auto other: neighbor_face) { //Point_3, or Segment_3, or Triangle_3, or std::vector < Point_3 >
                            CGAL::cpp11::result_of<K2::Intersect_3(K2::Triangle_3, K2::Triangle_3)>::type
                                    res_tt = intersection(tri_this, other);
                            if (res_tt) {
                                if (const K2::Segment_3 *s = boost::get<K2::Segment_3>(&*res_tt)) {
                                    vs_tmp.push_back(*s);
                                } else if (const K2::Triangle_3 *t = boost::get<K2::Triangle_3>(&*res_tt)) {
                                    vs_tmp.emplace_back(t->vertex(0), t->vertex(1));
                                    vs_tmp.emplace_back(t->vertex(1), t->vertex(2));
                                    vs_tmp.emplace_back(t->vertex(2), t->vertex(0));
                                } else if (std::vector<K2::Point_3> *vs = boost::get<std::vector<K2::Point_3 >>(&*res_tt)) {
                                    sort_by_polar_order(*vs, tri_this.supporting_plane().orthogonal_vector());
                                    for (int k = 0; k < vs->size(); k++) {
                                        vs_tmp.emplace_back(vs->operator[](k), vs->operator[]((k + 1) % vs->size()));
                                        // cerr << "run iiiiiiiiiiiiiiit" << endl;
                                    }
                                }
                            }
                        }
                        vector<K2::Segment_3> vs;
                        for (auto se: vs_tmp) {
                            if (!segment_in_line(se, e0) && !segment_in_line(se, e1) && !segment_in_line(se, e2)) {
                                vs.push_back(se);
                            }
                        }
                        vector<vector<K2::Point_3> > res = CGAL_CDT_NEW({coverage_field_list[i].bound_face_vertex_exact[coverage_field_list[i].bound_face_id[j][0]],
                                                                     coverage_field_list[i].bound_face_vertex_exact[coverage_field_list[i].bound_face_id[j][1]],
                                                                     coverage_field_list[i].bound_face_vertex_exact[coverage_field_list[i].bound_face_id[j][2]]}, vs, tri_this);
                        //cout << res.size() <<endl;
                        for (auto each_tri: res) {
                            bool patch_flag = false;
                            K2::Point_3 center = CGAL::centroid(K2::Triangle_3(each_tri[0], each_tri[1], each_tri[2]));
                            for (auto j: neighbor_field) {
                                if (tri_this.supporting_plane().oriented_side(coverage_field_list[i].center) !=
                                    tri_this.supporting_plane().oriented_side(coverage_field_list[j].center) &&
                                    tri_this.supporting_plane().oriented_side(coverage_field_list[i].center)+
                                    tri_this.supporting_plane().oriented_side(coverage_field_list[j].center) ==0 &&
                                        coverage_field_list[j].in_or_on_field(center)) {
                                    patch_flag = true;
                                    break;
                                }
                            }
                            if(!patch_flag){
                                flag = true;
                                break;
                            }
                        }
                        coverage_field_list[i].bound_face_useful[j]=flag;
                        //cout <<i<<" "<<j <<" "<< (flag?"yes":"no") << endl;
                    }
                    else{
                        coverage_field_list[i].bound_face_useful[j]=false;
                    }

                }
            }
        },i);
    }

    for(int i=0;i<thread_num;i++)
        one_ring_select_thread_pool[i]->join();

//    int f12id = 1;
//    for(int i=0;i<mesh->FaceSize();i++) {
//        for (int j = 0; j < coverage_field_list[i].bound_face_id.size(); j++) {
//            K::Triangle_3 tri_this(
//                    coverage_field_list[i].bound_face_vertex_inexact[coverage_field_list[i].bound_face_id[j][0]],
//                    coverage_field_list[i].bound_face_vertex_inexact[coverage_field_list[i].bound_face_id[j][1]],
//                    coverage_field_list[i].bound_face_vertex_inexact[coverage_field_list[i].bound_face_id[j][2]]);
//            cout<<i<<" "<<j<<" "<<coverage_field_list[i].bound_face_useful[j] <<endl;
//            if (coverage_field_list[i].bound_face_useful[j] > 0) {
//                fprintf(file12, "v %lf %lf %lf\n", tri_this.vertex(0).x(), tri_this.vertex(0).y(),
//                        tri_this.vertex(0).z());
//                fprintf(file12, "v %lf %lf %lf\n", tri_this.vertex(1).x(), tri_this.vertex(1).y(),
//                        tri_this.vertex(1).z());
//                fprintf(file12, "v %lf %lf %lf\n", tri_this.vertex(2).x(), tri_this.vertex(2).y(),
//                        tri_this.vertex(2).z());
//                fprintf(file12, "f %d %d %d\n", f12id, f12id + 1, f12id + 2);
//                f12id+=3;
//            }
//        }
//    }


    for(auto  each_container_face : container_grid_face){
        each_grid_face_list.push_back({(size_t)each_container_face[0],(size_t)each_container_face[1],(size_t)each_container_face[2]});
        each_grid_face_list.push_back({(size_t)each_container_face[2],(size_t)each_container_face[3],(size_t)each_container_face[0]});
    }

    cout <<"build end "<< endl;

    cgal_polygon = make_shared<CGALPolygon>(mesh);

    std::list < K2::Triangle_3> triangles;

    mesh->initBBox();

    double max_move = 0;
    double min_move = 1e10;
    for (auto i: mesh->allfaces()) {
        max_move = max(max_move, abs(i.second.move_dist));
        min_move = min(min_move, abs(i.second.move_dist));
    }

    stx = mesh->BBoxMin.x() - max_move;
    sty = mesh->BBoxMin.y() - max_move;
    stz = mesh->BBoxMin.z() - max_move;

    printf("GL %lf\n", grid_len);

    unordered_map <grid, GridVertex, grid_hash, grid_equal> frame_grid_mp;

    vector <vector<int>> bfs_dir = {{0,  0,  1},
                                    {0,  1,  0},
                                    {1,  0,  0},
                                    {0,  0,  -1},
                                    {0,  -1, 0},
                                    {-1, 0,  0}};
    std::function < vector<grid>(grid) > get_neighbor = [&](grid g) {
        vector <grid> ret;
        for (auto i: bfs_dir) {
            int xx, yy, zz;
            xx = g.x + i[0];
            yy = g.y + i[1];
            zz = g.z + i[2];
            if (xx >= 0 && yy >= 0 && zz >= 0)
                ret.push_back({xx, yy, zz});
        }
        return ret;
    };
    int fsize = mesh->FaceSize();

    cout <<"bfs start \n" << endl;
    std::mutex bfs_mutex;
    std::vector <std::shared_ptr<std::thread> > bfs_thread_pool(thread_num);
    for(int i=0;i<thread_num;i++) {
        bfs_thread_pool[i] = make_shared<std::thread>([&](int now_id) {
            for (int face_id = 0; face_id < fsize; face_id++) {
                if(face_id % thread_num !=  now_id)continue;
                if (face_id % 1000 == 0)
                    printf("%d/%d\n", face_id, fsize);
                auto fh = make_pair(MeshKernel::iGameFaceHandle(face_id),
                                    mesh->faces(MeshKernel::iGameFaceHandle(face_id)));
                MeshKernel::iGameVertex center = (mesh->fast_iGameVertex[fh.second.vh(0)] +
                                                  mesh->fast_iGameVertex[fh.second.vh(1)] +
                                                  mesh->fast_iGameVertex[fh.second.vh(2)]) / 3;

                double move_limit = 0;
                for(int j=0;j<3;j++){
                    for(int k=0;k<field_move_vertices[j].size();k++){
                        move_limit = max(move_limit,(Point_K_to_iGameVertex(field_move_vertices[fh.second.vh(j)][k]) - mesh->fast_iGameVertex[fh.second.vh(j)]).norm());
                    }
                }

                grid now = vertex_to_grid(center);
                vector <MeshKernel::iGameFaceHandle> face_and_neighbor;
                face_and_neighbor.push_back(fh.first);
                for (auto i: mesh->NeighborFh(fh.first)) {
                    face_and_neighbor.push_back(i);
                }

                queue <grid> q;
                unordered_set <grid,grid_hash,grid_equal> is_visit;
                vector <grid> center_neighbor = get_neighbor(now);
                for (auto bfs_start_node: center_neighbor) {
                    is_visit.insert(bfs_start_node);
                    q.push(bfs_start_node);
                }
                while (!q.empty()) {
                    now = q.front();
                    q.pop();
                    //cout << now.x <<" "<< now.y <<""
                    std::unique_lock<std::mutex>lock1(bfs_mutex,std::defer_lock);
                    lock1.lock();
                    auto iter = frame_grid_mp.find(now);
                    if (iter == frame_grid_mp.end()) {
                        iter = frame_grid_mp.insert(make_pair(now, GridVertex())).first;
                    }
                    iter->second.field_list.push_back(fh.first);
                    lock1.unlock();
                    vector <grid> neighbor = get_neighbor(now);
                    for (auto j: neighbor) {
                        if (!is_visit.count(j)) {
                            double dist = cgal_vertex_triangle_dist(fh.second, getGridVertex(j, 0), mesh);

                            if (dist <  (grid_len*1.74/2+move_limit) ) { //TODO : zheli youhua cheng pianyi juli de shiji jisuan
                                q.push(j);
                                is_visit.insert(j);
                            }
                        }
                    }
                }
            }
        },i);
    }

    for(int i=0;i<thread_num;i++)
        bfs_thread_pool[i]->join();

    cout <<"bfs end \n" << endl;

    auto frame_grid_mp_bk =  frame_grid_mp;
    for(auto i : frame_grid_mp_bk) {
        for(auto j : container_grid_dir){
            int nx = i.first.x+j[0];
            int ny = i.first.y+j[1];
            int nz = i.first.z+j[2];
            //if(nx>=0 && ny>=0 && nz>=0){
                for(auto z : i.second.field_list)
                    frame_grid_mp[{nx,ny,nz}].field_list.push_back(z);
           // }
        }
    }


    std::vector <std::shared_ptr<std::thread> > each_grid_thread_pool(thread_num);
    for(int i=0;i<thread_num;i++) { //todo 存在偶发性多线程异常
        each_grid_thread_pool[i] = make_shared<std::thread>([&](int now_id) {
            int each_grid_cnt =-1;
            for (auto each_grid = frame_grid_mp.begin(); each_grid != frame_grid_mp.end(); each_grid++) { // todo 逻辑修改
                each_grid_cnt++;
                if(each_grid_cnt % thread_num != now_id)continue;
                if(each_grid_cnt % (thread_num) == now_id)
                    printf("thread num :%d each_grid_cnt %d/%d\n",now_id,each_grid_cnt,(int)frame_grid_mp.size());

                MeshKernel::iGameVertex grid_vertex = getGridVertex(each_grid->first, 0);

                sort(each_grid->second.field_list.begin(), each_grid->second.field_list.end());
                each_grid->second.field_list.resize(
                        unique(each_grid->second.field_list.begin(), each_grid->second.field_list.end()) -
                        each_grid->second.field_list.begin());


                K2::Point_3 small = iGameVertex_to_Point_K2(getGridVertex(each_grid->first, 0));
                K2::Point_3 big = iGameVertex_to_Point_K2(getGridVertex(each_grid->first, 7));

                std::vector<K2::Point_3> ps;
                for(int i=0;i<8;i++){
                    ps.push_back(getGridK2Vertex(small,big,i));
                }
                std::list<K2::Triangle_3 >  frame_faces_list;

                for(auto i : each_grid_face_list){
                    frame_faces_list.emplace_back(ps[i[0]],ps[i[1]],ps[i[2]]);
                }
                Tree frame_aabb_tree(frame_faces_list.begin(),frame_faces_list.end());
                CGAL::Polyhedron_3<K2> frame_poly;
                PMP::polygon_soup_to_polygon_mesh(ps, each_grid_face_list, frame_poly, CGAL::parameters::all_default());

                for(int j = 0; j <  each_grid->second.field_list.size() ; j++){
                    int field_id = each_grid->second.field_list[j];
                    for(int k=0;k<coverage_field_list[field_id].bound_face_id.size();k++){

                        int v0_id = coverage_field_list[field_id].bound_face_id[k][0];
                        int v1_id = coverage_field_list[field_id].bound_face_id[k][1];
                        int v2_id = coverage_field_list[field_id].bound_face_id[k][2];

                        K2::Triangle_3 this_tri(coverage_field_list[field_id].bound_face_vertex_exact[v0_id],
                                                coverage_field_list[field_id].bound_face_vertex_exact[v1_id],
                                                coverage_field_list[field_id].bound_face_vertex_exact[v2_id]
                                                );
                        if(!coverage_field_list[field_id].bound_face_useful[k])continue;
                        if(check_triangle_through_grid(small,big,frame_poly,this_tri)){

                            each_grid->second.field_face_though_list[field_id].push_back(k);
                            each_grid->second.face_hash_id_map[this_tri.id()] = {field_id,k};
                            each_grid->second.build_aabb_tree_triangle_list.push_back(this_tri);

                        }

                    }
                }
                Tree this_grid_aabb_tree( each_grid->second.build_aabb_tree_triangle_list.begin(),
                                          each_grid->second.build_aabb_tree_triangle_list.end());

                //X: 这里写查询 因为aabb树指针的bug，所以反着做，用每个格子去查询每个面

                for(auto item : each_grid->second.field_face_though_list){
                    int field_id = item.first;
                   for(int bound_face_id : item.second){
                       int v0_id = coverage_field_list[field_id].bound_face_id[bound_face_id][0];
                       int v1_id = coverage_field_list[field_id].bound_face_id[bound_face_id][1];
                       int v2_id = coverage_field_list[field_id].bound_face_id[bound_face_id][2];
                        K2::Triangle_3 this_tri(coverage_field_list[field_id].bound_face_vertex_exact[v0_id],
                                                coverage_field_list[field_id].bound_face_vertex_exact[v1_id],
                                                coverage_field_list[field_id].bound_face_vertex_exact[v2_id]
                        );

                       std::list< Tree::Intersection_and_primitive_id<K2::Triangle_3>::Type> intersections;

                       this_grid_aabb_tree.all_intersections(this_tri,std::back_inserter(intersections));

                       unique_lock<mutex> lock(mutex);
                       for(auto item : intersections) { //todo 这里改成批量插入
                           if(const K2::Segment_3 * s = boost::get<K2::Segment_3>(&(item.first))){

                               map<size_t,pair<int,int> >::iterator iter = each_grid->second.face_hash_id_map.find(item.second->id());
                               if(iter->second.first != field_id ) {
                                   coverage_field_list[field_id].bound_face_cutting_segment[bound_face_id].push_back(
                                           *s);
                               }
                           }
                           else if(const K2::Point_3 * p = boost::get<K2::Point_3>(&(item.first))){
                               map<size_t,pair<int,int> >::iterator iter = each_grid->second.face_hash_id_map.find(item.second->id());
                               if(iter->second.first != field_id ) {
                                   coverage_field_list[field_id].bound_face_cutting_point[bound_face_id].push_back(*p);
                               }
                           }
                           else if(const std::vector<K2::Point_3> * v = boost::get<std::vector<K2::Point_3> >(&(item.first)) ){
                               map<size_t,pair<int,int> >::iterator iter = each_grid->second.face_hash_id_map.find(item.second->id());
                               if(iter->second.first != field_id ){
                                   for(int j=0;j<v->size();j++){
                                       coverage_field_list[field_id].bound_face_cutting_point[bound_face_id].push_back(v->at(j));
                                   }
                               }
                           }
                       }
                       coverage_field_list[field_id].bound_face_cross_field_list[bound_face_id].push_back(each_grid->first);
                   }
                }


              }
            cout << "each_grid_cnt end thread num:"<< now_id<<endl;
        }, i);
    }

    for(int i=0;i<thread_num;i++)
        each_grid_thread_pool[i]->join();



    // 这里修正每个面

    // 应该先做点的序号合并


    for (auto each_grid= frame_grid_mp.begin(); each_grid != frame_grid_mp.end(); each_grid++) {
        //if(!(each_grid->first.x == 26 && each_grid->first.y == 28 && each_grid->first.z == 9  ))continue;
        //if(!(each_grid->first.x == 19 && each_grid->first.y == 19 && each_grid->first.z == 1 ))continue;
        //if(!(each_grid->first.x == 2 && each_grid->first.y == 3 && each_grid->first.z == 0 ))continue;
        auto small  = getGridVertex(each_grid->first,0);
        auto big  = getGridVertex(each_grid->first,7);
        static int f3_id = 1;
        for (int ii = 0; ii < 7; ii++) {
            for (int jj = 0; jj < DirectedGridEdge[ii].size(); jj++) {
                int from = ii;
                int to = DirectedGridEdge[ii][jj];
                MeshKernel::iGameVertex fv = getGridiGameVertex(small, big, from);
                MeshKernel::iGameVertex tv = getGridiGameVertex(small, big, to);
                fprintf(file11, "v %lf %lf %lf\n", fv.x(), fv.y(), fv.z());
                fprintf(file11, "v %lf %lf %lf\n", tv.x(), tv.y(), tv.z());
                fprintf(file11, "l %d %d\n", f3_id, f3_id + 1);
                f3_id += 2;
            }
        }
    }

//    {
//        FILE *file20 = fopen( (input_filename + "_20.obj").c_str(), "w");
//        FILE *file21 = fopen( (input_filename + "_21.obj").c_str(), "w");
//        FILE *file22 = fopen( (input_filename + "_22.obj").c_str(), "w");
//        for(int i=0;i< coverage_field_list[292].bound_face_cutting_point.size();i++){ // 每一个面
//            for(int j=0;j<coverage_field_list[292].bound_face_cutting_point[i].size();j++){
//                fprintf(file20,"v %lf %lf %lf\n",CGAL::to_double(coverage_field_list[292].bound_face_cutting_point[i][j].x()),
//                        CGAL::to_double(coverage_field_list[292].bound_face_cutting_point[i][j].y()),
//                        CGAL::to_double(coverage_field_list[292].bound_face_cutting_point[i][j].z())
//                        );
//            }
//        }
//        int f21id = 1;
//        for(int i=0;i< coverage_field_list[292].bound_face_cutting_segment.size();i++){ // 每一个面
//            for(int j=0;j<coverage_field_list[292].bound_face_cutting_segment[i].size();j++){
//                fprintf(file21,"v %lf %lf %lf\n",CGAL::to_double(coverage_field_list[292].bound_face_cutting_segment[i][j].vertex(0).x()),
//                        CGAL::to_double(coverage_field_list[292].bound_face_cutting_segment[i][j].vertex(0).y()),
//                        CGAL::to_double(coverage_field_list[292].bound_face_cutting_segment[i][j].vertex(0).z())
//                );
//                fprintf(file21,"v %lf %lf %lf\n",CGAL::to_double(coverage_field_list[292].bound_face_cutting_segment[i][j].vertex(1).x()),
//                        CGAL::to_double(coverage_field_list[292].bound_face_cutting_segment[i][j].vertex(1).y()),
//                        CGAL::to_double(coverage_field_list[292].bound_face_cutting_segment[i][j].vertex(1).z())
//                );
//                fprintf(file21,"l %d %d\n",f21id,f21id+1);
//                f21id += 2;
//            }
//        }
//        int f22id = 1;
//        for(int i=0;i< coverage_field_list[292].bound_face_id.size();i++){ // 每一个面
//            K::Point_3 v0 = coverage_field_list[292].bound_face_vertex_inexact[coverage_field_list[292].bound_face_id[i][0]];
//            K::Point_3 v1 = coverage_field_list[292].bound_face_vertex_inexact[coverage_field_list[292].bound_face_id[i][1]];
//            K::Point_3 v2 = coverage_field_list[292].bound_face_vertex_inexact[coverage_field_list[292].bound_face_id[i][2]];
//            fprintf(file22,"v %lf %lf %lf\n",v0.x(),v0.y(),v0.z());
//            fprintf(file22,"v %lf %lf %lf\n",v1.x(),v1.y(),v1.z());
//            fprintf(file22,"v %lf %lf %lf\n",v2.x(),v2.y(),v2.z());
//            fprintf(file22,"f %d %d %d\n",f22id,f22id+1,f22id+2);
//            f22id += 3;
//        }
//    }
//
//    exit(0);
//    for (int field_id = 0; field_id < fsize; field_id++) {
//
//
//        coverage_field_list[field_id].field_id = 1;
//        coverage_field_list[field_id].do_cdt();
//        cout <<"field_id:" << field_id <<" "<<coverage_field_list[field_id].bound_face_id.size() << endl;
//        int tt0 = 0;
//        int tt1 = 0;
//        int tt2 = 0;
//
//        for(int i=0;i< coverage_field_list[field_id].bound_face_cutting_point.size();i++){
//            tt0 += coverage_field_list[field_id].bound_face_cutting_point[i].size();
//        }
//        for(int i=0;i< coverage_field_list[field_id].bound_face_cutting_segment.size();i++){
//            tt1 += coverage_field_list[field_id].bound_face_cutting_segment[i].size();
//        }
//        for(int i=0;i< coverage_field_list[field_id].bound_face_cross_field_list.size();i++){
//            tt2 += coverage_field_list[field_id].bound_face_cross_field_list[i].size();
//        }
//        cout <<  tt0 <<" "<< tt1 <<" "<<tt2<< endl;
//
//        coverage_field_list[field_id].renumber();
//
//    }
//
//    exit(0);

    std::vector <std::shared_ptr<std::thread> > field_vertex_numbering_thread_pool(thread_num);
    std::atomic<int>global_vertex_id_sum;
    for(int i=0;i<thread_num;i++) {
        field_vertex_numbering_thread_pool[i] = make_shared<std::thread>([&](int now_id) {
            for (int field_id = 0; field_id < fsize; field_id++) {
                if (field_id % thread_num != now_id)continue;
                if(field_id %(10*thread_num) ==0)
                    cout << "start do cdt "<< field_id <<"/"<<fsize << endl;
                coverage_field_list[field_id].field_id = field_id;
                coverage_field_list[field_id].do_cdt();
                coverage_field_list[field_id].renumber();
                global_vertex_id_sum+= coverage_field_list[field_id].renumber_bound_face_vertex.size();
            }
        },i);
    }
    for(int i=0;i<thread_num;i++)
        field_vertex_numbering_thread_pool[i]->join();

    //exit(0);
    global_vertex_list.resize(global_vertex_id_sum);

    int global_cnt = 0;
    for (int field_id = 0; field_id < fsize; field_id++) {
        for (int i = 0; i < coverage_field_list[field_id].renumber_bound_face_vertex.size(); i++) {

            coverage_field_list[field_id].renumber_bound_face_vertex_global_id[i] = global_cnt + i;
            global_vertex_list[global_cnt + i] =  coverage_field_list[field_id].renumber_bound_face_vertex[i];
        }
        global_cnt += coverage_field_list[field_id].renumber_bound_face_vertex.size();
    }
//    cout <<"start kd tree join"<<endl;
//    std::vector<K::Point_3> kd_tree_points;
//    map<unsigned long long,int> global_kd_tree_mp;
//    for(int i=0;i<global_vertex_list.size();i++){
//        K::Point_3 p = PointK2_Point(global_vertex_list[i]);
//        unsigned long long id = unique_hash_value(p);
//        global_kd_tree_mp[id] = i;
//        kd_tree_points.push_back(p);
//    }
//    DSUMultiThread dsu_multi_thread((int)kd_tree_points.size());
//    std::vector <std::shared_ptr<std::thread> > global_dsu_thread_pool(thread_num);
//    dsu_multi_thread.run();
//    for(int i=0;i<thread_num;i++) {
//        global_dsu_thread_pool[i] = make_shared<std::thread>([&](int now_id) {
//            Kd_tree global_kd_tree(kd_tree_points.begin(),kd_tree_points.end());
//            for(int i=0;i<kd_tree_points.size();i++){
//                if (i % thread_num != now_id)continue;
//                std::vector<K::Point_3> result;
//                Fuzzy_circle fs(kd_tree_points[i], myeps);
//                global_kd_tree.search(std::back_inserter(result), fs);
//                for (const Point& p : result) {
//                    dsu_multi_thread.join(i,global_kd_tree_mp[unique_hash_value(p)]);
//                    //dsu.join(CGAL::hash_value(kd_tree_points[i]),hash_value(p));
//                }
//            }
//        },i);
//    }
//    for(int i=0;i<thread_num;i++)
//        global_dsu_thread_pool[i]->join();
//
//    cout <<"wait dsu_multi_thread stop"<< endl;
//    dsu_multi_thread.stop();
//    cout <<"wait dsu_multi_thread.join_thread->join()"<< endl;
//    dsu_multi_thread.join_thread->join();
//
//   // exit(0);
//    cout <<"global_vertex_id_sum:" << global_vertex_id_sum <<endl;
//    for (int field_id = 0; field_id < fsize; field_id++) {
//        for (int i = 0; i < coverage_field_list[field_id].renumber_bound_face_vertex_global_id.size(); i++) {
//
//            //coverage_field_list[field_id].renumber_bound_face_vertex_global_id[i] = dsu_multi_thread.find_root(coverage_field_list[field_id].renumber_bound_face_vertex_global_id[i]);
//        }
//    }
    int global_face_cnt = 0;
    for (int field_id = 0; field_id < fsize; field_id++) {
        for (int i = 0; i < coverage_field_list[field_id].renumber_bound_face_id.size(); i++) {
            int tv0 =  coverage_field_list[field_id].renumber_bound_face_vertex_global_id[coverage_field_list[field_id].renumber_bound_face_id[i][0]];
            int tv1 =  coverage_field_list[field_id].renumber_bound_face_vertex_global_id[coverage_field_list[field_id].renumber_bound_face_id[i][1]];
            int tv2 =  coverage_field_list[field_id].renumber_bound_face_vertex_global_id[coverage_field_list[field_id].renumber_bound_face_id[i][2]];
            //if(set<int>{tv0,tv1,tv2}.size() < 3)continue;
            global_face_cnt++;
        }
    }
    //exit(0);
    global_face_list.resize(global_face_cnt);
    global_face_cnt = 0;
    for (int field_id = 0; field_id < fsize; field_id++) {
        for (int i = 0; i < coverage_field_list[field_id].renumber_bound_face_id.size(); i++) {
            int tv0 =  coverage_field_list[field_id].renumber_bound_face_vertex_global_id[coverage_field_list[field_id].renumber_bound_face_id[i][0]];
            int tv1 =  coverage_field_list[field_id].renumber_bound_face_vertex_global_id[coverage_field_list[field_id].renumber_bound_face_id[i][1]];
            int tv2 =  coverage_field_list[field_id].renumber_bound_face_vertex_global_id[coverage_field_list[field_id].renumber_bound_face_id[i][2]];

            int useful = coverage_field_list[field_id].renumber_bound_face_useful[i];
            if(set<int>{tv0,tv1,tv2}.size() < 3)continue;
            for(int j =0 ; j < coverage_field_list[field_id].renumber_bound_face_cross_field_list[i].size();j++){
                grid cross_grid = coverage_field_list[field_id].renumber_bound_face_cross_field_list[i][j];
                frame_grid_mp[cross_grid].global_face_list.push_back(global_face_cnt);
            }
            coverage_field_list[field_id].renumber_bound_face_global_id[i] = global_face_cnt;
            global_face_list[global_face_cnt].field_id = field_id;
            global_face_list[global_face_cnt].idx0 = tv0;
            global_face_list[global_face_cnt].idx1 = tv1;
            global_face_list[global_face_cnt].idx2 = tv2;
            global_face_list[global_face_cnt].useful = useful;
            global_face_list[global_face_cnt].center = CGAL::centroid(K2::Triangle_3(global_vertex_list[tv0],global_vertex_list[tv1],global_vertex_list[tv2]));
            global_face_cnt++;

        }
    }


    //exit(0);
    // 接下去就是求交

    cout << " 接下去时求交"<<endl;
    std::vector <std::shared_ptr<std::thread> > face_generate_ray_detect_thread_pool(thread_num);
    for(int i=0;i<thread_num;i++) {
        face_generate_ray_detect_thread_pool[i] = make_shared<std::thread>([&](int now_id) {
            int each_grid_cnt = -1;
            for (auto each_grid = frame_grid_mp.begin(); each_grid != frame_grid_mp.end(); each_grid++) {
                each_grid_cnt++;
                if (each_grid_cnt % thread_num != now_id)continue; //todo 这里需要开启
                if (each_grid_cnt % (thread_num * 5) == now_id)
                    printf("face_generate_ray_detect_thread_pool %d/%d\n", each_grid_cnt, (int) frame_grid_mp.size());

                unordered_map<unsigned long long , pair<int,int> > triangle_mapping;
                std::list<K2::Triangle_3>l;
                K2::Point_3 small = iGameVertex_to_Point_K2(getGridVertex(each_grid->first, 0));
                K2::Point_3 big = iGameVertex_to_Point_K2(getGridVertex(each_grid->first, 7));
                std::vector<K2::Point_3> ps;
                for(int i=0;i<8;i++){
                    ps.push_back(getGridK2Vertex(small,big,i));
                }

                CGAL::Polyhedron_3<K2> frame_poly;
                PMP::polygon_soup_to_polygon_mesh(ps, each_grid_face_list, frame_poly, CGAL::parameters::all_default());

                std::function<bool(K2::Point_3)> vertex_in_frame = [&](K2::Point_3 v){
                    return small.x() <= v.x() && v.x() <= big.x() &&
                           small.y() <= v.y() && v.y() <= big.y() &&
                           small.z() <= v.z() && v.z() <= big.z();
                };

                for(int i=0;i<each_grid->second.field_list.size();i++){
                    int field_id = each_grid->second.field_list[i];
                    bool useful = false;
                    if(coverage_field_list[field_id].in_field(midpoint(K2::Segment_3(ps[0],ps[7])))){
                        useful = true;
                    }
                    if(!useful)
                    {
                        if(vertex_in_frame(coverage_field_list[field_id].center)){
                            useful = true;
                        }
                    }
                    if(!useful)
                        useful = CGAL::Polygon_mesh_processing::do_intersect(*coverage_field_list[field_id].poly,frame_poly);
                    if(useful){
                        for(int j=0;j<coverage_field_list[field_id].renumber_bound_face_global_id.size();j++){
                            int fid_global = coverage_field_list[field_id].renumber_bound_face_global_id[j];
                            K2::Point_3 v0 = global_vertex_list[global_face_list[fid_global].idx0];
                            K2::Point_3 v1 = global_vertex_list[global_face_list[fid_global].idx1];
                            K2::Point_3 v2 = global_vertex_list[global_face_list[fid_global].idx2];
                            K2::Triangle_3 tri(v0,v1,v2);
                            triangle_mapping[tri.id()] = {field_id,fid_global};
                            l.push_back(tri);
                        }
                    }
                }
                Tree aabb_tree(l.begin(),l.end());//todo 明天把这里打开debug一下看看？？

                for(int i=0;i<each_grid->second.global_face_list.size();i++){
                    int global_face_id = each_grid->second.global_face_list[i];
                    if(global_face_list[global_face_id].useful == false)continue;
                    K2::Point_3 v0 = global_vertex_list[global_face_list[global_face_id].idx0];
                    K2::Point_3 v1 = global_vertex_list[global_face_list[global_face_id].idx1];
                    K2::Point_3 v2 = global_vertex_list[global_face_list[global_face_id].idx2];

                    K2::Vector_3 ray_vec = K2::Triangle_3(v0,v1,v2).supporting_plane().orthogonal_vector();

                    K2::Ray_3 ray(global_face_list[global_face_id].center,-ray_vec);
                    std::list< Tree::Intersection_and_primitive_id<K2::Ray_3>::Type> intersections;
                    aabb_tree.all_intersections(ray,std::back_inserter(intersections));
                    for(auto item : intersections) {
                        pair<int,int> belong = triangle_mapping[item.second->id()];
                        int which_field = belong.first;
                        int which_id = belong.second;
                        if(which_field == global_face_list[global_face_id].field_id)continue;
                        if(const K2::Point_3* p = boost::get<K2::Point_3>(&(item.first))){
                            if(*p != global_face_list[global_face_id].center)
                                global_face_list[global_face_id].ray_detect_map[which_field].emplace_back(*p,which_id);
//                            else
//                                global_face_list[global_face_id].special_face_id.insert(which_field);
                        }
                        else{
                            global_face_list[global_face_id].special_field_id.insert(which_field);
                        }
                    }
                    // 参考3164行的写法写
                    // 注意事项：1. 建立射线的时候要所有面
                    // 注意事项：2。可能产生重复
                }

            }
        },i);
    }
    for(int i=0;i<thread_num;i++)
        face_generate_ray_detect_thread_pool[i]->join();


    cout <<"start generate"<<endl;

    std::vector <std::shared_ptr<std::thread> > global_face_final_generate_thread_pool(thread_num);
    atomic<int> flag = 0;
    for(int i=0;i<thread_num;i++) {
        global_face_final_generate_thread_pool[i] = make_shared<std::thread>([&](int now_id) {
            std::list<K2::Triangle_3>origin_face_list;
            for(int i=0;i<mesh->FaceSize();i++){
                K2::Point_3 v0(mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(0)].x(),
                               mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(0)].y(),
                               mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(0)].z());

                K2::Point_3 v1(mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(1)].x(),
                               mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(1)].y(),
                               mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(1)].z());

                K2::Point_3 v2(mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(2)].x(),
                               mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(2)].y(),
                               mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(2)].z());
                origin_face_list.emplace_back(v0,v1,v2);
            }

            Tree origin_face_tree(origin_face_list.begin(),origin_face_list.end());

            for (int i = 0; i < global_face_list.size(); i++) {
                if(i %(thread_num*10) ==0)
                    cout << "global_face_list:" << i <<"/"<<global_face_list.size()<< endl;
                if (i % thread_num != now_id)continue;
                if(global_face_list[i].useful == false )continue;
                for(std::unordered_map<int,vector<pair<K2::Point_3,int> > >::iterator it = global_face_list[i].ray_detect_map.begin();
                it != global_face_list[i].ray_detect_map.end(); it++
                ){////field 交点 对应的面
                    unordered_set<int>se;
                    std::vector<K2::Point_3> intersection_v;
                    for(int j=0;j<it->second.size();j++){
                        if(se.count(it->second[j].second))continue;
                        se.insert(it->second[j].second);
                        intersection_v.push_back(it->second[j].first);
                    }
                    sort(intersection_v.begin(),intersection_v.end(),[&](const K2::Point_3 &a,const K2::Point_3 &b){
                        if(a.x() != b.x()){
                            return a.x() < b.x();
                        }
                        else if(a.y() != b.y()){
                            return a.y() < b.y();
                        }
                        return a.z() < b.z();
                    });
                    int old_size = intersection_v.size();
                    intersection_v.resize(unique(intersection_v.begin(),intersection_v.end())-intersection_v.begin());
                    if(old_size != intersection_v.size()){
                        global_face_list[i].special_field_id.insert(it->first);
                        continue;
                    }
                    if(intersection_v.size()%2 == 1) {

                        global_face_list[i].useful = -200;
                        break;
                    }
                    // 1 是自己  // 若全无交点，可想而知一定在外，若有交点，那不一定在外，则正方向被覆盖，然后需要探测负方向，和自己
                    // 这里加入自己射线正负方向的寻找
                    // 2 出现多次 加入
                }


                if(  global_face_list[i].useful>0 ){
                    K2::Point_3 v0 = global_vertex_list[global_face_list[i].idx0];
                    K2::Point_3 v1 = global_vertex_list[global_face_list[i].idx1];
                    K2::Point_3 v2 = global_vertex_list[global_face_list[i].idx2];
                    if (running_mode == 2) {
                        if (!cgal_polygon->inMesh(centroid(K2::Triangle_3(v0, v1, v2)))) {
                            global_face_list[i].useful = -300;
                        }
                    }
                    else{
                        if (!cgal_polygon->outMesh(centroid(K2::Triangle_3(v0, v1, v2)))) {
                            global_face_list[i].useful = -300;
                        }
                    }

                    if(sqrt(CGAL::to_double(origin_face_tree.squared_distance(centroid(K2::Triangle_3(v0,v1,v2))))) < CGAL::Epeck::FT(myeps)){
                        global_face_list[i].useful = -300;
                    }
                }
            }
        },i);
    }
    for(int i=0;i<thread_num;i++)
        global_face_final_generate_thread_pool[i]->join();

    std::vector <std::shared_ptr<std::thread> > special_case_detect(thread_num);
    for(int i=0;i<thread_num;i++) {
        special_case_detect[i] = make_shared<std::thread>([&](int now_id) {
            for (int i = 0; i < global_face_list.size(); i++) {
                if (i % thread_num != now_id)continue;
                bool inner_flag = true;
                if(global_face_list[i].useful > 0 && !global_face_list[i].special_field_id.empty() )[[unlikely]]{
                    cout <<"meeting rare case, do special method"<<std::endl;
                    K2::Triangle_3 tri(global_vertex_list[global_face_list[i].idx0],
                                       global_vertex_list[global_face_list[i].idx1],
                                       global_vertex_list[global_face_list[i].idx2]);
                    for(auto field_id : global_face_list[i].special_field_id){
                        if(in_single_coverage_field_ray_detect(field_id,tri)){
                            inner_flag = false;
                            break;
                        }
                    }
                    if(!inner_flag){
                        global_face_list[i].useful = -200;
                    }
                }
            }
        },i);
    }
    for(int i=0;i<thread_num;i++)
        special_case_detect[i]->join();


    for(int i=0;i<global_vertex_list.size();i++){
        fprintf(file6,"v %lf %lf %lf\n",CGAL::to_double(global_vertex_list[i].x()),
                CGAL::to_double(global_vertex_list[i].y()),
                CGAL::to_double(global_vertex_list[i].z()));
    }

    double sum_avg_edge = 0;

    for (int i = 0; i < global_face_list.size(); i++) {
        cout <<i <<" "<< global_face_list[i].useful << endl;
        if(global_face_list[i].useful>0) {
            if(running_mode == 1) {
                sum_avg_edge += sqrt(CGAL::to_double((CGAL::squared_distance(global_vertex_list[global_face_list[i].idx0] , global_vertex_list[global_face_list[i].idx1]))));
                sum_avg_edge += sqrt(CGAL::to_double((CGAL::squared_distance(global_vertex_list[global_face_list[i].idx2] , global_vertex_list[global_face_list[i].idx1]))));
                sum_avg_edge += sqrt(CGAL::to_double((CGAL::squared_distance(global_vertex_list[global_face_list[i].idx0] , global_vertex_list[global_face_list[i].idx2]))));

                fprintf(file6, "f %d %d %d\n", global_face_list[i].idx0 + 1, global_face_list[i].idx2 + 1,
                        global_face_list[i].idx1 + 1);
            }
            else {
                fprintf(file6, "f %d %d %d\n", global_face_list[i].idx0 + 1, global_face_list[i].idx1 + 1,
                        global_face_list[i].idx2 + 1);
            }
        }
    }
    fclose(file6);
    sum_avg_edge /=(global_face_list.size()*3*20);


    string new_name = (input_filename + "_tmp.obj");
    string cmd;
//    for(int i=0;i<4;i++)new_name.pop_back();
    if(0) {
        cmd =
                ("../fTetWild/build/FloatTetwild_bin -i" + (input_filename + "_tmp.obj")+
                (" && mv " + new_name + "__sf.obj " + input_filename + "_final_result.obj");
    }
    else{
        cmd = ("../TetWild/build/Tetwild " + (input_filename + "_tmp.obj") ) +
                     (" && mv "+new_name+"__sf.obj " + input_filename+"_final_result.obj" );

    }

    //system(("mv "+new_name+"__sf.obj " + input_filename+"_result.obj" ).c_str());
    //Remeshing().run((input_filename + "_result.obj").c_str());
    cout << cmd << endl;
    system(cmd.c_str());
    return 0;
}
// 1 2 3 4 5 6
// 5 1 2 3 7 8
//




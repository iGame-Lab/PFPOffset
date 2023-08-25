//
// Created by rainbowwing on 2023/8/25.
//

#ifndef THICKEN2_GRID_H
#define THICKEN2_GRID_H

double stx, sty, stz;
double grid_len = 1e10;

struct grid {
    int x;
    int y;
    int z;
    grid(){
        x=-1;
        y=-1;
        z=-1;
    }
    grid(int x,int y,int z){
        this->x=x;
        this->y=y;
        this->z=z;
    }
    friend bool operator<(const grid &a, const grid &b) {
        if (a.x != b.x)
            return a.x < b.x;
        if (a.y != b.y)
            return a.y < b.y;
        return a.z < b.z;
    }

    friend bool operator==(const grid &a, const grid &b) {
        return a.x == b.x && a.y == b.y && a.z == b.z;
    }

    bool valid(){
        return x>=0 && y>=0 && z>=0;
    }
};

std::hash<int>int_hash;
struct grid_hash{
    size_t operator () (grid x) const {
        return int_hash(x.x) ^ int_hash(x.y<<6) ^ int_hash(x.z<<12);
    }};
struct grid_equal
{
    bool operator() (grid a,  grid b) const {
        return a.x == b.x  &&  a.y == b.y &&  a.z == b.z;
    }
};



vector <vector<int>> GridVertexDir = {{0, 0, 0}, //0
                                      {1, 0, 0}, //1
                                      {0, 1, 0}, //2
                                      {0, 0, 1}, //3
                                      {1, 1, 0}, //4
                                      {1, 0, 1}, //5
                                      {0, 1, 1}, //6
                                      {1, 1, 1}}; //7


vector <vector<int> > container_grid_face = {{0,3,5,1},
                                             {0,2,6,3},
                                             {0,1,4,2},
                                             {2,4,7,6},
                                             {3,6,7,5},
                                             {1,5,7,4}
};


MeshKernel::iGameVertex getGridVertex(grid g, int k) {
    double gsx = stx + g.x * grid_len;
    double gsy = sty + g.y * grid_len;
    double gsz = stz + g.z * grid_len;
    return MeshKernel::iGameVertex(gsx + GridVertexDir[k][0] * grid_len,
                                   gsy + GridVertexDir[k][1] * grid_len,
                                   gsz + GridVertexDir[k][2] * grid_len);
}

vector <vector<int>> DirectedGridEdge = {{1, 2, 3},
                                         { 4, 5},
                                         { 4, 6},
                                         {5, 6},
                                         { 7},
                                         { 7},
                                         { 7}};

MeshKernel::iGameVertex getGridiGameVertex(const MeshKernel::iGameVertex& small, const MeshKernel::iGameVertex& big, int k) {

    double x = small.x();
    double y = small.y();
    double z = small.z();
    if(GridVertexDir[k][0] > 0 )
        x = big.x();
    if(GridVertexDir[k][1] > 0 )
        y = big.y();
    if(GridVertexDir[k][2] > 0 )
        z = big.z();
    return MeshKernel::iGameVertex(x,y,z);
}

K2::Point_3 getGridK2Vertex(const K2::Point_3& small, const K2::Point_3& big, int k) {

    CGAL::Epeck::FT x = small.x();
    CGAL::Epeck::FT y = small.y();
    CGAL::Epeck::FT z = small.z();
    if(GridVertexDir[k][0] > 0 )
        x = big.x();
    if(GridVertexDir[k][1] > 0 )
        y = big.y();
    if(GridVertexDir[k][2] > 0 )
        z = big.z();
    return K2::Point_3(x,y,z);
}

grid vertex_to_grid(MeshKernel::iGameVertex v) {
    int x = int((v.x() - stx) / grid_len + myeps); // 先不管
    int y = int((v.y() - sty) / grid_len + myeps); // 先不管
    int z = int((v.z() - stz) / grid_len + myeps); // 先不管
    return grid{x, y, z};
}



struct GridVertex {
    int grid_type;
    vector <MeshKernel::iGameFaceHandle> field_list;
    vector<MeshKernel::iGameVertex>build_v;
    vector<int>build_v_global;
    vector<vector<size_t> >build_face;
    vector<K2::Triangle_3> generate_face_list;
    map<int,vector<int> >  field_face_though_list;
    map<size_t,pair<int,int> > face_hash_id_map;
    vector<int> global_face_list;
    std::list<K2::Triangle_3 > build_aabb_tree_triangle_list;
    GridVertex() {
        grid_type = -1;
    }
};

vector<vector<size_t> >face_type_012{{0,1,2}};
std::vector<std::vector<std::size_t> > each_grid_face_list;

bool check_vertex_in_grid(K2::Point_3 small, K2::Point_3 big, K2::Point_3 v){
    return small.x() <= v.x() && v.x() <= big.x() &&
           small.y() <= v.y() && v.y() <= big.y() &&
           small.z() <= v.z() && v.z() <= big.z();
};
bool check_triangle_through_grid(K2::Point_3 small, K2::Point_3 big,CGAL::Polyhedron_3<K2>& frame_poly , K2::Triangle_3 tri) {
    if(check_vertex_in_grid(small,big,centroid(tri))){
        return true;
    }
    CGAL::Polyhedron_3<K2> this_face;
    vector<K2::Point_3> vs{tri.vertex(0),tri.vertex(1),tri.vertex(2)};
    PMP::polygon_soup_to_polygon_mesh(vs, face_type_012, this_face, CGAL::parameters::all_default());
    return CGAL::Polygon_mesh_processing::do_intersect(frame_poly,this_face);
};


vector<vector<int> > container_grid_dir{{-1,-1,-1},{-1,-1,0},{-1,-1,1},
                                        {-1,0,-1},{-1,0,0},{-1,0,1},
                                        {-1,1,-1},{-1,1,0},{-1,1,1},
                                        {0,-1,-1},{0,-1,0},{0,-1,1},
                                        {0,0,-1},{0,0,1},
                                        {0,1,-1},{0,1,0},{0,1,1},
                                        {1,-1,-1},{1,-1,0},{1,-1,1},
                                        {1,0,-1},{1,0,0},{1,0,1},
                                        {1,1,-1},{1,1,0},{1,1,1},};
#endif //THICKEN2_GRID_H

#include"Mesh.h"
#include"string"
#include <queue>

// Mesh 定义
namespace MeshKernel {

    void Mesh::initBBox() {
        double bbox_min_x = 99999999, bbox_min_y = 99999999, bbox_min_z = 99999999;
        double bbox_max_x = -99999999, bbox_max_y = -99999999, bbox_max_z = -99999999;
        for (auto& vp : vertices_) {
            bbox_min_x = std::min(bbox_min_x, vp.second.x());
            bbox_min_y = std::min(bbox_min_y, vp.second.y());
            bbox_min_z = std::min(bbox_min_z, vp.second.z());
            bbox_max_x = std::max(bbox_max_x, vp.second.x());
            bbox_max_y = std::max(bbox_max_y, vp.second.y());
            bbox_max_z = std::max(bbox_max_z, vp.second.z());
        }
        BBoxMin = iGameVertex(bbox_min_x, bbox_min_y, bbox_min_z);
        BBoxMax = iGameVertex(bbox_max_x, bbox_max_y, bbox_max_z);
        printf("volume mesh: BBox: (%.3f, %.3f, %.3f) --> (%.3f, %.3f, %.3f)\n",
               BBoxMin.x(), BBoxMin.y(), BBoxMin.z(), BBoxMax.x(), BBoxMax.y(), BBoxMax.z());
    }
    iGameFaceHandle Mesh::AddFace(const std::vector<iGameVertexHandle>& _vhs,double face_move ) {
        //std::cout << " 现在加的面中的点为 : " << std::endl;
        //for (auto& vh : _vhs) {
        //	std::cout << vh << " : ";
        //}
        //std::cout << std::endl;
        std::vector<iGameEdgeHandle> ehs(_vhs.size());
        // 得到该面边的handle
        for (int i = 0; i < _vhs.size(); ++i) {
            if (i == 0) {
                ehs[i] = AddEdge(_vhs[_vhs.size() - 1], _vhs[i]);
            }
            else {
                ehs[i] = AddEdge(_vhs[i], _vhs[i - 1]);
            }
        }
        iGameFace f(_vhs, ehs);
        f.move_dist=face_move;
        // 如果该面已经存在，则返回面的handle
        if (Face2Fh_.count(f)) {
            //std::cout << "该面已经存在 序号为 : " << Face2Fh_[f] << "以及面中的各个边的序号为 : " << std::endl;
            iGameFaceHandle fh = Face2Fh_[f];
            //for (int i = 0; i < ehs.size(); i++) {
            //	std::cout << ehs[i] << " ";
            //}
            //for (int i = 0; i < 4; i++) {
            //	std::cout << faces_[fh].eh(i) << " ";
            //}
            //std::cout << std::endl;
            // 仍然更新 handle 和 面直接的映射
            faces_[fh] = f;											// 建立handle与该面之间的映射
            return Face2Fh_[f];
        }
            // 否则
        else {
            iGameFaceHandle fh = GenFaceHandle();						// 生成一个新的handle
            //std::cout << " 新的面的序号为 : " << fh << " 和面中的各个边的序号为 : " << std::endl;
            //for (int i = 0; i < ehs.size(); i++) {
            //	std::cout << ehs[i] << " ";
            //}
            faces_[fh] = f;											// 建立handle与该面之间的映射
            Face2Fh_[f] = fh;										// 建立该面与handle之间的映射
            AddFace2Neighbor(fh);									// 将该面添加至面所包含的点和边的相邻面中
            //for (int i = 0; i < 4; i++) {
            //	std::cout << faces_[fh].eh(i) << " ";
            //}
            //std::cout << std::endl;
            return Face2Fh_[f];										// 返回该面的handle
        }
    }
    iGameFaceHandle Mesh::AddFace(const std::vector<iGameVertexHandle>& _vhs) {
        //std::cout << " 现在加的面中的点为 : " << std::endl;
        //for (auto& vh : _vhs) {
        //	std::cout << vh << " : ";
        //}
        //std::cout << std::endl;
        std::vector<iGameEdgeHandle> ehs(_vhs.size());
        // 得到该面边的handle
        for (int i = 0; i < _vhs.size(); ++i) {
            if (i == 0) {
                ehs[i] = AddEdge(_vhs[_vhs.size() - 1], _vhs[i]);
            }
            else {
                ehs[i] = AddEdge(_vhs[i], _vhs[i - 1]);
            }
        }
        iGameFace f(_vhs, ehs);
        // 如果该面已经存在，则返回面的handle
        if (Face2Fh_.count(f)) {
            //std::cout << "该面已经存在 序号为 : " << Face2Fh_[f] << "以及面中的各个边的序号为 : " << std::endl;
            iGameFaceHandle fh = Face2Fh_[f];
            //for (int i = 0; i < ehs.size(); i++) {
            //	std::cout << ehs[i] << " ";
            //}
            //for (int i = 0; i < 4; i++) {
            //	std::cout << faces_[fh].eh(i) << " ";
            //}
            //std::cout << std::endl;
            // 仍然更新 handle 和 面直接的映射
            faces_[fh] = f;											// 建立handle与该面之间的映射
            return Face2Fh_[f];
        }
            // 否则
        else {
            iGameFaceHandle fh = GenFaceHandle();						// 生成一个新的handle
            //std::cout << " 新的面的序号为 : " << fh << " 和面中的各个边的序号为 : " << std::endl;
            //for (int i = 0; i < ehs.size(); i++) {
            //	std::cout << ehs[i] << " ";
            //}
            faces_[fh] = f;											// 建立handle与该面之间的映射
            Face2Fh_[f] = fh;										// 建立该面与handle之间的映射
            AddFace2Neighbor(fh);									// 将该面添加至面所包含的点和边的相邻面中
            //for (int i = 0; i < 4; i++) {
            //	std::cout << faces_[fh].eh(i) << " ";
            //}
            //std::cout << std::endl;
            return Face2Fh_[f];										// 返回该面的handle
        }
    }
    iGameEdgeHandle Mesh::AddEdge(const iGameVertexHandle& vh1, const iGameVertexHandle& vh2) {
        iGameEdge e(vh1, vh2);
        //std::cout << "边的两个端点的序号为 : " << vh1 << " : " << vh2 << std::endl;
        // 如果该边已经存在，则返回该边的handle
        if (Edge2Eh_.count(e)) {
            //std::cout << "此时该边已经存在 且边的序号为 : " << Edge2Eh_[e] << std::endl;
            return Edge2Eh_[e];
        }
            // 否则
        else {
            iGameEdgeHandle eh = GenEdgeHandle();						// 生成一个新的handle
            //std::cout << "生成边的序号为 : " << eh << std::endl;
            edges_[eh] = e;											// 建立handle与该边之间的映射
            Edge2Eh_[e] = eh;										// 建立该边与handle之间的映射
            AddEdge2Neighbor(eh);									// 将该边记录为其两个顶点的邻接边
            return Edge2Eh_[e];										// 返回该边的handle
        }
    }
    iGameVertexHandle Mesh::AddVertex(const iGameVertex& _v) {
        if (Vertex2Vh_.count(_v)) return Vertex2Vh_[_v];
        else {
            iGameVertexHandle vh = GenVertexHandle();
            vertices_[vh] = _v;
            Vertex2Vh_[_v] = vh;
            return Vertex2Vh_[_v];
        }
    }

    iGameVertexHandle Mesh::DeleteVertex(iGameVertexHandle _vh) {
        if (!vertices_.count(_vh)) return iGameVertexHandle(-1);        // 如果该点不存在，则返回-1
        else {
            // 删除相邻边元素（删除边时自然删除了相邻面）
            auto ve = NeighborEhOfVertex_[_vh];
            for (iGameEdgeHandle eh : ve) {
                DeleteEdge(eh);
            }
            // 删除顶点元素和以其为根据的所有邻接关系
            Vertex2Vh_.erase(vertices_[_vh]);
            vertices_.erase(_vh);
            NeighborEhOfVertex_.erase(_vh);
            NeighborFhOfVertex_.erase(_vh);
            return _vh;
        }

    }
    iGameEdgeHandle Mesh::DeleteEdge(iGameEdgeHandle _eh) {
        if (!edges_.count(_eh)) return iGameEdgeHandle(-1);             // 如果该边不存在，则返回-1
        else {
            // 删除邻接关系
            iGameEdge e(edges_[_eh]);
            for (int i = 0; i < 2; ++i) {
                iGameVertexHandle ev = e.vh(i);
                NeighborEhOfVertex_[ev].erase(_eh);                // 删除点邻接边
            }
            // 删除相邻面元素
            auto ef = NeighborFhOfEdge_[_eh];
            for (iGameFaceHandle fh : ef) {
                DeleteFace(fh);
            }
            // 删除边元素和以其为根据的所有邻接关系
            Edge2Eh_.erase(edges_[_eh]);
            edges_.erase(_eh);
            NeighborFhOfEdge_.erase(_eh);
            return _eh;
        }
    }
    iGameFaceHandle Mesh::DeleteFace(iGameFaceHandle _fh) {
        if (!faces_.count(_fh)) return iGameFaceHandle(-1);             // 如果该面不存在，则返回-1
        else {                                                     // 如果该面存在，则返回删除的这个面的handle
            // 删除邻接关系
            iGameFace f(faces_[_fh]);
            for (int i = 0; i < f.size(); ++i) {
                iGameVertexHandle fv = f.vh(i);
                iGameEdgeHandle fe = f.eh(i);
                NeighborFhOfVertex_[fv].erase(_fh);               // 删除点邻接面
                NeighborFhOfEdge_[fe].erase(_fh);                 // 删除边邻接面
            }
            // 删除面元素
            Face2Fh_.erase(faces_[_fh]);
            faces_.erase(_fh);
            return _fh;
        }
    }

    void Mesh::AddFace2Neighbor(const iGameFaceHandle& _fh)
    {
        iGameFace f = faces_[_fh];
        size_t n = f.size();
        for (int i = 0; i < n; ++i) {
            NeighborFhOfVertex_[f.vh(i)].insert(_fh);
        }
        for (int i = 0; i < n; ++i) {
            NeighborFhOfEdge_[f.eh(i)].insert(_fh);
        }
    }

    void Mesh::AddEdge2Neighbor(const iGameEdgeHandle& _eh)
    {
        iGameEdge e = edges_[_eh];
        NeighborEhOfVertex_[e.vh1()].insert(_eh);
        NeighborEhOfVertex_[e.vh2()].insert(_eh);

    }

    void Mesh::DeleteFace2Neighbor(const iGameFaceHandle& _fh) {
        iGameFace f = faces_[_fh];
        size_t n = f.size();
        for (int i = 0; i < n; ++i) {
            NeighborFhOfVertex_[f.vh(i)].erase(_fh);
        }
        for (int i = 0; i < n; ++i) {
            NeighborFhOfEdge_[f.eh(i)].erase(_fh);
        }
    }
    void Mesh::DeleteEdge2Neighbor(const iGameEdgeHandle& _eh) {
        iGameEdge e = edges_[_eh];
        NeighborEhOfVertex_[e.vh1()].erase(_eh);
        NeighborEhOfVertex_[e.vh2()].erase(_eh);
    }

    Mesh& Mesh::operator=(const Mesh& _surfacemesh)
    {
        vertices_ = _surfacemesh.vertices_;
        edges_ = _surfacemesh.edges_;
        faces_ = _surfacemesh.faces_;
        Vertex2Vh_ = _surfacemesh.Vertex2Vh_;
        Edge2Eh_ = _surfacemesh.Edge2Eh_;
        Face2Fh_ = _surfacemesh.Face2Fh_;
        NeighborEhOfVertex_ = _surfacemesh.NeighborEhOfVertex_;
        NeighborFhOfVertex_ = _surfacemesh.NeighborFhOfVertex_;
        NeighborFhOfEdge_ = _surfacemesh.NeighborFhOfEdge_;
        VertexHandleID_ = _surfacemesh.VertexHandleID_;
        EdgeHandleID_ = _surfacemesh.EdgeHandleID_;
        FaceHandleID_ = _surfacemesh.FaceHandleID_;
        fast_iGameVertex = _surfacemesh.fast_iGameVertex;
        fast_iGameEdge =  _surfacemesh.fast_iGameEdge;
        fast_iGameFace=  _surfacemesh.fast_iGameFace;
        FastNeighborEhOfVertex_ = _surfacemesh.FastNeighborEhOfVertex_;
        FastNeighborFhOfVertex_ = _surfacemesh.FastNeighborFhOfVertex_;
        FastNeighborFhOfEdge_ = _surfacemesh.FastNeighborFhOfEdge_;

        is_plane_vertex = _surfacemesh.is_plane_vertex;

        return *this;

    }

    /*=========================读写元素===============================*/
    // 读取ID为i的顶点
    iGameVertex& Mesh::vertices(iGameVertexHandle _vh) {
        assert(vertices_.count(_vh));
        return vertices_[_vh];
    }
    const iGameVertex Mesh::vertices(iGameVertexHandle _vh) const {
        assert(vertices_.count(_vh));
        return vertices_.find(_vh)->second;                // unordered_map 的 [] 操作符不是常量成员函数，无法对常量函数使用
    }
    void Mesh::build_fast(){
        fast_iGameVertex.clear();
        fast_iGameEdge.clear();
        fast_iGameFace.clear();
        fast_iGameVertex.resize(VertexSize()+10);
        fast_iGameEdge.resize(EdgeSize()+10);
        fast_iGameFace.resize(FaceSize()+10);
        /*
        std::unordered_map<iGameVertexHandle, std::unordered_set<iGameEdgeHandle> > NeighborEhOfVertex_;          //点的邻接边
		std::unordered_map<iGameVertexHandle, std::unordered_set<iGameFaceHandle> > NeighborFhOfVertex_;          //点的邻接面
		std::unordered_map<iGameEdgeHandle, std::unordered_set<iGameFaceHandle> > NeighborFhOfEdge_;
         */
        FastNeighborEhOfVertex_.clear();
        FastNeighborEhOfVertex_.resize(NeighborEhOfVertex_.size()+10);
        FastNeighborFhOfVertex_.clear();
        FastNeighborFhOfVertex_.resize(NeighborFhOfVertex_.size()+10);
        FastNeighborFhOfEdge_.clear();
        FastNeighborFhOfEdge_.resize(NeighborFhOfEdge_.size()+10);

        for(auto i : NeighborFhOfEdge_){
            FastNeighborFhOfEdge_[i.first] = i.second;
        }
        for(auto i : NeighborFhOfVertex_){
            FastNeighborFhOfVertex_[i.first] = i.second;
        }
        for(auto i : NeighborEhOfVertex_){
            FastNeighborEhOfVertex_[i.first] = i.second;
        }


        for(auto i : allvertices()){
            fast_iGameVertex[i.first]=i.second;
        }
        for(auto i : alledges()){
            fast_iGameEdge[i.first]=i.second;
        }
        for(auto i : allfaces()){
            fast_iGameFace[i.first]=i.second;
        }
    }

//    void Mesh::init_plane_point(){
//        is_plane_vertex.resize(VertexSize());
//        for(auto i : allfaces()){
//            for(auto n1 : NeighborFh(i.first))
//                for(auto n2 : NeighborFh(i.first)){
//                    if(n1!= n2){
//                       auto normal1 = ((fast_iGameVertex[fast_iGameFace[n1].vh(1)] - fast_iGameVertex[fast_iGameFace[n1].vh(0)])
//                               % (fast_iGameVertex[fast_iGameFace[n1].vh(2)] - fast_iGameVertex[fast_iGameFace[n1].vh(0)])).normalize();
//
//                       auto normal2 = ((fast_iGameVertex[fast_iGameFace[n2].vh(1)] - fast_iGameVertex[fast_iGameFace[n2].vh(0)])
//                               % (fast_iGameVertex[fast_iGameFace[n2].vh(2)] - fast_iGameVertex[fast_iGameFace[n2].vh(0)])).normalize();
//
//                       if(normal1 * normal2 <0.8){
//                           is_plane_vertex[i.first] = true;
//                       }
//                    }
//                }
//        }
//    }

    // 读取ID为i的边
    iGameEdge& Mesh::edges(iGameEdgeHandle _eh) {
        assert(edges_.count(_eh));
        return edges_[_eh];
    }
    bool Mesh::edges_vaild(iGameEdgeHandle _eh) {
        return edges_.count(_eh);
    }
    const iGameEdge& Mesh::edges(iGameEdgeHandle _eh) const {
        assert(edges_.count(_eh));
        return edges_.find(_eh)->second;
    }
    // 读取ID为i的面
    iGameFace& Mesh::faces(iGameFaceHandle _fh) {
        assert(faces_.count(_fh));
        return faces_[_fh];
    }
    const iGameFace Mesh::faces(iGameFaceHandle _fh) const {
        assert(faces_.count(_fh));
        return faces_.find(_fh)->second;
    }
    bool Mesh::faces_vaild(iGameFaceHandle _fh) {
        return  faces_.count(_fh);
    }
    /*====================根据元素得到对应ID=========================*/
    const iGameVertexHandle Mesh::vertexhandle(iGameVertex _vertex) const {
        if (Vertex2Vh_.find(_vertex) != Vertex2Vh_.end()) return Vertex2Vh_.find(_vertex)->second;
        else return iGameVertexHandle(-1);
    }
    const iGameEdgeHandle Mesh::edgehandle(iGameEdge& _edge) const {
        if (Edge2Eh_.find(_edge) != Edge2Eh_.end()) return Edge2Eh_.find(_edge)->second;
        else return iGameEdgeHandle(-1);
    }
    const iGameFaceHandle Mesh::facehandle(iGameFace& _face) const {
        if (Face2Fh_.find(_face) != Face2Fh_.end()) return Face2Fh_.find(_face)->second;
        else return iGameFaceHandle(-1);
    }

    /*======================得到邻接关系============================*/
    // 顶点的邻接点
    // 先找邻接边，再找邻接点
    std::unordered_set<iGameVertexHandle> Mesh::NeighborVh(iGameVertexHandle _vh) {
        std::unordered_set<iGameVertexHandle> neighborvh;
        auto neighboreh = NeighborEh(_vh);
        // 存在邻接边是前提
        if (neighboreh.size()) {
            for (iGameEdgeHandle eh : neighboreh) {
                if (edges_[eh].vh1() != _vh) neighborvh.insert(edges_[eh].vh1());
                if (edges_[eh].vh2() != _vh) neighborvh.insert(edges_[eh].vh2());
            }
        }
        return neighborvh;
    }
    // 顶点的邻接边
    std::unordered_set<iGameEdgeHandle>& Mesh::NeighborEh(iGameVertexHandle _vh) {
        if (NeighborEhOfVertex_.count(_vh)) return NeighborEhOfVertex_[_vh];
        else return empty_ehs;               // 返回一个空的集合
    }
    // 顶点的邻接面
    std::unordered_set<iGameFaceHandle>& Mesh::NeighborFh(iGameVertexHandle _vh) {
        if (NeighborFhOfVertex_.count(_vh)) return NeighborFhOfVertex_[_vh];
        else return empty_fhs;               // 返回一个空的集合
    }
    // 顶点在边上的另一个顶点
    MeshKernel::iGameVertexHandle Mesh::NeighborVhFromEdge(iGameVertexHandle _vh, iGameEdgeHandle _eh) {
        iGameVertexHandle vh = (edges_[_eh].vh(0) == _vh) ? edges_[_eh].vh(1) : edges_[_eh].vh(0);
        // assert(vh == edges_[_eh].vh(0) || vh == edges_[_eh].vh(1));
        return vh;
    }

    // 边的邻接边
    // 两个顶点的所有邻接边去除当前边
    std::unordered_set<iGameEdgeHandle> Mesh::NeighborEh(iGameEdgeHandle _eh) {
        assert(edges_.count(_eh));                     // 保证该边handle存在
        std::unordered_set<iGameEdgeHandle> neighboreh;     // 保存输出的结果
        int k = 0;                                     // 遍历两个顶点
        while (k < 2) {
            iGameVertexHandle vh = edges_[_eh].vh(k);
            auto vhneighboreh = NeighborEh(vh);        // 得到点的邻接边
            for (iGameEdgeHandle eh : vhneighboreh) {
                if (eh != _eh) neighboreh.insert(eh);
            }
            ++k;
        }

        return neighboreh;
    }
    // 边的邻接面
    std::unordered_set<iGameFaceHandle>& Mesh::NeighborFh(iGameEdgeHandle _eh) {
        if (NeighborFhOfEdge_.count(_eh)) return NeighborFhOfEdge_[_eh];
        else return empty_fhs;               // 返回一个空的集合
    }
    // 面的邻接面
    // 邻接面：有一条相同边
    std::unordered_set<iGameFaceHandle> Mesh::NeighborFh(iGameFaceHandle _fh) {
        assert(faces_.count(_fh));                     // 保证该边handle存在
        std::unordered_set<iGameFaceHandle> neigborface;
        int k = 0;                                     // 遍历两个顶点
        size_t facesize = faces_[_fh].size();
        while (k < facesize) {
            iGameEdgeHandle eh = faces_[_fh].eh(k);
            auto ehneighborfh = NeighborFh(eh);        // 得到点的邻接边
            for (iGameFaceHandle fh : ehneighborfh) {
                if (fh != _fh) neigborface.insert(fh);
            }
            ++k;
        }
        return neigborface;
    }
    // 邻接面2：有一个公共顶点
    std::unordered_set<iGameFaceHandle> Mesh::Neighbor2Fh(iGameFaceHandle _fh) {
        assert(faces_.count(_fh));                     // 保证该边handle存在
        std::unordered_set<iGameFaceHandle> neigborface;
        auto v_indices = faces_[_fh].getVertexHandle();
        for (auto& v_idx : v_indices) {
            auto adjF = NeighborFh(v_idx);
            for (iGameFaceHandle fh : adjF) {
                if (fh != _fh) neigborface.insert(fh);
            }
        }
        return neigborface;
    }
}



// SurfaceMesh 定义
namespace MeshKernel {
    void SurfaceMesh::InitMesh(const std::vector<iGameVertex>& _vertices,
                               const std::vector<std::vector<iGameVertexHandle>>& _elements,
                               const std::vector<double >face_move) {
        std::vector<iGameVertexHandle>new_handle(_vertices.size());
        for (int i=0;i<_vertices.size();i++) {
            auto v =  _vertices[i];
            auto vh = AddVertex(iGameVertex(v.x(), v.y(), v.z()));
            new_handle[i]=vh;
        }
        for(int i=0;i<_elements.size();i++){
            auto f= _elements[i];
            if(f.size()<3)continue;
//            auto vh0 = AddVertex(iGameVertex(_vertices[f[0]].x(), v.y(), v.z()))
            auto tmp = AddFace({new_handle[f[0]],new_handle[f[1]],new_handle[f[2]]},face_move[i]);

        }

    }
    SurfaceMesh& SurfaceMesh::operator=(const SurfaceMesh& _surfacemesh) {
        if (this != &_surfacemesh) {
            Mesh::operator=(_surfacemesh);
        }
        return *this;
    }

    bool SurfaceMesh::isOnBoundary(iGameEdgeHandle eh) {
        auto fcnt = NeighborFh(eh).size();
        return fcnt == 1;
    }

    bool SurfaceMesh::isOnBoundary(iGameVertexHandle vh) {
        for (auto eh : NeighborEh(vh)) {
            if (isOnBoundary(eh)) {
                return true;
            }
        }
        return false;
    }

    bool SurfaceMesh::isOnBoundary(iGameFaceHandle fh) {
        auto face = faces(fh);
        for (auto eh : face.getEdgeHandle()) {
            if (isOnBoundary(eh)) {
                return true;
            }
        }
        return false;
    }


    bool SurfaceMesh::isTriangleMesh() {
        for (auto& fp : faces_) {
            if (fp.second.getVertexHandle().size() != 3) return false;
        }
        return true;
    }

    void SurfaceMesh::updateAllHandles() {
        int vcnt = VertexSize(), fcnt = FaceSize();
        std::vector<MeshKernel::iGameVertex> newVertices;
        std::vector<std::vector<MeshKernel::iGameVertexHandle>> newFaces;
        std::vector<double>v;
        std::unordered_map<int, int> mp;// old id to new id
        int idx = 0;
        for (auto& fp : allfaces()) {
            auto vhs = fp.second.getVertexHandle();
            for (auto& vh : vhs) {
                if (!mp.count(vh)) {
                    mp[vh] = idx++;
                    newVertices.push_back(vertices_[vh]);
                }
                vh = iGameVertexHandle(mp[vh]);
            }
            newFaces.push_back(vhs);
            v.push_back(fp.second.move_dist);
        }

        *this = MeshKernel::SurfaceMesh(newVertices, newFaces,v);
    }

    bool SurfaceMesh::isConnected(iGameVertexHandle vh1, iGameVertexHandle vh2) {
        for (auto vh : NeighborVh(vh1)) {
            if (vh == vh2) return true;
        }
        return false;
    }

    bool SurfaceMesh::isConnected(iGameEdgeHandle eh1, iGameEdgeHandle eh2) {
        if (!isValid(eh1) || !isValid(eh2)) return false;
        auto e1 = edges_[eh1];
        auto e2 = edges_[eh2];
        auto vh1 = e1.vh1(), vh2 = e1.vh2();
        auto vh3 = e2.vh1(), vh4 = e2.vh2();
        return (vh1 == vh3 || vh1 == vh4 || vh2 == vh3 || vh2 == vh4);
    }




    void SurfaceMesh::genAllFacesNormal() {

    }

    void SurfaceMesh::genAllVerticesNormal() {

    }

    void SurfaceMesh::genAllEdgesLength() {
        for (auto& ep : this->edges_) {

        }
    }


    iGameEdgeHandle SurfaceMesh::getEdgeHandle(iGameVertexHandle vh1, iGameVertexHandle vh2) {
        iGameEdgeHandle ret(-1);
        int vh_sum = vh1 + vh2;
        if (!isValid(vh1) || !isValid(vh2)) return ret;
        for (auto& eh : NeighborEh(vh1)) {
            auto& e = edges_[eh];
            if (vh_sum = e.vh1() + e.vh2()) {
                ret = eh;
                break;
            }
        }
        return ret;
    }

    size_t SurfaceMesh::getBoundaryVerticesCount() {
        size_t cnt = 0;
        for (auto& vp : vertices_) {
            if (isOnBoundary(vp.first))
                cnt++;
        }
        return cnt;
    }
    bool SurfaceMesh::isClosure() {

        for(auto i : this->alledges()){
            if(NeighborFh(i.first).size() == 1)
                return false;
        }
        return true;
    }




    double SurfaceMesh::getLength(iGameEdgeHandle eh) {
        assert(edges_.count(eh));
        auto& edge = edges_[eh];
        auto& v1 = vertices_[edge.vh1()];
        auto& v2 = vertices_[edge.vh2()];
        return (v1 - v2).norm();
    }


}




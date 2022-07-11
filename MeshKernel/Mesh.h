#pragma once
#include "Cell.h"


namespace MeshKernel {
	// 网格的概念
	class Mesh {
	public:
		iGameVertex BBoxMin, BBoxMax;
		iGameEdgeHandle getEdgeHandle(iGameVertexHandle vh1, iGameVertexHandle vh2) {// return the edge whose vex is vh1 and vh2
			iGameEdgeHandle ret(-1);
			for (auto& eh : NeighborEh(vh1)) {
				auto& e = edges_[eh];
				auto adj = e.vh1() + e.vh2() - vh1;
				if (vh2 == adj) {
					ret = eh;
					break;
				}
			}
			return ret;
		}

		void initBBox();
		inline bool isValid(iGameVertexHandle _vh) { return vertices_.count(_vh); }
		inline bool isValid(iGameEdgeHandle _eh) { return edges_.count(_eh); }
		inline bool isValid(iGameFaceHandle _fh) { return faces_.count(_fh); }
		inline bool isSingular(iGameVertexHandle vh) {
			return NeighborEhOfVertex_[vh].empty();
		}


		/*=========================读写元素===============================*/
			// 读取ID为i的顶点
		iGameVertex& vertices(iGameVertexHandle _vh);
		const iGameVertex vertices(iGameVertexHandle _vh) const;        // unordered_map 的 [] 操作符不是常量成员函数，无法对常量函数使用
		// 读取ID为i的边
		iGameEdge& edges(iGameEdgeHandle _eh);
		const iGameEdge& edges(iGameEdgeHandle _eh) const;
        bool edges_vaild(iGameEdgeHandle _eh);
		// 读取ID为i的面
		iGameFace& faces(iGameFaceHandle _fh);
		const iGameFace faces(iGameFaceHandle _fh) const;

        bool faces_vaild(iGameFaceHandle _fh);

		size_t vsize() const { return vertices_.size(); }
		size_t esize() const { return edges_.size(); }
		size_t fsize() const { return faces_.size(); }
		const std::unordered_map<iGameVertexHandle, iGameVertex>& allvertices() const { return vertices_; }
		const std::unordered_map<iGameEdgeHandle, iGameEdge>& alledges() const { return edges_; }
		const std::unordered_map<iGameFaceHandle, iGameFace>& allfaces() const { return faces_; }
        iGameFaceHandle AddFace(const std::vector<iGameVertexHandle>& _vhs,double face_move );
        /*====================根据元素得到对应ID=========================*/
		const iGameVertexHandle vertexhandle(iGameVertex _vertex) const;
		const iGameEdgeHandle edgehandle(iGameEdge& _edge) const;
		const iGameFaceHandle facehandle(iGameFace& _face) const;
		/*======================得到邻接关系============================*/
		// 顶点的邻接点
		std::unordered_set<iGameVertexHandle> NeighborVh(iGameVertexHandle _vh);
		// 顶点的邻接边
		std::unordered_set<iGameEdgeHandle>& NeighborEh(iGameVertexHandle _vh);
		// 顶点的邻接面
		std::unordered_set<iGameFaceHandle>& NeighborFh(iGameVertexHandle _vh);
		// 顶点在边上的另一个顶点
		iGameVertexHandle NeighborVhFromEdge(iGameVertexHandle _vh, iGameEdgeHandle _eh);
		// 边的邻接边
		std::unordered_set<iGameEdgeHandle> NeighborEh(iGameEdgeHandle _eh);
		// 边的邻接面
		std::unordered_set<iGameFaceHandle>& NeighborFh(iGameEdgeHandle _eh);
		// 面的邻接面
		std::unordered_set<iGameFaceHandle> NeighborFh(iGameFaceHandle _fh);// share common edge
		std::unordered_set<iGameFaceHandle> Neighbor2Fh(iGameFaceHandle _fh);// share common vertex
		/*=========================添加元素=============================*/
		iGameVertexHandle AddVertex(const iGameVertex& _v);
		iGameEdgeHandle AddEdge(const iGameVertexHandle& _vh1, const iGameVertexHandle& _vh2);
		iGameFaceHandle AddFace(const std::vector<iGameVertexHandle>& _vhs);

		/*=========================删除元素=============================*/
		// 删除低级元素时删除一切邻接的高级元素
		// 删除高级元素时保留仍被其他元素使用的低级元素
		iGameVertexHandle DeleteVertex(iGameVertexHandle _vh);
		iGameEdgeHandle DeleteEdge(iGameEdgeHandle _eh);
		iGameFaceHandle DeleteFace(iGameFaceHandle _fh);


		/*=========================生成唯一ID========================*/
		iGameVertexHandle GenVertexHandle() { return (iGameVertexHandle)VertexHandleID_++; }
		iGameEdgeHandle GenEdgeHandle() { return (iGameEdgeHandle)EdgeHandleID_++; }
		iGameFaceHandle GenFaceHandle() { return (iGameFaceHandle)FaceHandleID_++; }

		size_t VertexSize() { return vertices_.size(); }
		size_t EdgeSize() { return edges_.size(); }
		size_t FaceSize() { return faces_.size(); }

        void build_fast();

        //void init_plane_point();

        std::vector<iGameVertex>fast_iGameVertex;
        std::vector<iGameEdge>fast_iGameEdge;
        std::vector<iGameFace>fast_iGameFace;
        std::vector<bool>is_plane_vertex;


        std::vector<std::unordered_set<iGameEdgeHandle> > FastNeighborEhOfVertex_;          //点的邻接边
        std::vector<std::unordered_set<iGameFaceHandle> > FastNeighborFhOfVertex_;          //点的邻接面
        std::vector<std::unordered_set<iGameFaceHandle> > FastNeighborFhOfEdge_;

    protected:
		/*============下一个可使用的ID=========*/
		int VertexHandleID_ = 0;
		int EdgeHandleID_ = 0;
		int FaceHandleID_ = 0;
		std::unordered_set<iGameVertexHandle> empty_vhs;
		std::unordered_set<iGameEdgeHandle> empty_ehs;
		std::unordered_set<iGameFaceHandle> empty_fhs;

		/*=============修改邻接关系============*/
		// 将该面添加至面所包含的点和边的相邻面中
		void AddFace2Neighbor(const iGameFaceHandle& _fh);
		// 将该边添加至边所包含的点相邻边中
		void AddEdge2Neighbor(const iGameEdgeHandle& _eh);
		// Todo: Delete Neighbor
		void DeleteFace2Neighbor(const iGameFaceHandle& _fh);
		void DeleteEdge2Neighbor(const iGameEdgeHandle& _eh);
	protected:
		// handle到元素的对应
		std::unordered_map<iGameVertexHandle, iGameVertex> vertices_;
		std::unordered_map<iGameEdgeHandle, iGameEdge> edges_;
		std::unordered_map<iGameFaceHandle, iGameFace> faces_;




		// 元素到handle的对应
		std::unordered_map<iGameVertex, iGameVertexHandle> Vertex2Vh_;
		std::unordered_map<iGameEdge, iGameEdgeHandle> Edge2Eh_;
		std::unordered_map<iGameFace, iGameFaceHandle> Face2Fh_;

		// 邻接关系
		std::unordered_map<iGameVertexHandle, std::unordered_set<iGameEdgeHandle> > NeighborEhOfVertex_;          //点的邻接边
		std::unordered_map<iGameVertexHandle, std::unordered_set<iGameFaceHandle> > NeighborFhOfVertex_;          //点的邻接面
		std::unordered_map<iGameEdgeHandle, std::unordered_set<iGameFaceHandle> > NeighborFhOfEdge_;              //边的邻接面
	protected:
		virtual void InitMesh(const std::vector<iGameVertex>& _vertices,
			const std::vector<std::vector<iGameVertexHandle>>& _elements,std::vector<double>) = 0;                              // 必须重写
		Mesh& operator=(const Mesh& _mesh);


    };

	// 曲面网格
	class SurfaceMesh : public Mesh {
	public:
		/*==========================构造函数==============================*/
		SurfaceMesh() {};
		SurfaceMesh(const std::vector<iGameVertex>& _vertices, const std::vector<std::vector<iGameVertexHandle>>& _faces
        , const std::vector<double >face_move) {
			InitMesh(_vertices, _faces,face_move);
		}
		//SurfaceMesh(const SurfaceMesh& _surfacemesh);
		/*=============初始化网格=============*/
		virtual void InitMesh(const std::vector<iGameVertex>& _vertices,
			const std::vector<std::vector<iGameVertexHandle>>& _elements,
        const std::vector<double >face_move) ;
		SurfaceMesh& operator=(const SurfaceMesh& _surfacemesh);

		bool isOnBoundary(iGameEdgeHandle);
		bool isOnBoundary(iGameVertexHandle);
		bool isOnBoundary(iGameFaceHandle);

		void genNormal(iGameFaceHandle);
		void genNormal(iGameVertexHandle);
		void genAllFacesNormal();
		void genAllVerticesNormal();// 自然会生成所有面的法向量

		bool isClosure();
        void RevertMeshNormal();

		void genAllEdgesLength();
		void genLength(iGameEdgeHandle);
		double getLength(iGameEdgeHandle);
		void addNoise(double);

		iGameEdgeHandle getEdgeHandle(iGameVertexHandle vh1, iGameVertexHandle vh2);


		void updateAllHandles();

		bool isConnected(iGameVertexHandle, iGameVertexHandle);
		bool isConnected(iGameEdgeHandle, iGameEdgeHandle);

		bool isTriangleMesh();
		size_t getBoundaryVerticesCount();
	};


}
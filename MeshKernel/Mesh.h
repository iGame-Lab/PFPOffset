#pragma once
#include "Face.h"


namespace MeshKernel {
	// 网格的概念


	// 曲面网格
	class SurfaceMesh {
    public:
        MeshKernel::iGameVertex BBoxMax;
        MeshKernel::iGameVertex BBoxMin;

        std::vector<std::unordered_set<iGameFaceHandle> > FastNeighborFhOfVertex_;
        std::vector<std::unordered_set<iGameFaceHandle> > FastNeighborFhOfFace_;

        std::vector<iGameVertex>fast_iGameVertex;
       // std::vector<iGameEdge>fast_iGameEdge;
        std::vector<iGameFace>fast_iGameFace;

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
        void initBBox();
        int FaceSize(){
            return fast_iGameFace.size();
        }
        int VertexSize(){
            return fast_iGameVertex.size();
        }

	};


}
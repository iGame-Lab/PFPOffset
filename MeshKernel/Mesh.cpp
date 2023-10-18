#include"Mesh.h"
#include"string"
#include <queue>
#include <set>

// Mesh 定义




// SurfaceMesh 定义
namespace MeshKernel {
    void SurfaceMesh::initBBox() {
        double bbox_min_x = 99999999, bbox_min_y = 99999999, bbox_min_z = 99999999;
        double bbox_max_x = -99999999, bbox_max_y = -99999999, bbox_max_z = -99999999;
        for (auto &vp: fast_iGameVertex) {
            bbox_min_x = std::min(bbox_min_x, vp.x());
            bbox_min_y = std::min(bbox_min_y, vp.y());
            bbox_min_z = std::min(bbox_min_z, vp.z());
            bbox_max_x = std::max(bbox_max_x, vp.x());
            bbox_max_y = std::max(bbox_max_y, vp.y());
            bbox_max_z = std::max(bbox_max_z, vp.z());
        }
        BBoxMin = iGameVertex(bbox_min_x, bbox_min_y, bbox_min_z);
        BBoxMax = iGameVertex(bbox_max_x, bbox_max_y, bbox_max_z);
        printf("volume mesh: BBox: (%.3f, %.3f, %.3f) --> (%.3f, %.3f, %.3f)\n",
               BBoxMin.x(), BBoxMin.y(), BBoxMin.z(), BBoxMax.x(), BBoxMax.y(), BBoxMax.z());
    }


    void SurfaceMesh::InitMesh(const std::vector<iGameVertex>& _vertices,
                               const std::vector<std::vector<iGameVertexHandle>>& _elements,
                               const std::vector<double >face_move) {
        for(int i=0;i<_vertices.size();i++){
            this->fast_iGameVertex.push_back(_vertices[i]);
            this->FastNeighborFhOfVertex_.emplace_back();
        }
        for(int i=0;i<_elements.size();i++){
            this->fast_iGameFace.emplace_back(_elements[i], face_move[i]);
            this->FastNeighborFhOfFace_.emplace_back();
            for(int j=0;j<_elements[i].size();j++)
                this->FastNeighborFhOfVertex_[_elements[i][j].idx()].insert(MeshKernel::iGameFaceHandle(i));

        }
        for(int i =0;i< this->FastNeighborFhOfVertex_.size();i++){
            for(iGameFaceHandle j : this->FastNeighborFhOfVertex_[i]){
                std::set<int>se;
                se.insert(this->fast_iGameFace[j].vh(0));
                se.insert(this->fast_iGameFace[j].vh(1));
                se.insert(this->fast_iGameFace[j].vh(2));
                for(iGameFaceHandle k: this->FastNeighborFhOfVertex_[i])
                {
                    if(j==k)continue;
                    int cnt = se.count(this->fast_iGameFace[k].vh(0)) +
                            se.count(this->fast_iGameFace[k].vh(1)) +
                            se.count(this->fast_iGameFace[k].vh(2));
                    if(cnt == 2){
                        this->FastNeighborFhOfFace_[j].insert(k);
                    }
                }
            }
        }
    }
}




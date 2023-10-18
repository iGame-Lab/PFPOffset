#pragma once
#include"Edge.h"
#include<algorithm>
namespace MeshKernel {
	// 网格中的面片
	class iGameFace {
	public:
		/*=========================构造函数=============================*/

		iGameFace(const std::vector<iGameVertexHandle>& _vertices) {

			vertices_.assign(_vertices.begin(), _vertices.end());
            move_dist = -1;
		}
        iGameFace(const std::vector<iGameVertexHandle>& _vertices,double move_dist) {

            vertices_.assign(_vertices.begin(), _vertices.end());
            this->move_dist = move_dist;
        }
		iGameFace(const iGameFace& _f) {
			*this = _f;
		}
		/*=======================基本操作符==========================*/
		inline iGameFace& operator=(const iGameFace& _f) {
			vertices_.assign(_f.vertices_.begin(), _f.vertices_.end());
            move_dist=_f.move_dist;
			return *this;
		}

		/*=======================读写元素==========================*/
		const iGameVertexHandle vh(int k) const { return vertices_[k]; }            // 常量对象只支持读
		iGameVertexHandle& vh(int k) { return vertices_[k]; }		               // 非常量对象支持读写



		/*================顶点handle排列==========================*/


        double move_dist;

	private:
		// 保存该面所有边的handle
		// 保存该面所有点的handle
		std::vector<iGameVertexHandle> vertices_;
		// 保存该面的边数
		// 法向量
	};
}

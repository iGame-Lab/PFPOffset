#pragma once
#include <vector>
#include <iostream>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <cassert>
#include "Kernel.h"
#include "Handle.h"
#include <cmath>

namespace MeshKernel {
	// 网格中的顶点
	class iGameVertex {
	public:
		iGameVertex() :position(3, 0.0) {}
		iGameVertex(const iGameVertex& _v) :position(_v.position) {}
		iGameVertex(double x, double y, double z) :position{ x,y,z } {}

		/*===========读写单一分量元素=================*/
		// 常量成员只能读元素
		const double x() const { return position[0]; }
		const double y() const { return position[1]; }
		const double z() const { return position[2]; }

		// 普通成员可以读写元素
		double& x() { return position[0]; }
		double& y() { return position[1]; }
		double& z() { return position[2]; }

		/*============基本运算=================*/
		inline iGameVertex operator+(const iGameVertex& _rhs) const {                             //加法
			return iGameVertex(x() + _rhs.x(), y() + _rhs.y(), z() + _rhs.z());
		}
		inline iGameVertex operator-(const iGameVertex& _rhs) const {                             //减法
			return iGameVertex(x() - _rhs.x(), y() - _rhs.y(), z() - _rhs.z());
		}
		inline iGameVertex operator*(double k) const {                                        //数乘
			return iGameVertex(x() * k, y() * k, z() * k);
		}
		inline iGameVertex operator/(double k) const {                                        //数除
			return iGameVertex(x() / k, y() / k, z() / k);
		}
		inline double operator*(const iGameVertex& _rhs) const {                              //点乘
			return double(x() * _rhs.x() + y() * _rhs.y() + z() * _rhs.z());
		}
		inline iGameVertex operator%(const iGameVertex& _rhs) const {                              //叉乘
			return iGameVertex(y() * _rhs.z() - z() * _rhs.y(),
				z() * _rhs.x() - x() * _rhs.z(),
				x() * _rhs.y() - y() * _rhs.x());
		}
		inline iGameVertex& operator+=(const iGameVertex& _rhs) {
			position[0] += _rhs.x();
			position[1] += _rhs.y();
			position[2] += _rhs.z();
			return *this;
		}
		inline iGameVertex& operator-=(const iGameVertex& _rhs) {
			position[0] -= _rhs.x();
			position[1] -= _rhs.y();
			position[2] -= _rhs.z();
			return *this;
		}
		inline double operator[](int i) {
			assert(i >= 0 && i < position.size());
			return position[i];
		}
		inline iGameVertex& operator*=(double k) {
			position[0] *= k;
			position[1] *= k;
			position[2] *= k;
			return *this;
		}
		inline iGameVertex& operator/=(double k) {
			assert(k != 0.f);
			position[0] /= k;
			position[1] /= k;
			position[2] /= k;
			return *this;
		}
		inline double norm() const {                                                      //模长
			return sqrt(x() * x() + y() * y() + z() * z());
		}
		inline double norm2() const {                                                      //模长方
			return x() * x() + y() * y() + z() * z();
		}
		inline iGameVertex normalize() const {                                                 //单位化
			auto m = this->norm();
			return iGameVertex(x() / m, y() / m, z() / m);
		}

		inline double dist(const iGameVertex& _rhs) const {                                   //计算俩个向量的距离
			return (*this - _rhs).norm();
		}
		inline void setPosition(double x, double y, double z) {
			position.resize(3);
			position[0] = x;
			position[1] = y;
			position[2] = z;
		}
		inline void setPosition(std::vector<double> Pos) {
			assert(Pos.size() == 3);
			position = Pos;
		}
		inline void setX(double x) {
			position[0] = x;
		}
		inline void setY(double y) {
			position[1] = y;
		}
		inline void setZ(double z) {
			position[2] = z;
		}
		inline void setNormal(double x, double y, double z) {
			normal.resize(3);
			normal[0] = x;
			normal[1] = y;
			normal[2] = z;
		}
		inline void setNormal(std::vector<double> N) {
			assert(N.size() == 3);
			normal = N;
		}
		inline double getNormalX() {
			if (normal.size() == 3) return normal[0];
			return 0.f;
		}
		inline double getNormalY() {
			if (normal.size() == 3) return normal[1];
			return 0.f;
		}
		inline double getNormalZ() {
			if (normal.size() == 3) return normal[2];
			return 0.f;
		}
		inline void setColor(double x, double y, double z) {
			color.resize(3);
			color[0] = x;
			color[1] = y;
			color[2] = z;
		}
		inline void setColor(std::vector<double> C) {
			assert(C.size() == 3);
			color = C;
		}
		inline double getColorR() {
			if (color.size() == 3) return color[0];
			return 0.f;
		}
		inline double getColorG() {
			if (color.size() == 3) return color[1];
			return 0.f;
		}
		inline double getColorB() {
			if (color.size() == 3) return color[2];
			return 0.f;
		}

		/*===========比较运算=================*/
		inline bool operator==(const iGameVertex& _rhs) const {
			return (x() == _rhs.x() && y() == _rhs.y() && z() == _rhs.z());
		}
		/*===========赋值运算=================*/
		//返回值为引用，支持连等
		iGameVertex& operator=(const iGameVertex& _v) {
			position = _v.position;
            return *this;
		}
	private:
		std::vector<double> position;
		std::vector<double> normal;
		std::vector<double> color;
	};
}

/*======================特化V3f的哈希映射============================*/
//冲突时会调用相等函数
namespace std
{
	template<> struct hash<MeshKernel::iGameVertex>
	{
		size_t operator()(const MeshKernel::iGameVertex& v)const
		{
			return hash<double>()(v.x()) ^ hash<double>()(v.y()) ^ hash<double>()(v.z());
		}
	};
}
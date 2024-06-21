#include "Cartesian.h"
namespace CGAL
{
	template <typename T>
	static Vector_3<T> cross_product(Vector_3<T> vt1, Vector_3<T> vt2)
	{
		return Vector_3<T>(vt1.m_tY * vt2.m_tZ - vt1.m_tZ * vt2.m_tY, vt1.m_tZ * vt2.m_tX - vt1.m_tX * vt2.m_tZ, vt1.m_tX * vt2.m_tY - vt1.m_tY * vt2.m_tX);
	}
	template <typename T>
	static Point_3<T> barycenter(Point_3<T> pt1, T _tW1, Point_3<T> pt2, T _tW2)
	{
		T _tSum = _tW1 + _tW2;
		_tW1 /= _tSum;
		_tW2 /= _tSum;
		return Point_3<T>(pt1.m_tX * _tW1 + pt2.m_tX * _tW2, pt1.m_tY * _tW1 + pt2.m_tY * _tW2, pt1.m_tZ * _tW1 + pt2.m_tZ * _tW2);
	}
	template <typename T>
	static Point_3<T> barycenter(Point_3<T> pt1, T _tW1, Point_3<T> pt2, T _tW2, Point_3<T> pt3, T _tW3)
	{
		T _tSum = _tW1 + _tW2 + _tW3;
		_tW1 /= _tSum;
		_tW2 /= _tSum;
		_tW3 /= _tSum;
		return Point_3<T>(pt1.m_tX * _tW1 + pt2.m_tX * _tW2 + pt3.m_tX * _tW3, pt1.m_tY * _tW1 + pt2.m_tY * _tW2 + pt3.m_tY * _tW3, pt1.m_tZ * _tW1 + pt2.m_tZ * _tW2 + pt3.m_tZ *_tW3);
	}
	template <typename T>
	static Point_3<T> barycenter(Point_3<T> pt1, T _tW1, Point_3<T> pt2, T _tW2, Point_3<T> pt3, T _tW3, Point_3<T> pt4, T _tW4)
	{
		T _tSum = _tW1 + _tW2 + _tW3 + _tW4;
		_tW1 /= _tSum;
		_tW2 /= _tSum;
		_tW3 /= _tSum;
		_tW4 /= _tSum;
		return Point_3<T>(pt1.m_tX * _tW1 + pt2.m_tX * _tW2 + pt3.m_tX * _tW3 + pt4.m_tX * _tW4, pt1.m_tY * _tW1 + pt2.m_tY * _tW2 + pt3.m_tY * _tW3 + pt4.m_tY * _tW4, pt1.m_tZ * _tW1 + pt2.m_tZ * _tW2 + pt3.m_tZ *_tW3 + pt4.m_tZ * _tW4);
	}
	template <typename T>
	static T squared_length(Vector_3<T> vtVector)
	{
		return vtVector.m_tX * vtVector.m_tX + vtVector.m_tY * vtVector.m_tY + vtVector.m_tZ * vtVector.m_tZ;
	}
	template <typename T>
	static T squared_distance(Point_3<T> pt1, Point_3<T> pt2)
	{
		T _tX, _tY, _tZ;
		_tX = pt2.m_tX - pt1.m_tX;
		_tY = pt2.m_tY - pt1.m_tY;
		_tZ = pt2.m_tZ - pt1.m_tZ;
		return _tX * _tX + _tY * _tY + _tZ * _tZ;
	}
}
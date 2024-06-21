#include <ostream>
#define ORIGINAL Point_3(0, 0, 0)
namespace CGAL
{
	template <typename T>
	struct Coord_3
	{
		union
		{
			struct
			{
				T m_tX, m_tY, m_tZ;
			};
			struct
			{
				T m_lptCoord[3];
			};
		};
		Coord_3(T &_tX, T &_tY, T &_tZ)
		{
			m_tX = _tX;
			m_tY = _tY;
			m_tZ = _tZ;
		};
		T& x(){ return m_tX; };
		T& y(){ return m_tY; };
		T& z(){ return m_tZ; };
		bool operator==(Coord_3 coord)
		{
			return (m_tX == coord.m_tX && m_tY == coord.m_tY && m_tZ == coord.m_tZ);
		}
	};
	template <typename T>
	static std::ostream& operator<<(std::ostream &ostr, Coord_3<T> &c3)
	{
		ostr << c3.m_tX << '\t' << c3.m_tY << '\t' << c3.m_tZ << "\t\n";
		return ostr;
	}
	template <typename T>
	struct Point_3: public Coord_3<T>
	{
		Point_3(T _tX = 0, T _tY = 0, T _tZ = 0) : Coord_3<T>(_tX, _tY, _tZ){};
	};
	template <typename T>
	struct Vector_3 : public Coord_3 < T >
	{
		Vector_3(T _tX = 0, T _tY = 0, T _tZ = 0) : Coord_3<T>(_tX, _tY, _tZ){};
	};
	template <typename T>
	struct Cartesian
	{
		typedef T FT;
		typedef Point_3<T> Point_3;
		typedef Vector_3<T> Vector_3;
	};
	template <typename T>
	static Point_3<T> operator+(Point_3<T> ptPoint, Vector_3<T> vtTranslate)
	{
		return Point_3<T>(ptPoint.m_tX + vtTranslate.m_tX, ptPoint.m_tY + vtTranslate.m_tY, ptPoint.m_tZ + vtTranslate.m_tZ);
	}
	template <typename T>
	static Point_3<T> operator+(Vector_3<T> vtTranslate, Point_3<T> ptPoint)
	{
		return Point_3<T>(vtTranslate.m_tX + ptPoint.m_tX, vtTranslate.m_tY + ptPoint.m_tY, vtTranslate.m_tZ + ptPoint.m_tZ);
	}
	template <typename T>
	static Vector_3<T> operator+(Vector_3<T> vt1, Vector_3<T> vt2)
	{
		return Vector_3<T>(vt1.m_tX + vt2.m_tX, vt1.m_tY + vt2.m_tY, vt1.m_tZ + vt2.m_tZ);
	}
	template <typename T>
	static Point_3<T> operator-(Point_3<T> ptPoint, Vector_3<T> vtTranslate)
	{
		return Point_3<T>(ptPoint.m_tX - vtTranslate.m_tX, ptPoint.m_tY - vtTranslate.m_tY, ptPoint.m_tZ - vtTranslate.m_tZ);
	}
	template <typename T>
	static Vector_3<T> operator-(Point_3<T> pt1, Point_3<T> pt2)
	{
		return Vector_3<T>(pt1.m_tX - pt2.m_tX, pt1.m_tY - pt2.m_tY, pt1.m_tZ - pt2.m_tZ);
	}
	template <typename T>
	static Vector_3<T> operator-(Vector_3<T> vt1, Vector_3<T> vt2)
	{
		return Vector_3<T>(vt1.m_tX - vt2.m_tX, vt1.m_tY - vt2.m_tY, vt1.m_tZ - vt2.m_tZ);
	}
	template <typename T>
	static Vector_3<T> operator-(Vector_3<T> vtVector)
	{
		return Vector_3<T>(-vtVector.m_tX, -vtVector.m_tY, -vtVector.m_tZ);
	}
	template <typename T>
	static Vector_3<T> operator*(Vector_3<T> vtVector, T _tScale)
	{
		return Vector_3<T>(vtVector.m_tX * _tScale, vtVector.m_tY * _tScale, vtVector.m_tZ * _tScale);
	}
	template <typename T>
	static Vector_3<T> operator/(Vector_3<T> vtVector, T _tScale)
	{
		return Vector_3<T>(vtVector.m_tX / _tScale, vtVector.m_tY / _tScale, vtVector.m_tZ / _tScale);
	}
	template <typename T>
	static Vector_3<T> operator*(T _tScale, Vector_3<T> vtVector)
	{
		return Vector_3<T>(_tScale * vtVector.m_tX, _tScale * vtVector.m_tY, _tScale * vtVector.m_tZ);
	}
	template <typename T>
	static T operator*(Vector_3<T> vt1, Vector_3<T> vt2)
	{
		return vt1.m_tX * vt2.m_tX + vt1.m_tY * vt2.m_tY + vt1.m_tZ * vt2.m_tZ;
	}
};

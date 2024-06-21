#include "CTracer.h"
#include <vector>
//included by "CDistanceHeap.h"
#ifndef CPOLYHEDRONMODIFIER_H
#define CPOLYHEDRONMODIFIER_H

struct CIntegerParameter
{
    CMeshBase::FT lpftPara[2];
};

struct CSamplePoint
{
    int iIndex, iBlockId;
    CMeshBase::FT lpftPara[2];
    CMeshBase::Point_3 ptCoordinate;
    std::vector<int> m_vecAdjacent;
    ~CSamplePoint()
    {
        if(!m_vecAdjacent.empty())
        {
            m_vecAdjacent.clear();
        }
    }
};
struct CSampleBlock
{
    int iPointStartIndex, iPointEndIndex;
    CMeshBase::CPolyhedron::Halfedge_handle hhBase;
};
struct CFacetLinkage
{
    int lpiVertices[4];
};

class CSamplePolyhedronModifier//:public CGAL::Modifier_base<CMeshBase::CSamplerPolyhedron::HalfedgeDS>
{
public:
    CSamplePolyhedronModifier(std::vector<CSamplePoint> &vecSamplePoints, std::vector<CFacetLinkage> &vecFacetLinkages, CMeshBase::CSamplerPolyhedron &mshRemeshedSurface);
	void GenerateMesh();
    //void operator()(CMeshBase::CSamplerPolyhedron::HalfedgeDS &hds);
protected:
    std::vector<CSamplePoint> &m_vecSamplePoints;
    std::vector<CFacetLinkage> &m_vecFacetLinkages;
	CMeshBase::CSamplerPolyhedron &m_mshRemeshedSurface;
    int m_nNumBoundaries;
private:
};

class CTriPolyhedronModifier//:public CGAL::Modifier_base<CMeshBase::CPolyhedron::HalfedgeDS>
{
public:
	CTriPolyhedronModifier(std::vector<double> &vecCoordinates, std::vector<int> &vecIndex, CMeshBase::CPolyhedron &mshSurface);
	void GenerateMesh();
    //void operator()(CMeshBase::CPolyhedron::HalfedgeDS &hds);
protected:
    std::vector<double> &m_vecCoordinates;
    std::vector<int> &m_vecIndex;
	CMeshBase::CPolyhedron &m_mshSurface;
//    int m_nNumHalfEdges;
};

#endif // CPOLYHEDRONMODIFIER_H

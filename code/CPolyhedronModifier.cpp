#include "CPolyhedronModifier.h"

CSamplePolyhedronModifier::CSamplePolyhedronModifier(std::vector<CSamplePoint> &vecSamplePoints, std::vector<CFacetLinkage> &vecFacetLinkages, CMeshBase::CSamplerPolyhedron &mshRemeshedSurface) : m_vecSamplePoints(vecSamplePoints), m_vecFacetLinkages(vecFacetLinkages), m_mshRemeshedSurface(mshRemeshedSurface){}

void CSamplePolyhedronModifier::GenerateMesh()
{
	std::vector<CSamplePoint>::iterator iterPoint;
	std::vector<CFacetLinkage>::iterator iterFacet;
	for (iterPoint = m_vecSamplePoints.begin(); iterPoint != m_vecSamplePoints.end(); ++iterPoint)
	{
		m_mshRemeshedSurface.append_vertex(iterPoint->ptCoordinate.x(), iterPoint->ptCoordinate.y(), iterPoint->ptCoordinate.z());
	}
	m_mshRemeshedSurface.begin_facets();
	for (iterFacet = m_vecFacetLinkages.begin(); iterFacet != m_vecFacetLinkages.end(); ++iterFacet)
	{
		m_mshRemeshedSurface.append_facet(4, iterFacet->lpiVertices[0], iterFacet->lpiVertices[1], iterFacet->lpiVertices[2], iterFacet->lpiVertices[3]);
	}
	m_mshRemeshedSurface.end_facets();
}

CTriPolyhedronModifier::CTriPolyhedronModifier(std::vector<double> &vecCoordinates, std::vector<int> &vecIndex, CMeshBase::CPolyhedron &mshSurface): m_vecCoordinates(vecCoordinates), m_vecIndex(vecIndex), m_mshSurface(mshSurface){}

void CTriPolyhedronModifier::GenerateMesh()
{
	std::vector<double>::iterator iterCoord;
	std::vector<int>::iterator iterIndex;
	double dblX, dblY, dblZ;
	int iV1, iV2, iV3;
	iterCoord = m_vecCoordinates.begin();
	while (iterCoord != m_vecCoordinates.end())
	{
		dblX = *iterCoord;
		++iterCoord;
		dblY = *iterCoord;
		++iterCoord;
		dblZ = *iterCoord;
		++iterCoord;
		m_mshSurface.append_vertex(dblX, dblY, dblZ);
	}
	m_mshSurface.begin_facets();
	iterIndex = m_vecIndex.begin();
	while (iterIndex != m_vecIndex.end())
	{
		iV1 = *iterIndex;
		++iterIndex;
		iV2 = *iterIndex;
		++iterIndex;
		iV3 = *iterIndex;
		++iterIndex;
		m_mshSurface.append_facet(3, iV1, iV2, iV3);
	}
	m_mshSurface.end_facets();
}

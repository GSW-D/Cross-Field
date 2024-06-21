#include "CQuadMeshEvaluator.h"
CQuadMeshEvaluator::CQuadMeshEvaluator(CMeshBase::CSamplerPolyhedron &mshSurface): m_mshSurface(mshSurface)
{
	m_lpAspectRatioDistribution = NULL;
	m_lpAreaDistribution = NULL;
}

CQuadMeshEvaluator::~CQuadMeshEvaluator()
{
	if (m_lpAspectRatioDistribution != NULL)
	{
		delete[]m_lpAspectRatioDistribution;
	}
	if (m_lpAreaDistribution != NULL)
	{
		delete[]m_lpAreaDistribution;
	}
}


void CQuadMeshEvaluator::m_fnJumpSampling()
{
	std::queue<CMeshBase::CSamplerPolyhedron::Vertex_handle> quvhQueue;
	CMeshBase::CSamplerPolyhedron::Vertex_iterator vi;
	CMeshBase::CSamplerPolyhedron::Vertex_handle vhSeed, vhNew;;
	CMeshBase::CSamplerPolyhedron::Halfedge_handle hh0, hh;
	for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
	{
		if (vi->iMark == 1)
		{
			quvhQueue.push(&(*vi));
		}
	}
	while (!quvhQueue.empty())
	{
		vhSeed = quvhQueue.front();
		quvhQueue.pop();
		hh0 = vhSeed->halfedge()->opposite();
		hh = hh0;
		do
		{
			vhNew = NULL;
			if (hh->vertex()->iBorderId == -1 && hh->vertex()->degree() == 4)
			{
				vhNew = hh->next()->opposite()->next()->vertex();
				if (vhNew->iMark == 0)
				{
					vhNew->iMark = 1;
					quvhQueue.push(vhNew);
				}
			}
			hh = hh->opposite()->next();
		} while (hh != hh0);
	}
}

void CQuadMeshEvaluator::m_fnCountAngle()
{
	CMeshBase::CSamplerPolyhedron::Facet_iterator fi;
	CMeshBase::CSamplerPolyhedron::Halfedge_handle hh0, hh;
	CMeshBase::FT ftDot, ftCross;
	int i, n;
	for (i = 0; i < 18; ++i)
	{
		m_lpdblAngleDistribution[i] = 0;
	}
	n = 0;
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
	{
		hh0 = fi->halfedge();
		hh = hh0;
		do
		{
			hh->vtVector = hh->vertex()->point() - hh->prev()->vertex()->point();
			hh = hh->next();
		} while (hh != hh0);
		do
		{
			ftDot = hh->vtVector * hh->next()->vtVector;
			ftCross = sqrt(CGAL::squared_length(CGAL::cross_product(hh->vtVector, hh->next()->vtVector)));
			hh->ftAngle = atan2(ftCross, -ftDot) * 180.0 / ftPi;
			i = int((hh->ftAngle - 45) / 5.0);
			if (i < 0)
			{
				i = 0;
			}
			if (i > 17)
			{
				i = 17;
			}
			m_lpdblAngleDistribution[i] += 1.0;
			++n;
			hh = hh->next();
		} while (hh != hh0);
	}
	for (i = 0; i < 18; ++i)
	{
		m_lpdblAngleDistribution[i] /= n;
	}
}

void CQuadMeshEvaluator::m_fnCountAspectRatio(double dblUnitWidth, int nGropus)
{
	CMeshBase::CSamplerPolyhedron::Facet_iterator fi;
	CMeshBase::CSamplerPolyhedron::Halfedge_handle hh0, hh;
	int i, n;
	if (m_lpAspectRatioDistribution != NULL)
	{
		delete[]m_lpAspectRatioDistribution;
	}
	m_lpAspectRatioDistribution = new double[nGropus];
	for (i = 0; i < nGropus; ++i)
	{
		m_lpAspectRatioDistribution[i] = 0;
	}
	n = 0;
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
	{
		hh0 = fi->halfedge();
		hh = hh0;
		do
		{
			hh->vtVector = hh->vertex()->point() - hh->prev()->vertex()->point();
			hh->ftSqLen = CGAL::squared_length(hh->vtVector);
			hh->ftLen = sqrt(hh->ftSqLen);
			hh = hh->next();
		} while (hh != hh0);
		do
		{
			hh->ftLogLenRatio = log(hh->next()->ftLen / hh->ftLen);
			i = round(hh->ftLogLenRatio / dblUnitWidth + (nGropus - 1) * 0.5);
			if (i < 0)
			{
				i = 0;
			}
			if (i >= nGropus)
			{
				i = nGropus - 1;
			}
			m_lpAspectRatioDistribution[i] += 1.0;
			++n;
			hh = hh->next();
		} while (hh != hh0);
	}
	for (i = 0; i < nGropus; ++i)
	{
		m_lpAspectRatioDistribution[i] /= n;
	}
}

void CQuadMeshEvaluator::m_fnCountFacetArea(double dblUnitWidth, int nGropus)
{
	CMeshBase::CSamplerPolyhedron::Facet_iterator fi;
	CMeshBase::CSamplerPolyhedron::Halfedge_handle hh0, hh;
	CMeshBase::Vector_3 vtNorm;
	CMeshBase::FT ftAreaAve, ftSumArea;
	int i, n;
	if (m_lpAreaDistribution != NULL)
	{
		delete[]m_lpAreaDistribution;
	}
	m_lpAreaDistribution = new double[nGropus];
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
	{
		vtNorm = CMeshBase::Vector_3(0, 0, 0);
		hh0 = fi->halfedge();
		hh = hh0;
		do
		{
			vtNorm = vtNorm + CGAL::cross_product(hh->vtVector, hh->next()->vtVector);
			hh = hh->next()->next();
		} while (hh != hh0);
		fi->ftArea = sqrt(CGAL::squared_length(vtNorm)) / 2;
		fi->ftLogArea = log(CGAL::squared_length(vtNorm) / 4) / 2;
	}
	n = 0;
	ftAreaAve = 0;
	ftSumArea = 0;
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
	{
		ftAreaAve += fi->ftLogArea;
		ftSumArea += fi->ftArea;
		++n;
	}
	ftAreaAve /= n;
	ftSumArea /= n;
	for (i = 0; i < nGropus; ++i)
	{
		m_lpAreaDistribution[i] = 0;
	}
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
	{
		i = round((fi->ftLogArea - ftAreaAve) / dblUnitWidth + (nGropus - 1) * 0.5);
		if (i < 0)
		{
			i = 0;
		}
		if (i >= nGropus)
		{
			i = nGropus - 1;
		}
		m_lpAreaDistribution[i] += 1.0;
	}
	for (i = 0; i < nGropus; ++i)
	{
		m_lpAreaDistribution[i] /= n;
	}
	m_ftAreaDistortion = 0;
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
	{
		m_ftAreaDistortion += (fi->ftArea * fi->ftArea / ftSumArea + ftSumArea) / 2;
	}
	m_ftAreaDistortion /= (ftSumArea * n);
}

double CQuadMeshEvaluator::m_fnCountAngularDistortion()
{
	CMeshBase::CSamplerPolyhedron::Facet_iterator fi;
	CMeshBase::CSamplerPolyhedron::Halfedge_handle hh0, hh;
	CMeshBase::Vector_3 vtNorm, vtCross;
	CMeshBase::FT ftSumArea, ftSumDistortion, ftLocalDistortion, ftAveDistortion, ftLocalArea;
	int nEdges;
	ftSumDistortion = 0;
	ftSumArea = 0;
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
	{
		vtNorm = CMeshBase::Vector_3(0, 0, 0);
		ftAveDistortion = 0;
		nEdges = 0;
		hh0 = fi->halfedge();
		hh = hh0;
		do
		{
			vtCross = CGAL::cross_product(hh->vtVector, hh->next()->vtVector);
			ftLocalDistortion = hh->vtVector * hh->vtVector + hh->next()->vtVector * hh->next()->vtVector;
			ftLocalDistortion = ftLocalDistortion * ftLocalDistortion / (vtCross * vtCross) - 2;
			ftAveDistortion = ftAveDistortion + ftLocalDistortion;
			vtNorm = vtNorm + vtCross;
			++nEdges;
			hh = hh->next();
		} while (hh != hh0);
		ftAveDistortion /= nEdges;
		vtNorm = vtNorm / CMeshBase::FT(nEdges);
		ftLocalArea = sqrt(vtNorm * vtNorm);
		//if (ftAveDistortion > 0 && ftAveDistortion < 10)
		{
			ftSumDistortion += ftAveDistortion * ftLocalArea;
			ftSumArea += ftLocalArea;
		}
	}
	ftAveDistortion = ftSumDistortion / (2 * ftSumArea);
	return ftAveDistortion;
}

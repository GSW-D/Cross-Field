#include "CMeshInfoInitializer.h"


void CMeshInfoInitializer::m_fnNormalizeMesh()
{
	CMeshBase::FT ftAveX, ftAveY, ftAveZ, ftAveRadius;
	CMeshBase::Point_3 ptCenter;
	CMeshBase::CPolyhedron::Vertex_iterator vi;
	ftAveX = ftAveY = ftAveZ = 0;
	for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
	{
		ftAveX += vi->point().x();
		ftAveY += vi->point().y();
		ftAveZ += vi->point().z();
	}
	ftAveX /= m_mshSurface.size_of_vertices();
	ftAveY /= m_mshSurface.size_of_vertices();
	ftAveZ /= m_mshSurface.size_of_vertices();
	ptCenter = CMeshBase::Point_3(ftAveX, ftAveY, ftAveZ);
	ftAveRadius = 0;
	for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
	{
		ftAveRadius += CGAL::squared_distance(vi->point(), ptCenter);
	}
	ftAveRadius = sqrt(ftAveRadius / m_mshSurface.size_of_vertices());
	for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
	{
		vi->point() = CMeshBase::Point_3((vi->point().x() - ftAveX) / ftAveRadius, (vi->point().y() - ftAveY) / ftAveRadius, (vi->point().z() - ftAveZ) / ftAveRadius);
	}
}

void CMeshInfoInitializer::m_fnCountNewMeshSize()
{
	CMeshBase::CSamplerPolyhedron::Edge_iterator ei;
	CMeshBase::Vector_3 vtVector;
	m_ftNewMeanRadius = 0;
	for (ei = m_mshRemeshedSurface.edges_begin(); ei != m_mshRemeshedSurface.edges_end(); ++ei)
	{
		vtVector = ei->vertex()->point() - ei->prev()->vertex()->point();
		m_ftNewMeanRadius += sqrt(vtVector * vtVector);
	}
	m_ftNewMeanRadius /= (m_mshRemeshedSurface.size_of_edges() * 4);
}

void CMeshInfoInitializer::m_fnMarkNewMeshBorder()
{
	CMeshBase::CSamplerPolyhedron::Vertex_iterator vi;
	CMeshBase::CSamplerPolyhedron::Halfedge_iterator hi;
	CMeshBase::CSamplerPolyhedron::Halfedge_handle hh0, hh;
	unsigned int iIndex;
	for (vi = m_mshRemeshedSurface.vertices_begin(); vi != m_mshRemeshedSurface.vertices_end(); ++vi)
	{
		vi->iBorderId = -1;
		vi->iMark = 0;
	}
	iIndex = 0;
	for (hi = m_mshRemeshedSurface.halfedges_begin(); hi != m_mshRemeshedSurface.halfedges_end(); ++hi)
	{
		if (hi->is_border() && hi->vertex()->iBorderId == -1)
		{
			hh0 = &(*hi);
			hh = hh0;
			do
			{
				hh->vertex()->iBorderId = iIndex;
				hh = hh->next();
			} while (hh != hh0);
			++iIndex;
		}
	}
}

void CMeshInfoInitializer::m_fnMarkNewMeshCrevasse()
{
	CMeshBase::CSamplerPolyhedron::Edge_iterator ei;
	CMeshBase::CSamplerPolyhedron::Vertex_iterator vi;
	CMeshBase::CSamplerPolyhedron::Halfedge_handle hh0, hh1, hh;
	for (ei = m_mshRemeshedSurface.edges_begin(); ei != m_mshRemeshedSurface.edges_end(); ++ei)
	{
		ei->bCrevasse = false;
		ei->opposite()->bCrevasse = false;
	}
	for (vi = m_mshRemeshedSurface.vertices_begin(); vi != m_mshRemeshedSurface.vertices_end(); ++vi)
	{
		if (vi->iBorderId == -1 && vi->degree() != 4)
		{
			hh0 = vi->halfedge();
			hh1 = hh0;
			do
			{
				hh = hh1;
				do
				{
					hh->bCrevasse = true;
					hh->opposite()->bCrevasse = true;
					hh = hh->prev()->opposite()->prev();
				} while (!hh->bCrevasse && hh->vertex()->degree() == 4 && hh->vertex()->iBorderId == -1);
				hh1 = hh1->next()->opposite();
			} while (hh1 != hh0);
		}
	}
}

void CMeshInfoInitializer::m_fnInitVector()
{
	CMeshBase::CPolyhedron::Halfedge_iterator hi;
	for (hi = m_mshSurface.halfedges_begin(); hi != m_mshSurface.halfedges_end(); ++hi)
	{
		hi->vtVector = hi->vertex()->point() - hi->prev()->vertex()->point();
	}
}

void CMeshInfoInitializer::m_fnInitDot()
{
	CMeshBase::CPolyhedron::Halfedge_iterator hi;
	for (hi = m_mshSurface.halfedges_begin(); hi != m_mshSurface.halfedges_end(); ++hi)
	{
		hi->ftDot = hi->vtVector * hi->next()->vtVector;
	}
}

void CMeshInfoInitializer::m_fnInitFacetNormal()
{
	CMeshBase::CPolyhedron::Facet_iterator fi;
	m_ftGeometricDimension = 0;
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
	{
		fi->vtNorm = CGAL::cross_product(fi->halfedge()->vtVector, fi->halfedge()->next()->vtVector);
		fi->ftTwArea = sqrt(fi->vtNorm * fi->vtNorm);
		m_ftGeometricDimension += fi->ftTwArea;
	}
	m_ftGeometricDimension = sqrt(m_ftGeometricDimension / 2.0);
}


void CMeshInfoInitializer::m_fnInitLength()
{
	CMeshBase::CPolyhedron::Edge_iterator ei;
	for (ei = m_mshSurface.edges_begin(); ei != m_mshSurface.edges_end(); ++ei)
	{
		ei->ftSqLen = ei->vtVector * ei->vtVector;
		ei->ftLen = sqrt(ei->ftSqLen);
		ei->opposite()->ftSqLen = ei->ftSqLen;
		ei->opposite()->ftLen = ei->ftLen;
	}
}

void CMeshInfoInitializer::m_fnInitAngle()
{
	CMeshBase::CPolyhedron::Halfedge_iterator hi;
	CMeshBase::Vector_3 vtNorm;
	for (hi = m_mshSurface.halfedges_begin(); hi != m_mshSurface.halfedges_end(); ++hi)
	{
		if (hi->is_border())
		{
			vtNorm = CGAL::cross_product(hi->vtVector, hi->next()->vtVector);
			hi->ftAngle = atan2(sqrt(vtNorm * vtNorm), hi->ftDot);
		}
		else
		{
			if (hi->prev()->ftLen == 0.0)
			{
				if (hi->ftLen == 0.0)
				{
					hi->ftAngle = ftPi * 2.0 / 3.0;
				}
				else
				{
					hi->ftAngle = ftPi;
				}
			}
			else
			{
				hi->ftAngle = atan2(hi->facet()->ftTwArea, hi->ftDot);
			}
		}
	}
}



void CMeshInfoInitializer::m_fnInitHELocalInd()
{
	CMeshBase::CPolyhedron::Halfedge_iterator hi;
	CMeshBase::CPolyhedron::Facet_iterator fi;
	CMeshBase::CPolyhedron::Halfedge_handle hh0, hh;
	int iLocalInd;
	for (hi = m_mshSurface.halfedges_begin(); hi != m_mshSurface.halfedges_end(); ++hi)
	{
		hi->iLocalInd = 0;
	}
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
	{
		iLocalInd = 0;
		hh0 = fi->halfedge();
		hh = hh0;
		do
		{
			hh->iLocalInd = iLocalInd;
			++iLocalInd;
			hh = hh->next();
		} while (hh != hh0);
	}
}

void CMeshInfoInitializer::m_fnMarkBorder()
{
	CMeshBase::CPolyhedron::Vertex_iterator vi;
	CMeshBase::CPolyhedron::Halfedge_iterator hi;
	CMeshBase::CPolyhedron::Facet_iterator fi;
	CMeshBase::CPolyhedron::Halfedge_handle hh;
	int iBordId;
	for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
	{
		vi->iBordId = -1;
	}
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
	{
		fi->bDirFixed = false;
	}
	iBordId = 0;
	for (hi = m_mshSurface.halfedges_begin(); hi != m_mshSurface.halfedges_end(); ++hi)
	{
		if (hi->is_border() && hi->vertex()->iBordId == -1)
		{
			hh = &(*hi);
			do
			{
				hh->opposite()->facet()->bDirFixed = true;//bBordFacet = true;
				hh->vertex()->iBordId = iBordId;
				hh = hh->next();
			} while (hh != &(*hi));
			++iBordId;
		}
	}
}

void CMeshInfoInitializer::m_fnInitAngleDefect()
{
	CMeshBase::CPolyhedron::Vertex_iterator vi;
	CMeshBase::CPolyhedron::Halfedge_handle hh0, hh;
	int iInnerAdjacent;
//	int nIdealFacets;
	for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
	{
		if (vi->iBordId == -1)
		{
			vi->ftAngleDefect = ftTwoPi;
			hh0 = vi->halfedge();
			hh = hh0;
			do
			{
				vi->ftAngleDefect -= (ftPi - hh->ftAngle);
				hh = hh->next()->opposite();
			} while (hh != hh0);
			vi->ftGaussCurvature = vi->ftAngleDefect / vi->ftVoronoiArea;
		}
		else
		{
			vi->ftAngleDefect = ftPi;
			hh0 = vi->halfedge();
			hh = hh0;
			do
			{
				if (!hh->is_border())
				{
					vi->ftAngleDefect -= (ftPi - hh->ftAngle);
				}
				hh = hh->next()->opposite();
			} while (hh != hh0);
			vi->ftGaussCurvature = 0;
		}
	}
	for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
	{
		if (vi->iBordId != -1)
		{
			iInnerAdjacent = 0;
			hh0 = vi->halfedge()->opposite();
			hh = hh0;
			do
			{
				if (hh->vertex()->iBordId == -1)
				{
					vi->ftGaussCurvature += hh->vertex()->ftGaussCurvature;
					++iInnerAdjacent;
				}
				hh = hh->prev()->opposite();
			} while (hh != hh0);
			vi->ftGaussCurvature /= iInnerAdjacent;
		}
	}
}

void CMeshInfoInitializer::m_fnInitLaplaceCoef()
{
	CMeshBase::CPolyhedron::Edge_iterator ei;
	CMeshBase::CPolyhedron::Halfedge_handle hh0, hh;
	CMeshBase::CPolyhedron::Vertex_iterator vi;
	for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
	{
		vi->ftVoronoiArea = 0;
	}
	for (ei = m_mshSurface.edges_begin(); ei != m_mshSurface.edges_end(); ++ei)
	{
		ei->ftLaplaceCoef = 0;
		ei->opposite()->ftLaplaceCoef = 0;
		hh0 = &(*ei);
		hh = hh0;
		do
		{
			if (!hh->is_border())
			{
				ei->ftLaplaceCoef -= hh->next()->ftDot / hh->facet()->ftTwArea;
			}
			hh = hh->opposite();
		} while (hh != hh0);
		if (ei->ftLaplaceCoef < 1.0/65536.0)
		{
			ei->ftLaplaceCoef = 1.0/65536.0;
		}
		ei->opposite()->ftLaplaceCoef = ei->ftLaplaceCoef;
		ei->vertex()->ftVoronoiArea += ei->ftSqLen * ei->ftLaplaceCoef;
		ei->prev()->vertex()->ftVoronoiArea += ei->ftSqLen * ei->ftLaplaceCoef;
	}
}


void CMeshInfoInitializer::m_fnInscribedCircle()
{
	CMeshBase::CPolyhedron::Facet_iterator fi;
	CMeshBase::CPolyhedron::Halfedge_handle hh0, hh;
	CMeshBase::FT ftPerimeter, ftSumX, ftSumY, ftSumZ;
	m_ftMeanRadius = 0;
	m_ftCurlFieldDimension = 0;
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
	{
		ftPerimeter = 0;
		ftSumX = 0;
		ftSumY = 0;
		ftSumZ = 0;
		hh0 = fi->halfedge();
		hh = hh0;
		do
		{
			ftSumX += hh->vertex()->point().x() * hh->prev()->ftLen;
			ftSumY += hh->vertex()->point().y() * hh->prev()->ftLen;
			ftSumZ += hh->vertex()->point().z() * hh->prev()->ftLen;
			ftPerimeter += hh->ftLen;
			hh = hh->next();
		} while (hh != hh0);
		ftSumX /= ftPerimeter;
		ftSumY /= ftPerimeter;
		ftSumZ /= ftPerimeter;
		fi->ptIncent = CMeshBase::Point_3(ftSumX, ftSumY, ftSumZ);
		fi->ftRadius = fi->ftTwArea / ftPerimeter;
		m_ftMeanRadius += fi->ftRadius;
		if (fi->halfedge()->ftLen == 0)
		{
			fi->vtAuxiliaryAxis = hh->next()->vtVector / hh->next()->ftLen;
			fi->vtPrincipalAxis = CGAL::cross_product(fi->vtAuxiliaryAxis, fi->vtNorm) / fi->ftTwArea;
		}
		else
		{
			fi->vtPrincipalAxis = fi->halfedge()->vtVector * (fi->ftRadius / fi->halfedge()->ftLen);
			fi->vtAuxiliaryAxis = CGAL::cross_product(fi->vtNorm, fi->vtPrincipalAxis) / fi->ftTwArea;
		}
		m_ftCurlFieldDimension += fi->ftTwArea;
	}
	m_ftCurlFieldDimension = sqrt(1 / m_ftCurlFieldDimension);
	m_ftMeanRadius /= m_mshSurface.size_of_facets();
}

void CMeshInfoInitializer::m_fnInitPolarAxisDif()
{
	CMeshBase::CPolyhedron::Edge_iterator ei;
	CMeshBase::CPolyhedron::Halfedge_handle hh;
	for (ei = m_mshSurface.edges_begin(); ei != m_mshSurface.edges_end(); ++ei)
	{
		if (ei->is_border_edge())
		{
			hh = &(*ei);
			if (hh->is_border())
			{
				hh = hh->opposite();
			}
			hh->ftPolarAxisDif = ftPi;
			if (hh->iLocalInd == 1)
			{
				hh->ftPolarAxisDif += hh->prev()->ftAngle;
			}
			else if (hh->iLocalInd == 2)
			{
				hh->ftPolarAxisDif -= hh->ftAngle;
			}
			hh->ftPolarAxisDif -= CMeshBase::FT(int(ei->ftPolarAxisDif / ftTwoPi + 15.5) - 15) * ftTwoPi;
			hh->opposite()->ftPolarAxisDif = ftTwoPi - hh->ftPolarAxisDif;
			hh->opposite()->ftPolarAxisDif -= CMeshBase::FT(int(ei->opposite()->ftPolarAxisDif / ftTwoPi + 15.5) - 15) * ftTwoPi;

		}
		else
		{
			ei->ftPolarAxisDif = ftPi;
			if (ei->opposite()->iLocalInd == 1)
			{
				ei->ftPolarAxisDif -= ei->opposite()->prev()->ftAngle;
			}
			else if (ei->opposite()->iLocalInd == 2)
			{
				ei->ftPolarAxisDif += ei->opposite()->ftAngle;
			}
			if (ei->iLocalInd == 1)
			{
				ei->ftPolarAxisDif += ei->prev()->ftAngle;
			}
			else if (ei->iLocalInd == 2)
			{
				ei->ftPolarAxisDif -= ei->ftAngle;
			}
			ei->ftPolarAxisDif -= CMeshBase::FT(int(ei->ftPolarAxisDif / ftTwoPi + 15.5) - 15) * ftTwoPi;
			ei->opposite()->ftPolarAxisDif = ftTwoPi - ei->ftPolarAxisDif;
			ei->opposite()->ftPolarAxisDif -= CMeshBase::FT(int(ei->opposite()->ftPolarAxisDif / ftTwoPi + 15.5) - 15) * ftTwoPi;
		}
	}
}


void CMeshInfoInitializer::m_fnFacetPrincipleDirection(CMeshBase::CPolyhedron::Facet_handle fh)
{
	CMeshBase::CPolyhedron::Halfedge_handle hh0, hh, hh1, hh2;
	int i, j, k, n;
	CMeshBase::FT ftMinCurDif = 4.0, ftMaxErr = 4.0;
	CMeshBase::FT lpftWeightedCurvature[3];
	std::vector<CMeshBase::FT> lpvecftQuadraticItem[3], vecftZ;
	std::vector<CMeshBase::FT>::iterator iterftI, iterftJ, iterftSqX, iterftXY, iterftSqY, iterftZ;
	CMeshBase::Vector_3 vtSample, lpvtRow[3], lpvtCross[3], vtProjection, ftDot;
	CMeshBase::FT ftXCur, ftYCur, ftZCur, ftVolume, ftErr, ftSumSqErr, ftDirU, ftDirV, ftTestScale, ftTestErr;
	CMeshBase::FT lpftSymmetryMatrix[9], lpftQuadraticCoef[3], lpftTestCoef[3];
	CMeshBase::FT ftEquationCoef0, ftEquationCoef1, ftSquareRootedDiscriminant, ftEigenValue0, ftEigenValue1;
	hh0 = fh->halfedge();
	hh = hh0;
	n = 0;
	do
	{
		hh1 = hh->next()->opposite();
		hh2 = hh->opposite()->prev();
		while (hh1 != hh2)
		{
			vtSample = hh1->facet()->ptIncent - fh->ptIncent;
			ftXCur = vtSample * fh->vtPrincipalAxis / fh->ftRadius;
			ftYCur = vtSample * fh->vtAuxiliaryAxis / fh->ftRadius;
			ftZCur = vtSample * fh->vtNorm / fh->ftTwArea;
			lpvecftQuadraticItem[0].push_back((ftXCur * ftXCur));
			lpvecftQuadraticItem[1].push_back((ftXCur * ftYCur));
			lpvecftQuadraticItem[2].push_back((ftYCur * ftYCur));
			vecftZ.push_back(ftZCur);
			++n;
			hh1 = hh1->next()->opposite();
		}
		hh = hh->next();
	} while (hh != hh0);
	for (i = 0; i < 3; ++i)
	{
		for (j = i; j < 3; ++j)
		{
			iterftI = lpvecftQuadraticItem[i].begin();
			iterftJ = lpvecftQuadraticItem[j].begin();
			lpftSymmetryMatrix[i * 3 + j] = 0;
			for (k = 0; k < n; ++k)
			{
				lpftSymmetryMatrix[i * 3 + j] += (*iterftI) * (*iterftJ);
				++iterftI;
				++iterftJ;
			}
		}
		iterftI = lpvecftQuadraticItem[i].begin();
		iterftJ = vecftZ.begin();
		lpftWeightedCurvature[i] = 0;
		for (k = 0; k < n; ++k)
		{
			lpftWeightedCurvature[i] += (*iterftI) * (*iterftJ);
			++iterftI;
			++iterftJ;
		}
	}
	lpftSymmetryMatrix[3] = lpftSymmetryMatrix[1];
	lpftSymmetryMatrix[6] = lpftSymmetryMatrix[2];
	lpftSymmetryMatrix[7] = lpftSymmetryMatrix[5];
	for (i = 0; i < 3; ++i)
	{
		lpvtRow[i] = CMeshBase::Vector_3(lpftSymmetryMatrix[i * 3], lpftSymmetryMatrix[i * 3 + 1], lpftSymmetryMatrix[i * 3 + 2]);
	}
	for (i = 0; i < 3; ++i)
	{
		j = (i + 1) % 3;
		k = (i + 2) % 3;
		lpvtCross[i] = CGAL::cross_product(lpvtRow[j], lpvtRow[k]);
	}
	ftVolume = lpvtCross[0] * lpvtRow[0];
	vtProjection = CMeshBase::Vector_3(lpftWeightedCurvature[0], lpftWeightedCurvature[1], lpftWeightedCurvature[2]) / ftVolume;
	for (i = 0; i < 3; ++i)
	{
		lpftQuadraticCoef[i] = lpvtCross[i] * vtProjection;
	}
	iterftSqX = lpvecftQuadraticItem[0].begin();
	iterftXY = lpvecftQuadraticItem[1].begin();
	iterftSqY = lpvecftQuadraticItem[2].begin();
	iterftZ = vecftZ.begin();
	ftSumSqErr = 0;
	for (k = 0; k < n; ++k)
	{
		ftErr = (*iterftZ) - (*iterftSqX) * lpftQuadraticCoef[0] - (*iterftXY) * lpftQuadraticCoef[1] - (*iterftSqY) * lpftQuadraticCoef[2];
		ftSumSqErr += ftErr * ftErr;
		++iterftSqX;
		++iterftXY;
		++iterftSqY;
		++iterftZ;
	}
	lpftQuadraticCoef[1] /= 2;
	ftEquationCoef0 = lpftQuadraticCoef[0] * lpftQuadraticCoef[2] - lpftQuadraticCoef[1] * lpftQuadraticCoef[1];
	ftEquationCoef1 = (lpftQuadraticCoef[0] + lpftQuadraticCoef[2]) / 2;
	ftSquareRootedDiscriminant = sqrt(ftEquationCoef1 * ftEquationCoef1 - ftEquationCoef0);
	ftEigenValue0 = ftEquationCoef1 + ftSquareRootedDiscriminant;
	ftEigenValue1 = ftEquationCoef1 - ftSquareRootedDiscriminant;
	if (fabs(ftEigenValue0) > fabs(ftEigenValue1) * ftMinCurDif)
	{
		if (lpftQuadraticCoef[0] > 0)
		{
			ftDirU = ftEigenValue0 - (lpftQuadraticCoef[0] - lpftQuadraticCoef[1]);
			ftDirV = -ftEigenValue0 - (lpftQuadraticCoef[1] - lpftQuadraticCoef[2]);
		}
		else
		{
			ftDirU = ftEigenValue0 - (lpftQuadraticCoef[0] + lpftQuadraticCoef[1]);
			ftDirV = ftEigenValue0 - (lpftQuadraticCoef[1] + lpftQuadraticCoef[2]);
		}
		ftTestScale = ftEigenValue1 / (ftDirU * ftDirU + ftDirV * ftDirV);
		lpftTestCoef[0] = ftDirU * ftDirU * ftTestScale;
		lpftTestCoef[1] = 2 * ftDirU * ftDirV * ftTestScale;
		lpftTestCoef[2] = ftDirV * ftDirV * ftTestScale;
		ftTestErr = 0;
		iterftSqX = lpvecftQuadraticItem[0].begin();
		iterftXY = lpvecftQuadraticItem[1].begin();
		iterftSqY = lpvecftQuadraticItem[2].begin();
		iterftZ = vecftZ.begin();
		for (k = 0; k < n; ++k)
		{
			ftErr = *iterftZ - (*iterftSqX) * (lpftTestCoef[0]) - (*iterftXY) * (lpftTestCoef[1]) - (*iterftSqY) * (lpftTestCoef[2]);
			ftTestErr += ftErr * ftErr;
			++iterftSqX;
			++iterftXY;
			++iterftSqY;
			++iterftZ;
		}
		if (ftSumSqErr * ftMaxErr < ftTestErr)
		{
			fh->ftCurvatureChart = atan2(ftDirV, ftDirU);
			fh->iDistOrd = 0;
		}
	}
	else if (fabs(ftEigenValue1) > fabs(ftEigenValue0) * ftMinCurDif)
	{
		if (lpftQuadraticCoef[1] < 0)
		{
			ftDirU = ftEigenValue1 - (lpftQuadraticCoef[0] - lpftQuadraticCoef[1]);
			ftDirV = -ftEigenValue1 - (lpftQuadraticCoef[1] - lpftQuadraticCoef[2]);
		}
		else
		{
			ftDirU = ftEigenValue1 - (lpftQuadraticCoef[0] + lpftQuadraticCoef[1]);
			ftDirV = ftEigenValue1 - (lpftQuadraticCoef[1] + lpftQuadraticCoef[2]);
		}
		ftTestScale = ftEigenValue0 / (ftDirU * ftDirU + ftDirV * ftDirV);
		lpftTestCoef[0] = ftDirU * ftDirU * ftTestScale;
		lpftTestCoef[1] = 2 * ftDirU * ftDirV * ftTestScale;
		lpftTestCoef[2] = ftDirV * ftDirV * ftTestScale;
		ftTestErr = 0;
		iterftSqX = lpvecftQuadraticItem[0].begin();
		iterftXY = lpvecftQuadraticItem[1].begin();
		iterftSqY = lpvecftQuadraticItem[2].begin();
		iterftZ = vecftZ.begin();
		for (k = 0; k < n; ++k)
		{
			ftErr = *iterftZ - (*iterftSqX) * (lpftTestCoef[0]) - (*iterftXY) * (lpftTestCoef[1]) - (*iterftSqY) * (lpftTestCoef[2]);
			ftTestErr += ftErr * ftErr;
			++iterftSqX;
			++iterftXY;
			++iterftSqY;
			++iterftZ;
		}
		if (ftSumSqErr * ftMaxErr < ftTestErr)
		{
			fh->ftCurvatureChart = atan2(ftDirV, ftDirU);
			fh->iDistOrd = 0;
		}
	}
	for (i = 0; i < 3; ++i)
	{
		lpvecftQuadraticItem[i].clear();
	}
	vecftZ.clear();
}


void CMeshInfoInitializer::m_fnInitMainDirection()
{
	CMeshBase::CPolyhedron::Facet_iterator fi;
	CMeshBase::CPolyhedron::Halfedge_handle hh0, hh;
	bool bNotBorder;
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
	{
		bNotBorder = true;
		hh0 = fi->halfedge();
		hh = hh0;
		do
		{
			if (hh->opposite()->is_border())
			{
				bNotBorder = false;
			}
			hh = hh->next();
		} while (hh != hh0);
		if (bNotBorder)
		{
			//m_fnFacetPrincipleDirection(&(*fi));
			fi->ftChartDir = fi->ftCurvatureChart;
		}
	}
}


void CMeshInfoInitializer::m_fnInitOtherDefaultInfo()
{
	CMeshBase::CPolyhedron::Vertex_iterator vi;
	CMeshBase::CPolyhedron::Edge_iterator ei;
	CMeshBase::CPolyhedron::Facet_iterator fi;
	int iIndex;
	iIndex = 0;
	for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
	{
		vi->ftDist = DBL_MAX;
		vi->iDegreeBias = 0;
		vi->nNumPorts = 0;
		vi->iHeapPos = -1;
		vi->bDistConf = false;
		vi->iIndex = iIndex;
		vi->hhCome = NULL;
		vi->iSelected = 0;
		++iIndex;
	}
	for (ei = m_mshSurface.edges_begin(); ei != m_mshSurface.edges_end(); ++ei)
	{
		ei->bFreeEdge = false;
		ei->bInStack = false;
		ei->hhCrevassePrev = NULL;
		ei->hhCrevasseNext = NULL;
		if (!m_bProtectGlobalPara)
		{
			ei->lpftGlobalPara = NULL;
		}
		ei->bSelected = false;
		ei->ftForceBarrier = 0;
		ei->bSharp = false;
		ei->opposite()->bFreeEdge = false;
		ei->opposite()->bInStack = false;
		ei->opposite()->hhCrevassePrev = NULL;
		ei->opposite()->hhCrevasseNext = NULL;
		if (!m_bProtectGlobalPara)
		{
			ei->opposite()->lpftGlobalPara = NULL;
		}
		ei->opposite()->bSelected = false;
		ei->opposite()->ftForceBarrier = 0;
		ei->opposite()->bSharp = false;
	}
	iIndex = 0;
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
	{
		fi->bDirFixed = false;
		fi->iIndex = iIndex;
		fi->bUnif = false;
		fi->bSelected = false;
		++iIndex;
	}
}


void CMeshInfoInitializer::m_fnFillGeometricInfo()
{
	m_fnInitVector();
	m_fnInitDot();
	m_fnInitFacetNormal();
	m_fnInitLength();
	m_fnInitAngle();
	m_fnInitHELocalInd();
	m_fnMarkBorder();
	m_fnInitLaplaceCoef();
	m_fnInitAngleDefect();
	m_fnInitPolarAxisDif();
	m_fnInscribedCircle();
	m_fnInitMainDirection();
	m_fnInitOtherDefaultInfo();
	m_nNumExpectedFacets = 1000;
}

void CMeshInfoInitializer::m_fnDetectSeamType()
{
    CMeshBase::CPolyhedron::Edge_iterator ei;
    CMeshBase::FT ftAngleDif;
	CMeshBase::CPolyhedron::Halfedge_handle hh;
    for(ei = m_mshSurface.edges_begin(); ei != m_mshSurface.edges_end(); ++ei)
    {
        if(ei->is_border_edge())
        {
			hh = &(*ei);
			if (hh->is_border())
			{
				hh = hh->opposite();
			}
			ftAngleDif = -hh->facet()->ftChartDir + hh->ftPolarAxisDif;
			ftAngleDif -= CMeshBase::FT(int(ftAngleDif / ftTwoPi + 15.5) - 15) * ftTwoPi;
			hh->iSeamType = (int(ftAngleDif / ftHalfPi + 15.5) - 13) % 4 - 2;
			hh->opposite()->iSeamType = (2 - hh->iSeamType) % 4 - 2;
        }
        else
        {
            ftAngleDif = ei->opposite()->facet()->ftChartDir - ei->facet()->ftChartDir + ei->ftPolarAxisDif;
            ftAngleDif = ftAngleDif - CMeshBase::FT(int(ftAngleDif / ftTwoPi + 15.5) - 15) * ftTwoPi;
            ei->iSeamType = (int((ftAngleDif / ftHalfPi) + 15.5) - 13)%4 - 2;
            ei->opposite()->iSeamType = (2 - ei->iSeamType) % 4 - 2;
        }
    }
}

void CMeshInfoInitializer::m_fnCorrectSeams()
{
	;
}

void CMeshInfoInitializer::m_fnMarkSingularity()
{
	CMeshBase::CPolyhedron::Vertex_iterator vi;
	CMeshBase::CPolyhedron::Halfedge_handle hh0, hh;
	CMeshBase::FT  ftCurrentRotation;
	int nExpectedDegree;
	for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
	{
		if (vi->iBordId == -1)
		{
			hh0 = vi->halfedge();
			hh = hh0;
			vi->ftAngleBias = 0;
			do
			{
				ftCurrentRotation = hh->opposite()->facet()->ftChartDir - hh->facet()->ftChartDir - hh->iSeamType * ftHalfPi + hh->ftPolarAxisDif;
				ftCurrentRotation = ftCurrentRotation - CMeshBase::FT(int(ftCurrentRotation / ftTwoPi + 15.5) - 15) * ftTwoPi;
				vi->ftAngleBias += ftCurrentRotation;
				hh = hh->next()->opposite();
			} while (hh != hh0);
			vi->iDegreeBias = int(vi->ftAngleBias / ftHalfPi + 15.5) - 15;
			nExpectedDegree = 4 - int((vi->ftAngleBias + vi->ftAngleDefect) / ftHalfPi + 15.5) + 15;
			if (nExpectedDegree == 4)
			{
				vi->nNumPorts = 0;
			}
			else
			{
				vi->nNumPorts = nExpectedDegree;
			}
		}
		else
		{
			hh0 = vi->halfedge();
			while (!hh0->opposite()->is_border())
			{
				hh0 = hh0->next()->opposite();
			}
			vi->ftAngleBias = 0;
			hh = hh0->opposite()->prev()->opposite()->prev();
			while (hh != hh0)
			{
				ftCurrentRotation = hh->opposite()->facet()->ftChartDir - hh->facet()->ftChartDir - hh->iSeamType * ftHalfPi + hh->ftPolarAxisDif;
				ftCurrentRotation = ftCurrentRotation - CMeshBase::FT(int(ftCurrentRotation / ftTwoPi + 15.5) - 15) * ftTwoPi;
				vi->ftAngleBias += ftCurrentRotation;
				hh = hh->opposite()->prev();
			}
			vi->iDegreeBias = int(vi->ftAngleBias / ftHalfPi + 15.5) - 15;
			nExpectedDegree = 2 - int((vi->ftAngleBias + vi->ftAngleDefect) / ftHalfPi + 15.5) + 15;
			if (nExpectedDegree == 2)
			{
				vi->nNumPorts = 0;
			}
			else
			{
				vi->nNumPorts = nExpectedDegree;
			}
		}
	}
}

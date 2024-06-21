#include "CMeshInfoRepresentor.h"
void CMeshInfoRepresentor::m_fnRepresentGaussCurvature()
{
	CMeshBase::FT ftMax, ftMin;
	CMeshBase::CPolyhedron::Vertex_iterator vi;
	ftMax = -DBL_MAX;
	ftMin = DBL_MAX;
	for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
	{
		if (vi->ftGaussCurvature < ftMin)
		{
			ftMin = vi->ftGaussCurvature;
		}
		if (vi->ftGaussCurvature > ftMax)
		{
			ftMax = vi->ftGaussCurvature;
		}
	}
	for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
	{
		m_fnColorMap(vi->ftGaussCurvature, ftMax, ftMin, vi->lpdblColor);
	}
}
void CMeshInfoRepresentor::m_fnRepresentConformalCoefficient()
{
	CMeshBase::FT ftMax, ftMin, ftLen;
	CMeshBase::CPolyhedron::Vertex_iterator vi;
	ftMax = -DBL_MAX;
	ftMin = DBL_MAX;
	for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
	{
		ftLen = vi->ftDensity;
		if (ftLen < ftMin)
		{
			ftMin = ftLen;
		}
		if (ftLen > ftMax)
		{
			ftMax = ftLen;
		}
	}
	for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
	{
		ftLen = vi->ftDensity;
		m_fnColorMap(ftLen, ftMax, ftMin, vi->lpdblColor);
	}
}

void CMeshInfoRepresentor::m_fnRepresentAngleError(int iDir, CMeshBase::FT*lpftErr)
{
	CMeshBase::CPolyhedron::Facet_iterator fi;
	if (iDir == 0 || iDir == 1)
	{
		for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
		{
			m_fnColorMap(fi->lpftAngleErr[iDir], lpftErr[iDir * 3 + 2], 0, fi->lpdblColor);
		}
	}
	else if (iDir == 2)
	{
		for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
		{
			m_fnColorMap((fi->lpftAngleErr[0] + fi->lpftAngleErr[1]) / 2, lpftErr[iDir * 3 + 2], 0, fi->lpdblColor);
		}
	}
}


void CMeshInfoRepresentor::m_fnRepresentParaCurl(int iDir, CMeshBase::FT* lpftCurl)
{
	CMeshBase::CPolyhedron::Facet_iterator fi;
	if (iDir == 0 || iDir == 1)
	{
		for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
		{
			m_fnColorMap(fi->lpftCurl[iDir], lpftCurl[iDir * 3 + 2], 0, fi->lpdblColor);
		}
	}
	else if (iDir == 2)
	{
		for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
		{
			m_fnColorMap(sqrt((fi->lpftCurl[0] * fi->lpftCurl[0] + fi->lpftCurl[1] * fi->lpftCurl[1])) / 2, lpftCurl[iDir * 3 + 2], 0, fi->lpdblColor);
		}
	}
}

void CMeshInfoRepresentor::m_fnRepresentLocalPara()
{
	CMeshBase::CPolyhedron::Facet_iterator fi;
	int i;
	CMeshBase::CPolyhedron::Halfedge_handle hh0, hh;
	CMeshBase::Vector_3 lpvtRotatedPara[2], lpvtLocalPara[2];
	CMeshBase::FT ftMaxCoef, ftSqScl;
	ftMaxCoef = 0;
	m_dblMeanLocalParaScl = 0;
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
	{
		for (i = 0; i < 2; ++i)
		{
			lpvtLocalPara[i] = CMeshBase::Vector_3(0, 0, 0);
		}
		hh0 = fi->halfedge();
		hh = hh0;
		do
		{
			for (i = 0; i < 2; ++i)
			{
				lpvtRotatedPara[i] = (hh->next()->lpftGlobalPara[i]) * hh->vtVector + (hh->prev()->lpftGlobalPara[i]) * hh->next()->vtVector + (hh->lpftGlobalPara[i]) * hh->prev()->vtVector;
				lpvtLocalPara[i] = lpvtLocalPara[i] + lpvtRotatedPara[i];
			}
			hh = hh->next();
		} while (hh != hh0);
		for (i = 0; i < 2; ++i)
		{
			fi->lpftLocalParaDir[i] = atan2((-lpvtLocalPara[i] * fi->vtPrincipalAxis), (lpvtLocalPara[i] * fi->vtAuxiliaryAxis));
			ftSqScl = lpvtLocalPara[i] * lpvtLocalPara[i];
			fi->lpftLocalParaScl[i] = sqrt(ftSqScl) / fi->ftRadius;
			m_dblMeanLocalParaScl += ftSqScl;
			if (fi->lpftLocalParaScl[i] > ftMaxCoef)
			{
				ftMaxCoef = fi->lpftLocalParaScl[i];
			}
		}
	}
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
	{
		for (i = 0; i < 2; ++i)
		{
			fi->lpftLocalParaScl[i] /= ftMaxCoef;
		}
	}
	m_dblMeanLocalParaScl = sqrt((m_dblMeanLocalParaScl / m_mshSurface.size_of_facets())) / 3.0;
}

void CMeshInfoRepresentor::m_fnColorMap(double dblValue, double dblMaxValue, double dblMinValue, double *lpdblColor)
{
	int i;
	double lpdblOffSet[3] = { -0.25, 0, 0.25 };
	double dblUniformizedValue = (dblValue * 2 - dblMaxValue - dblMinValue) / ((dblMaxValue - dblMinValue) * 2.0);
	double dblOffsetValue;
	for (i = 0; i < 3; ++i)
	{
		dblOffsetValue = dblUniformizedValue + lpdblOffSet[i];
		lpdblColor[i] = (fabs(dblOffsetValue + 0.375) - fabs(dblOffsetValue + 0.125) - fabs(dblOffsetValue - 0.125) + fabs(dblOffsetValue - 0.375)) * 2;
		if (lpdblColor[i] > 1)
		{
			lpdblColor[i] = 1;
		}
		else if (lpdblColor[i] < 0)
		{
			lpdblColor[i] = 0;
		}
	}
	if (dblValue <= dblMinValue)
	{
		lpdblColor[0] = 0.0;
		lpdblColor[1] = 0.0;
		lpdblColor[2] = 0.5;
	}
	if (dblValue >= dblMaxValue)
	{
		lpdblColor[0] = 0.5;
		lpdblColor[1] = 0.0;
		lpdblColor[2] = 0.0;
	}
}

void CMeshInfoRepresentor::m_fnRepresentGlobalPara()
{
	CMeshBase::CPolyhedron::Halfedge_iterator hi;
	CMeshBase::FT lpftParaExtrema[4];
	int i;
	for (i = 0; i < 2; ++i)
	{
		lpftParaExtrema[i] = DBL_MAX;
		lpftParaExtrema[i + 2] = -DBL_MAX;
	}
	for (hi = m_mshSurface.halfedges_begin(); hi != m_mshSurface.halfedges_end(); ++hi)
	{
		if (!hi->is_border())
		{
			for (i = 0; i < 2; ++i)
			{
				if (hi->lpftGlobalPara[i] < lpftParaExtrema[i])
				{
					lpftParaExtrema[i] = hi->lpftGlobalPara[i];
				}
				if (hi->lpftGlobalPara[i] > lpftParaExtrema[i + 2])
				{
					lpftParaExtrema[i + 2] = hi->lpftGlobalPara[i];
				}
			}
		}
	}
	for (hi = m_mshSurface.halfedges_begin(); hi != m_mshSurface.halfedges_end(); ++hi)
	{
		for (i = 0; i < 2; ++i)
		{
			m_fnColorMap(hi->lpftGlobalPara[i], lpftParaExtrema[i + 2], lpftParaExtrema[i], (hi->lpdblColor + 3 + i * 3));
		}
	}
}

void CMeshInfoRepresentor::m_fnRepresentElectricField()
{
	CMeshBase::CPolyhedron::Facet_iterator fi;
	CMeshBase::FT ftLen;
	m_ftMaxElectricField = 0;
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
	{
		fi->lpftElectricFieldSum[0] = fi->vtNablaPhi * fi->vtPrincipalAxis / fi->ftRadius;
		fi->lpftElectricFieldSum[1] = fi->vtNablaPhi * fi->vtAuxiliaryAxis / fi->ftRadius;
		fi->lpftElectricFieldDif[0] = fi->vtDifNabla * fi->vtPrincipalAxis / fi->ftRadius;
		fi->lpftElectricFieldDif[1] = fi->vtDifNabla * fi->vtAuxiliaryAxis / fi->ftRadius;
		ftLen = fi->vtNablaPhi * fi->vtNablaPhi;
		if (ftLen > m_ftMaxElectricField)
		{
			m_ftMaxElectricField = ftLen;
		}
	}
	m_ftMaxElectricField = sqrt(m_ftMaxElectricField);
}

void CMeshInfoRepresentor::m_fnReportSingularities(int &iPositive, int &iNegative)
{
	iPositive = 0;
	iNegative = 0;
	CMeshBase::CPolyhedron::Vertex_iterator vi;
	for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
	{
		if (vi->iBordId == -1)
		{
			if (vi->nNumPorts != 0)
			{
				if (vi->nNumPorts < 4)
				{
					iPositive += (4 - vi->nNumPorts);
				}
				else
				{
					iNegative += (vi->nNumPorts - 4);
				}
			}
		}
	}
}
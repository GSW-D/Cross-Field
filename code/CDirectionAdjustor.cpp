#include "CDirectionAdjustor.h"
CDirectionAdjustor::CDirectionAdjustor(CMeshBase::CPolyhedron &mshSurface): m_mshSurface(mshSurface)
{
    //ctor
}


void CDirectionAdjustor::m_fnGenerateMoment()
{
    CMeshBase::CPolyhedron::Edge_iterator ei;
    CMeshBase::CPolyhedron::Facet_iterator fi;
    CMeshBase::FT ftAngleDif;
    for(ei = m_mshSurface.edges_begin(); ei != m_mshSurface.edges_end(); ++ei)
    {
        if(ei->is_border_edge())
        {
            ei->ftMom = 0;
            ei->opposite()->ftMom = 0;
        }
        else
        {
            ftAngleDif = ei->opposite()->facet()->ftChartDir - ei->facet()->ftChartDir + ei->ftPolarAxisDif;
            ei->ftMom = ftAngleDif - ei->iSeamType * ftHalfPi;
            ei->ftMom = ei->ftMom - CMeshBase::FT(int(ei->ftMom / ftTwoPi + 15.5) - 15) * ftTwoPi;
            ei->opposite()->ftMom = -ei->ftMom;
        }
    }
}



void CDirectionAdjustor::m_fnInitDirectionMatrix(CMeshBase::CLDLTSv &LDLTS)
{
	CMeshBase::CPolyhedron::Facet_iterator fi;
	std::queue<CMeshBase::CPolyhedron::Facet_handle> qufh;
	CMeshBase::CPolyhedron::Facet_handle fh0, fh;
	CMeshBase::CPolyhedron::Halfedge_handle hh0, hh;
	MatElem<CMeshBase::FT> meDiag, meTri;
	int iEquationIndex;
	bool bAddConstraint;
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
	{
		fi->iEquationIndex = -1;
	}
	bAddConstraint = true;
	iEquationIndex = 0;
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
	{
		if (fi->bDirFixed)
		{
			bAddConstraint = false;
		}
		else if (fi->iEquationIndex == -1)
		{
			fi->iEquationIndex = iEquationIndex;
			qufh.push(&(*fi));
			while (!qufh.empty())
			{
				fh0 = qufh.front();
				qufh.pop();
				hh0 = fh0->halfedge();
				hh = hh0;
				do
				{
					if (hh->ftLaplaceCoef <= 0.0 && !hh->opposite()->is_border())
					{
						fh = hh->opposite()->facet();
						if (fh->iEquationIndex == -1)
						{
							fh->ftChartDir -= hh->ftMom;
							hh->opposite()->next()->ftMom += hh->ftMom;
							hh->opposite()->next()->opposite()->ftMom = -hh->opposite()->next()->ftMom;
							hh->opposite()->prev()->ftMom += hh->ftMom;
							hh->opposite()->prev()->opposite()->ftMom = -hh->opposite()->prev()->ftMom;
							hh->ftMom = 0;
							hh->opposite()->ftMom = 0;
							fh->iEquationIndex = iEquationIndex;
							qufh.push(fh);
						}
					}
					hh = hh->next();
				} while (hh != hh0);
			}
			++iEquationIndex;
		}
	}
	if (bAddConstraint)
	{
		--iEquationIndex;
		for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
		{
			if (fi->iEquationIndex == iEquationIndex)
			{
				fi->iEquationIndex = -1;
			}
		}
	}
	LDLTS.init(iEquationIndex);
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
	{
		if (fi->iEquationIndex != -1)
		{
			meDiag.row = fi->iEquationIndex;
			meDiag.col = fi->iEquationIndex;
			meDiag.data = 0;
			hh0 = fi->halfedge();
			hh = hh0;
			do
			{
				if (!hh->opposite()->is_border() && hh->ftLaplaceCoef > 0.0)
				{
					meDiag.data += 1.0 / hh->ftLaplaceCoef;
					if (hh->opposite()->facet()->iEquationIndex != -1)
					{
						meTri.row = fi->iEquationIndex;
						meTri.col = hh->opposite()->facet()->iEquationIndex;
						meTri.data = -1.0 / hh->ftLaplaceCoef;
						LDLTS.add_data(meTri.row, meTri.col, meTri.data);
					}
				}
				hh = hh->next();
			} while (hh != hh0);
			LDLTS.add_data(meDiag.row, meDiag.col, meDiag.data);
		}
	}
}

void CDirectionAdjustor::m_fnFillMoment(CMeshBase::CLDLTSv &LDLTS)
{
	CMeshBase::CPolyhedron::Facet_iterator fi;
	CMeshBase::CPolyhedron::Halfedge_handle hh0, hh;
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
	{
		if (fi->iEquationIndex != -1)
		{
			LDLTS.m_b[fi->iEquationIndex] = 0;
			hh0 = fi->halfedge();
			hh = hh0;
			do
			{
				if (!hh->opposite()->is_border() && hh->ftLaplaceCoef > 0.0)
				{
					LDLTS.m_b[fi->iEquationIndex] += hh->ftMom / hh->ftLaplaceCoef;
				}
				hh = hh->next();
			} while (hh != hh0);
		}
	}
}


void CDirectionAdjustor::m_fnAdjustDirection(CMeshBase::CLDLTSv &LDLTS)
{
	CMeshBase::CPolyhedron::Facet_iterator fi;
	int i;
	CMeshBase::FT ftAve;
	bool bAddConstraint;
	bAddConstraint = true;
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
	{
		if (fi->bDirFixed)
		{
			bAddConstraint = false;
		}
	}
	ftAve = 0;
	if (bAddConstraint)
	{
		for (i = 0; i < LDLTS.m_dim; ++i)
		{
			ftAve += LDLTS.m_x[i];
		}
		ftAve /= LDLTS.m_dim;
	}
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
	{
		if (fi->iEquationIndex != -1)
		{
			fi->ftChartDir += (LDLTS.m_x[fi->iEquationIndex] - ftAve);
			fi->ftChartDir = fi->ftChartDir - CMeshBase::FT(int(fi->ftChartDir / ftTwoPi + 15.5) - 15) * ftTwoPi;
		}
		else if (bAddConstraint)
		{
			fi->ftChartDir -= ftAve;
		}
	}
}

CMeshBase::FT CDirectionAdjustor::m_fnMaxCurlFlux(CMeshBase::CPolyhedron::Facet_handle &fh)
{
	CMeshBase::CPolyhedron::Facet_iterator fi;
	CMeshBase::FT ftCurlFlux, ftMaxCurlFlux, ftSumWeight;
	CMeshBase::CPolyhedron::Halfedge_handle hh0, hh;
	ftMaxCurlFlux = 0;
	fh = NULL;
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
	{
		if (fi->bDirFixed)
		{
			ftCurlFlux = 0;
			ftSumWeight = 0;
			hh0 = fi->halfedge();
			hh = hh0;
			do
			{
				ftCurlFlux += (hh->ftMom / hh->ftLaplaceCoef);
				ftSumWeight += 1.0 / hh->ftLaplaceCoef;
				hh = hh->next();
			} while (hh != hh0);
			ftCurlFlux /= ftSumWeight;
			if (fabs(ftCurlFlux) > fabs(ftMaxCurlFlux))
			{
				ftMaxCurlFlux = ftCurlFlux;
				fh = &(*fi);
			}
		}
	}
	return ftMaxCurlFlux;
}

void CDirectionAdjustor::m_fnAdjustCurlFlux(CMeshBase::CPolyhedron::Facet_handle fh, int iSign)
{
	CMeshBase::CPolyhedron::Halfedge_handle hh0, hh;
	if (fh != NULL && iSign != 0)
	{
		hh0 = fh->halfedge();
		hh = hh0;
		do
		{
			if (iSign > 0)
			{
				hh->iSeamType = (hh->iSeamType + 3) % 4 - 2;
				hh->opposite()->iSeamType = (6 - hh->iSeamType) % 4 - 2;
			}
			else
			{
				hh->iSeamType = (hh->iSeamType + 5) % 4 - 2;
				hh->opposite()->iSeamType = (6 - hh->iSeamType) % 4 - 2;
			}
			hh = hh->next();
		} while (hh != hh0);
	}
}

void CDirectionAdjustor::m_fnInitDensityMatrix(CMeshBase::CLDLTSv &LDLTS)
{
	CMeshBase::CPolyhedron::Vertex_iterator vi;
	CMeshBase::CPolyhedron::Halfedge_handle hh0, hh;
	MatElem<CMeshBase::FT> meDiag, meTri;
	int iEquationIndex;
	iEquationIndex = 0;
	for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
	{
		if (vi->degree() < 3)
		{
			vi->iEquationIndex = -1;
		}
		else
		{
			vi->iEquationIndex = iEquationIndex;
			++iEquationIndex;
		}
	}
	--vi;
	vi->iEquationIndex = -1;
	vi->ftDensity = 0;
	--iEquationIndex;
	LDLTS.init(iEquationIndex);
	for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
	{
		if (vi->iEquationIndex != -1)
		{
			meDiag.row = vi->iEquationIndex;
			meDiag.col = vi->iEquationIndex;
			meDiag.data = 0;
			hh0 = vi->halfedge();
			hh = hh0;
			do
			{
				if (!hh->is_border_edge())
				{
					meDiag.data += hh->ftLaplaceCoef;
					if (hh->prev()->vertex()->iEquationIndex != -1)
					{
						meTri.row = vi->iEquationIndex;
						meTri.col = hh->prev()->vertex()->iEquationIndex;
						meTri.data = -hh->ftLaplaceCoef;
						LDLTS.add_data(meTri.row, meTri.col, meTri.data);
					}
				}
				hh = hh->next()->opposite();
			} while (hh != hh0);
			LDLTS.add_data(meDiag.row, meDiag.col, meDiag.data);
		}
	}
}

void CDirectionAdjustor::m_fnFillPoisson(CMeshBase::CLDLTSv &LDLTS)
{
	CMeshBase::CPolyhedron::Vertex_iterator vi;
	for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
	{
		if (vi->iEquationIndex != -1)
		{
			if (vi->iBordId == -1)
			{
				if (vi->nNumPorts == 0)
				{
					LDLTS.m_b[vi->iEquationIndex] = -vi->ftAngleDefect;
				}
				else
				{
					LDLTS.m_b[vi->iEquationIndex] = -vi->ftAngleDefect + (4 - int(vi->nNumPorts)) * ftHalfPi;
				}
			}
			else
			{
				if (vi->nNumPorts == 0)
				{
					LDLTS.m_b[vi->iEquationIndex] = -vi->ftAngleDefect;
				}
				else
				{
					LDLTS.m_b[vi->iEquationIndex] = -vi->ftAngleDefect + (2 - int(vi->nNumPorts)) * ftHalfPi;
				}
			}
		}
	}
}

void CDirectionAdjustor::m_fnAssignDensity(CMeshBase::CLDLTSv &LDLTS)
{
	CMeshBase::CPolyhedron::Vertex_iterator vi;
	for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
	{
		if (vi->iEquationIndex != -1)
		{
			vi->ftDensity = LDLTS.m_x[vi->iEquationIndex];
		}
	}
}

void CDirectionAdjustor::m_fnModifyDensity(CMeshBase::FT ftChargeRadius, double &dblBottom, double &dblTop)
{
	CMeshBase::CPolyhedron::Vertex_iterator vi;
	std::queue<CMeshBase::CPolyhedron::Vertex_handle> quvh;
	CMeshBase::CPolyhedron::Halfedge_handle hh0, hh;
	CMeshBase::FT ftNum, ftDen, ftElectricQuantity, ftArea, ftAreaRadius, ftZeroLevel;
	CMeshBase::CPolyhedron::Vertex_handle vhSeed;
	CMeshBase::CPolyhedron::Facet_iterator fi;
	for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
	{
		if (vi->iBordId != -1 && vi->nNumPorts != 0)
		{
			ftNum = 0;
			ftDen = 0;
			hh0 = vi->halfedge()->opposite();
			hh = hh0;
			do
			{
				ftNum += hh->vertex()->ftDensity;
				ftDen += 1.0;
				hh = hh->opposite()->next();
			} while (hh != hh0);
			vi->ftViewDensity = ftNum / ftDen;
		}
		else
		{
			vi->ftViewDensity = vi->ftDensity;
		}
	}
	for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
	{
		ftElectricQuantity = vi->iDegreeBias * ftHalfPi - vi->ftAngleDefect;
		if (vi->iBordId == -1 && (ftElectricQuantity > ftHalfPi / 2 || ftElectricQuantity < -ftHalfPi / 2))
		{
			ftNum = 0;
			ftDen = 0;
			ftArea = 0;
			hh0 = vi->halfedge();
			hh = hh0;
			do
			{
				ftNum += hh->ftLaplaceCoef * hh->prev()->vertex()->ftDensity;
				ftDen += hh->ftLaplaceCoef;
				if (!hh->is_border())
				{
					ftArea += hh->facet()->ftTwArea;
				}
				hh = hh->next()->opposite();
			} while (hh != hh0);
			if (vi->iBordId == -1)
			{
				ftAreaRadius = sqrt(ftArea / (2 * ftPi - vi->ftAngleDefect));
				vi->ftViewDensity = ftNum / ftDen + ftElectricQuantity * (log(ftAreaRadius) - log(ftChargeRadius)) / ((ftTwoPi - vi->ftAngleDefect) * 2);
			}
			else
			{
				ftAreaRadius = sqrt(ftArea / (ftPi - vi->ftAngleDefect));
				vi->ftViewDensity = ftNum / ftDen + ftElectricQuantity * (log(ftAreaRadius) - log(ftChargeRadius)) / ((ftPi - vi->ftAngleDefect) * 2);
			}
			if (ftElectricQuantity > 0)
			{
				quvh.push(&(*vi));
				while (!quvh.empty())
				{
					vhSeed = quvh.front();
					quvh.pop();
					hh0 = vhSeed->halfedge()->opposite();
					hh = hh0;
					do
					{
						if (hh->vertex()->ftViewDensity > vi->ftViewDensity)
						{
							hh->vertex()->ftViewDensity = vi->ftViewDensity;
							quvh.push(hh->vertex());
						}
						hh = hh->opposite()->next();
					} while (hh != hh0);
				}
			}
			else
			{
				quvh.push(&(*vi));
				while (!quvh.empty())
				{
					vhSeed = quvh.front();
					quvh.pop();
					hh0 = vhSeed->halfedge()->opposite();
					hh = hh0;
					do
					{
						if (hh->vertex()->ftViewDensity < vi->ftViewDensity)
						{
							hh->vertex()->ftViewDensity = vi->ftViewDensity;
							quvh.push(hh->vertex());
						}
						hh = hh->opposite()->next();
					} while (hh != hh0);
				}
			}
		}
	}
	ftNum = 0;
	ftDen = 0;
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
	{
		hh0 = fi->halfedge();
		ftNum += (hh0->vertex()->ftViewDensity + hh0->next()->vertex()->ftViewDensity + hh0->prev()->vertex()->ftViewDensity) * fi->ftTwArea / 3;
		ftDen += fi->ftTwArea;
	}
	ftZeroLevel = ftNum / ftDen;
	dblBottom = DBL_MAX;
	dblTop = -DBL_MAX;
	for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
	{
		vi->ftViewDensity -= ftZeroLevel;
		if (vi->ftViewDensity > dblTop)
		{
			dblTop = vi->ftViewDensity;
		}
		if (vi->ftViewDensity < dblBottom)
		{
			dblBottom = vi->ftViewDensity;
		}
	}
}

void CDirectionAdjustor::m_fnMoveSingularity()
{
	CMeshBase::CPolyhedron::Vertex_iterator vi;
	CMeshBase::CPolyhedron::Vertex_handle lpvh[2];
	std::vector<CMeshBase::CPolyhedron::Halfedge_handle> vecvh;
	CMeshBase::CPolyhedron::Halfedge_handle hh0, hh, hhPath;
	lpvh[0] = NULL;
	lpvh[1] = NULL;
	for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
	{
		if (vi->iSelected == 1)
		{
			lpvh[0] = &(*vi);
			vi->iSelected = 0;
		}
		if (vi->iSelected == 2)
		{
			lpvh[1] = &(*vi);
			vi->iSelected = 0;
		}
	}
	if (lpvh[0] != NULL && lpvh[1] != NULL)
	{
		hhPath = NULL;
		hh0 = lpvh[0]->halfedge();
		hh = hh0;
		do
		{
			if (hh->bFreeEdge)
			{
				hhPath = hh;
			}
			hh = hh->next()->opposite();
		} while (hh != hh0);
		while (hhPath != NULL)
		{
			vecvh.push_back(hhPath);
			hh0 = hhPath->opposite();
			hhPath = NULL;
			hh = hh0->next()->opposite();
			while (hh != hh0)
			{
				if (hh->bFreeEdge)
				{
					hhPath = hh;
				}
				hh = hh->next()->opposite();
			}
		}
		if (vecvh.size() > 0 && (*(vecvh.rbegin()))->prev()->vertex() == lpvh[1])
		{
			while (vecvh.size() > 0)
			{
				hhPath = *(vecvh.rbegin());
				vecvh.pop_back();
				hhPath->iSeamType = (hhPath->iSeamType + 3) % 4 - 2;
				hhPath->opposite()->iSeamType = (6 - hhPath->iSeamType) % 4 - 2;
			}
		}
	}
}

void CDirectionAdjustor::m_fnMoveSingularity(CMeshBase::CPolyhedron::Vertex_handle vhMin, CMeshBase::CPolyhedron::Vertex_handle vhMax)
{
	CMeshBase::CPolyhedron::Vertex_iterator vi;
	CMeshBase::CPolyhedron::Halfedge_handle hh0, hh, hhPath;
	CMeshBase::CPolyhedron::Vertex_handle vh;
	std::vector<CMeshBase::CPolyhedron::Vertex_handle> vecVH;
	CMeshBase::CPolyhedron::Edge_iterator ei;
	int iHeapPos;
	CMeshBase::FT ftDist, ftMinDist;
	ftMinDist = DBL_MAX;
	for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
	{
		vi->iSourceId = -1;
		vi->hhCome = NULL;
		vi->ftDist = -1;
		vi->bDistConf = false;
		vi->iHeapPos = -1;
	}
	vhMin->iSourceId = 0;
	vhMin->ftDist = 0;
	vhMin->iHeapPos = 0;
	vecVH.push_back(vhMin);
	vhMax->iSourceId = 1;
	vhMax->ftDist = 0;
	vhMax->iHeapPos = 1;
	vecVH.push_back(vhMax);
	while (!vecVH.empty())
	{
		vh = vecVH[0];
		vecVH[0] = *(vecVH.rbegin());
		vecVH[0]->iHeapPos = 0;
		vecVH.pop_back();
		m_fnSift(vecVH, 0);
		vh->iHeapPos = -1;
		vh->bDistConf = true;
		hh0 = vh->halfedge()->opposite();
		hh = hh0;
		do
		{
			if (!hh->vertex()->bDistConf && hh->vertex()->iBordId == -1)
			{
				ftDist = vh->ftDist + hh->ftLen;
				if (hh->vertex()->iHeapPos == -1)
				{
					iHeapPos = vecVH.size();
					hh->vertex()->iSourceId = vh->iSourceId;
					hh->vertex()->hhCome = hh;
					hh->vertex()->ftDist = ftDist;
					hh->vertex()->iHeapPos = iHeapPos;
					vecVH.push_back(hh->vertex());
					m_fnCheck(vecVH, iHeapPos);
				}
				else
				{
					if (ftDist < hh->vertex()->ftDist)
					{
						hh->vertex()->iSourceId = vh->iSourceId;
						hh->vertex()->hhCome = hh;
						hh->vertex()->ftDist = ftDist;
						m_fnCheck(vecVH, hh->vertex()->iHeapPos);
					}
				}
			}
			hh = hh->prev()->opposite();
		} while (hh != hh0);
	}
	ftMinDist = DBL_MAX;
	hhPath = NULL;
	for (ei = m_mshSurface.edges_begin(); ei != m_mshSurface.edges_end(); ++ei)
	{
		if (ei->prev()->vertex()->iSourceId == 0 && ei->vertex()->iSourceId == 1)
		{
			hh0 = &(*ei);
		}
		else if (ei->prev()->vertex()->iSourceId == 1 && ei->vertex()->iSourceId == 0)
		{
			hh0 = ei->opposite();
		}
		else
		{
			hh0 = NULL;
		}
		if (hh0 != NULL)
		{
			ftDist = 0;
			hh = hh0;
			while (hh->vertex()->hhCome != NULL)
			{
				hh = hh->vertex()->hhCome->opposite();
				ftDist -= hh->ftMom / hh->ftLaplaceCoef;
			}
			hh = hh0;
			while (hh != NULL)
			{
				ftDist -= hh->ftMom / hh->ftLaplaceCoef;
				hh = hh->opposite()->vertex()->hhCome;
			}
			if (ftDist < ftMinDist)
			{
				ftMinDist = ftDist;
				hhPath = hh0;
			}
		}
	}
	if (hhPath != NULL)
	{
		if (hhPath->vertex()->iSourceId == 0)
		{
			hhPath = hhPath->opposite();
		}
		hhPath->iSeamType = (hhPath->iSeamType + 3) % 4 - 2;
		hhPath->bFreeEdge = true;
		hh = hhPath->vertex()->hhCome;
		while (hh != NULL)
		{
			hh->iSeamType = (hh->iSeamType + 5) % 4 - 2;
			hh->opposite()->iSeamType = (6 - hh->iSeamType) % 4 - 2;
			hh->bFreeEdge = true;
			hh->opposite()->bFreeEdge = true;
			hh = hh->opposite()->vertex()->hhCome;
		}
		hhPath->opposite()->iSeamType = (6 - hhPath->iSeamType) % 4 - 2;
		hhPath->opposite()->bFreeEdge = true;
		hh = hhPath->opposite()->vertex()->hhCome;
		while (hh != NULL)
		{
			hh->iSeamType = (hh->iSeamType + 3) % 4 - 2;
			hh->opposite()->iSeamType = (6 - hh->iSeamType) % 4 - 2;
			hh->bFreeEdge = true;
			hh->opposite()->bFreeEdge = true;
			hh = hh->opposite()->vertex()->hhCome;
		}
	}
}

void CDirectionAdjustor::m_fnFindPoles(CMeshBase::FT ftChargeRadius, CMeshBase::CPolyhedron::Vertex_handle &vhMin, CMeshBase::CPolyhedron::Vertex_handle &vhMax)
{
	CMeshBase::CPolyhedron::Vertex_iterator vi;
	CMeshBase::CPolyhedron::Halfedge_handle hh0, hh;
	std::vector<CMeshBase::CPolyhedron::Vertex_handle> vecvhPeaks, vecvhValleys;
	std::vector<CMeshBase::FT> vecftPeakDensities, vecftValleyDensities;
	CMeshBase::FT ftArea, ftAreaRadius, ftElectricQuantity, ftNum, ftDen;
	int iMax, iMin, i;
	unsigned int iTag;
	for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
	{
		if (vi->iBordId == -1 || vi->ftAngleBias - vi->iDegreeBias * ftHalfPi < -ftHalfPi / 2 || vi->ftAngleBias - vi->iDegreeBias * ftHalfPi > ftHalfPi / 2)
		{
			iTag = 0;
			hh0 = vi->halfedge()->opposite();
			hh = hh0;
			do
			{
				if (hh->vertex()->iBordId == -1)
				{
					if (iTag == 0)
					{
						if (vi->ftDensity > hh->vertex()->ftDensity)
						{
							iTag = 1;
						}
						else if (vi->ftDensity < hh->vertex()->ftDensity)
						{
							iTag = 2;
						}
					}
					else if (iTag == 1)
					{
						if (vi->ftDensity < hh->vertex()->ftDensity)
						{
							iTag = 3;
						}
					}
					else if (iTag == 2)
					{
						if (vi->ftDensity > hh->vertex()->ftDensity)
						{
							iTag = 3;
						}
					}
				}
				hh = hh->opposite()->next();
			} while (hh != hh0 && iTag != 3);
			if (iTag == 1)
			{
				vecvhPeaks.push_back(&(*vi));
				vecftPeakDensities.push_back(vi->ftDensity);
			}
			else if (iTag == 2)
			{
				vecvhValleys.push_back(&(*vi));
				vecftValleyDensities.push_back(vi->ftDensity);
			}
		}
	}
	for (i = 0; i < vecvhPeaks.size(); ++i)
	{
		ftElectricQuantity = vecvhPeaks[i]->iDegreeBias * ftHalfPi - vecvhPeaks[i]->ftAngleDefect;
		if (ftElectricQuantity > ftHalfPi / 2)
		{
			ftNum = 0;
			ftDen = 0;
			ftArea = 0;
			hh0 = vecvhPeaks[i]->halfedge();
			hh = hh0;
			do
			{
				ftNum += hh->ftLaplaceCoef * hh->prev()->vertex()->ftDensity;
				ftDen += hh->ftLaplaceCoef;
				if (!hh->is_border())
				{
					ftArea += hh->facet()->ftTwArea;
				}
				hh = hh->next()->opposite();
			} while (hh != hh0);
			if (vecvhPeaks[i]->iBordId == -1)
			{
				ftAreaRadius = sqrt(ftArea / (2 * ftPi - vecvhPeaks[i]->ftAngleDefect));
				vecftPeakDensities[i] = ftNum / ftDen + ftElectricQuantity * (log(ftAreaRadius) - log(ftChargeRadius)) / ((ftTwoPi - vecvhPeaks[i]->ftAngleDefect) * 2);
			}
			else
			{
				ftAreaRadius = sqrt(ftArea / (ftPi - vecvhPeaks[i]->ftAngleDefect));
				vecftPeakDensities[i] = ftNum / ftDen + ftElectricQuantity * (log(ftAreaRadius) - log(ftChargeRadius)) / ((ftPi - vecvhPeaks[i]->ftAngleDefect) * 2);
			}
		}
	}
	for (i = 0; i < vecvhValleys.size(); ++i)
	{
		ftElectricQuantity = vecvhValleys[i]->iDegreeBias * ftHalfPi - vecvhValleys[i]->ftAngleDefect;
		if (ftElectricQuantity < -ftHalfPi / 2)
		{
			ftNum = 0;
			ftDen = 0;
			ftArea = 0;
			hh0 = vecvhValleys[i]->halfedge();
			hh = hh0;
			do
			{
				ftNum += hh->ftLaplaceCoef * hh->prev()->vertex()->ftDensity;
				ftDen += hh->ftLaplaceCoef;
				if (!hh->is_border())
				{
					ftArea += hh->facet()->ftTwArea;
				}
				hh = hh->next()->opposite();
			} while (hh != hh0);
			if (vecvhValleys[i]->iBordId == -1)
			{
				ftAreaRadius = sqrt(ftArea / (2 * ftPi - vecvhValleys[i]->ftAngleDefect));
				vecftValleyDensities[i] = ftNum / ftDen + ftElectricQuantity * (log(ftAreaRadius) - log(ftChargeRadius)) / ((ftTwoPi - vecvhValleys[i]->ftAngleDefect) * 2);
			}
			else
			{
				ftAreaRadius = sqrt(ftArea / (ftPi - vecvhValleys[i]->ftAngleDefect));
				vecftValleyDensities[i] = ftNum / ftDen + ftElectricQuantity * (log(ftAreaRadius) - log(ftChargeRadius)) / ((ftPi - vecvhValleys[i]->ftAngleDefect) * 2);
			}
		}
	}
	iMax = 0;
	for (i = 1; i < vecvhPeaks.size(); ++i)
	{
		if (vecftPeakDensities[i] > vecftPeakDensities[iMax])
		{
			iMax = i;
		}
	}
	iMin = 0;
	for (i = 1; i < vecvhValleys.size(); ++i)
	{
		if (vecftValleyDensities[i] < vecftValleyDensities[iMin])
		{
			iMin = i;
		}
	}
	if (iMax < vecvhPeaks.size() && iMin < vecvhValleys.size())
	{
		vhMax = vecvhPeaks[iMax];
		vhMin = vecvhValleys[iMin];
	}
	else
	{
		vhMax = NULL;
		vhMin = NULL;
	}
}

void CDirectionAdjustor::m_fnLinkPoles(CMeshBase::CPolyhedron::Vertex_handle &vhMin, CMeshBase::CPolyhedron::Vertex_handle &vhMax)
{
	CMeshBase::CPolyhedron::Vertex_iterator vi;
	std::vector<CMeshBase::CPolyhedron::Vertex_handle> vecVH;
	CMeshBase::CPolyhedron::Vertex_handle vh;
	CMeshBase::CPolyhedron::Halfedge_handle hh0, hh, hhPath;
	CMeshBase::CPolyhedron::Edge_iterator ei;
	CMeshBase::FT ftDist, ftMinDist;
	int iHeapPos;
	for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
	{
		vi->iSourceId = -1;
		vi->hhCome = NULL;
		vi->ftDist = -1;
		vi->bDistConf = false;
		vi->iHeapPos = -1;
	}
	vhMin->iSourceId = 0;
	vhMin->ftDist = 0;
	vhMin->iHeapPos = 0;
	vecVH.push_back(vhMin);
	vhMax->iSourceId = 1;
	vhMax->ftDist = 0;
	vhMax->iHeapPos = 1;
	vecVH.push_back(vhMax);
	while (!vecVH.empty())
	{
		vh = vecVH[0];
		vecVH[0] = *(vecVH.rbegin());
		vecVH[0]->iHeapPos = 0;
		vecVH.pop_back();
		m_fnSift(vecVH, 0);
		vh->iHeapPos = -1;
		vh->bDistConf = true;
		hh0 = vh->halfedge()->opposite();
		hh = hh0;
		do
		{
			if (!hh->vertex()->bDistConf && hh->vertex()->iBordId == -1)
			{
				ftDist = vh->ftDist + hh->ftLen;
				if (hh->vertex()->iHeapPos == -1)
				{
					iHeapPos = vecVH.size();
					hh->vertex()->iSourceId = vh->iSourceId;
					hh->vertex()->hhCome = hh;
					hh->vertex()->ftDist = ftDist;
					hh->vertex()->iHeapPos = iHeapPos;
					vecVH.push_back(hh->vertex());
					m_fnCheck(vecVH, iHeapPos);
				}
				else
				{
					if (ftDist < hh->vertex()->ftDist)
					{
						hh->vertex()->iSourceId = vh->iSourceId;
						hh->vertex()->hhCome = hh;
						hh->vertex()->ftDist = ftDist;
						m_fnCheck(vecVH, hh->vertex()->iHeapPos);
					}
				}
			}
			hh = hh->prev()->opposite();
		} while (hh != hh0);
	}
	ftMinDist = DBL_MAX;
	hhPath = NULL;
	for (ei = m_mshSurface.edges_begin(); ei != m_mshSurface.edges_end(); ++ei)
	{
		if (ei->prev()->vertex()->iSourceId == 0 && ei->vertex()->iSourceId == 1)
		{
			hh0 = &(*ei);
		}
		else if (ei->prev()->vertex()->iSourceId == 1 && ei->vertex()->iSourceId == 0)
		{
			hh0 = ei->opposite();
		}
		else
		{
			hh0 = NULL;
		}
		if (hh0 != NULL)
		{
			ftDist = hh0->ftLen;
			hh = hh0;
			while (hh->vertex()->hhCome != NULL)
			{
				hh = hh->vertex()->hhCome->opposite();
				ftDist += hh->ftLen;
			}
			hh = hh0;
			while (hh != NULL)
			{
				ftDist += hh->ftLen;
				hh = hh->opposite()->vertex()->hhCome;
			}
			if (ftDist < ftMinDist)
			{
				ftMinDist = ftDist;
				hhPath = hh0;
			}
		}
	}
	if (hhPath != NULL)
	{
		if (hhPath->vertex()->iSourceId == 0)
		{
			hhPath = hhPath->opposite();
		}
		hhPath->bFreeEdge = true;
		hh = hhPath->vertex()->hhCome;
		while (hh != NULL)
		{
			hh->bFreeEdge = true;
			hh->opposite()->bFreeEdge = true;
			hh = hh->opposite()->vertex()->hhCome;
		}
		hhPath->opposite()->bFreeEdge = true;
		hh = hhPath->opposite()->vertex()->hhCome;
		while (hh != NULL)
		{
			hh->bFreeEdge = true;
			hh->opposite()->bFreeEdge = true;
			hh = hh->opposite()->vertex()->hhCome;
		}
	}
}



CMeshBase::FT CDirectionAdjustor::m_fnEnergy(CMeshBase::FT ftChargeRadius)
{
	CMeshBase::FT ftEnergy, ftElectricQuality, ftArea, ftAreaRadius, ftDif;
	CMeshBase::CPolyhedron::Edge_iterator ei;
	CMeshBase::CPolyhedron::Vertex_iterator vi;
	CMeshBase::CPolyhedron::Halfedge_handle hh0, hh;
	CMeshBase::Vector_3 vtRotatedField;
	ftEnergy = 0;
	for (ei = m_mshSurface.edges_begin(); ei != m_mshSurface.edges_end(); ++ei)
	{
		if (!ei->is_border_edge())
		{
			ftDif = ei->vertex()->ftDensity - ei->prev()->vertex()->ftDensity;
			ftEnergy += ftDif * ftDif * ei->ftLaplaceCoef;
		}
	}
	for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
	{
		if (vi->iDegreeBias != 0)
		{
			ftElectricQuality = vi->iDegreeBias * ftHalfPi - vi->ftAngleDefect;
			if (ftElectricQuality < -ftHalfPi / 2 || ftElectricQuality > ftHalfPi / 2)
			{
				ftArea = 0;
				hh0 = vi->halfedge();
				hh = hh0;
				do
				{
					if (!hh->is_border())
					{
						vtRotatedField = hh->vertex()->ftDensity * hh->prev()->vtVector + hh->next()->vertex()->ftDensity * hh->vtVector + hh->prev()->vertex()->ftDensity * hh->next()->vtVector;
						ftEnergy -= vtRotatedField * vtRotatedField / hh->facet()->ftTwArea;
						ftArea += hh->facet()->ftTwArea;
					}
					hh = hh->next()->opposite();
				} while (hh != hh0);
				if (vi->iBordId == -1)
				{
					ftAreaRadius = sqrt(ftArea / (2 * ftPi - vi->ftAngleDefect));
					ftEnergy += (log(ftAreaRadius) - log(ftChargeRadius)) * ftElectricQuality * ftElectricQuality / (2 * (ftTwoPi - vi->ftAngleDefect));
				}
				else
				{
					ftAreaRadius = sqrt(ftArea / (ftPi - vi->ftAngleDefect));
					ftEnergy += (log(ftAreaRadius) - log(ftChargeRadius)) * ftElectricQuality * ftElectricQuality / (2 * (ftPi - vi->ftAngleDefect));
				}
			}
		}
	}
	return ftEnergy;
}

CMeshBase::FT CDirectionAdjustor::m_fnSmoothEnergy(CMeshBase::FT ftChargeRadius)
{
	CMeshBase::FT ftEnergy, ftElectricQuality, ftArea, ftAreaRadius, ftAngleDif;
	CMeshBase::CPolyhedron::Edge_iterator ei;
	CMeshBase::CPolyhedron::Vertex_iterator vi;
	CMeshBase::CPolyhedron::Halfedge_handle hh0, hh;
	CMeshBase::Vector_3 vtRotatedField;
	ftEnergy = 0;
	for (ei = m_mshSurface.edges_begin(); ei != m_mshSurface.edges_end(); ++ei)
	{
		if (!ei->is_border_edge())
		{
			ftAngleDif = ei->opposite()->facet()->ftChartDir - ei->facet()->ftChartDir + ei->ftPolarAxisDif - ei->iSeamType * ftHalfPi;
			ftAngleDif = ftAngleDif - CMeshBase::FT(int(ftAngleDif / ftTwoPi + 15.5) - 15) * ftTwoPi;
			ftEnergy += ftAngleDif * ftAngleDif * ei->ftLaplaceCoef;
		}
	}
	for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
	{
		if (vi->iDegreeBias != 0)
		{
			ftElectricQuality = vi->iDegreeBias * ftHalfPi - vi->ftAngleDefect;
			if (ftElectricQuality < -ftHalfPi / 2 || ftElectricQuality > ftHalfPi / 2)
			{
				ftArea = 0;
				hh0 = vi->halfedge();
				hh = hh0;
				do
				{
					if (!hh->is_border())
					{
						vtRotatedField = hh->vertex()->ftDensity * hh->prev()->vtVector + hh->next()->vertex()->ftDensity * hh->vtVector + hh->prev()->vertex()->ftDensity * hh->next()->vtVector;
						ftEnergy -= vtRotatedField * vtRotatedField / hh->facet()->ftTwArea;
						ftArea += hh->facet()->ftTwArea;
					}
					hh = hh->next()->opposite();
				} while (hh != hh0);
				if (vi->iBordId == -1)
				{
					ftAreaRadius = sqrt(ftArea / (2 * ftPi - vi->ftAngleDefect));
					ftEnergy += (log(ftAreaRadius) - log(ftChargeRadius)) * ftElectricQuality * ftElectricQuality / (2 * (ftTwoPi - vi->ftAngleDefect));
				}
				else
				{
					ftAreaRadius = sqrt(ftArea / (ftPi - vi->ftAngleDefect));
					ftEnergy += (log(ftAreaRadius) - log(ftChargeRadius)) * ftElectricQuality * ftElectricQuality / (2 * (ftPi - vi->ftAngleDefect));
				}
			}
		}
	}
	return ftEnergy;
}

CMeshBase::FT CDirectionAdjustor::m_fnCurlEnergy()
{
	CMeshBase::CPolyhedron::Facet_iterator fi;
	CMeshBase::FT ftEnergy;
	ftEnergy = 0;
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
	{
		ftEnergy += (fi->vtDifNabla * fi->vtDifNabla) * fi->ftTwArea;
	}
	return ftEnergy;
}

void CDirectionAdjustor::m_fnMarkCurlPath()
{
	CMeshBase::CPolyhedron::Vertex_iterator vi;
	CMeshBase::CPolyhedron::Edge_iterator ei;
	CMeshBase::CPolyhedron::Facet_iterator fi;
	CMeshBase::FT ftMaxDensity, ftDensity;
	CMeshBase::CPolyhedron::Halfedge_handle hhSeed, hhStart, hh0, hh1, hh2, hhClose;
	CMeshBase::CPolyhedron::Facet_handle fhStart;
	std::queue<CMeshBase::CPolyhedron::Halfedge_handle> lpquhh[2];
	int i;
	for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
	{
		vi->bDistConf = 0;
		vi->hhCome = NULL;
	}
	for (ei = m_mshSurface.edges_begin(); ei != m_mshSurface.edges_end(); ++ei)
	{
		ei->bFreeEdge = false;
		ei->opposite()->bFreeEdge = false;
		if (ei->is_border_edge() || ei->ftLaplaceCoef == 0.0)
		{
			ei->ftCurl = 0;
			ei->opposite()->ftCurl = 0;
		}
		else
		{
			ei->ftCurl = ei->ftMom / ei->ftLaplaceCoef - (ei->vertex()->ftDensity - ei->prev()->vertex()->ftDensity);
			ei->opposite()->ftCurl = -ei->ftCurl;
		}
	}
	ftMaxDensity = 0;
	//for (ei = m_mshSurface.edges_begin(); ei != m_mshSurface.edges_end(); ++ei)
	//{
	//	if (!ei->is_border_edge())
	//	{
	//		ftDensity = ei->ftCurl / ei->ftLen;
	//		if (ftDensity > ftMaxDensity)
	//		{
	//			ftMaxDensity = ftDensity;
	//			hhStart = &(*ei);
	//		}
	//		if (ftDensity < -ftMaxDensity)
	//		{
	//			ftMaxDensity = -ftDensity;
	//			hhStart = ei->opposite();
	//		}
	//	}
	//}
	fhStart = NULL;
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
	{
		ftDensity = fi->vtDifNabla * fi->vtDifNabla;
		if (ftDensity > ftMaxDensity)
		{
			ftMaxDensity = ftDensity;
			fhStart = &(*fi);
		}
	}
	hhStart = NULL;
	ftMaxDensity = 0;
	hh0 = fhStart->halfedge();
	hh1 = hh0;
	do
	{
		if (!hh1->is_border_edge() && hh1->ftLaplaceCoef > 0)
		{
			ftDensity = hh1->ftCurl / hh1->ftLen;
			if (ftDensity > ftMaxDensity)
			{
				ftMaxDensity = ftDensity;
				hhStart = hh1;
			}
			if (ftDensity < -ftMaxDensity)
			{
				ftMaxDensity = -ftDensity;
				hhStart = hh1->opposite();
			}
		}
		hh1 = hh1->next();
	} while (hh1 != hh0);
	hhStart->vertex()->bDistConf = 1;
	hhStart->opposite()->vertex()->bDistConf = 2;
	hhClose = NULL;
	lpquhh[0].push(hhStart);
	lpquhh[1].push(hhStart->opposite());
	i = 0;
	while (!(lpquhh[0].empty() && lpquhh[1].empty()))
	{
		if (!lpquhh[i].empty())
		{
			hhSeed = lpquhh[i].front();
			lpquhh[i].pop();
			ftMaxDensity = 0;
			hh0 = hhSeed->opposite();
			hh1 = NULL;
			hh2 = hhSeed->next();
			while (hh2 != hh0)
			{
				ftDensity = (hh2->ftCurl / hh2->ftLen) * (1 - i * 2);
				if (ftDensity > ftMaxDensity)
				{
					hh1 = hh2;
					ftMaxDensity = ftDensity;
				}
				hh2 = hh2->opposite()->next();
			}
			if (hh1 != NULL)
			{
				hh2 = hh1->opposite()->next();
				while (!(hh1 == hh0 || hh2 == hh0))
				{
					if (hh2 == hh0 || ((hh2->ftCurl / hh2->ftLen - hh1->ftCurl / hh1->ftLen) * (1 - i * 2) < 0))
					{
						if ((hh1->ftCurl / hh1->ftLen) * (1 - i * 2) > 0)
						{
							if (hh1->vertex()->bDistConf == 0)
							{
								hh1->vertex()->hhCome = hh1;
								hh1->vertex()->bDistConf = i + 1;
								lpquhh[i].push(hh1);
							}
							else if (hh1->vertex()->bDistConf != i + 1)
							{
								hhClose = hh1;
								while (!lpquhh[0].empty())
								{
									lpquhh[0].pop();
								}
								while (!lpquhh[1].empty())
								{
									lpquhh[1].pop();
								}
							}
						}
						hh1 = hh1->prev()->opposite();
					}
					else if (hh1 == hh0 || ((hh1->ftCurl / hh1->ftLen - hh2->ftCurl / hh2->ftLen) * (1 - i * 2) < 0))
					{
						if ((hh2->ftCurl / hh2->ftLen) * (1 - i * 2) > 0)
						{
							if (hh2->vertex()->bDistConf == 0)
							{
								hh2->vertex()->hhCome = hh2;
								hh2->vertex()->bDistConf = i + 1;
								lpquhh[i].push(hh2);
							}
							else if (hh2->vertex()->bDistConf != i + 1)
							{
								hhClose = hh2;
								while (!lpquhh[0].empty())
								{
									lpquhh[0].pop();
								}
								while (!lpquhh[1].empty())
								{
									lpquhh[1].pop();
								}
							}
						}
						hh2 = hh2->opposite()->next();
					}
				}
			}
		}
		i = (i + 1) % 2;
	}
	if (hhClose != NULL)
	{
		hhStart->bFreeEdge = true;
		hhStart->opposite()->bFreeEdge = true;
		hh1 = hhClose;
		do
		{
			hh1->bFreeEdge = true;
			hh1->opposite()->bFreeEdge = true;
			hh1 = hh1->prev()->vertex()->hhCome;
		} while (hh1 != NULL);
		hh2 = hhClose->vertex()->hhCome;
		while (hh2 != NULL)
		{
			hh2->bFreeEdge = true;
			hh2->opposite()->bFreeEdge = true;
			hh2 = hh2->prev()->vertex()->hhCome;
		}
	}
}

bool CDirectionAdjustor::m_fnTestCurlPath()
{
	CMeshBase::CPolyhedron::Facet_iterator fi;
	std::queue<CMeshBase::CPolyhedron::Facet_handle> qufh;
	CMeshBase::CPolyhedron::Facet_handle fh0;
	CMeshBase::CPolyhedron::Halfedge_handle hh0, hh;
	CMeshBase::CPolyhedron::Edge_iterator ei;
	bool bCurl;
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
	{
		fi->iTempIndex = -1;
	}
	fh0 = &(*m_mshSurface.facets_begin());
	fh0->iTempIndex = 0;
	qufh.push(fh0);
	while (!qufh.empty())
	{
		fh0 = qufh.front();
		qufh.pop();
		hh0 = fh0->halfedge();
		hh = hh0;
		do
		{
			if (!(hh->opposite()->is_border() || hh->bFreeEdge) && hh->opposite()->facet()->iTempIndex == -1)
			{
				hh->opposite()->facet()->iTempIndex = 0;
				qufh.push(hh->opposite()->facet());
			}
			hh = hh->next();
		} while (hh != hh0);
	}
	bCurl = true;
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_begin(); ++fi)
	{
		if (fi->iTempIndex == -1)
		{
			bCurl = false;
			break;
		}
	}
	if (!bCurl)
	{
		for (ei = m_mshSurface.edges_begin(); ei != m_mshSurface.edges_end(); ++ei)
		{
			ei->bFreeEdge = false;
			ei->opposite()->bFreeEdge = false;
		}
	}
	return bCurl;
}

void CDirectionAdjustor::m_fnReleaseCurl()
{
	CMeshBase::CPolyhedron::Edge_iterator ei;
	for (ei = m_mshSurface.edges_begin(); ei != m_mshSurface.edges_end(); ++ei)
	{
		if (ei->bFreeEdge)
		{
			if (ei->ftCurl > 0)
			{
				ei->iSeamType = (ei->iSeamType + 3) % 4 - 2;
				ei->opposite()->iSeamType = (6 - ei->iSeamType) % 4 - 2;
			}
			else
			{
				ei->iSeamType = (ei->iSeamType + 5) % 4 - 2;
				ei->opposite()->iSeamType = (6 - ei->iSeamType) % 4 - 2;
			}
			ei->bFreeEdge = false;
			ei->opposite()->bFreeEdge = false;
		}
	}
}

void CDirectionAdjustor::m_fnGenerateNablaField()
{
	CMeshBase::CPolyhedron::Facet_iterator fi;
	CMeshBase::CPolyhedron::Halfedge_handle hh;
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
	{
		hh = fi->halfedge();
		if (!(hh->opposite()->is_border() || hh->next()->opposite()->is_border() || hh->prev()->opposite()->is_border()))
		{
			if (hh->next()->ftLaplaceCoef <= hh->ftLaplaceCoef && hh->next()->ftLaplaceCoef <= hh->prev()->ftLaplaceCoef)
			{
				hh = hh->next();
			}
			else if (hh->prev()->ftLaplaceCoef <= hh->ftLaplaceCoef && hh->prev()->ftLaplaceCoef <= hh->next()->ftLaplaceCoef)
			{
				hh = hh->prev();
			}
			if (hh->prev()->ftLaplaceCoef == 0 || hh->next()->ftLaplaceCoef == 0)
			{
				hh = NULL;
			}
		}
		else
		{
			if (!(hh->opposite()->is_border() || hh->next()->opposite()->is_border()))
			{
				hh = hh->prev();
			}
			else if (!(hh->opposite()->is_border() || hh->prev()->opposite()->is_border()))
			{
				hh = hh->next();
			}
			if (hh->opposite()->is_border())
			{
				hh = NULL;
			}
		}
		if (hh != NULL)
		{
			fi->vtNablaTheta = ((hh->prev()->ftMom / hh->prev()->ftLaplaceCoef) * hh->next()->vtVector - (hh->next()->ftMom / hh->next()->ftLaplaceCoef) * hh->prev()->vtVector) / fi->ftTwArea;
			fi->vtNablaTheta = CGAL::cross_product(fi->vtNorm, fi->vtNablaTheta) / fi->ftTwArea;
			fi->vtNablaPhi = (hh->vertex()->ftDensity * hh->prev()->vtVector + hh->next()->vertex()->ftDensity * hh->vtVector + hh->prev()->vertex()->ftDensity * hh->next()->vtVector) / fi->ftTwArea;
			fi->vtNablaPhi = CGAL::cross_product(fi->vtNorm, fi->vtNablaPhi) / fi->ftTwArea;
			fi->vtDifNabla = fi->vtNablaTheta - fi->vtNablaPhi;
		}
		else
		{
			fi->vtNablaPhi = CMeshBase::Vector_3(0, 0, 0);
			fi->vtDifNabla = CMeshBase::Vector_3(0, 0, 0);
		}
	}
}

void CDirectionAdjustor::m_fnGenerateForces()
{
	CMeshBase::CPolyhedron::Vertex_iterator vi;
	CMeshBase::CPolyhedron::Halfedge_handle hh0, hh;
	CMeshBase::Vector_3 vtOuterNorm;
	for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
	{
		vi->vtThetaForce = CMeshBase::Vector_3(0, 0, 0);
		vi->vtPhiForce = CMeshBase::Vector_3(0, 0, 0);
		if (vi->iDegreeBias != 0)
		{
			hh0 = vi->halfedge();
			hh = hh0;
			do
			{
				if (!hh->is_border())
				{
					vtOuterNorm = CGAL::cross_product(hh->prev()->vtVector, hh->facet()->vtNorm) / hh->facet()->ftTwArea;
					vi->vtPhiForce = vi->vtPhiForce + (vtOuterNorm * hh->facet()->vtNablaPhi) * hh->facet()->vtNablaPhi - (hh->facet()->vtNablaPhi * hh->facet()->vtNablaPhi / 2) * vtOuterNorm;
					vi->vtThetaForce = vi->vtThetaForce + (vtOuterNorm * hh->facet()->vtNablaTheta) * hh->facet()->vtNablaTheta - (hh->facet()->vtNablaTheta * hh->facet()->vtNablaTheta / 2) * vtOuterNorm;
				}
				hh = hh->next()->opposite();
			} while (hh != hh0);
			vi->vtCurlForce = vi->vtThetaForce - vi->vtPhiForce;
		}
	}
}

void CDirectionAdjustor::m_fnGenerateParaForce()
{
	CMeshBase::CPolyhedron::Vertex_iterator vi;
	CMeshBase::CPolyhedron::Halfedge_handle hh0, hh;
	CMeshBase::Vector_3 vtOuterNorm, vtNablaPara, vtLocalForce;
	int i;
	for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
	{
		vi->vtCurlForce = CMeshBase::Vector_3(0, 0, 0);
		if (vi->iDegreeBias != 0)
		{
			hh0 = vi->halfedge();
			hh = hh0;
			do
			{
				if (!hh->is_border())
				{
					vtOuterNorm = CGAL::cross_product(hh->prev()->vtVector, hh->facet()->vtNorm) / hh->facet()->ftTwArea;
					for (i = 0; i < 2; ++i)
					{
						vtNablaPara = hh->vtVector * hh->next()->lpftGlobalPara[i] + hh->next()->vtVector * hh->prev()->lpftGlobalPara[i] + hh->prev()->vtVector * hh->lpftGlobalPara[i];
						//vtNablaPara = CGAL::cross_product(hh->facet()->vtNorm, vtNablaPara) / hh->facet()->ftTwArea;
						vi->vtCurlForce = vi->vtCurlForce + (vtOuterNorm * vtNablaPara) * vtNablaPara - (vtNablaPara * vtNablaPara / 2) * vtOuterNorm;
					}
				}
				hh = hh->next()->opposite();
			} while (hh != hh0);
		}
	}
}


int CDirectionAdjustor::m_fnAdjustSingularitiesTheta()
{
	CMeshBase::CPolyhedron::Vertex_iterator vi;
	CMeshBase::CPolyhedron::Edge_iterator ei;
	CMeshBase::CPolyhedron::Halfedge_handle hh0, hh, hhPath;
	CMeshBase::FT ftMinDot, ftDot;
	int nAdjustment, iMerge;
	for (ei = m_mshSurface.edges_begin(); ei != m_mshSurface.edges_end(); ++ei)
	{
		ei->ftForce = 0;
		ei->opposite()->ftForce = 0;
	}
	for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
	{
		if (vi->iDegreeBias != 0)
		{
			hh0 = vi->halfedge();
			hh = hh0;
			ftMinDot = 0;
			hhPath = NULL;
			do
			{
				if (hh->prev()->vertex()->iDegreeBias * vi->iDegreeBias < 0)
				{
					ftDot = hh->vtVector * vi->vtThetaForce / hh->ftLen;
					if (ftDot < ftMinDot && ftDot < 0.0)
					{
						ftMinDot = ftDot;
						hhPath = hh;
					}
				}
				hh = hh->next()->opposite();
			} while (hh != hh0);
			if (hhPath == NULL)
			{
				do
				{
					ftDot = hh->vtVector * vi->vtThetaForce / hh->ftLen;
					if (ftDot < ftMinDot)
					{
						ftMinDot = ftDot;
						hhPath = hh;
					}
					hh = hh->next()->opposite();
				} while (hh != hh0);
			}
			if (hhPath != NULL)
			{
				if (vi->iDegreeBias < 0)
				{
					hhPath->ftForce += ftMinDot;
					hhPath->opposite()->ftForce -= ftMinDot;
				}
				else
				{
					hhPath->ftForce -= ftMinDot;
					hhPath->opposite()->ftForce += ftMinDot;
				}
			}
		}
	}
	for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
	{
		vi->iHeapPos = -1;
	}
	nAdjustment = 0;
	for (ei = m_mshSurface.edges_begin(); ei != m_mshSurface.edges_end(); ++ei)
	{
		if (ei->ftForce != 0.0)
		{
			if (fabs(ei->ftForce) > fabs(ei->ftForceBarrier))
			{
				++nAdjustment;
			}
		}
	}
	nAdjustment = 0;
	iMerge = 0;
	for (ei = m_mshSurface.edges_begin(); ei != m_mshSurface.edges_end(); ++ei)
	{
		if (ei->ftForce != 0.0)
		{
			if (fabs(ei->ftForce) > fabs(ei->ftForceBarrier))
			{
				if (ei->prev()->vertex()->iDegreeBias == 0 || ei->vertex()->iDegreeBias == 0)
				{
					ei->ftForceBarrier = ei->ftForce;
					ei->opposite()->ftForceBarrier = ei->opposite()->ftForce;
				}
				else
				{
					ei->ftForceBarrier = 0;
					ei->opposite()->ftForceBarrier = 0;
					++iMerge;
				}
				if (ei->is_border_edge())
				{
					if (ei->is_border())
					{
						hh0 = ei->opposite();
					}
					else
					{
						hh0 = &(*ei);
					}
					if (hh0->ftForce > 0)
					{
						hh0->next()->iSeamType = (hh0->next()->iSeamType + 5) % 4 - 2;
						hh0->next()->opposite()->iSeamType = (6 - hh0->next()->iSeamType) % 4 - 2;
						hh0->prev()->iSeamType = (hh0->prev()->iSeamType + 5) % 4 - 2;
						hh0->prev()->opposite()->iSeamType = (6 - hh0->prev()->iSeamType) % 4 - 2;
					}
					else
					{
						hh0->next()->iSeamType = (hh0->next()->iSeamType + 3) % 4 - 2;
						hh0->next()->opposite()->iSeamType = (6 - hh0->next()->iSeamType) % 4 - 2;
						hh0->prev()->iSeamType = (hh0->prev()->iSeamType + 3) % 4 - 2;
						hh0->prev()->opposite()->iSeamType = (6 - hh0->prev()->iSeamType) % 4 - 2;
					}
				}
				else
				{
					if (ei->ftForce > 0)
					{
						ei->iSeamType = (ei->iSeamType + 3) % 4 - 2;
						ei->opposite()->iSeamType = (6 - ei->iSeamType) % 4 - 2;
					}
					else
					{
						ei->iSeamType = (ei->iSeamType + 5) % 4 - 2;
						ei->opposite()->iSeamType = (6 - ei->iSeamType) % 4 - 2;
					}
				}
				++nAdjustment;
			}
		}
	}
	if (iMerge > 1)
	{
		for (ei = m_mshSurface.edges_begin(); ei != m_mshSurface.edges_end(); ++ei)
		{
			ei->ftForceBarrier = 0;
			ei->opposite()->ftForce = 0;
		}
	}
	return nAdjustment;
}

int CDirectionAdjustor::m_fnAdjustSingularitiesPhi()
{
	CMeshBase::CPolyhedron::Vertex_iterator vi;
	CMeshBase::CPolyhedron::Edge_iterator ei;
	CMeshBase::CPolyhedron::Halfedge_handle hh0, hh, hhPath;
	CMeshBase::FT ftMinDot, ftDot;
	int nAdjustment, iMerge;
	for (ei = m_mshSurface.edges_begin(); ei != m_mshSurface.edges_end(); ++ei)
	{
		ei->ftForce = 0;
		ei->opposite()->ftForce = 0;
	}
	for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
	{
		if (vi->iDegreeBias != 0)
		{
			hh0 = vi->halfedge();
			hh = hh0;
			ftMinDot = 0;
			hhPath = NULL;
			do
			{
				if (hh->prev()->vertex()->iDegreeBias * vi->iDegreeBias < 0)
				{
					ftDot = hh->vtVector * vi->vtPhiForce / hh->ftLen;
					if (ftDot < ftMinDot && ftDot < 0.0)
					{
						ftMinDot = ftDot;
						hhPath = hh;
					}
				}
				hh = hh->next()->opposite();
			} while (hh != hh0);
			if (hhPath == NULL)
			{
				do
				{
					ftDot = hh->vtVector * vi->vtPhiForce / hh->ftLen;
					if (ftDot < ftMinDot)
					{
						ftMinDot = ftDot;
						hhPath = hh;
					}
					hh = hh->next()->opposite();
				} while (hh != hh0);
			}
			if (hhPath != NULL)
			{
				if (vi->iDegreeBias < 0)
				{
					hhPath->ftForce += ftMinDot;
					hhPath->opposite()->ftForce -= ftMinDot;
				}
				else
				{
					hhPath->ftForce -= ftMinDot;
					hhPath->opposite()->ftForce += ftMinDot;
				}
			}
		}
	}
	for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
	{
		vi->iHeapPos = -1;
	}
	nAdjustment = 0;
	for (ei = m_mshSurface.edges_begin(); ei != m_mshSurface.edges_end(); ++ei)
	{
		if (ei->ftForce != 0.0)
		{
			if (fabs(ei->ftForce) > fabs(ei->ftForceBarrier))
			{
				++nAdjustment;
			}
		}
	}
	nAdjustment = 0;
	iMerge = 0;
	for (ei = m_mshSurface.edges_begin(); ei != m_mshSurface.edges_end(); ++ei)
	{
		if (ei->ftForce != 0.0)
		{
			if (fabs(ei->ftForce) > fabs(ei->ftForceBarrier))
			{
				if (ei->prev()->vertex()->iDegreeBias == 0 || ei->vertex()->iDegreeBias == 0)
				{
					ei->ftForceBarrier = ei->ftForce;
					ei->opposite()->ftForceBarrier = ei->opposite()->ftForce;
				}
				else
				{
					ei->ftForceBarrier = 0;
					ei->opposite()->ftForceBarrier = 0;
					++iMerge;
				}
				if (ei->is_border_edge())
				{
					if (ei->is_border())
					{
						hh0 = ei->opposite();
					}
					else
					{
						hh0 = &(*ei);
					}
					if (hh0->ftForce > 0)
					{
						hh0->next()->iSeamType = (hh0->next()->iSeamType + 5) % 4 - 2;
						hh0->next()->opposite()->iSeamType = (6 - hh0->next()->iSeamType) % 4 - 2;
						hh0->prev()->iSeamType = (hh0->prev()->iSeamType + 5) % 4 - 2;
						hh0->prev()->opposite()->iSeamType = (6 - hh0->prev()->iSeamType) % 4 - 2;
					}
					else
					{
						hh0->next()->iSeamType = (hh0->next()->iSeamType + 3) % 4 - 2;
						hh0->next()->opposite()->iSeamType = (6 - hh0->next()->iSeamType) % 4 - 2;
						hh0->prev()->iSeamType = (hh0->prev()->iSeamType + 3) % 4 - 2;
						hh0->prev()->opposite()->iSeamType = (6 - hh0->prev()->iSeamType) % 4 - 2;
					}
				}
				else
				{
					if (ei->ftForce > 0)
					{
						ei->iSeamType = (ei->iSeamType + 3) % 4 - 2;
						ei->opposite()->iSeamType = (6 - ei->iSeamType) % 4 - 2;
					}
					else
					{
						ei->iSeamType = (ei->iSeamType + 5) % 4 - 2;
						ei->opposite()->iSeamType = (6 - ei->iSeamType) % 4 - 2;
					}
				}
				++nAdjustment;
			}
		}
	}
	if (iMerge > 1)
	{
		for (ei = m_mshSurface.edges_begin(); ei != m_mshSurface.edges_end(); ++ei)
		{
			ei->ftForceBarrier = 0;
			ei->opposite()->ftForce = 0;
		}
	}
	return nAdjustment;
}

int CDirectionAdjustor::m_fnAdjustSingularitiesPara()
{
	CMeshBase::CPolyhedron::Vertex_iterator vi;
	CMeshBase::CPolyhedron::Edge_iterator ei;
	CMeshBase::CPolyhedron::Halfedge_handle hh0, hh, hhPath;
	CMeshBase::FT ftMinDot, ftDot;
	int nAdjustment;
	for (ei = m_mshSurface.edges_begin(); ei != m_mshSurface.edges_end(); ++ei)
	{
		ei->ftForce = 0;
		ei->opposite()->ftForce = 0;
	}
	for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
	{
		if (vi->iDegreeBias != 0)
		{
			hh0 = vi->halfedge();
			hh = hh0;
			ftMinDot = 0;
			hhPath = NULL;
			do
			{
				ftDot = hh->vtVector * vi->vtCurlForce / hh->ftLen;
				if (ftDot < ftMinDot)
				{
					ftMinDot = ftDot;
					hhPath = hh;
				}
				hh = hh->next()->opposite();
			} while (hh != hh0);
			if (hhPath != NULL)
			{
				if (vi->iDegreeBias < 0)
				{
					hhPath->ftForce += ftMinDot;
					hhPath->opposite()->ftForce -= ftMinDot;
				}
				else
				{
					hhPath->ftForce -= ftMinDot;
					hhPath->opposite()->ftForce += ftMinDot;
				}
			}
		}
	}
	for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
	{
		vi->iHeapPos = -1;
	}
	nAdjustment = 0;
	for (ei = m_mshSurface.edges_begin(); ei != m_mshSurface.edges_end(); ++ei)
	{
		if (ei->ftForce != 0.0)
		{
			if (fabs(ei->ftForce) > fabs(ei->ftForceBarrier))
			{
				++nAdjustment;
			}
		}
	}
	nAdjustment = 0;
	for (ei = m_mshSurface.edges_begin(); ei != m_mshSurface.edges_end(); ++ei)
	{
		if (ei->ftForce != 0.0)
		{
			if (fabs(ei->ftForce) > fabs(ei->ftForceBarrier))
			{
				ei->ftForceBarrier = ei->ftForce;
				ei->opposite()->ftForceBarrier = ei->opposite()->ftForce;
				if (ei->is_border_edge())
				{
					if (ei->is_border())
					{
						hh0 = ei->opposite();
					}
					else
					{
						hh0 = &(*ei);
					}
					if (hh0->ftForce > 0)
					{
						hh0->next()->iSeamType = (hh0->next()->iSeamType + 5) % 4 - 2;
						hh0->next()->opposite()->iSeamType = (6 - hh0->next()->iSeamType) % 4 - 2;
						hh0->prev()->iSeamType = (hh0->prev()->iSeamType + 5) % 4 - 2;
						hh0->prev()->opposite()->iSeamType = (6 - hh0->prev()->iSeamType) % 4 - 2;
					}
					else
					{
						hh0->next()->iSeamType = (hh0->next()->iSeamType + 3) % 4 - 2;
						hh0->next()->opposite()->iSeamType = (6 - hh0->next()->iSeamType) % 4 - 2;
						hh0->prev()->iSeamType = (hh0->prev()->iSeamType + 3) % 4 - 2;
						hh0->prev()->opposite()->iSeamType = (6 - hh0->prev()->iSeamType) % 4 - 2;
					}
				}
				else
				{
					if (ei->ftForce > 0)
					{
						ei->iSeamType = (ei->iSeamType + 3) % 4 - 2;
						ei->opposite()->iSeamType = (6 - ei->iSeamType) % 4 - 2;
					}
					else
					{
						ei->iSeamType = (ei->iSeamType + 5) % 4 - 2;
						ei->opposite()->iSeamType = (6 - ei->iSeamType) % 4 - 2;
					}
				}
				++nAdjustment;
			}
		}
	}
	return nAdjustment;
}


int CDirectionAdjustor::m_fnConstrainedDecurl()
{
	CMeshBase::CPolyhedron::Vertex_iterator vi;
	CMeshBase::CPolyhedron::Halfedge_handle hh0, hh, hhPath;
	CMeshBase::CPolyhedron::Vertex_handle vh;
	CMeshBase::FT ftMinDot, ftDot, ftWeight, ftMaxWeight;
	int nAdjustment;
	ftMaxWeight = 0.0;
	vh = NULL;
	nAdjustment = 0;
	for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
	{
		if (vi->iDegreeBias != 0)
		{
			if (vi->vtCurlForce * vi->vtThetaForce >= 0)
			{
				ftWeight = vi->vtCurlForce * vi->vtCurlForce;
				if (ftWeight > ftMaxWeight)
				{
					ftMaxWeight = ftWeight;
					vh = &(*vi);
				}
			}
		}
	}
	if (vh == NULL)
	{
		ftMaxWeight = 0.0;
		for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
		{
			if (vi->iDegreeBias != 0)
			{
				ftWeight = -vi->vtCurlForce * vi->vtCurlForce / (vi->vtCurlForce * vi->vtThetaForce);
				if (ftWeight > ftMaxWeight)
				{
					ftMaxWeight = ftWeight;
					vh = &(*vi);
				}
			}
		}
	}
	if (vh != NULL)
	{
		hh0 = vh->halfedge();
		hh = hh0;
		ftMinDot = 0;
		hhPath = NULL;
		do
		{
			if (!hh->is_border_edge())
			{
				ftDot = hh->vtVector * vh->vtCurlForce / hh->ftLen;
				if (ftDot < ftMinDot)
				{
					ftMinDot = ftDot;
					hhPath = hh;
				}
			}
			hh = hh->next()->opposite();
		} while (hh != hh0);
		if (hhPath != NULL)
		{
			if (vh->iDegreeBias < 0)
			{
				hhPath->iSeamType = (hhPath->iSeamType + 5) % 4 - 2;
				hhPath->opposite()->iSeamType = (6 - hhPath->iSeamType) % 4 - 2;
			}
			else
			{
				hhPath->iSeamType = (hhPath->iSeamType + 3) % 4 - 2;
				hhPath->opposite()->iSeamType = (6 - hhPath->iSeamType) % 4 - 2;
			}
			++nAdjustment;
		}
	}
	return nAdjustment;
}

int CDirectionAdjustor::m_fnCurlAdjust()
{
	CMeshBase::CPolyhedron::Vertex_iterator vi;
	CMeshBase::CPolyhedron::Halfedge_handle hh0, hh, hhPath;
	CMeshBase::CPolyhedron::Vertex_handle vh;
	CMeshBase::FT ftMinDot, ftDot, ftWeight, ftMaxWeight;
	int nAdjustment;
	ftMaxWeight = 0.0;
	vh = NULL;
	nAdjustment = 0;
	for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
	{
		if (vi->iDegreeBias != 0)
		{
			ftWeight = vi->vtCurlForce * vi->vtCurlForce;
			if (ftWeight > ftMaxWeight)
			{
				ftMaxWeight = ftWeight;
				vh = &(*vi);
			}
		}
	}
	if (vh != NULL)
	{
		hh0 = vh->halfedge();
		hh = hh0;
		ftMinDot = 0;
		hhPath = NULL;
		do
		{
			if (!hh->is_border_edge())
			{
				ftDot = hh->vtVector * vh->vtCurlForce / hh->ftLen;
				if (ftDot < ftMinDot)
				{
					ftMinDot = ftDot;
					hhPath = hh;
				}
			}
			hh = hh->next()->opposite();
		} while (hh != hh0);
		if (hhPath != NULL)
		{
			if (vh->iDegreeBias < 0)
			{
				hhPath->iSeamType = (hhPath->iSeamType + 5) % 4 - 2;
				hhPath->opposite()->iSeamType = (6 - hhPath->iSeamType) % 4 - 2;
			}
			else
			{
				hhPath->iSeamType = (hhPath->iSeamType + 3) % 4 - 2;
				hhPath->opposite()->iSeamType = (6 - hhPath->iSeamType) % 4 - 2;
			}
			++nAdjustment;
		}
	}
	return nAdjustment;
}


int CDirectionAdjustor::m_fnConstrainedOptimize()
{
	CMeshBase::CPolyhedron::Vertex_iterator vi;
	CMeshBase::CPolyhedron::Edge_iterator ei;
	CMeshBase::CPolyhedron::Halfedge_handle hh0, hh, hhPath;
	CMeshBase::FT ftMinDot, ftDot;
	int nAdjustment, iMerge;
	for (ei = m_mshSurface.edges_begin(); ei != m_mshSurface.edges_end(); ++ei)
	{
		ei->ftForce = 0;
		ei->opposite()->ftForce = 0;
	}
	for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
	{
		if (vi->iDegreeBias != 0 && vi->vtThetaForce * vi->vtCurlForce >= 0.0)
		{
			hh0 = vi->halfedge();
			hh = hh0;
			ftMinDot = 0;
			hhPath = NULL;
			do
			{
				if (hh->prev()->vertex()->iDegreeBias * vi->iDegreeBias < 0)
				{
					ftDot = hh->vtVector * vi->vtThetaForce / hh->ftLen;
					if (ftDot < ftMinDot && ftDot < 0.0)
					{
						ftMinDot = ftDot;
						hhPath = hh;
					}
				}
				hh = hh->next()->opposite();
			} while (hh != hh0);
			if (hhPath == NULL)
			{
				do
				{
					ftDot = hh->vtVector * vi->vtThetaForce / hh->ftLen;
					if (ftDot < ftMinDot)
					{
						ftMinDot = ftDot;
						hhPath = hh;
					}
					hh = hh->next()->opposite();
				} while (hh != hh0);
			}
			if (hhPath != NULL)
			{
				if (vi->iDegreeBias < 0)
				{
					hhPath->ftForce += ftMinDot;
					hhPath->opposite()->ftForce -= ftMinDot;
				}
				else
				{
					hhPath->ftForce -= ftMinDot;
					hhPath->opposite()->ftForce += ftMinDot;
				}
			}
		}
	}
	for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
	{
		vi->iHeapPos = -1;
	}
	nAdjustment = 0;
	for (ei = m_mshSurface.edges_begin(); ei != m_mshSurface.edges_end(); ++ei)
	{
		if (ei->ftForce != 0.0)
		{
			if (fabs(ei->ftForce) > fabs(ei->ftForceBarrier))
			{
				++nAdjustment;
			}
		}
	}
	nAdjustment = 0;
	iMerge = 0;
	for (ei = m_mshSurface.edges_begin(); ei != m_mshSurface.edges_end(); ++ei)
	{
		if (ei->ftForce != 0.0)
		{
			if (fabs(ei->ftForce) > fabs(ei->ftForceBarrier))
			{
				if (ei->prev()->vertex()->iDegreeBias == 0 || ei->vertex()->iDegreeBias == 0)
				{
					ei->ftForceBarrier = ei->ftForce;
					ei->opposite()->ftForceBarrier = ei->opposite()->ftForce;
				}
				else
				{
					ei->ftForceBarrier = 0;
					ei->opposite()->ftForceBarrier = 0;
					++iMerge;
				}
				if (ei->is_border_edge())
				{
					if (ei->is_border())
					{
						hh0 = ei->opposite();
					}
					else
					{
						hh0 = &(*ei);
					}
					if (hh0->ftForce > 0)
					{
						hh0->next()->iSeamType = (hh0->next()->iSeamType + 5) % 4 - 2;
						hh0->next()->opposite()->iSeamType = (6 - hh0->next()->iSeamType) % 4 - 2;
						hh0->prev()->iSeamType = (hh0->prev()->iSeamType + 5) % 4 - 2;
						hh0->prev()->opposite()->iSeamType = (6 - hh0->prev()->iSeamType) % 4 - 2;
					}
					else
					{
						hh0->next()->iSeamType = (hh0->next()->iSeamType + 3) % 4 - 2;
						hh0->next()->opposite()->iSeamType = (6 - hh0->next()->iSeamType) % 4 - 2;
						hh0->prev()->iSeamType = (hh0->prev()->iSeamType + 3) % 4 - 2;
						hh0->prev()->opposite()->iSeamType = (6 - hh0->prev()->iSeamType) % 4 - 2;
					}
				}
				else
				{
					if (ei->ftForce > 0)
					{
						ei->iSeamType = (ei->iSeamType + 3) % 4 - 2;
						ei->opposite()->iSeamType = (6 - ei->iSeamType) % 4 - 2;
					}
					else
					{
						ei->iSeamType = (ei->iSeamType + 5) % 4 - 2;
						ei->opposite()->iSeamType = (6 - ei->iSeamType) % 4 - 2;
					}
				}
				++nAdjustment;
			}
		}
	}
	if (iMerge > 1)
	{
		for (ei = m_mshSurface.edges_begin(); ei != m_mshSurface.edges_end(); ++ei)
		{
			ei->ftForceBarrier = 0;
			ei->opposite()->ftForce = 0;
		}
	}
	return nAdjustment;
}

int CDirectionAdjustor::m_fnMergeAdjacentSingularities()
{
	CMeshBase::CPolyhedron::Edge_iterator ei;
	for (ei = m_mshSurface.edges_begin(); ei != m_mshSurface.edges_end(); ++ei)
	{
		if (ei->vertex()->nNumPorts != 0 && ei->prev()->vertex()->nNumPorts != 0)
		{
			;
		}
	}
	return 0;
}

int CDirectionAdjustor::m_fnNumberOfSingularities()
{
	CMeshBase::CPolyhedron::Vertex_iterator vi;
	int nNum;
	nNum = 0;
	for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
	{
		if (vi->iBordId == -1 && vi->nNumPorts != 0)
		{
			++nNum;
		}
	}
	return nNum;
}

CMeshBase::FT CDirectionAdjustor::m_fnDiscreteCurl()
{
	CMeshBase::CPolyhedron::Edge_iterator ei;
	CMeshBase::CPolyhedron::Halfedge_handle hh0, hh1;
	CMeshBase::CPolyhedron::Vector_3 lpvtChart[4];
	CMeshBase::FT ftAve, ftMax, lpftErr[2];
	int nNum, i;
	ftAve = 0;
	ftMax = 0;
	nNum = 0;
	for (ei = m_mshSurface.edges_begin(); ei != m_mshSurface.edges_end(); ++ei)
	{
		if (!ei->is_border_edge())
		{
			hh0 = &(*ei);
			hh1 = hh0->opposite();
			lpvtChart[0] = (cos(hh0->facet()->ftChartDir) * hh0->facet()->vtPrincipalAxis + sin(hh0->facet()->ftChartDir) * hh0->facet()->vtAuxiliaryAxis) / hh0->facet()->ftRadius;
			lpvtChart[1] = (-sin(hh0->facet()->ftChartDir) * hh0->facet()->vtPrincipalAxis + cos(hh0->facet()->ftChartDir) * hh0->facet()->vtAuxiliaryAxis) / hh0->facet()->ftRadius;
			lpvtChart[2] = (cos(hh1->facet()->ftChartDir + hh1->iSeamType * ftHalfPi) * hh1->facet()->vtPrincipalAxis + sin(hh1->facet()->ftChartDir + hh1->iSeamType * ftHalfPi) * hh1->facet()->vtAuxiliaryAxis) / hh1->facet()->ftRadius;
			lpvtChart[3] = (-sin(hh1->facet()->ftChartDir + hh1->iSeamType * ftHalfPi) * hh1->facet()->vtPrincipalAxis + cos(hh1->facet()->ftChartDir + hh1->iSeamType * ftHalfPi) * hh1->facet()->vtAuxiliaryAxis) / hh1->facet()->ftRadius;
			lpftErr[0] = fabs((lpvtChart[0] * hh0->vtVector - lpvtChart[2] * hh0->vtVector) / hh0->ftLen);
			lpftErr[1] = fabs((lpvtChart[3] * hh0->vtVector - lpvtChart[1] * hh0->vtVector) / hh0->ftLen);
			for (i = 0; i < 2; ++i)
			{
				ftAve += lpftErr[i];
				if (lpftErr[i] > ftMax)
				{
					ftMax = lpftErr[i];
				}
			}
			++nNum;
		}
	}
	ftAve /= (nNum * 2);
	return ftAve;
}
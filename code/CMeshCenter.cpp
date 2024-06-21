#include "CMeshCenter.h"

CMeshCenter::CMeshCenter()
{
    //ctor
    //m_psProcessStatus = PS_EMPTY;
    m_iDensityMaximum = -1;
    m_iDensityMinimum = -1;
	m_bExtremaSelecting = false;
	m_lpiGlobalIndex = NULL;
	m_lpftGlobalPara = NULL;
	m_nGlobalParameterSize = 0;
	m_hhDirEdge = NULL;
	m_fhFacetFix = NULL;
}

CMeshCenter::~CMeshCenter()
{
}




void CMeshCenter::m_fnClearBarrier()
{
    CMeshBase::CPolyhedron::Edge_iterator ei;
	for (ei = m_mshSurface.edges_begin(); ei != m_mshSurface.edges_end(); ++ei)
	{
		ei->ftForceBarrier = 0;
		ei->opposite()->ftForceBarrier = 0;
		ei->iTempIndex = 0;
		ei->opposite()->iTempIndex = 0;
	}
}

void CMeshCenter::m_fnBackupField()
{
	CMeshBase::CPolyhedron::Edge_iterator ei;
	CMeshBase::CPolyhedron::Facet_iterator fi;
	for (ei = m_mshSurface.edges_begin(); ei != m_mshSurface.edges_end(); ++ei)
	{
		ei->iSeamTypeCopy = ei->iSeamType;
		ei->opposite()->iSeamTypeCopy = ei->opposite()->iSeamType;
	}
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
	{
		fi->ftChartDirCopy = fi->ftChartDir;
	}
}

void CMeshCenter::m_fnRestoreField()
{
	CMeshBase::CPolyhedron::Edge_iterator ei;
	CMeshBase::CPolyhedron::Facet_iterator fi;
	for (ei = m_mshSurface.edges_begin(); ei != m_mshSurface.edges_end(); ++ei)
	{
		ei->iSeamType = ei->iSeamTypeCopy;
		ei->opposite()->iSeamType = ei->opposite()->iSeamTypeCopy;
	}
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
	{
		fi->ftChartDir = fi->ftChartDirCopy;
	}
}


void CMeshCenter::m_fnClearCrevasses()
{
    std::vector<CCrevasse>::iterator iterCrevasse;
    if(m_vecCrevasses.size() > 0)
    {
        for(iterCrevasse = m_vecCrevasses.begin(); iterCrevasse != m_vecCrevasses.end(); ++iterCrevasse)
        {
            if(iterCrevasse->m_vechhElements.size() > 0)
            {
                iterCrevasse->m_vechhElements.clear();
            }
        }
        m_vecCrevasses.clear();
    }
}

void CMeshCenter::m_fnClearGlobalPara()
{
	CMeshBase::CPolyhedron::Halfedge_iterator hi;
	if (m_lpftGlobalPara != NULL)
	{
		delete[]m_lpftGlobalPara;
		m_lpftGlobalPara = NULL;
	}
	if (m_lpiGlobalIndex != NULL)
	{
		delete[]m_lpiGlobalIndex;
	}
	for (hi = m_mshSurface.halfedges_begin(); hi != m_mshSurface.halfedges_end(); ++hi)
	{
		hi->lpftGlobalPara = NULL;
		hi->lpiGlobalIndex = NULL;
	}
}


void CMeshCenter::m_fnClearFreeEdge()
{
    CMeshBase::CPolyhedron::Halfedge_iterator hi;
    for(hi = m_mshSurface.halfedges_begin(); hi != m_mshSurface.halfedges_end(); ++hi)
    {
        hi->bFreeEdge = false;
    }
}

void CMeshCenter::m_fnClearFixedDirection()
{
	CMeshBase::CPolyhedron::Facet_iterator fi;
	CMeshBase::CPolyhedron::Halfedge_handle hh0, hh;
	bool bClear;
	//m_lstFixedDirections.clear();
	//m_lstFixedFacets.clear();
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
	{
		if (fi->bDirFixed)
		{
			bClear = true;
			hh0 = fi->halfedge();
			hh = hh0;
			do
			{
				if (hh->opposite()->is_border())
				{
					bClear = false;
				}
				hh = hh->next();
			} while (hh != hh0);
			if (bClear)
			{
				//m_lstFixedFacets.push_back(&(*fi));
				//m_lstFixedDirections.push_back(fi->ftChartDir);
				fi->bDirFixed = false;
			}
		}
	}
}

void CMeshCenter::m_fnDelaunayFlip()
{
	CMeshBase::Kernel::FT ftSquaredNormLen2, ftSquaredNormLen4, ftCross1, ftCross2, ftCross3, ftCross4, ftDot1, ftDot2, ftDot3, ftDot4;
	CMeshBase::Vector_3 vtNormal2, vtNormal4;
	bool bComplete, bSplit, bEdgeFlip, bMeshSplit;
	CMeshBase::CPolyhedron::Vertex_handle vh1, vh2;
	CMeshBase::CPolyhedron::Halfedge_handle hh0, hh1, hhCandidate;
	CMeshBase::CPolyhedron::Edge_iterator ei;
	std::vector<CMeshBase::CPolyhedron::Halfedge_handle> vecEdgeStack;
	bComplete = false;
	bEdgeFlip = false;
	bMeshSplit = false;
	while (!bComplete)
	{
		bComplete = true;
		for (ei = m_mshSurface.edges_begin(); ei != m_mshSurface.edges_end(); ++ei)
		{
			if (!(ei->is_border_edge()))
			{
				ftDot1 = -ei->next()->ftDot;
				ftDot2 = -(ei->prev()->vtVector) * (ei->opposite()->next()->vtVector);
				ftDot3 = -ei->opposite()->next()->ftDot;
				ftDot4 = -(ei->opposite()->prev()->vtVector) * (ei->next()->vtVector);
				vtNormal2 = CGAL::cross_product(ei->prev()->vtVector, ei->opposite()->next()->vtVector);
				vtNormal4 = CGAL::cross_product(ei->opposite()->prev()->vtVector, ei->next()->vtVector);
				ftSquaredNormLen2 = vtNormal2 * vtNormal2;
				ftSquaredNormLen4 = vtNormal4 * vtNormal4;
				ftCross1 = ei->facet()->ftTwArea;
				ftCross2 = sqrt(ftSquaredNormLen2);
				ftCross3 = ei->opposite()->facet()->ftTwArea;
				ftCross4 = sqrt(ftSquaredNormLen4);
				if (atan2(ftCross1, ftDot1) + atan2(ftCross3, ftDot3) > ftPi)//(atan2(ftCross1, ftDot1) + atan2(ftCross3, ftDot3) > atan2(ftCross2, ftDot2) + atan2(ftCross4, ftDot4))
				{
					bEdgeFlip = true;
					if (ei->bInStack == false)
					{
						ei->bInStack = true;
						ei->opposite()->bInStack = true;
						vecEdgeStack.push_back(&(*ei));
					}
				}
			}
		}
		while (vecEdgeStack.size() != 0)
		{
			hhCandidate = *vecEdgeStack.rbegin();
			vecEdgeStack.pop_back();
			hhCandidate->bInStack = false;
			hhCandidate->opposite()->bInStack = false;
			ftDot1 = -hhCandidate->next()->ftDot;
			ftDot2 = -(hhCandidate->prev()->vtVector) * (hhCandidate->opposite()->next()->vtVector);
			ftDot3 = -hhCandidate->opposite()->next()->ftDot;
			ftDot4 = -(hhCandidate->opposite()->prev()->vtVector) * (hhCandidate->next()->vtVector);
			vtNormal2 = CGAL::cross_product(hhCandidate->prev()->vtVector, hhCandidate->opposite()->next()->vtVector);
			vtNormal4 = CGAL::cross_product(hhCandidate->opposite()->prev()->vtVector, hhCandidate->next()->vtVector);
			ftSquaredNormLen2 = (vtNormal2) * (vtNormal2);
			ftSquaredNormLen4 = (vtNormal4) * (vtNormal4);
			ftCross1 = hhCandidate->facet()->ftTwArea;
			ftCross2 = sqrt(ftSquaredNormLen2);
			ftCross3 = hhCandidate->opposite()->facet()->ftTwArea;
			ftCross4 = sqrt(ftSquaredNormLen4);
			if (atan2(ftCross1, ftDot1) + atan2(ftCross3, ftDot3) > ftPi)//(atan2(ftCross1, ftDot1) + atan2(ftCross3, ftDot3) > atan2(ftCross2, ftDot2) + atan2(ftCross4, ftDot4))
			{
				bSplit = false;
				vh1 = hhCandidate->next()->vertex();
				vh2 = hhCandidate->opposite()->next()->vertex();
				hh0 = vh1->halfedge()->opposite();
				hh1 = hh0;
				do
				{
					if (hh1->vertex() == vh2)
					{
						bSplit = true;
						break;
					}
					hh1 = hh1->opposite()->next();
				} while (hh1 != hh0);
				if (bSplit)
				{
				}
				else
				{
					m_mshSurface.flip_edge(&(*hhCandidate));
					hhCandidate->vtVector = hhCandidate->vertex()->point() - hhCandidate->opposite()->vertex()->point();
					hhCandidate->opposite()->vtVector = -hhCandidate->vtVector;
					hhCandidate->ftSqLen = hhCandidate->vtVector * hhCandidate->vtVector;
					hhCandidate->opposite()->ftSqLen = hhCandidate->ftSqLen;
					hhCandidate->facet()->vtNorm = vtNormal2;
					hhCandidate->opposite()->facet()->vtNorm = vtNormal4;
					hhCandidate->facet()->ftTwArea = ftSquaredNormLen2;
					hhCandidate->facet()->ftTwArea = ftCross2;
					hhCandidate->opposite()->facet()->ftTwArea = ftSquaredNormLen4;
					hhCandidate->opposite()->facet()->ftTwArea = ftCross4;
					hhCandidate->next()->ftDot = -ftDot2;
					hhCandidate->opposite()->next()->ftDot = -ftDot4;
					hhCandidate->ftDot = (hhCandidate->vtVector) * (hhCandidate->next()->vtVector);
					hhCandidate->prev()->ftDot = (hhCandidate->prev()->vtVector) * (hhCandidate->vtVector);
					hhCandidate->opposite()->ftDot = (hhCandidate->opposite()->vtVector) * (hhCandidate->opposite()->next()->vtVector);
					hhCandidate->opposite()->prev()->ftDot = (hhCandidate->opposite()->prev()->vtVector) * (hhCandidate->opposite()->vtVector);
					if (!hhCandidate->prev()->bInStack && !hhCandidate->prev()->is_border_edge())
					{
						hhCandidate->prev()->bInStack = true;
						hhCandidate->prev()->opposite()->bInStack = true;
						vecEdgeStack.push_back(hhCandidate->prev());
					}
					if (!hhCandidate->next()->bInStack && !hhCandidate->next()->is_border_edge())
					{
						hhCandidate->next()->bInStack = true;
						hhCandidate->next()->opposite()->bInStack = true;
						vecEdgeStack.push_back(hhCandidate->next());
					}
					if (!hhCandidate->opposite()->prev()->bInStack && !hhCandidate->opposite()->prev()->is_border_edge())
					{
						hhCandidate->opposite()->prev()->bInStack = true;
						hhCandidate->opposite()->prev()->opposite()->bInStack = true;
						vecEdgeStack.push_back(hhCandidate->opposite()->prev());
					}
					if (!hhCandidate->opposite()->next()->bInStack && !hhCandidate->opposite()->next()->is_border_edge())
					{
						hhCandidate->opposite()->next()->bInStack = true;
						hhCandidate->opposite()->next()->opposite()->bInStack = true;
						vecEdgeStack.push_back(hhCandidate->opposite()->next());
					}
				}
			}
		}
	}
	return;
}

void CMeshCenter::m_fnGenerateMoment()
{
    CMeshBase::CPolyhedron::Edge_iterator ei;
//    CMeshBase::CPolyhedron::Halfedge_handle hh0, hh;
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
            ftAngleDif = ei->opposite()->facet()->ftChartDir - ei->facet()->ftChartDir + ei->ftPolarAxisDif - ei->iSeamType * ftHalfPi;
            ftAngleDif = ftAngleDif - CMeshBase::FT(int(ftAngleDif / ftTwoPi + 15.5) - 15) * ftTwoPi;
            ei->ftMom = ftAngleDif;
            ei->opposite()->ftMom = -ei->ftMom;
        }
    }
}


void CMeshCenter::m_fnSolveParaScale()
{
	CMeshBase::CPolyhedron::Facet_iterator fi;
	CMeshBase::CPolyhedron::Edge_iterator ei;
	MatElem<CMeshBase::FT> me;
	CLeastSquareSolver<CMeshBase::FT> LSS;
	CMeshBase::Vector_3 lpvtChart[4];
	CMeshBase::CPolyhedron::Facet_handle fh;
	int iIndex, i, iSeamType;
	iIndex = 0;
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
	{
		fi->iEquationIndex = iIndex;
		++iIndex;
	}
	--fi;
	fi->iEquationIndex = -1;
	fi->ftLocalParaScl = 1.0;
	iIndex = 0;
	for (ei = m_mshSurface.edges_begin(); ei != m_mshSurface.edges_end(); ++ei)
	{
		if (!ei->is_border_edge())
		{
			iIndex += 2;
		}
	}
	LSS.init(iIndex, m_mshSurface.size_of_facets() - 1);
	for (i = 0; i < iIndex; ++i)
	{
		LSS.m_b[i] = 0;
	}
	me.row = 0;
	for (ei = m_mshSurface.edges_begin(); ei != m_mshSurface.edges_end(); ++ei)
	{
		if (!ei->is_border_edge())
		{
			fh = ei->facet();
			lpvtChart[0] = (fh->vtPrincipalAxis * cos(fh->ftChartDir) + fh->vtAuxiliaryAxis * sin(fh->ftChartDir)) / fh->ftRadius;
			lpvtChart[1] = (fh->vtPrincipalAxis * (-sin(fh->ftChartDir)) + fh->vtAuxiliaryAxis * cos(fh->ftChartDir)) / fh->ftRadius;
			fh = ei->opposite()->facet();
			iSeamType = ei->iSeamType;
			lpvtChart[2] = (fh->vtPrincipalAxis * cos(fh->ftChartDir - iSeamType * ftHalfPi) + fh->vtAuxiliaryAxis * sin(fh->ftChartDir - iSeamType * ftHalfPi)) / fh->ftRadius;
			lpvtChart[3] = (fh->vtPrincipalAxis * (-sin(fh->ftChartDir - iSeamType * ftHalfPi)) + fh->vtAuxiliaryAxis * cos(fh->ftChartDir - iSeamType * ftHalfPi)) / fh->ftRadius;
			for (i = 0; i < 2; ++i)
			{
				if (ei->facet()->iEquationIndex == -1)
				{
					LSS.m_b[me.row] = -ei->facet()->ftLocalParaScl * lpvtChart[i] * ei->vtVector / ei->ftLen;
				}
				else
				{
					me.col = ei->facet()->iEquationIndex;
					me.data = lpvtChart[i] * ei->vtVector / ei->ftLen;
					LSS.add_coef(me.row, me.col, me.data);
				}
				if (ei->opposite()->facet()->iEquationIndex == -1)
				{
					LSS.m_b[me.row] = ei->opposite()->facet()->ftLocalParaScl * lpvtChart[i + 2] * ei->vtVector / ei->ftLen;
				}
				else
				{
					me.col = ei->opposite()->facet()->iEquationIndex;
					me.data = -lpvtChart[i + 2] * ei->vtVector / ei->ftLen;
					LSS.add_coef(me.row, me.col, me.data);
				}
				++me.row;
			}
		}
	}
	LSS.decompose();
	LSS.solve();
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
	{
		if (fi->iEquationIndex != -1)
		{
			fi->ftLocalParaScl = LSS.m_x[fi->iEquationIndex];
		}
	}
	LSS.clear();
}

void CMeshCenter::m_fnAdjustParaScale()
{
	CMeshBase::CPolyhedron::Facet_iterator fi;
	CMeshBase::CPolyhedron::Edge_iterator ei;
	MatElem<CMeshBase::FT> me;
	CLeastSquareSolver<CMeshBase::FT> LSS;
	CMeshBase::Vector_3 lpvtChart[4];
	CMeshBase::CPolyhedron::Facet_handle fh;
	int iIndex, i, iSeamType;
	CMeshBase::FT ftScale;
	iIndex = 0;
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
	{
		fi->iEquationIndex = iIndex;
		++iIndex;
	}
	--fi;
	fi->iEquationIndex = -1;
	fi->ftLocalParaScl = 1.0;
	iIndex = 0;
	for (ei = m_mshSurface.edges_begin(); ei != m_mshSurface.edges_end(); ++ei)
	{
		if (!ei->is_border_edge())
		{
			iIndex += 2;
		}
	}
	LSS.init(iIndex, m_mshSurface.size_of_facets() - 1);
	for (i = 0; i < iIndex; ++i)
	{
		LSS.m_b[i] = 0;
	}
	me.row = 0;
	for (ei = m_mshSurface.edges_begin(); ei != m_mshSurface.edges_end(); ++ei)
	{
		if (!ei->is_border_edge())
		{
			fh = ei->facet();
			ftScale = fh->ftLocalParaScl;
			lpvtChart[0] = (fh->vtPrincipalAxis * cos(fh->ftChartDir) + fh->vtAuxiliaryAxis * sin(fh->ftChartDir)) / fh->ftRadius;
			lpvtChart[1] = (fh->vtPrincipalAxis * (-sin(fh->ftChartDir)) + fh->vtAuxiliaryAxis * cos(fh->ftChartDir)) / fh->ftRadius;
			fh = ei->opposite()->facet();
			ftScale = fh->ftLocalParaScl;
			ftScale /= 2;
			iSeamType = ei->iSeamType;
			lpvtChart[2] = (fh->vtPrincipalAxis * cos(fh->ftChartDir - iSeamType * ftHalfPi) + fh->vtAuxiliaryAxis * sin(fh->ftChartDir - iSeamType * ftHalfPi)) / fh->ftRadius;
			lpvtChart[3] = (fh->vtPrincipalAxis * (-sin(fh->ftChartDir - iSeamType * ftHalfPi)) + fh->vtAuxiliaryAxis * cos(fh->ftChartDir - iSeamType * ftHalfPi)) / fh->ftRadius;
			for (i = 0; i < 2; ++i)
			{
				if (ei->facet()->iEquationIndex == -1)
				{
					LSS.m_b[me.row] = -ei->facet()->ftLocalParaScl * lpvtChart[i] * ei->vtVector / (ei->ftLen * ftScale);
				}
				else
				{
					me.col = ei->facet()->iEquationIndex;
					me.data = lpvtChart[i] * ei->vtVector / (ei->ftLen * ftScale);
					LSS.add_coef(me.row, me.col, me.data);
				}
				if (ei->opposite()->facet()->iEquationIndex == -1)
				{
					LSS.m_b[me.row] = ei->opposite()->facet()->ftLocalParaScl * lpvtChart[i + 2] * ei->vtVector / (ei->ftLen * ftScale);
				}
				else
				{
					me.col = ei->opposite()->facet()->iEquationIndex;
					me.data = -lpvtChart[i + 2] * ei->vtVector / (ei->ftLen * ftScale);
					LSS.add_coef(me.row, me.col, me.data);
				}
				++me.row;
			}
		}
	}
	LSS.decompose();
	LSS.solve();
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
	{
		if (fi->iEquationIndex != -1)
		{
			fi->ftLocalParaScl = LSS.m_x[fi->iEquationIndex];
		}
	}
	LSS.clear();
}

CMeshBase::FT CMeshCenter::m_fnMIQEnergy()
{
	CMeshBase::CPolyhedron::Edge_iterator ei;
	CMeshBase::FT ftAngleDif, ftEnergy;
	ftEnergy = 0.0;
	for (ei = m_mshSurface.edges_begin(); ei != m_mshSurface.edges_end(); ++ei)
	{
		if (ei->is_border_edge())
		{
			ei->ftMom = 0;
			ei->opposite()->ftMom = 0;
		}
		else
		{
			ftAngleDif = ei->opposite()->facet()->ftChartDir - ei->facet()->ftChartDir + ei->ftPolarAxisDif - ei->iSeamType * ftHalfPi;
			ftAngleDif = ftAngleDif - CMeshBase::FT(int(ftAngleDif / ftTwoPi + 15.5) - 15) * ftTwoPi;
			ftEnergy = ftAngleDif * ftAngleDif;
		}
	}
	return ftEnergy;
}


void CMeshCenter::m_fnEstimateAngleError(CMeshBase::FT *lpftErr)
{
	CMeshBase::CPolyhedron::Facet_iterator fi;
	CMeshBase::CPolyhedron::Halfedge_handle hh0, hh;
	CMeshBase::Vector_3 lpvtGrad[2], lpvtChart[4], lpvtPartial[2];
	CMeshBase::FT lpftTmpErr[2], ftMatchErr, ftCorrectMatch, ftMeanErr;
	int i, iDir;
	lpftErr[0] = lpftErr[3] = lpftErr[6] = ftHalfPi;
	lpftErr[1] = lpftErr[4] = lpftErr[7] = 0;
	lpftErr[2] = lpftErr[5] = lpftErr[8] = 0;
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
	{
		lpvtGrad[1] = lpvtGrad[0] = CMeshBase::Vector_3(0, 0, 0);
		hh0 = fi->halfedge();
		hh = hh0;
		do
		{
			lpvtGrad[0] = lpvtGrad[0] + hh->vtVector * hh->next()->lpftGlobalPara[0];
			lpvtGrad[1] = lpvtGrad[1] + hh->vtVector * hh->next()->lpftGlobalPara[1];
			hh = hh->next();
		} while (hh != hh0);
		lpvtPartial[0] = (lpvtGrad[1] * lpvtGrad[1]) * lpvtGrad[0] - (lpvtGrad[0] * lpvtGrad[1]) * lpvtGrad[1];
		lpvtPartial[1] = (lpvtGrad[0] * lpvtGrad[0]) * lpvtGrad[1] - (lpvtGrad[0] * lpvtGrad[1]) * lpvtGrad[0];
		lpvtChart[0] = cos(fi->ftChartDir) * fi->vtPrincipalAxis + sin(fi->ftChartDir) * fi->vtAuxiliaryAxis;
		lpvtChart[1] = cos(fi->ftFrameDir) * fi->vtPrincipalAxis + sin(fi->ftFrameDir) * fi->vtAuxiliaryAxis;
		lpvtChart[2] = -lpvtChart[0];
		lpvtChart[3] = -lpvtChart[1];
		ftCorrectMatch = ftPi * ftPi * 2;
		for (i = 0; i < 4; ++i)
		{
			lpftTmpErr[0] = atan2(sqrt(CGAL::squared_length(CGAL::cross_product(lpvtGrad[0], lpvtChart[i]))), (lpvtGrad[0] * lpvtChart[i]));
			lpftTmpErr[1] = atan2(sqrt(CGAL::squared_length(CGAL::cross_product(lpvtGrad[1], lpvtChart[(i + 1) % 4]))), (lpvtGrad[1] * lpvtChart[(i + 1) % 4]));
			ftMatchErr = lpftTmpErr[0] * lpftTmpErr[0] + lpftTmpErr[1] * lpftTmpErr[1];
			if (ftMatchErr < ftCorrectMatch)
			{
				iDir = i;
				ftCorrectMatch = ftMatchErr;
				fi->lpftAngleErr[0] = fabs(lpftTmpErr[0]);
				fi->lpftAngleErr[1] = fabs(lpftTmpErr[1]);
			}
		}
		for (i = 0; i < 4; ++i)
		{
			lpftTmpErr[0] = atan2(sqrt(CGAL::squared_length(CGAL::cross_product(lpvtGrad[0], lpvtChart[i]))), (lpvtGrad[0] * lpvtChart[i]));
			lpftTmpErr[1] = atan2(sqrt(CGAL::squared_length(CGAL::cross_product(lpvtGrad[1], lpvtChart[(i + 3) % 4]))), (lpvtGrad[1] * lpvtChart[(i + 3) % 4]));
			ftMatchErr = lpftTmpErr[0] * lpftTmpErr[0] + lpftTmpErr[1] * lpftTmpErr[1];
			if (ftMatchErr < ftCorrectMatch)
			{
				iDir = i;
				ftCorrectMatch = ftMatchErr;
				fi->lpftAngleErr[0] = fabs(lpftTmpErr[0]);
				fi->lpftAngleErr[1] = fabs(lpftTmpErr[1]);
			}
		}
		for (i = 0; i < 4; ++i)
		{
			lpftTmpErr[0] = atan2(sqrt(CGAL::squared_length(CGAL::cross_product(lpvtPartial[0], lpvtChart[i]))), (lpvtPartial[0] * lpvtChart[i]));
			lpftTmpErr[1] = atan2(sqrt(CGAL::squared_length(CGAL::cross_product(lpvtPartial[1], lpvtChart[(i + 1) % 4]))), (lpvtPartial[1] * lpvtChart[(i + 1) % 4]));
			ftMatchErr = lpftTmpErr[0] * lpftTmpErr[0] + lpftTmpErr[1] * lpftTmpErr[1];
			if (ftMatchErr < ftCorrectMatch)
			{
				iDir = i;
				ftCorrectMatch = ftMatchErr;
				fi->lpftAngleErr[0] = fabs(lpftTmpErr[0]);
				fi->lpftAngleErr[1] = fabs(lpftTmpErr[1]);
			}
		}
		for (i = 0; i < 4; ++i)
		{
			lpftTmpErr[0] = atan2(sqrt(CGAL::squared_length(CGAL::cross_product(lpvtPartial[0], lpvtChart[i]))), (lpvtPartial[0] * lpvtChart[i]));
			lpftTmpErr[1] = atan2(sqrt(CGAL::squared_length(CGAL::cross_product(lpvtPartial[1], lpvtChart[(i + 3) % 4]))), (lpvtPartial[1] * lpvtChart[(i + 3) % 4]));
			ftMatchErr = lpftTmpErr[0] * lpftTmpErr[0] + lpftTmpErr[1] * lpftTmpErr[1];
			if (ftMatchErr < ftCorrectMatch)
			{
				iDir = i;
				ftCorrectMatch = ftMatchErr;
				fi->lpftAngleErr[0] = fabs(lpftTmpErr[0]);
				fi->lpftAngleErr[1] = fabs(lpftTmpErr[1]);
			}
		}
		ftMeanErr = (fi->lpftAngleErr[0] + fi->lpftAngleErr[1]) / 2;
		if (fi->lpftAngleErr[0] < lpftErr[0])
		{
			lpftErr[0] = fi->lpftAngleErr[0];
		}
		if (fi->lpftAngleErr[0] > lpftErr[2])
		{
			lpftErr[2] = fi->lpftAngleErr[0];
		}
		if (fi->lpftAngleErr[1] < lpftErr[3])
		{
			lpftErr[3] = fi->lpftAngleErr[1];
		}
		if (fi->lpftAngleErr[1] > lpftErr[5])
		{
			lpftErr[5] = fi->lpftAngleErr[1];
		}
		if (ftMeanErr < lpftErr[6])
		{
			lpftErr[6] = ftMeanErr;
		}
		if (ftMeanErr > lpftErr[8])
		{
			lpftErr[8] = ftMeanErr;
		}
		lpftErr[1] += fi->lpftAngleErr[0];
		lpftErr[4] += fi->lpftAngleErr[1];
		lpftErr[7] += ftMeanErr;
	}
	lpftErr[1] /= m_mshSurface.size_of_facets();
	lpftErr[4] /= m_mshSurface.size_of_facets();
	lpftErr[7] /= m_mshSurface.size_of_facets();
}

CMeshBase::FT CMeshCenter::m_fnAlignEnergy()
{
	CMeshBase::CPolyhedron::Facet_iterator fi;
	CMeshBase::CPolyhedron::Halfedge_handle hh0, hh;
	CMeshBase::Vector_3 lpvtGrad[2], lpvtChart[4], lpvtPartial[2];
	CMeshBase::FT lpftTmpErr[2], ftMatchErr, ftCorrectMatch, ftAlignEnergy, ftSin, ftSq1, ftSq2;
	int i, iDir;
	ftAlignEnergy = 0;
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
	{
		lpvtGrad[1] = lpvtGrad[0] = CMeshBase::Vector_3(0, 0, 0);
		hh0 = fi->halfedge();
		hh = hh0;
		do
		{
			lpvtGrad[0] = lpvtGrad[0] + hh->vtVector * hh->next()->lpftGlobalPara[0];
			lpvtGrad[1] = lpvtGrad[1] + hh->vtVector * hh->next()->lpftGlobalPara[1];
			hh = hh->next();
		} while (hh != hh0);
		lpvtPartial[0] = (lpvtGrad[1] * lpvtGrad[1]) * lpvtGrad[0] - (lpvtGrad[0] * lpvtGrad[1]) * lpvtGrad[1];
		lpvtPartial[1] = (lpvtGrad[0] * lpvtGrad[0]) * lpvtGrad[1] - (lpvtGrad[0] * lpvtGrad[1]) * lpvtGrad[0];
		lpvtChart[0] = cos(fi->ftChartDir) * fi->vtPrincipalAxis + sin(fi->ftChartDir) * fi->vtAuxiliaryAxis;
		lpvtChart[1] = cos(fi->ftFrameDir) * fi->vtPrincipalAxis + sin(fi->ftFrameDir) * fi->vtAuxiliaryAxis;
		lpvtChart[2] = -lpvtChart[0];
		lpvtChart[3] = -lpvtChart[1];
		ftCorrectMatch = ftPi * ftPi * 2;
		for (i = 0; i < 4; ++i)
		{
			lpftTmpErr[0] = atan2(sqrt(CGAL::squared_length(CGAL::cross_product(lpvtGrad[0], lpvtChart[i]))), (lpvtGrad[0] * lpvtChart[i]));
			lpftTmpErr[1] = atan2(sqrt(CGAL::squared_length(CGAL::cross_product(lpvtGrad[1], lpvtChart[(i + 1) % 4]))), (lpvtGrad[1] * lpvtChart[(i + 1) % 4]));
			ftMatchErr = lpftTmpErr[0] * lpftTmpErr[0] + lpftTmpErr[1] * lpftTmpErr[1];
			if (ftMatchErr < ftCorrectMatch)
			{
				iDir = i;
				ftCorrectMatch = ftMatchErr;
				fi->lpftAngleErr[0] = fabs(lpftTmpErr[0]);
				fi->lpftAngleErr[1] = fabs(lpftTmpErr[1]);
				ftSin = sin(fi->lpftAngleErr[0]);
				ftSq1 = ftSin * ftSin;
				ftSin = sin(fi->lpftAngleErr[1]);
				ftSq2 = ftSin * ftSin;
			}
		}
		for (i = 0; i < 4; ++i)
		{
			lpftTmpErr[0] = atan2(sqrt(CGAL::squared_length(CGAL::cross_product(lpvtGrad[0], lpvtChart[i]))), (lpvtGrad[0] * lpvtChart[i]));
			lpftTmpErr[1] = atan2(sqrt(CGAL::squared_length(CGAL::cross_product(lpvtGrad[1], lpvtChart[(i + 3) % 4]))), (lpvtGrad[1] * lpvtChart[(i + 3) % 4]));
			ftMatchErr = lpftTmpErr[0] * lpftTmpErr[0] + lpftTmpErr[1] * lpftTmpErr[1];
			if (ftMatchErr < ftCorrectMatch)
			{
				iDir = i;
				ftCorrectMatch = ftMatchErr;
				fi->lpftAngleErr[0] = fabs(lpftTmpErr[0]);
				fi->lpftAngleErr[1] = fabs(lpftTmpErr[1]);
				ftSin = sin(fi->lpftAngleErr[0]);
				ftSq1 = ftSin * ftSin;
				ftSin = sin(fi->lpftAngleErr[1]);
				ftSq2 = ftSin * ftSin;
			}
		}
		for (i = 0; i < 4; ++i)
		{
			lpftTmpErr[0] = atan2(sqrt(CGAL::squared_length(CGAL::cross_product(lpvtPartial[0], lpvtChart[i]))), (lpvtPartial[0] * lpvtChart[i]));
			lpftTmpErr[1] = atan2(sqrt(CGAL::squared_length(CGAL::cross_product(lpvtPartial[1], lpvtChart[(i + 1) % 4]))), (lpvtPartial[1] * lpvtChart[(i + 1) % 4]));
			ftMatchErr = lpftTmpErr[0] * lpftTmpErr[0] + lpftTmpErr[1] * lpftTmpErr[1];
			if (ftMatchErr < ftCorrectMatch)
			{
				iDir = i;
				ftCorrectMatch = ftMatchErr;
				fi->lpftAngleErr[0] = fabs(lpftTmpErr[0]);
				fi->lpftAngleErr[1] = fabs(lpftTmpErr[1]);
				ftSin = sin(fi->lpftAngleErr[0]);
				ftSq1 = ftSin * ftSin;
				ftSin = sin(fi->lpftAngleErr[1]);
				ftSq2 = ftSin * ftSin;
			}
		}
		for (i = 0; i < 4; ++i)
		{
			lpftTmpErr[0] = atan2(sqrt(CGAL::squared_length(CGAL::cross_product(lpvtPartial[0], lpvtChart[i]))), (lpvtPartial[0] * lpvtChart[i]));
			lpftTmpErr[1] = atan2(sqrt(CGAL::squared_length(CGAL::cross_product(lpvtPartial[1], lpvtChart[(i + 3) % 4]))), (lpvtPartial[1] * lpvtChart[(i + 3) % 4]));
			ftMatchErr = lpftTmpErr[0] * lpftTmpErr[0] + lpftTmpErr[1] * lpftTmpErr[1];
			if (ftMatchErr < ftCorrectMatch)
			{
				iDir = i;
				ftCorrectMatch = ftMatchErr;
				fi->lpftAngleErr[0] = fabs(lpftTmpErr[0]);
				fi->lpftAngleErr[1] = fabs(lpftTmpErr[1]);
				ftSin = sin(fi->lpftAngleErr[0]);
				ftSq1 = ftSin * ftSin;
				ftSin = sin(fi->lpftAngleErr[1]);
				ftSq2 = ftSin * ftSin;
			}
		}
		ftAlignEnergy += (ftSq1 + ftSq2);
	}
	return ftAlignEnergy;
}

void CMeshCenter::m_fnEvaluatePara(CMeshBase::FT &ftGrids, CMeshBase::FT &ftAreaDistortion, CMeshBase::FT &ftShapeDistortion)
{
	CMeshBase::CPolyhedron::Facet_iterator fi;
	CMeshBase::CPolyhedron::Halfedge_handle hh0, hh;
	CMeshBase::FT ftArea, ftLocalDistortion, *lpftLocalParaArea, ftAreaRatio, lpftSqGradLen[2], ftGradDot, ftTrace, ftCross;
	CMeshBase::Vector_3 lpvtPara[3], lpvtParaGrad[2];
	int i, j, k;
	ftGrids = 0;
	ftArea = 0;
	ftShapeDistortion = 0;
	lpftLocalParaArea = new CMeshBase::FT[m_mshSurface.size_of_facets()];
	i = 0;
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
	{
		ftArea += fi->ftTwArea;
		hh0 = fi->halfedge();
		hh = hh0;
		for (j = 0; j < 3; ++j)
		{
			lpvtPara[j] = CMeshBase::Vector_3(hh->lpftGlobalPara[0], hh->lpftGlobalPara[1], 0);
			hh = hh->next();
		}
		lpftLocalParaArea[i] = (lpvtPara[1].x() - lpvtPara[0].x()) * (lpvtPara[2].y() - lpvtPara[0].y()) - (lpvtPara[1].y() - lpvtPara[0].y()) * (lpvtPara[2].x() - lpvtPara[0].x());
		if (fabs(lpftLocalParaArea[i]) > 1.0e-8)
		{
			ftGrids += lpftLocalParaArea[i];
		}
		else
		{
			lpftLocalParaArea[i] = 0;
		}
		for (j = 0; j < 2; ++j)
		{
			lpvtParaGrad[j] = CMeshBase::Vector_3(0, 0, 0);
			hh = hh0;
			for (k = 0; k < 3; ++k)
			{
				lpvtParaGrad[j] = lpvtParaGrad[j] + hh->lpftGlobalPara[j] * hh->prev()->vtVector;
				hh = hh->next();
			}
		}
		lpftSqGradLen[0] = lpvtParaGrad[0] * lpvtParaGrad[0];
		lpftSqGradLen[1] = lpvtParaGrad[1] * lpvtParaGrad[1];
		ftGradDot = lpvtParaGrad[0] * lpvtParaGrad[1];
		ftTrace = lpftSqGradLen[0] + lpftSqGradLen[1];
		ftCross = lpftSqGradLen[0] * lpftSqGradLen[1] - ftGradDot * ftGradDot;
		if (ftCross > 0)
		{
			ftLocalDistortion = (ftTrace * ftTrace) / (2 * ftCross) - 1;
		}
		else
		{
			ftLocalDistortion = 0;
		}
		if (ftLocalDistortion < 2)
		{
			ftShapeDistortion += (ftLocalDistortion * fi->ftTwArea);
		}
		++i;
	}
	ftShapeDistortion /= ftArea;
	i = 0;
	ftAreaRatio = ftArea / ftGrids;
	ftAreaDistortion = 0;
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
	{
		if (fabs(fi->ftTwArea) >= 1e-8 && fabs(lpftLocalParaArea[i]) > 1e-8)
		{
			ftLocalDistortion = (lpftLocalParaArea[i] * ftAreaRatio / (fi->ftTwArea) + (fi->ftTwArea) / (lpftLocalParaArea[i] * ftAreaRatio));
			if (ftLocalDistortion > 10)
			{
				ftLocalDistortion = 10;
			}
			ftAreaDistortion += fabs(ftLocalDistortion * fi->ftTwArea);
		}
		++i;
	}
	ftAreaDistortion /= fabs(ftArea * 2.0);
	ftGrids /= 2.0;
	delete[]lpftLocalParaArea;
}

void CMeshCenter::m_fnExtendPath(CMeshBase::CPolyhedron::Halfedge_handle hhStart)
{
	CMeshBase::CPolyhedron::Halfedge_handle hh0, hh1, hh2;
	hh0 = hhStart;
	bool bTerminal;
	bTerminal = false;
	while (!bTerminal)
	{
		if (hh0->vertex()->iDegreeBias != 0)
		{
			bTerminal = true;
		}
		else
		{
			hh1 = hh0->next()->opposite();
			while (hh1 != hh0)
			{
				if (hh1->bFreeEdge)
				{
					break;
				}
				hh1 = hh1->next()->opposite();
			}
			if (hh1 == hh0)
			{
				bTerminal = true;
			}
			else
			{
				hh2 = hh0->opposite()->prev();
				while (hh2 != hh0)
				{
					if (hh2->bFreeEdge)
					{
						break;
					}
					hh2 = hh2->opposite()->prev();
				}
				if (hh1 == hh2 && !hh1->bSelected)
				{
					hh1->bSelected = true;
					hh1->opposite()->bSelected = true;
					hh0 = hh1->opposite();
				}
				else
				{
					bTerminal = true;
				}
			}
		}
	}
}

void CMeshCenter::m_fnDetectLocalPath(CMeshBase::CPolyhedron::Vertex_handle vh0, CMeshBase::CPolyhedron::Vertex_handle vh1)
{
	CMeshBase::CPolyhedron::Vertex_iterator vi;
	CMeshBase::CPolyhedron::Halfedge_handle hh0, hh, hhPath;
	CMeshBase::CPolyhedron::Vertex_handle vh;
	std::vector<CMeshBase::CPolyhedron::Vertex_handle> vecVH;
	//int iHeapPos;
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
	vh0->iSourceId = 0;
	vh0->ftDist = 0;
	vh0->iHeapPos = 0;
	vecVH.push_back(vh0);
	vh1->iSourceId = 1;
	vh1->ftDist = 0;
	vh1->iHeapPos = 1;
	vecVH.push_back(vh1);
	hhPath = NULL;
	while (!vecVH.empty())
	{
		vh = vecVH[0];
		vecVH[0] = *(vecVH.rbegin());
		vecVH[0]->iHeapPos = 0;
		vecVH.pop_back();
		m_fnSift(vecVH, 0);
		vh->iHeapPos = -1;
		vh->bDistConf = true;
		if (vh->iBordId == -1)
		{
			hh0 = vh->halfedge()->opposite();
			hh = hh0;
			do
			{
				if (!hh->bFreeEdge)
				{
					if (!hh->vertex()->bDistConf)
					{
						ftDist = vh->ftDist + hh->ftLen;
						if (hh->vertex()->iHeapPos == -1)
						{
							//iHeapPos = vecVH.size();
							hh->vertex()->iSourceId = vh->iSourceId;
							hh->vertex()->hhCome = hh;
							hh->vertex()->ftDist = ftDist;
							hh->vertex()->iHeapPos = vecVH.size();
							vecVH.push_back(hh->vertex());
							m_fnCheck(vecVH, hh->vertex()->iHeapPos);
						}
						else if (ftDist < hh->vertex()->ftDist)
						{
							hh->vertex()->iSourceId = vh->iSourceId;
							hh->vertex()->hhCome = hh;
							hh->vertex()->ftDist = ftDist;
							m_fnCheck(vecVH, hh->vertex()->iHeapPos);
						}
					}
					else if (hh->vertex()->iSourceId != -1 && hh->vertex()->iSourceId != vh->iSourceId)
					{
						ftDist = vh->ftDist + hh->vertex()->ftDist + hh->ftLen;
						if (ftDist < ftMinDist)
						{
							ftMinDist = ftDist;
							hhPath = hh;
						}
					}
				}
				hh = hh->prev()->opposite();
			} while (hh != hh0);
		}
	}
	if (hhPath != NULL)
	{
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

//void CMeshCenter::m_fnDetectPath()
//{
//    CMeshBase::CPolyhedron::Vertex_iterator vi;
//    CMeshBase::CPolyhedron::Edge_iterator ei;
//    CMeshBase::CPolyhedron::Halfedge_handle hh0, hh;
//    CMeshBase::CPolyhedron::Vertex_handle vh;
//    std::vector<CMeshBase::CPolyhedron::Vertex_handle> vecVH;
//    std::vector<CMeshBase::CPolyhedron::Halfedge_handle> vecHH;
//    std::vector<CMeshBase::CPolyhedron::Halfedge_handle>::iterator iterHH, iterMinHH;
//    std::vector<CMeshBase::FT> vecftDist;
//    std::vector<CMeshBase::FT>::iterator iterftDist, iterftMinDist;
//    int iHeapPos;
//    CMeshBase::FT ftDist;
//    bool bNew;
//    for(vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
//    {
//        vi->hhSource = NULL;
//        vi->hhCome = NULL;
//        vi->ftDist = -1;
//        vi->bDistConf = false;
//        vi->iHeapPos = -1;
//    }
//    for(vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
//    {
//        if(vi->iDegreeBias != 0)//(vi->nExpectedDegree != 6)
//        {
//            iHeapPos = vecVH.size();
//            //vi->vhSource = &(*vi);
//            vi->ftDist = 0;
//			vi->ftComeAngle = 0.0;
//            vi->iHeapPos = iHeapPos;
//            vecVH.push_back(&(*vi));
//        }
//    }
//    if(!vecVH.empty())
//    {
//        while(!vecVH.empty())
//        {
//            vh = vecVH[0];
//            vecVH[0] = *(vecVH.rbegin());
//            vecVH[0]->iHeapPos = 0;
//            vecVH.pop_back();
//            m_fnSift(vecVH, 0);
//            vh->iHeapPos = -1;
//            vh->bDistConf = true;
//            hh0 = vh->halfedge()->opposite();
//            hh = hh0;
//            do
//            {
//                if(!hh->vertex()->bDistConf)
//                {
//                    ftDist = vh->ftDist + hh->ftLen;
//                    if(hh->vertex()->iHeapPos == -1)
//                    {
//                        iHeapPos = vecVH.size();
//                        hh->vertex()->vhSource = vh->vhSource;
//                        hh->vertex()->hhCome = hh;
//                        hh->vertex()->ftDist = ftDist;
//                        hh->vertex()->iHeapPos = iHeapPos;
//                        vecVH.push_back(hh->vertex());
//                        m_fnCheck(vecVH, iHeapPos);
//                    }
//                    else
//                    {
//                        if(ftDist < hh->vertex()->ftDist)
//                        {
//                            hh->vertex()->vhSource = vh->vhSource;
//                            hh->vertex()->hhCome = hh;
//                            hh->vertex()->ftDist = ftDist;
//                            m_fnCheck(vecVH, hh->vertex()->iHeapPos);
//                        }
//                    }
//                }
//                hh = hh->prev()->opposite();
//            }while(hh != hh0);
//        }
//        for(ei = m_mshSurface.edges_begin(); ei != m_mshSurface.edges_end(); ++ei)
//        {
//            if(ei->vertex()->vhSource != NULL && ei->opposite()->vertex()->vhSource != NULL)
//            {
//                if(ei->vertex()->vhSource->iDegreeBias < 0 && ei->opposite()->vertex()->vhSource->iDegreeBias > 0)
//                {
//                    bNew = true;
//                    ftDist = (ei->vertex()->ftDist + ei->opposite()->vertex()->ftDist + ei->ftLen);
//                    if(!vecHH.empty())
//                    {
//                        iterftDist = vecftDist.begin();
//                        for(iterHH = vecHH.begin(); iterHH != vecHH.end(); ++iterHH)
//                        {
//                            if((*iterHH)->vertex()->vhSource == ei->vertex()->vhSource && (*iterHH)->opposite()->vertex()->vhSource == ei->opposite()->vertex()->vhSource)
//                            {
//                                if(ftDist < (*iterftDist))
//                                {
//                                    *iterHH = &(*ei);
//                                    *iterftDist = ftDist;
//                                }
//                                bNew = false;
//                                break;
//                            }
//                            ++iterftDist;
//                        }
//                    }
//                    if(bNew)
//                    {
//                        vecHH.push_back(&(*ei));
//                        vecftDist.push_back(ftDist);
//                    }
//                }
//                else if(ei->vertex()->vhSource->iDegreeBias > 0 && ei->opposite()->vertex()->vhSource->iDegreeBias < 0)//(ei->vertex()->vhSource->nExpectedDegree > 6 && ei->opposite()->vertex()->vhSource->nExpectedDegree < 6)
//                {
//                    bNew = true;
//                    ftDist = ei->vertex()->ftDist + ei->opposite()->vertex()->ftDist + ei->ftLen;
//                    if(!vecHH.empty())
//                    {
//                        iterftDist = vecftDist.begin();
//                        for(iterHH = vecHH.begin(); iterHH != vecHH.end(); ++iterHH)
//                        {
//                            if((*iterHH)->vertex()->vhSource == ei->opposite()->vertex()->vhSource && (*iterHH)->opposite()->vertex()->vhSource == ei->vertex()->vhSource)
//                            {
//                                if(ftDist < *(iterftDist))
//                                {
//                                    *iterHH = ei->opposite();
//                                    *iterftDist = ftDist;
//                                }
//                                bNew = false;
//                                break;
//                            }
//                            ++iterftDist;
//                        }
//                    }
//                    if(bNew)
//                    {
//                        vecHH.push_back(ei->opposite());
//                        vecftDist.push_back(ftDist);
//                    }
//                }
//            }
//            ei->bFreeEdge = false;
//            ei->opposite()->bFreeEdge = false;
//        }
//        while(!vecHH.empty())
//        {
//            iterMinHH = vecHH.begin();
//            iterftMinDist = vecftDist.begin();
//            iterHH = iterMinHH;
//            iterftDist = iterftMinDist;
//            ++iterHH;
//            ++iterftDist;
//            while(iterHH != vecHH.end())
//            {
//                if(*iterftDist < *iterftMinDist)
//                {
//                    iterftMinDist = iterftDist;
//                    iterMinHH = iterHH;
//                }
//                ++iterHH;
//                ++iterftDist;
//            }
//            hh = *iterMinHH;
//            hh->bFreeEdge = true;
//            while(hh != NULL)
//            {
//                hh->bFreeEdge = true;
//                hh->opposite()->bFreeEdge = true;
//                hh = hh->opposite()->vertex()->hhCome;
//            }
//            hh = (*iterMinHH)->opposite();
//            hh->bFreeEdge = true;
//            while(hh != NULL)
//            {
//                hh->bFreeEdge = true;
//                hh->opposite()->bFreeEdge = true;
//                hh = hh->opposite()->vertex()->hhCome;
//            }
//            hh = (*iterMinHH);
//            *iterMinHH = *(vecHH.rbegin());
//            *iterftMinDist = *(vecftDist.rbegin());
//            vecHH.pop_back();
//            vecftDist.pop_back();
//            iterHH = vecHH.begin();
//            iterftDist = vecftDist.begin();
//            while(iterHH != vecHH.end())
//            {
//                if((*iterHH)->vertex()->vhSource == hh->vertex()->vhSource || (*iterHH)->opposite()->vertex()->vhSource == hh->opposite()->vertex()->vhSource)
//                {
//                    *iterHH = *(vecHH.rbegin());
//                    *iterftDist = *(vecftDist.rbegin());
//                    vecHH.pop_back();
//                    vecftDist.pop_back();
//                }
//                else
//                {
//                    ++iterHH;
//                    ++iterftDist;
//                }
//            }
//        }
//    }
//}

void CMeshCenter::m_fnClearSelectedPath()
{
	CMeshBase::CPolyhedron::Edge_iterator ei;
	for (ei = m_mshSurface.edges_begin(); ei != m_mshSurface.edges_end(); ++ei)
	{
		if (ei->bFreeEdge && ei->bSelected)
		{
			ei->bFreeEdge = false;
			ei->bSelected = false;
			ei->opposite()->bFreeEdge = false;
			ei->opposite()->bSelected = false;
		}
	}
}

void CMeshCenter::m_fnFreeBoundary()
{
	CMeshBase::CPolyhedron::Edge_iterator ei;
	for (ei = m_mshSurface.edges_begin(); ei != m_mshSurface.edges_end(); ++ei)
	{
		if (ei->is_border_edge())
		{
			ei->bFreeEdge = true;
			ei->opposite()->bFreeEdge = true;
		}
	}
}

bool CMeshCenter::m_fnIsAligned()
{
	std::vector<CCrevasse>::iterator iterCrevasse;
	m_bAligned = true;
	for (iterCrevasse = m_vecCrevasses.begin(); (iterCrevasse != m_vecCrevasses.end() && m_bAligned); ++iterCrevasse)
	{
		if (!iterCrevasse->is_border())
		{
			if (iterCrevasse->next()->next()->next()->next() != &(*iterCrevasse))
			{
				m_bAligned = false;
			}
		}
	}
	return m_bAligned;
}

//void CMeshCenter::m_fnCopyChartPara()
//{
//	CMeshBase::CPolyhedron::Edge_iterator ei;
//	for (ei = m_mshSurface.edges_begin(); ei != m_mshSurface.edges_end(); ++ei)
//	{
//		memcpy(ei->lpftCurlFreePara, ei->lpftChartPara, 2 * sizeof(CMeshBase::FT));
//		memcpy(ei->opposite()->lpftCurlFreePara, ei->opposite()->lpftChartPara, 2 * sizeof(CMeshBase::FT));
//	}
//}

void CMeshCenter::m_fnFillVertexCrossField()
{
	CMeshBase::CPolyhedron::Vertex_iterator vi;
	CMeshBase::CPolyhedron::Halfedge_handle hh0, hh;
	CMeshBase::Vector_3 lpvtVector[2];
	CMeshBase::FT ftAngle;
	CMeshBase::CPolyhedron::Facet_iterator fi;
	bool bBlank;
	int iAccSeamType;
	for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
	{
		lpvtVector[0] = lpvtVector[1] = CMeshBase::Vector_3(0, 0, 0);
		iAccSeamType = 0;
		if (vi->iBordId == -1)
		{
			hh0 = vi->halfedge();
			if (vi->nNumPorts == 0)
			{
				hh = hh0;
				do
				{
					ftAngle = hh->facet()->ftChartDir - iAccSeamType * ftHalfPi;
					lpvtVector[0] = lpvtVector[0] + hh->facet()->vtPrincipalAxis * cos(ftAngle) + hh->facet()->vtAuxiliaryAxis * sin(ftAngle);
					lpvtVector[1] = lpvtVector[1] - hh->facet()->vtPrincipalAxis * sin(ftAngle) + hh->facet()->vtAuxiliaryAxis * cos(ftAngle);
					iAccSeamType = (iAccSeamType + hh->next()->iSeamType + 6) % 4 - 2;
					hh = hh->next()->opposite();
				} while (hh != hh0);
				lpvtVector[0] = lpvtVector[0] / sqrt(CGAL::squared_length(lpvtVector[0]));
				lpvtVector[1] = lpvtVector[1] / sqrt(CGAL::squared_length(lpvtVector[1]));
				hh = hh0;
				do
				{
					hh->lpvtViewField[0] = lpvtVector[(iAccSeamType + 2) % 2] * CMeshBase::FT(((iAccSeamType + 2) % 4) / 2 * 2 - 1);
					hh->lpvtViewField[1] = lpvtVector[1 - (iAccSeamType + 2) % 2] * CMeshBase::FT(((iAccSeamType + 3) % 4) / 2 * 2 - 1);
					iAccSeamType = (iAccSeamType + hh->next()->iSeamType + 6) % 4 - 2;
					hh = hh->next()->opposite();
				} while (hh != hh0);
			}
			else
			{
				hh = hh0;
				do
				{
					hh->lpvtViewField[0] = lpvtVector[0];
					hh->lpvtViewField[1] = lpvtVector[1];
					hh = hh->next()->opposite();
				} while (hh != hh0);
			}
		}
		else
		{
			hh0 = vi->halfedge();
			while (!hh0->opposite()->is_border())
			{
				hh0 = hh0->opposite()->prev();
			}
			hh = hh0;
			while (!hh->is_border())
			{
				ftAngle = hh->facet()->ftChartDir - iAccSeamType * ftHalfPi;
				lpvtVector[0] = lpvtVector[0] + hh->facet()->vtPrincipalAxis * cos(ftAngle) + hh->facet()->vtAuxiliaryAxis * sin(ftAngle);
				lpvtVector[1] = lpvtVector[1] - hh->facet()->vtPrincipalAxis * sin(ftAngle) + hh->facet()->vtAuxiliaryAxis * cos(ftAngle);
				iAccSeamType = (iAccSeamType + hh->next()->iSeamType + 6) % 4 - 2;
				hh = hh->next()->opposite();
			}
			lpvtVector[0] = lpvtVector[0] / sqrt(CGAL::squared_length(lpvtVector[0]));
			lpvtVector[1] = lpvtVector[1] / sqrt(CGAL::squared_length(lpvtVector[1]));
			iAccSeamType = 0;
			hh = hh0;
			while (!hh->is_border())
			{
				hh->lpvtViewField[0] = lpvtVector[(iAccSeamType + 2) % 2] * CMeshBase::FT(((iAccSeamType + 2) % 4) / 2 * 2 - 1);
				hh->lpvtViewField[1] = lpvtVector[1 - (iAccSeamType + 2) % 2] * CMeshBase::FT(((iAccSeamType + 3) % 4) / 2 * 2 - 1);
				iAccSeamType = (iAccSeamType + hh->next()->iSeamType + 6) % 4 - 2;
				hh = hh->next()->opposite();
			}
		}
	}
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
	{
		bBlank = true;
		hh0 = fi->halfedge();
		hh = hh0;
		do
		{
			if (hh->vertex()->iBordId != -1 || hh->vertex()->nNumPorts == 0)
			{
				bBlank = false;
			}
			hh = hh->next();
		} while (hh != hh0);
		if (bBlank)
		{
			ftAngle = hh->facet()->ftChartDir;
			lpvtVector[0] = hh->facet()->vtPrincipalAxis * cos(ftAngle) + hh->facet()->vtAuxiliaryAxis * sin(ftAngle);
			lpvtVector[1] = - hh->facet()->vtPrincipalAxis * sin(ftAngle) + hh->facet()->vtAuxiliaryAxis * cos(ftAngle);
			lpvtVector[0] = lpvtVector[0] / sqrt(CGAL::squared_length(lpvtVector[0]));
			lpvtVector[1] = lpvtVector[1] / sqrt(CGAL::squared_length(lpvtVector[1]));
			hh = hh0;
			do
			{
				hh->lpvtViewField[0] = lpvtVector[0];
				hh->lpvtViewField[1] = lpvtVector[1];
				hh = hh->next();
			} while (hh != hh0);
		}
	}
}

void CMeshCenter::m_fnFillVertexFrameField()
{
	CMeshBase::CPolyhedron::Vertex_iterator vi;
	CMeshBase::CPolyhedron::Halfedge_handle hh0, hh;
	CMeshBase::Vector_3 lpvtVector[2];
	CMeshBase::FT lpftAngle[4];
	CMeshBase::CPolyhedron::Facet_iterator fi;
	bool bBlank;
	int iAccSeamType;
	for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
	{
		lpvtVector[0] = lpvtVector[1] = CMeshBase::Vector_3(0, 0, 0);
		iAccSeamType = 0;
		if (vi->iBordId == -1)
		{
			hh0 = vi->halfedge();
			if (vi->nNumPorts == 0)
			{
				hh = hh0;
				do
				{
					//ftAngle = hh->facet()->ftChartDir - iAccSeamType * ftHalfPi;
					lpftAngle[0] = hh->facet()->ftChartDir;
					lpftAngle[1] = hh->facet()->ftFrameDir;
					lpftAngle[2] = hh->facet()->ftChartDir + ftPi;
					lpftAngle[3] = hh->facet()->ftFrameDir + ftPi;
					lpvtVector[0] = lpvtVector[0] + hh->facet()->vtPrincipalAxis * cos(lpftAngle[(4 - iAccSeamType) % 4]) + hh->facet()->vtAuxiliaryAxis * sin(lpftAngle[(4 - iAccSeamType) % 4]);
					lpvtVector[1] = lpvtVector[1] + hh->facet()->vtPrincipalAxis * cos(lpftAngle[(5 - iAccSeamType) % 4]) + hh->facet()->vtAuxiliaryAxis * sin(lpftAngle[(5 - iAccSeamType) % 4]);
					iAccSeamType = (iAccSeamType + hh->next()->iSeamType + 6) % 4 - 2;
					hh = hh->next()->opposite();
				} while (hh != hh0);
				lpvtVector[0] = lpvtVector[0] / sqrt(CGAL::squared_length(lpvtVector[0]));
				lpvtVector[1] = lpvtVector[1] / sqrt(CGAL::squared_length(lpvtVector[1]));
				hh = hh0;
				do
				{
					hh->lpvtViewField[0] = lpvtVector[(iAccSeamType + 2) % 2] * CMeshBase::FT(((iAccSeamType + 2) % 4) / 2 * 2 - 1);
					hh->lpvtViewField[1] = lpvtVector[1 - (iAccSeamType + 2) % 2] * CMeshBase::FT(((iAccSeamType + 3) % 4) / 2 * 2 - 1);
					iAccSeamType = (iAccSeamType + hh->next()->iSeamType + 6) % 4 - 2;
					hh = hh->next()->opposite();
				} while (hh != hh0);
			}
			else
			{
				hh = hh0;
				do
				{
					hh->lpvtViewField[0] = lpvtVector[0];
					hh->lpvtViewField[1] = lpvtVector[1];
					hh = hh->next()->opposite();
				} while (hh != hh0);
			}
		}
		else
		{
			hh0 = vi->halfedge();
			while (!hh0->opposite()->is_border())
			{
				hh0 = hh0->opposite()->prev();
			}
			hh = hh0;
			while (!hh->is_border())
			{
				//ftAngle = hh->facet()->ftChartDir - iAccSeamType * ftHalfPi;
				lpftAngle[0] = hh->facet()->ftChartDir;
				lpftAngle[1] = hh->facet()->ftFrameDir;
				lpftAngle[2] = hh->facet()->ftChartDir + ftPi;
				lpftAngle[3] = hh->facet()->ftFrameDir + ftPi;
				lpvtVector[0] = lpvtVector[0] + hh->facet()->vtPrincipalAxis * cos(lpftAngle[(4 - iAccSeamType) % 4]) + hh->facet()->vtAuxiliaryAxis * sin(lpftAngle[(4 - iAccSeamType) % 4]);
				lpvtVector[1] = lpvtVector[1] + hh->facet()->vtPrincipalAxis * cos(lpftAngle[(5 - iAccSeamType) % 4]) + hh->facet()->vtAuxiliaryAxis * sin(lpftAngle[(5 - iAccSeamType) % 4]);
				iAccSeamType = (iAccSeamType + hh->next()->iSeamType + 6) % 4 - 2;
				hh = hh->next()->opposite();
			}
			lpvtVector[0] = lpvtVector[0] / sqrt(CGAL::squared_length(lpvtVector[0]));
			lpvtVector[1] = lpvtVector[1] / sqrt(CGAL::squared_length(lpvtVector[1]));
			iAccSeamType = 0;
			hh = hh0;
			while (!hh->is_border())
			{
				hh->lpvtViewField[0] = lpvtVector[(iAccSeamType + 2) % 2] * CMeshBase::FT(((iAccSeamType + 2) % 4) / 2 * 2 - 1);
				hh->lpvtViewField[1] = lpvtVector[1 - (iAccSeamType + 2) % 2] * CMeshBase::FT(((iAccSeamType + 3) % 4) / 2 * 2 - 1);
				iAccSeamType = (iAccSeamType + hh->next()->iSeamType + 6) % 4 - 2;
				hh = hh->next()->opposite();
			}
		}
	}
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
	{
		bBlank = true;
		hh0 = fi->halfedge();
		hh = hh0;
		do
		{
			if (hh->vertex()->iBordId != -1 || hh->vertex()->nNumPorts == 0)
			{
				bBlank = false;
			}
			hh = hh->next();
		} while (hh != hh0);
		if (bBlank)
		{
			lpftAngle[0] = hh->facet()->ftChartDir;
			lpftAngle[1] = hh->facet()->ftFrameDir;
			lpftAngle[2] = hh->facet()->ftChartDir + ftPi;
			lpftAngle[3] = hh->facet()->ftFrameDir + ftPi;
			lpvtVector[0] = lpvtVector[0] + hh->facet()->vtPrincipalAxis * cos(lpftAngle[(4 - iAccSeamType) % 4]) + hh->facet()->vtAuxiliaryAxis * sin(lpftAngle[(4 - iAccSeamType) % 4]);
			lpvtVector[1] = lpvtVector[1] + hh->facet()->vtPrincipalAxis * cos(lpftAngle[(5 - iAccSeamType) % 4]) + hh->facet()->vtAuxiliaryAxis * sin(lpftAngle[(5 - iAccSeamType) % 4]);
			lpvtVector[0] = lpvtVector[0] / sqrt(CGAL::squared_length(lpvtVector[0]));
			lpvtVector[1] = lpvtVector[1] / sqrt(CGAL::squared_length(lpvtVector[1]));
			hh = hh0;
			do
			{
				hh->lpvtViewField[0] = lpvtVector[0];
				hh->lpvtViewField[1] = lpvtVector[1];
				hh = hh->next();
			} while (hh != hh0);
		}
	}
}

void CMeshCenter::m_fnFillVertexElectricField()
{
	CMeshBase::CPolyhedron::Vertex_iterator vi;
	CMeshBase::CPolyhedron::Halfedge_handle hh0, hh;
	CMeshBase::Vector_3 vtVector;
	CMeshBase::CPolyhedron::Facet_iterator fi;
	int nNumFacets;
	bool bBlank;
	for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
	{
		vtVector = CMeshBase::Vector_3(0, 0, 0);
		hh0 = vi->halfedge();
		if (vi->iDegreeBias == 0)
		{
			nNumFacets = 0;
			hh = hh0;
			do
			{
				if (!hh->is_border())
				{
					vtVector = vtVector + hh->facet()->vtNablaPhi;
					++nNumFacets;
				}
				hh = hh->next()->opposite();
			} while (hh != hh0);
			vtVector = vtVector / CMeshBase::FT(nNumFacets);
			vi->ftFieldStrength = sqrt(CGAL::squared_length(vtVector)) ;
			if (vi->ftFieldStrength > 0.0)
			{
				vtVector = vtVector / vi->ftFieldStrength;
			}
		}
		else
		{
			vi->ftFieldStrength = 0;
			nNumFacets = 0;
			hh = hh0;
			do
			{
				if (!hh->is_border())
				{
					vi->ftFieldStrength += sqrt(CGAL::squared_length(hh->facet()->vtNablaPhi));
					++nNumFacets;
				}
				hh = hh->next()->opposite();
			} while (hh != hh0);
			vi->ftFieldStrength /= nNumFacets;
		}
		hh = hh0;
		do
		{
			if (!hh->is_border())
			{
				hh->lpvtViewField[0] = vtVector;
				hh->lpvtViewField[1] = CMeshBase::Vector_3(0, 0, 0);
			}
			hh = hh->next()->opposite();
		} while (hh != hh0);
	}
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
	{
		bBlank = true;
		hh0 = fi->halfedge();
		hh = hh0;
		do
		{
			if (hh->vertex()->iDegreeBias == 0)
			{
				bBlank = false;
			}
			hh = hh->next();
		} while (hh != hh0);
		if (bBlank)
		{
			vtVector = fi->vtNablaPhi / sqrt(CGAL::squared_length(fi->vtNablaPhi));
			hh = hh0;
			do
			{
				hh->lpvtViewField[0] = vtVector;
				hh = hh->next();
			} while (hh != hh0);
		}
	}
}

void CMeshCenter::m_fnFillVertexNablaTheta()
{
	CMeshBase::CPolyhedron::Vertex_iterator vi;
	CMeshBase::CPolyhedron::Halfedge_handle hh0, hh;
	CMeshBase::Vector_3 vtVector;
	CMeshBase::CPolyhedron::Facet_iterator fi;
	int nNumFacets;
	bool bBlank;
	for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
	{
		vtVector = CMeshBase::Vector_3(0, 0, 0);
		hh0 = vi->halfedge();
		if (vi->iDegreeBias == 0)
		{
			nNumFacets = 0;
			hh = hh0;
			do
			{
				if (!hh->is_border())
				{
					vtVector = vtVector + CGAL::cross_product(hh->facet()->vtNablaTheta, hh->facet()->vtNorm) / hh->facet()->ftTwArea;
					++nNumFacets;
				}
				hh = hh->next()->opposite();
			} while (hh != hh0);
			vtVector = vtVector / CMeshBase::FT(nNumFacets);
			vi->ftFieldStrength = sqrt(CGAL::squared_length(vtVector)) ;
			if (vi->ftFieldStrength > 0.0)
			{
				vtVector = vtVector / vi->ftFieldStrength;
			}
		}
		else
		{
			vi->ftFieldStrength = 0;
			nNumFacets = 0;
			hh = hh0;
			do
			{
				if (!hh->is_border())
				{
					vi->ftFieldStrength += sqrt(CGAL::squared_length(hh->facet()->vtNablaTheta));
					++nNumFacets;
				}
				hh = hh->next()->opposite();
			} while (hh != hh0);
			vi->ftFieldStrength /= nNumFacets;
		}
		hh = hh0;
		do
		{
			if (!hh->is_border())
			{
				hh->lpvtViewField[0] = vtVector;
				hh->lpvtViewField[1] = CMeshBase::Vector_3(0, 0, 0);
			}
			hh = hh->next()->opposite();
		} while (hh != hh0);
	}
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
	{
		bBlank = true;
		hh0 = fi->halfedge();
		hh = hh0;
		do
		{
			if (hh->vertex()->iDegreeBias == 0)
			{
				bBlank = false;
			}
			hh = hh->next();
		} while (hh != hh0);
		if (bBlank)
		{
			vtVector = CGAL::cross_product(fi->vtNablaTheta, fi->vtNorm) / (sqrt(CGAL::squared_length(fi->vtNablaTheta)) * fi->ftTwArea);
			hh = hh0;
			do
			{
				hh->lpvtViewField[0] = vtVector;
				hh = hh->next();
			} while (hh != hh0);
		}
	}
}

void CMeshCenter::m_fnFillVertexCurlField()
{
	CMeshBase::CPolyhedron::Vertex_iterator vi;
	CMeshBase::CPolyhedron::Halfedge_handle hh0, hh;
	CMeshBase::Vector_3 vtVector;
	int nNumFacets;
	CMeshBase::FT ftMinLen;
	for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
	{
		vtVector = CMeshBase::Vector_3(0, 0, 0);
		hh0 = vi->halfedge();
		nNumFacets = 0;
		hh = hh0;
		do
		{
			if (!(hh->is_border() || hh->opposite()->is_border() || hh->next()->opposite()->is_border() || hh->prev()->opposite()->is_border()))
			{
				vtVector = vtVector + hh->facet()->vtDifNabla;
				++nNumFacets;
			}
			hh = hh->next()->opposite();
		} while (hh != hh0);
		if (nNumFacets > 0)
		{
			vtVector = vtVector / CMeshBase::FT(nNumFacets);
		}
		vi->ftFieldStrength = sqrt(CGAL::squared_length(vtVector));
		if (vi->ftFieldStrength > 0.0)
		{
			vtVector = vtVector / vi->ftFieldStrength;
		}
		hh = hh0;
		do
		{
			if (!hh->is_border())
			{
				hh->lpvtViewField[0] = vtVector;
				hh->lpvtViewField[1] = CMeshBase::Vector_3(0, 0, 0);
			}
			hh = hh->next()->opposite();
		} while (hh != hh0);
	}
	for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
	{
		if (vi->iBordId != -1 && vi->ftFieldStrength == 0.0)
		{
			hh0 = vi->halfedge()->opposite();
			ftMinLen = DBL_MAX;
			hh = hh0;
			do
			{
				if (hh->vertex()->ftFieldStrength > 0.0 && hh->ftLen < ftMinLen)
				{
					vi->ftFieldStrength = hh->vertex()->ftFieldStrength;
					ftMinLen = hh->ftLen;
				}
				hh = hh->opposite()->next();
			} while (hh != hh0);
		}
	}
}

void CMeshCenter::m_fnDetectFieldMaximum()
{
	CMeshBase::CPolyhedron::Vertex_iterator vi;
	CMeshBase::CPolyhedron::Halfedge_handle hh0, hh;
	CMeshBase::FT ftCrossGradStrength, ftElectricStrength, ftCurlLength;
	CMeshBase::Vector_3 vtCrossGrad, vtElectric, vtCurlField;
	int nNumFacets;
	m_ftMaxCrossGrad = 0;
	m_ftMaxElectricField = 0;
	m_ftMaxCurlField = 0;
	for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
	{
		vtCrossGrad = CMeshBase::Vector_3(0, 0, 0);
		vtElectric = CMeshBase::Vector_3(0, 0, 0);
		hh0 = vi->halfedge();
		if (vi->iDegreeBias == 0)
		{
			nNumFacets = 0;
			hh = hh0;
			do
			{
				if (!hh->is_border())
				{
					vtCrossGrad = vtCrossGrad + hh->facet()->vtNablaTheta;
					vtElectric = vtElectric + hh->facet()->vtNablaPhi;
					++nNumFacets;
				}
				hh = hh->next()->opposite();
			} while (hh != hh0);
			vtCrossGrad = vtCrossGrad / CMeshBase::FT(nNumFacets);
			vtElectric = vtElectric / CMeshBase::FT(nNumFacets);
			ftCrossGradStrength = sqrt(CGAL::squared_length(vtCrossGrad));
			ftElectricStrength = sqrt(CGAL::squared_length(vtElectric));
		}
		else
		{
			vi->ftFieldStrength = 0;
			nNumFacets = 0;
			hh = hh0;
			do
			{
				if (!hh->is_border())
				{
					ftCrossGradStrength += sqrt(CGAL::squared_length(hh->facet()->vtNablaTheta));
					ftElectricStrength += sqrt(CGAL::squared_length(hh->facet()->vtNablaPhi));
					++nNumFacets;
				}
				hh = hh->next()->opposite();
			} while (hh != hh0);
			ftCrossGradStrength /= nNumFacets;
			ftElectricStrength /= nNumFacets;
		}
		if (m_ftMaxCrossGrad < ftCrossGradStrength)
		{
			m_ftMaxCrossGrad = ftCrossGradStrength;
		}
		if (m_ftMaxElectricField < ftElectricStrength)
		{
			m_ftMaxElectricField = ftElectricStrength;
		}
		vtCurlField = CMeshBase::Vector_3(0, 0, 0);
		hh0 = vi->halfedge();
		nNumFacets = 0;
		hh = hh0;
		do
		{
			if (!(hh->is_border() || hh->opposite()->is_border() || hh->next()->opposite()->is_border() || hh->prev()->opposite()->is_border()))
			{
				vtCurlField = vtCurlField + hh->facet()->vtDifNabla;
				++nNumFacets;
			}
			hh = hh->next()->opposite();
		} while (hh != hh0);
		if (nNumFacets > 0)
		{
			vtCurlField = vtCurlField / CMeshBase::FT(nNumFacets);
		}
		ftCurlLength = sqrt(CGAL::squared_length(vtCurlField));
		if (ftCurlLength > m_ftMaxCurlField)
		{
			m_ftMaxCurlField = ftCurlLength;
		}
	}
}

void CMeshCenter::m_fnFixDirection()
{
	CMeshBase::FT ftTheta, ftRot;
	if (m_fhFacetFix != NULL)
	{
		if (m_fhFacetFix->bDirFixed)
		{
			m_fhFacetFix->bDirFixed = false;
		}
		else if (m_hhDirEdge != NULL)
		{
			ftTheta = atan2(m_hhDirEdge->vtVector * m_fhFacetFix->vtAuxiliaryAxis, m_hhDirEdge->vtVector * m_fhFacetFix->vtPrincipalAxis);
			ftRot = ftTheta - m_fhFacetFix->ftChartDir;
			ftRot = ftRot - CMeshBase::FT(int(ftRot / ftHalfPi + 15.5) - 15) * ftHalfPi;
			m_fhFacetFix->ftChartDir += ftRot;
			m_fhFacetFix->bDirFixed = true;
		}
	}
}

void CMeshCenter::m_fnMultiFixDirection()
{
	CMeshBase::CPolyhedron::Edge_iterator ei;
	CMeshBase::CPolyhedron::Halfedge_handle hh0, hh;
	CMeshBase::FT ftTheta, ftRot;
	for (ei = m_mshSurface.edges_begin(); ei != m_mshSurface.edges_end(); ++ei)
	{
		if (ei->bSelected)
		{
			hh0 = &(*ei);
			hh = hh0;
			do
			{
				if (!hh->is_border())
				{
					ftTheta = atan2(hh->vtVector * hh->facet()->vtAuxiliaryAxis, hh->vtVector * hh->facet()->vtPrincipalAxis);
					ftRot = ftTheta - hh->facet()->ftChartDir;
					ftRot = ftRot - CMeshBase::FT(int(ftRot / ftHalfPi + 15.5) - 15) * ftHalfPi;
					hh->facet()->ftChartDir += ftRot;
					hh->facet()->bDirFixed = true;
				}
				hh = hh->opposite();
			} while (hh != hh0);
		}
	}
}

#include "CParameterGenerator.h"

CParameterGenerator::CParameterGenerator(CMeshBase::CPolyhedron &mshSurface, std::vector<CMeshCenter::CCrevasse>&vecCrevasses) : m_mshSurface(mshSurface), m_vecCrevasses(vecCrevasses)
{
    //ctor
}


void CParameterGenerator::m_fnEstimateParaCurl(CMeshBase::FT* lpftCurl)
{
	CMeshBase::CPolyhedron::Facet_iterator fi;
	CMeshBase::CPolyhedron::Halfedge_handle hh0, hh;
	CMeshBase::FT ftRoundDir, ftScale, ftMeanCurl;
	CMeshBase::Vector_3 lpvtParaGrad[2];
	lpftCurl[0] = lpftCurl[3] = lpftCurl[6] = DBL_MAX;
	lpftCurl[1] = lpftCurl[4] = lpftCurl[7] = 0;
	lpftCurl[2] = lpftCurl[5] = lpftCurl[8] = 0;
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
	{
		fi->lpftCurl[0] = fi->lpftCurl[1] = 0.0;
		hh0 = fi->halfedge();
		hh = hh0;
		do
		{
			ftRoundDir = fi->ftChartDir - hh->ftMom * (fi->ftTwArea / hh->ftDot) / hh->ftLaplaceCoef;
			ftScale = exp((hh->vertex()->ftDensity + hh->prev()->vertex()->ftDensity - hh->next()->vertex()->ftDensity * 2) / 12);
			lpvtParaGrad[0] = (fi->vtPrincipalAxis * cos(ftRoundDir) + fi->vtAuxiliaryAxis * sin(ftRoundDir)) * (ftScale / fi->ftRadius);
			lpvtParaGrad[1] = (-fi->vtPrincipalAxis * sin(ftRoundDir) + fi->vtAuxiliaryAxis * cos(ftRoundDir)) * (ftScale / fi->ftRadius);
			fi->lpftCurl[0] += lpvtParaGrad[0] * hh->vtVector;
			fi->lpftCurl[1] += lpvtParaGrad[1] * hh->vtVector;
			hh = hh->next();
		} while (hh != hh0);
		fi->lpftCurl[0] = fabs(fi->lpftCurl[0] / fi->ftTwArea);
		fi->lpftCurl[1] = fabs(fi->lpftCurl[1] / fi->ftTwArea);
		ftMeanCurl = sqrt(fi->lpftCurl[0] * fi->lpftCurl[0] + fi->lpftCurl[1] * fi->lpftCurl[1]);
		if (fi->lpftCurl[0] < lpftCurl[0])
		{
			lpftCurl[0] = fi->lpftCurl[0];
		}
		if (fi->lpftCurl[0] > lpftCurl[2])
		{
			lpftCurl[2] = fi->lpftCurl[0];
		}
		if (fi->lpftCurl[1] < lpftCurl[3])
		{
			lpftCurl[3] = fi->lpftCurl[1];
		}
		if (fi->lpftCurl[1] > lpftCurl[5])
		{
			lpftCurl[5] = fi->lpftCurl[1];
		}
		if (ftMeanCurl < lpftCurl[6])
		{
			lpftCurl[6] = ftMeanCurl;
		}
		if (ftMeanCurl > lpftCurl[8])
		{
			lpftCurl[8] = ftMeanCurl;
		}
		lpftCurl[1] += fi->lpftCurl[0];
		lpftCurl[4] += fi->lpftCurl[1];
		lpftCurl[7] += ftMeanCurl;
	}
	lpftCurl[1] /= m_mshSurface.size_of_facets();
	lpftCurl[4] /= m_mshSurface.size_of_facets();
	lpftCurl[7] /= m_mshSurface.size_of_facets();
}



void CParameterGenerator::m_fnMarkCrevasse()
{
	CMeshBase::CPolyhedron::Vertex_iterator vi;
	CMeshBase::CPolyhedron::Halfedge_handle hh0, hh1, hh2;
	for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
	{
		vi->nCrevasseSize = 0;
		hh0 = vi->halfedge();
		hh1 = hh0;
		do
		{
			if (hh1->bFreeEdge)
			{
				++(vi->nCrevasseSize);
			}
			hh1 = hh1->next()->opposite();
		} while (hh1 != hh0);
		if (vi->nCrevasseSize > 0)
		{
			do
			{
				if (hh1->bFreeEdge)
				{
					hh2 = hh1;
					do
					{
						hh2->hhCrevassePrev = hh1;
						hh2 = hh2->next()->opposite();
					} while (!hh2->bFreeEdge);
					hh2 = hh1->opposite();
					do
					{
						hh2->prev()->hhCrevasseNext = hh1->opposite();
						hh2 = hh2->prev()->opposite();
					} while (!hh2->bFreeEdge);
				}
				hh1 = hh1->next()->opposite();
			} while (hh1 != hh0);
		}
	}
}

void CParameterGenerator::m_fnExtractCrevasses()
{
	CMeshBase::CPolyhedron::Edge_iterator ei;
	CMeshCenter::CCrevasse newCrevasse;
	std::vector<CMeshCenter::CCrevasse>::iterator iterCrevasse;
	int iCrevasseId;
	CMeshBase::CPolyhedron::Halfedge_handle hh;
	std::vector<CMeshBase::CPolyhedron::Halfedge_handle>::iterator iterhhVisitor;
	for (ei = m_mshSurface.edges_begin(); ei != m_mshSurface.edges_end(); ++ei)
	{
		ei->iCrevasseId = -1;
		ei->opposite()->iCrevasseId = -1;
	}
	iCrevasseId = 0;
	for (ei = m_mshSurface.edges_begin(); ei != m_mshSurface.edges_end(); ++ei)
	{
		if (ei->bFreeEdge && ei->iCrevasseId == -1)
		{
			if (ei->is_border_edge())
			{
				if (ei->is_border())
				{
					hh = &(*ei);
				}
				else
				{
					hh = ei->opposite();
				}
			}
			else
			{
				hh = &(*ei);
			}
			while (hh->vertex()->nCrevasseSize == 2 && hh->vertex()->nNumPorts == 0)
			{
				hh = hh->hhCrevasseNext;
			}
			do
			{
				hh->iCrevasseId = iCrevasseId;
				newCrevasse.m_vechhElements.push_back(hh);
				hh = hh->prev()->hhCrevassePrev;
			} while (hh->vertex()->nCrevasseSize == 2 && hh->vertex()->nNumPorts == 0);
			memset(newCrevasse.lpftOffset, 0, 2 * sizeof(CMeshBase::FT));
			m_vecCrevasses.push_back(newCrevasse);
			++iCrevasseId;
			newCrevasse.m_vechhElements.clear();
			hh = hh->hhCrevasseNext->opposite();
			do
			{
				hh->iCrevasseId = iCrevasseId;
				newCrevasse.m_vechhElements.push_back(hh);
				hh = hh->prev()->hhCrevassePrev;
			} while (hh->vertex()->nCrevasseSize == 2 && hh->vertex()->nNumPorts == 0);
			m_vecCrevasses.push_back(newCrevasse);
			newCrevasse.m_vechhElements.clear();
			++iCrevasseId;
		}
	}
	for (iterCrevasse = m_vecCrevasses.begin(); iterCrevasse != m_vecCrevasses.end(); ++iterCrevasse)
	{
		iterCrevasse->m_lpPrev = &(m_vecCrevasses[iterCrevasse->last_element()->prev()->hhCrevassePrev->iCrevasseId]);
		iterCrevasse->m_lpNext = &(m_vecCrevasses[iterCrevasse->first_element()->hhCrevasseNext->iCrevasseId]);
		iterCrevasse->m_lpOpposite = &(m_vecCrevasses[iterCrevasse->first_element()->opposite()->iCrevasseId]);
	}
}



void CParameterGenerator::m_fnUniformizeChart()
{
	CMeshBase::CPolyhedron::Facet_handle fh0;
	std::queue<CMeshBase::CPolyhedron::Facet_handle> qufhQueue;
	CMeshBase::CPolyhedron::Halfedge_handle hh0, hh;
	CMeshBase::CPolyhedron::Facet_iterator fi;
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
	{
		fi->bUnif = false;
	}
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
	{
		fh0 = &(*fi);
		if (!fh0->bUnif)
		{
			fh0->bUnif = true;
			qufhQueue.push(fh0);
			while (!qufhQueue.empty())
			{
				fh0 = qufhQueue.front();
				qufhQueue.pop();
				hh0 = fh0->halfedge();
				hh = hh0;
				do
				{
					if (!(hh->bFreeEdge || hh->opposite()->facet()->bUnif))
					{
						hh->opposite()->facet()->ftChartDir -= hh->iSeamType * ftHalfPi;
						hh->opposite()->facet()->ftChartDir = hh->opposite()->facet()->ftChartDir - CMeshBase::FT(int(hh->opposite()->facet()->ftChartDir / ftTwoPi + 15.5) - 15) * ftTwoPi;
						hh->opposite()->facet()->ftFrameDir -= hh->iSeamType * ftHalfPi;
						hh->opposite()->facet()->ftFrameDir = hh->opposite()->facet()->ftFrameDir - CMeshBase::FT(int(hh->opposite()->facet()->ftFrameDir / ftTwoPi + 15.5) - 15) * ftTwoPi;
						hh->opposite()->next()->iSeamType = (hh->opposite()->next()->iSeamType + hh->iSeamType + 6) % 4 - 2;
						hh->opposite()->next()->opposite()->iSeamType = (6 - hh->opposite()->next()->iSeamType) % 4 - 2;
						hh->opposite()->prev()->iSeamType = (hh->opposite()->prev()->iSeamType + hh->iSeamType + 6) % 4 - 2;
						hh->opposite()->prev()->opposite()->iSeamType = (6 - hh->opposite()->prev()->iSeamType) % 4 - 2;
						hh->iSeamType = 0;
						hh->opposite()->iSeamType = 0;
						hh->opposite()->facet()->bUnif = true;
						qufhQueue.push(hh->opposite()->facet());
					}
					hh = hh->next();
				} while (hh != hh0);
			}
		}
	}
}

int CParameterGenerator::m_fnInitGlobalPara(CMeshBase::FT*&lpftGlobalParameter, int *&lpiGlobalIndex)
{
	int iIndex;
	CMeshBase::CPolyhedron::Vertex_iterator vi;
	CMeshBase::CPolyhedron::Halfedge_handle hh0, hhStart, hhEnd, hh;
	iIndex = 0;
	for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
	{
		if (vi->nCrevasseSize < 2)
		{
			++iIndex;
		}
		else
		{
			iIndex += vi->nCrevasseSize;
		}
	}
	lpftGlobalParameter = new CMeshBase::FT[iIndex * 2];
	lpiGlobalIndex = new int[iIndex * 2];
	memset(lpftGlobalParameter, 0, iIndex * 2 * sizeof(CMeshBase::FT));
	iIndex = 0;
	for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
	{
		if (vi->nCrevasseSize < 2)
		{
			hh0 = vi->halfedge();
			hh = hh0;
			do
			{
				hh->lpftGlobalPara = &(lpftGlobalParameter[iIndex * 2]);
				hh->lpiGlobalIndex = &(lpiGlobalIndex[iIndex * 2]);
				hh = hh->next()->opposite();
			} while (hh != hh0);
			++iIndex;
		}
		else
		{
			hh0 = vi->halfedge()->hhCrevassePrev;
			hhStart = hh0;
			do
			{
				hhEnd = hhStart->hhCrevasseNext->opposite();
				hh = hhStart;
				while (hh != hhEnd)
				{
					hh->lpftGlobalPara = &(lpftGlobalParameter[iIndex * 2]);
					hh->lpiGlobalIndex = &(lpiGlobalIndex[iIndex * 2]);
					hh = hh->next()->opposite();
				}
				++iIndex;
				hhStart = hhEnd;
			} while (hhStart != hh0);
		}
	}
	return iIndex;
}

void CParameterGenerator::m_fnSolveCrevasseParaRange()
{
	CMeshBase::CPolyhedron::Vertex_iterator vi;
	std::vector<CMeshCenter::CCrevasse>::iterator iterCrevasse;
	CMeshCenter::CCrevasse* lpCrevasse0, * lpCrevasse;
	CMeshBase::CPolyhedron::Halfedge_handle hh;
	std::vector<CMeshBase::CPolyhedron::Halfedge_handle>::iterator iterHH;
	CMeshBase::CMinNormSv MNS;
	MatElem<CMeshBase::FT> me;
	CMeshBase::FT ftChartDir, ftParaScale, *lpftX;
	CMeshBase::CPolyhedron::Vector_3 vtChart;
	int iEqIndex, i;
	iEqIndex = 0;
	CMeshBase::CPolyhedron::Edge_iterator ei;
	for (iterCrevasse = m_vecCrevasses.begin(); iterCrevasse != m_vecCrevasses.end(); ++iterCrevasse)
	{
		iterCrevasse->lpiEqInd[0] = -1;
		iterCrevasse->lpiEqInd[1] = -1;
	}
	for (iterCrevasse = m_vecCrevasses.begin(); iterCrevasse != m_vecCrevasses.end(); ++iterCrevasse)
	{
		if (!iterCrevasse->is_border())
		{
			if (iterCrevasse->opposite()->is_border())
			{
				iterCrevasse->lpiEqInd[(iterCrevasse->seamtype() + 3) % 2] = iEqIndex;
				++iEqIndex;
			}
			else
			{
				if (iterCrevasse->lpiEqInd[0] == -1 && iterCrevasse->opposite()->lpiEqInd[(iterCrevasse->seamtype() + 2) % 2] == -1)
				{
					iterCrevasse->lpiEqInd[0] = iEqIndex;
					++iEqIndex;
				}
				if (iterCrevasse->lpiEqInd[1] == -1 && iterCrevasse->opposite()->lpiEqInd[(iterCrevasse->seamtype() + 3) % 2] == -1)
				{
					iterCrevasse->lpiEqInd[1] = iEqIndex;
					++iEqIndex;
				}
			}
		}
	}
	lpftX = new CMeshBase::FT[iEqIndex];
	memset(lpftX, 0, iEqIndex * sizeof(CMeshBase::FT));
	MNS.init(2, iEqIndex);
	lpCrevasse0 = NULL;
	for (iterCrevasse = m_vecCrevasses.begin(); iterCrevasse != m_vecCrevasses.end() && lpCrevasse0 == NULL; ++iterCrevasse)
	{
		if (!iterCrevasse->is_border())
		{
			lpCrevasse0 = &(*iterCrevasse);
		}
	}
	lpCrevasse = lpCrevasse0;
	do
	{
		for (i = 0; i < 2; ++i)
		{
			if (lpCrevasse->lpiEqInd[i] != -1)
			{
				for (iterHH = lpCrevasse->m_vechhElements.begin(); iterHH != lpCrevasse->m_vechhElements.end(); ++iterHH)
				{
					hh = (*iterHH);
					if (hh->opposite()->is_border())
					{
						ftChartDir = hh->facet()->ftChartDir;
					}
					else
					{
						ftChartDir = hh->facet()->ftChartDir - hh->ftMom * (hh->opposite()->next()->ftDot / (hh->opposite()->facet()->ftTwArea * hh->ftLaplaceCoef));
					}
					vtChart = (hh->facet()->vtPrincipalAxis * cos(ftChartDir + ftHalfPi * (i + 1)) + hh->facet()->vtAuxiliaryAxis * sin(ftChartDir + ftHalfPi * (i + 1))) / hh->facet()->ftRadius;
					ftParaScale = exp(hh->vertex()->ftDensity + hh->prev()->vertex()->ftDensity);
					lpftX[lpCrevasse->lpiEqInd[i]] += ftParaScale * (vtChart * hh->vtVector);
				}
			}
		}
		lpCrevasse = lpCrevasse->next();
	} while (lpCrevasse != lpCrevasse0);
	lpCrevasse = lpCrevasse0;
	do
	{
		for (i = 0; i < 2; ++i)
		{
			me.row = i;
			if (lpCrevasse->lpiEqInd[i] != -1)
			{
				me.col = lpCrevasse->lpiEqInd[i];
				me.data = 1.0;
				MNS.add_coef(me.row, me.col, me.data);
			}
			if (lpCrevasse->opposite()->lpiEqInd[(i + lpCrevasse->seamtype() + 2) % 2] != -1)
			{
				me.col = lpCrevasse->opposite()->lpiEqInd[(i + lpCrevasse->seamtype() + 2) % 2];
				me.data = 1.0 * (((lpCrevasse->seamtype() - i + 5) % 4) / 2 * 2 - 1);
				MNS.add_coef(me.row, me.col, me.data);
			}
		}
		lpCrevasse = lpCrevasse->next();
	} while (lpCrevasse != lpCrevasse0);
	MNS.m_b[0] = 0.0;
	MNS.m_b[1] = 0.0;
	MNS.decompose();
	for (i = 0; i < MNS.m_mat_size; ++i)
	{
		MNS.m_b[MNS.m_ordered_mat[i].row] += (MNS.m_ordered_mat[i].data * lpftX[MNS.m_ordered_mat[i].col]);
	}
	MNS.solve();
	for (iterCrevasse = m_vecCrevasses.begin(); iterCrevasse != m_vecCrevasses.end(); ++iterCrevasse)
	{
		lpCrevasse = iterCrevasse->opposite();
		if (iterCrevasse->opposite()->is_border())
		{
			iterCrevasse->lpftParaRange[(iterCrevasse->seamtype() + 2) % 2] = 0.0;
			iterCrevasse->opposite()->lpftParaRange[0] = 0.0;
		}
		for (i = 0; i < 2; ++i)
		{
			if (iterCrevasse->lpiEqInd[i] != -1)
			{
				iterCrevasse->lpftParaRange[i] = lpftX[iterCrevasse->lpiEqInd[i]] - MNS.m_x[iterCrevasse->lpiEqInd[i]];
				iterCrevasse->opposite()->lpftParaRange[(iterCrevasse->seamtype() + i + 2) % 2] = iterCrevasse->lpftParaRange[i] * (1 - (iterCrevasse->seamtype() + 7 - i) % 4 / 2 * 2);
			}
		}
	}
	MNS.clear();
	delete[]lpftX;
}

void CParameterGenerator::m_fnAssignCrevassePara()
{
	std::vector<CMeshCenter::CCrevasse>::iterator iterCrevasse;
	CMeshCenter::CCrevasse* lpCrevasse0, * lpCrevasse;
	CMeshBase::FT lpftCrevassePara[2], lpftVertexPara[2], ftChartDir, ftParaScale;
	std::complex<CMeshBase::FT> ftcZ0, ftcZ1, ftcTrans;
	CMeshBase::CPolyhedron::Halfedge_handle hh;
	std::vector<CMeshBase::CPolyhedron::Halfedge_handle>::iterator iterHH;
	CMeshBase::Vector_3  vtChart;
	int i;
	lpCrevasse0 = NULL;
	for (iterCrevasse = m_vecCrevasses.begin(); iterCrevasse != m_vecCrevasses.end() && lpCrevasse0 == NULL; ++iterCrevasse)
	{
		if (!iterCrevasse->is_border())
		{
			lpCrevasse0 = &(*iterCrevasse);
		}
	}
	lpftCrevassePara[0] = lpftCrevassePara[1] = 0.0;
	lpCrevasse = lpCrevasse0;
	do
	{
		lpftVertexPara[0] = 0;
		lpftVertexPara[1] = 0;
		for (iterHH = lpCrevasse->m_vechhElements.begin(); iterHH != lpCrevasse->m_vechhElements.end(); ++iterHH)
		{
			hh = (*iterHH);
			hh->lpftGlobalPara[0] = lpftVertexPara[0];
			hh->lpftGlobalPara[1] = lpftVertexPara[1];
			for (i = 0; i < 2; ++i)
			{
				if (hh->opposite()->is_border())
				{
					ftChartDir = hh->facet()->ftChartDir;
				}
				else
				{
					ftChartDir = hh->facet()->ftChartDir - hh->ftMom * (hh->opposite()->next()->ftDot / (hh->opposite()->facet()->ftTwArea * hh->ftLaplaceCoef));
				}
				vtChart = (hh->facet()->vtPrincipalAxis * cos(ftChartDir + ftHalfPi * (i + 1)) + hh->facet()->vtAuxiliaryAxis * sin(ftChartDir + ftHalfPi * (i + 1))) / hh->facet()->ftRadius;
				//if ((hh->vertex()->iDegreeBias == 0) == (hh->prev()->vertex()->iDegreeBias == 0))
				{
					ftParaScale = exp(hh->vertex()->ftDensity + hh->prev()->vertex()->ftDensity);
				}
				//else if ((hh->vertex()->iDegreeBias == 0) && (hh->prev()->vertex()->iDegreeBias != 0))
				//{
				//	ftParaScale = exp(hh->vertex()->ftDensity * 2);
				//}
				//else if ((hh->vertex()->iDegreeBias != 0) && (hh->prev()->vertex()->iDegreeBias == 0))
				//{
				//	ftParaScale = exp(hh->prev()->vertex()->ftDensity * 2);
				//}
				lpftVertexPara[i] -= ftParaScale * (vtChart * hh->vtVector);
			}
		}
		ftcZ0.real(-lpftVertexPara[0]);
		ftcZ0.imag(-lpftVertexPara[1]);
		ftcZ1.real(lpCrevasse->lpftParaRange[0]);
		ftcZ1.imag(lpCrevasse->lpftParaRange[1]);
		ftcTrans = ftcZ1 / ftcZ0;
		for (iterHH = lpCrevasse->m_vechhElements.begin(); iterHH != lpCrevasse->m_vechhElements.end(); ++iterHH)
		{
			hh = (*iterHH);
			ftcZ0.real(hh->lpftGlobalPara[0]);
			ftcZ0.imag(hh->lpftGlobalPara[1]);
			ftcZ1 = ftcZ0 * ftcTrans;
			hh->lpftGlobalPara[0] = lpftCrevassePara[0] + ftcZ1.real();
			hh->lpftGlobalPara[1] = lpftCrevassePara[1] + ftcZ1.imag();
		}
		lpftCrevassePara[0] = lpftCrevassePara[0] - lpCrevasse->lpftParaRange[0];
		lpftCrevassePara[1] = lpftCrevassePara[1] - lpCrevasse->lpftParaRange[1];
		lpCrevasse = lpCrevasse->prev();
	} while (lpCrevasse != lpCrevasse0);
}

void CParameterGenerator::m_fnRescaleCrevassePara(int nFacets)
{
	std::vector<CMeshCenter::CCrevasse>::iterator iterCrevasse;
	CMeshCenter::CCrevasse* lpCrevasse0, * lpCrevasse;
	CMeshBase::CPolyhedron::Halfedge_handle hh;
	std::vector<CMeshBase::CPolyhedron::Halfedge_handle>::iterator iterHH;
	CMeshBase::FT ftArea, lpftLocalDif[2], ftCoef;
	lpCrevasse0 = NULL;
	for (iterCrevasse = m_vecCrevasses.begin(); iterCrevasse != m_vecCrevasses.end() && lpCrevasse0 == NULL; ++iterCrevasse)
	{
		if (!iterCrevasse->is_border())
		{
			lpCrevasse0 = &(*iterCrevasse);
		}
	}
	ftArea = 0;
	lpCrevasse = lpCrevasse0;
	do
	{
		for (iterHH = lpCrevasse->m_vechhElements.begin(); iterHH != lpCrevasse->m_vechhElements.end(); ++iterHH)
		{
			hh = *iterHH;
			lpftLocalDif[0] = hh->lpftGlobalPara[0] - hh->prev()->lpftGlobalPara[0];
			lpftLocalDif[1] = hh->lpftGlobalPara[1] - hh->prev()->lpftGlobalPara[1];
			ftArea += (hh->lpftGlobalPara[0] * lpftLocalDif[1] - hh->lpftGlobalPara[1] * lpftLocalDif[0]);
		}
		lpCrevasse = lpCrevasse->prev();
	} while (lpCrevasse != lpCrevasse0);
	ftCoef = sqrt(nFacets * 2 / fabs(ftArea));
	do
	{
		for (iterHH = lpCrevasse->m_vechhElements.begin(); iterHH != lpCrevasse->m_vechhElements.end(); ++iterHH)
		{
			hh = (*iterHH);
			hh->lpftGlobalPara[0] = hh->lpftGlobalPara[0] * ftCoef;
			hh->lpftGlobalPara[1] = hh->lpftGlobalPara[1] * ftCoef;
		}
		lpCrevasse = lpCrevasse->prev();
	} while (lpCrevasse != lpCrevasse0);
	do
	{
		lpCrevasse->lpftParaRange[0] = lpCrevasse->first_element()->lpftGlobalPara[0] - lpCrevasse->prev()->first_element()->lpftGlobalPara[0];
		lpCrevasse->lpftParaRange[1] = lpCrevasse->first_element()->lpftGlobalPara[1] - lpCrevasse->prev()->first_element()->lpftGlobalPara[1];
		lpCrevasse = lpCrevasse->prev();
	} while (lpCrevasse != lpCrevasse0);
}

void CParameterGenerator::m_fnSolveHarmonicPara()
{
	int iIndex;
	CMeshBase::CLDLTSv LDLTS;
	CMeshBase::CPolyhedron::Vertex_iterator vi;
	CMeshBase::CPolyhedron::Halfedge_handle hh0, hh;
	MatElem<CMeshBase::FT> meDiag, meTri;
	CMeshBase::FT *lpft_b[2];
	iIndex = 0;
	for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
	{
		if (vi->nCrevasseSize != 0)
		{
			vi->iEquationIndex = -1;
		}
		else
		{
			vi->iEquationIndex = iIndex;
			++iIndex;
		}
	}
	lpft_b[0] = new CMeshBase::FT[iIndex];
	lpft_b[1] = new CMeshBase::FT[iIndex];
	LDLTS.init(iIndex);
	for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
	{
		if (vi->iEquationIndex != -1)
		{
			meDiag.row = vi->iEquationIndex;
			meDiag.col = vi->iEquationIndex;
			meDiag.data = 0.0;
			lpft_b[0][meDiag.row] = 0.0;
			lpft_b[1][meDiag.row] = 0.0;
			hh0 = vi->halfedge();
			hh = hh0;
			do
			{
				if (hh->prev()->vertex()->iEquationIndex == -1)
				{
					lpft_b[0][meDiag.row] += hh->prev()->lpftGlobalPara[0] * hh->ftLaplaceCoef;
					lpft_b[1][meDiag.row] += hh->prev()->lpftGlobalPara[1] * hh->ftLaplaceCoef;
				}
				else
				{
					meTri.row = vi->iEquationIndex;
					meTri.col = hh->prev()->vertex()->iEquationIndex;
					meTri.data = -hh->ftLaplaceCoef;
					LDLTS.add_data(meTri.row, meTri.col, meTri.data);
				}
				meDiag.data += hh->ftLaplaceCoef;
				hh = hh->next()->opposite();
			} while (hh != hh0);
			LDLTS.add_data(meDiag.row, meDiag.col, meDiag.data);
		}
	}
	LDLTS.decompose();
	memcpy(LDLTS.m_b, lpft_b[0], iIndex * sizeof(CMeshBase::FT));
	LDLTS.solve();
	for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
	{
		if (vi->iEquationIndex != -1)
		{
			vi->halfedge()->lpftGlobalPara[0] = LDLTS.m_x[vi->iEquationIndex];
		}
	}
	memcpy(LDLTS.m_b, lpft_b[1], iIndex * sizeof(CMeshBase::FT));
	LDLTS.solve();
	for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
	{
		if (vi->iEquationIndex != -1)
		{
			vi->halfedge()->lpftGlobalPara[1] = LDLTS.m_x[vi->iEquationIndex];
		}
	}
	LDLTS.clear();
	delete[]lpft_b[1];
	delete[]lpft_b[0];
}

void CParameterGenerator::m_fnRescaleGlobal(int nFacets)
{
	CMeshBase::CPolyhedron::Facet_iterator fi;
	CMeshBase::CPolyhedron::Halfedge_handle hh0, hh;
	CMeshBase::CPolyhedron::Vertex_iterator vi;
	int i;
	CMeshBase::FT lpftLocalPara[4], ftParaArea, ftScale;
	ftParaArea = 0;
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
	{
		hh0 = fi->halfedge();
		for (i = 0; i < 2; ++i)
		{
			lpftLocalPara[i] = hh0->lpftGlobalPara[i] - hh0->prev()->lpftGlobalPara[i];
			lpftLocalPara[i + 2] = hh0->next()->lpftGlobalPara[i] - hh0->lpftGlobalPara[i];
		}
		ftParaArea += (lpftLocalPara[0] * lpftLocalPara[3] - lpftLocalPara[1] * lpftLocalPara[2]);
	}
	ftParaArea /= 2;
	ftScale = sqrt(nFacets / ftParaArea);
	for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
	{
		if (vi->iBordId == -1)
		{
			if (vi->nCrevasseSize < 2)
			{
				hh0 = vi->halfedge();
				hh0->lpftGlobalPara[0] *= ftScale;
				hh0->lpftGlobalPara[1] *= ftScale;
			}
			else
			{
				hh0 = vi->halfedge()->hhCrevassePrev;
				hh = hh0;
				do
				{
					hh->lpftGlobalPara[0] *= ftScale;
					hh->lpftGlobalPara[1] *= ftScale;
					hh = hh->hhCrevasseNext->opposite();
				} while (hh != hh0);
			}
		}
		else
		{
			hh0 = vi->halfedge();
			while (!hh0->is_border())
			{
				hh0 = hh0->next()->opposite();
			}
			hh = hh0->next()->opposite();
			while (hh != hh0)
			{
				hh->lpftGlobalPara[0] *= ftScale;
				hh->lpftGlobalPara[1] *= ftScale;
				hh = hh->hhCrevasseNext->opposite();
			}
		}
	}
}
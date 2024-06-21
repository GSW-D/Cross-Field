#include "CParameterAdjustor.h"

CParameterAdjustor::CParameterAdjustor(CMeshBase::CPolyhedron &mshSurface, std::vector<CMeshCenter::CCrevasse> &vecCrevasses): m_mshSurface(mshSurface), m_vecCrevasses(vecCrevasses)
{
	m_lpftRelativeOffset = NULL;
    //ctor
}

CParameterAdjustor::~CParameterAdjustor()
{
    //dtor
}

void CParameterAdjustor::m_fnIntegerizeCrevasseSpan()
{
	std::vector<CMeshCenter::CCrevasse>::iterator iterCrevasse;
	CMeshCenter::CCrevasse *lpCrevasse0, *lpCrevasse, *lpCrevasse1, *lpCrevasseAdjust;
	int lpiConstraint[2], i, j, iIntGroupInd, iHeadSize, iAdjustPara;
	CMeshBase::FT lpftSpanSum[2], lpftDisturb[2], ftAdjust, ftZeroDist, ftMinZeroDist, ftErr, ftMinErr, ftParaLen, ftBestAdjust, lpftLeastDisturb[2];
	lpCrevasse0 = NULL;
	for (iterCrevasse = m_vecCrevasses.begin(); iterCrevasse != m_vecCrevasses.end() && lpCrevasse0 == NULL; ++iterCrevasse)
	{
		if (!iterCrevasse->is_border())
		{
			lpCrevasse0 = &(*iterCrevasse);
		}
	}
	lpiConstraint[0] = lpiConstraint[1] = 1;
	lpCrevasse = lpCrevasse0;
	do
	{
		lpCrevasse->lpiIntHead[0] = lpCrevasse->lpiIntHead[1] = 1;
		lpCrevasse->lpiIntGroupInd[0] = lpCrevasse->lpiIntGroupInd[0] = -1;
		lpCrevasse->lpiIntGroupSize[0] = lpCrevasse->lpiIntGroupSize[1] = 0;
		lpCrevasse->lpftIntRange[0] = lpCrevasse->lpftIntRange[1] = 0.0;
		if (lpCrevasse->opposite()->is_border())
		{
			lpiConstraint[(lpCrevasse->seamtype() + 3) % 2] = 0;
			if (lpCrevasse->first_element()->vertex()->nNumPorts == 0)
			{
				lpCrevasse->lpiIntHead[(lpCrevasse->seamtype() + 3) % 2] = 0;
			}
		}
		if (lpCrevasse->next()->opposite()->is_border())
		{
			if (lpCrevasse->first_element()->vertex()->nNumPorts == 0)
			{
				lpCrevasse->lpiIntHead[(lpCrevasse->next()->seamtype() + 3) % 2] = 0;
			}
		}
		if (lpCrevasse->first_element()->vertex()->iBordId == -1 && lpCrevasse->first_element()->vertex()->nNumPorts == 0)
		{
			lpCrevasse->lpiIntHead[0] = lpCrevasse->lpiIntHead[1] = 0;
		}
		lpCrevasse = lpCrevasse->next();
	} while (lpCrevasse != lpCrevasse0);
	for (i = 0; i < 2; ++i)
	{
		iHeadSize = 0;
		lpCrevasse = lpCrevasse0;
		do
		{
			if (lpCrevasse->lpiIntHead[i] == 1)
			{
				++iHeadSize;
			}
			lpCrevasse = lpCrevasse->next();
		} while (lpCrevasse != lpCrevasse0);
		if (iHeadSize == 0)
		{
			lpCrevasse0->lpiIntHead[i] = 1;
		}
	}
	iIntGroupInd = 0;
	lpCrevasse = lpCrevasse0;
	do
	{
		for (i = 0; i < 2; ++i)
		{
			if (lpCrevasse->lpiIntHead[i] == 1 && lpCrevasse->lpiIntGroupSize[i] == 0)
			{
				lpCrevasse->lpftGroupRange[i] = 0.0;
				lpCrevasse1 = lpCrevasse;
				do
				{
					++lpCrevasse->lpiIntGroupSize[i];
					lpCrevasse->lpftGroupRange[i] += lpCrevasse1->lpftParaRange[i];
					lpCrevasse1 = lpCrevasse1->prev();
				} while (lpCrevasse1->lpiIntHead[i] == 0);
				lpCrevasse->lpftIntRange[i] = round(lpCrevasse->lpftGroupRange[i]);
				lpCrevasse1 = lpCrevasse;
				do
				{
					lpCrevasse1->lpiIntGroupSize[i] = lpCrevasse->lpiIntGroupSize[i];
					lpCrevasse1->lpiIntGroupInd[i] = iIntGroupInd;
					lpCrevasse1 = lpCrevasse1->prev();
				} while (lpCrevasse1->lpiIntHead[i] == 0);
				if (lpCrevasse->lpiIntGroupSize[i] == 1 && !lpCrevasse->opposite()->is_border())
				{
					lpCrevasse->opposite()->lpftGroupRange[(i + lpCrevasse->seamtype() + 2) % 2] = lpCrevasse->lpftGroupRange[i] * (1 - (lpCrevasse->seamtype() + 7 - i) % 4 / 2 * 2);
					lpCrevasse->opposite()->lpftIntRange[(i + lpCrevasse->seamtype() + 2) % 2] = lpCrevasse->lpftIntRange[i] * (1 - (lpCrevasse->seamtype() + 7 - i) % 4 / 2 * 2);
					lpCrevasse->opposite()->lpiIntGroupSize[(i + lpCrevasse->seamtype() + 2) % 2] = 1;
				}
				++iIntGroupInd;
			}
		}
		lpCrevasse = lpCrevasse->next();
	} while (lpCrevasse != lpCrevasse0);
	lpftSpanSum[0] = lpftSpanSum[1] = 0.0;
	lpCrevasse = lpCrevasse0;
	do
	{
		for (i = 0; i < 2; ++i)
		{
			if (lpCrevasse->lpiIntHead[i] == 1)
			{
				lpftSpanSum[i] += lpCrevasse->lpftIntRange[i];
			}
		}
		lpCrevasse = lpCrevasse->next();
	} while (lpCrevasse != lpCrevasse0);
	while (lpftSpanSum[0] > 0.5 || lpftSpanSum[0] < -0.5 || lpftSpanSum[1] > 0.5 || lpftSpanSum[1] < -0.5)
	{
		ftMinZeroDist = fabs(lpftSpanSum[0]) + fabs(lpftSpanSum[1]) - 1.0;
		ftMinErr = DBL_MAX;
		lpCrevasseAdjust = NULL;
		lpCrevasse = lpCrevasse0;
		do
		{
			lpCrevasse->lpiIntMark[0] = lpCrevasse->lpiIntMark[1] = 0;
			lpCrevasse = lpCrevasse->next();
		} while (lpCrevasse != lpCrevasse0);
		do
		{
			for (i = 0; i < 2; ++i)
			{
				if (lpCrevasse->lpiIntHead[i] == 1 && lpCrevasse->lpiIntMark[i] == 0)
				{
					if (lpCrevasse->lpiIntGroupSize[i] > 1)
					{
						lpftDisturb[i] = 1.0;
						lpftDisturb[1 - i] = 0.0;
					}
					else
					{
						switch (lpCrevasse->seamtype())
						{
						case -2:
							lpftDisturb[i] = 2.0;
							lpftDisturb[1 - i] = 0.0;
							break;
						case -1:
							lpftDisturb[i] = 1.0;
							lpftDisturb[1 - i] = 1.0 * (i * 2 - 1);
							break;
						case 0:
							lpftDisturb[i] = 0.0;
							lpftDisturb[1 - i] = 0.0;
							break;
						case 1:
							lpftDisturb[i] = 1.0;
							lpftDisturb[1 - i] = 1.0 * (1 - i * 2);
							break;
						}
					}
					if(fabs(lpftSpanSum[0] - lpftDisturb[0]) + fabs(lpftSpanSum[1] - lpftDisturb[1]) < fabs(lpftSpanSum[0] + lpftDisturb[0]) + fabs(lpftSpanSum[1] + lpftDisturb[1]))
					{
						ftAdjust = -1.0;
						ftZeroDist = fabs(lpftSpanSum[0] - lpftDisturb[0]) + fabs(lpftSpanSum[1] - lpftDisturb[1]);
						lpftDisturb[0] = -lpftDisturb[0];
						lpftDisturb[1] = -lpftDisturb[1];
					}
					else
					{
						ftAdjust = 1.0;
						ftZeroDist = fabs(lpftSpanSum[0] + lpftDisturb[0]) + fabs(lpftSpanSum[1] + lpftDisturb[1]);
					}
					if (ftZeroDist < ftMinZeroDist + 0.5)
					{
						ftMinZeroDist = ftZeroDist;
						ftParaLen = sqrt(lpCrevasse->lpftGroupRange[0] * lpCrevasse->lpftGroupRange[0] + lpCrevasse->lpftGroupRange[1] * lpCrevasse->lpftGroupRange[1]);
						ftErr = (fabs(lpCrevasse->lpftIntRange[i] - lpCrevasse->lpftGroupRange[i] + ftAdjust) - fabs(lpCrevasse->lpftIntRange[i] - lpCrevasse->lpftGroupRange[i])) / ftParaLen;
						if (ftErr < ftMinErr)
						{
							ftMinErr = ftErr;
							lpCrevasseAdjust = lpCrevasse;
							iAdjustPara = i;
							ftBestAdjust = ftAdjust;
							lpftLeastDisturb[0] = lpftDisturb[0];
							lpftLeastDisturb[1] = lpftDisturb[1];
						}
					}
					lpCrevasse->lpiIntMark[i] = 1;
					if (lpCrevasse->lpiIntGroupSize[i] > 1)
					{
						lpCrevasse->opposite()->lpiIntMark[(lpCrevasse->seamtype() + i + 2) % 2] = 1;
					}
				}
			}
			lpCrevasse = lpCrevasse->next();
		} while (lpCrevasse != lpCrevasse0);
		if (lpCrevasseAdjust != NULL)
		{
			lpCrevasseAdjust->lpftIntRange[iAdjustPara] += ftBestAdjust;
			if (lpCrevasseAdjust->lpiIntGroupSize[iAdjustPara] == 1)
			{
				lpCrevasseAdjust->opposite()->lpftIntRange[(lpCrevasseAdjust->seamtype() + iAdjustPara + 2) % 2] += ftBestAdjust * (1 - (lpCrevasseAdjust->seamtype() + 7 - iAdjustPara) % 4 / 2 * 2);
			}
			lpftSpanSum[0] += lpftLeastDisturb[0];
			lpftSpanSum[1] += lpftLeastDisturb[1];
		}
	}
	lpCrevasse = lpCrevasse0;
	do
	{
		for (i = 0; i < 2; ++i)
		{
			if (lpCrevasse->lpiIntHead[i] == 1 && lpCrevasse->lpiIntGroupSize[i] > 1)
			{
				lpftSpanSum[i] = lpCrevasse->lpftIntRange[i];
				lpCrevasse1 = lpCrevasse;
				for (j = 0; j < lpCrevasse->lpiIntGroupSize[i]; ++j)
				{
					if (!lpCrevasse1->opposite()->is_border())
					{
						lpCrevasse1->lpftIntRange[i] = lpCrevasse1->lpftParaRange[i];
						lpftSpanSum[i] -= lpCrevasse1->lpftIntRange[i];
					}
					lpCrevasse1 = lpCrevasse1->prev();
				}
				lpCrevasse1 = lpCrevasse;
				for (j = 0; j < lpCrevasse->lpiIntGroupSize[i]; ++j)
				{
					if (lpCrevasse1->opposite()->is_border())
					{
						lpCrevasse1->lpftIntRange[i] = lpftSpanSum[i];
					}
					lpCrevasse1 = lpCrevasse1->prev();
				}
			}
		}
		lpCrevasse = lpCrevasse->next();
	} while (lpCrevasse != lpCrevasse0);
	lpftSpanSum[0] = 0;
	lpftSpanSum[1] = 0;
	lpftDisturb[0] = 0;
	lpftDisturb[1] = 0;
	lpCrevasse = lpCrevasse0;
	do
	{
		lpftDisturb[0] += lpCrevasse->lpftParaRange[0];
		lpftDisturb[1] += lpCrevasse->lpftParaRange[1];
		lpftSpanSum[0] += lpCrevasse->lpftIntRange[0];
		lpftSpanSum[1] += lpCrevasse->lpftIntRange[1];
		lpCrevasse = lpCrevasse->next();
	} while (lpCrevasse != lpCrevasse0);
}

void CParameterAdjustor::m_fnAssignCrevasseIntPara()
{
	std::vector<CMeshCenter::CCrevasse>::iterator iterCrevasse;
	CMeshCenter::CCrevasse* lpCrevasse0, * lpCrevasse;
	CMeshBase::FT lpftCrevassePara[2];
	std::complex<CMeshBase::FT> ftcZ0, ftcZ1, ftcTrans;
	CMeshBase::CPolyhedron::Halfedge_handle hh;
	std::vector<CMeshBase::CPolyhedron::Halfedge_handle>::reverse_iterator riterHH;
	lpCrevasse0 = NULL;
	for (iterCrevasse = m_vecCrevasses.begin(); iterCrevasse != m_vecCrevasses.end() && lpCrevasse0 == NULL; ++iterCrevasse)
	{
		if (!iterCrevasse->is_border() && iterCrevasse->lpiIntHead[0] == 1 && iterCrevasse->lpiIntHead[1] == 1)
		{
			lpCrevasse0 = &(*iterCrevasse);
		}
	}
	lpftCrevassePara[0] = lpftCrevassePara[1] = 0.0;
	lpCrevasse = lpCrevasse0;
	do
	{
		ftcZ0.real(lpCrevasse->lpftParaRange[0]);
		ftcZ0.imag(lpCrevasse->lpftParaRange[1]);
		ftcZ1.real(lpCrevasse->lpftIntRange[0]);
		ftcZ1.imag(lpCrevasse->lpftIntRange[1]);
		ftcTrans = ftcZ1 / ftcZ0;
		for (riterHH = lpCrevasse->m_vechhElements.rbegin(); riterHH != lpCrevasse->m_vechhElements.rend(); ++riterHH)
		{
			hh = (*riterHH);
			ftcZ0.real(hh->lpftGlobalPara[0] - lpCrevasse->first_element()->lpftGlobalPara[0]);
			ftcZ0.imag(hh->lpftGlobalPara[1] - lpCrevasse->first_element()->lpftGlobalPara[1]);
			ftcZ1 = ftcZ0 * ftcTrans;
			hh->lpftGlobalPara[0] = lpftCrevassePara[0] + ftcZ1.real();
			hh->lpftGlobalPara[1] = lpftCrevassePara[1] + ftcZ1.imag();
		}
		lpftCrevassePara[0] = lpftCrevassePara[0] - lpCrevasse->lpftIntRange[0];
		lpftCrevassePara[1] = lpftCrevassePara[1] - lpCrevasse->lpftIntRange[1];
		lpCrevasse = lpCrevasse->prev();
	} while (lpCrevasse != lpCrevasse0);
}

void CParameterAdjustor::m_fnSolveHarmonicPara()
{
	int iIndex;
	CMeshBase::CLDLTSv LDLTS;
	CMeshBase::CPolyhedron::Vertex_iterator vi;
	CMeshBase::CPolyhedron::Halfedge_handle hh0, hh;
	MatElem<CMeshBase::FT> meDiag, meTri;
	CMeshBase::FT* lpft_b[2];
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

void CParameterAdjustor::m_fnExtractCrevasseOffset()
{
	std::vector<CMeshCenter::CCrevasse>::iterator iterCrevasse;
	for (iterCrevasse = m_vecCrevasses.begin(); iterCrevasse != m_vecCrevasses.end(); ++iterCrevasse)
	{
		if (iterCrevasse->is_border() || iterCrevasse->opposite()->is_border())
		{
			iterCrevasse->lpftOffset[0] = iterCrevasse->lpftOffset[1] = 0;
		}
		else
		{
			CParaConvertor::m_fnExtractOffset(iterCrevasse->seamtype(), iterCrevasse->first_element()->lpftGlobalPara, iterCrevasse->opposite()->prev()->first_element()->lpftGlobalPara, iterCrevasse->lpftOffset);
		}
	}
}
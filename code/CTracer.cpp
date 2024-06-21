#include "CTracer.h"

void CTracer::m_fnSift(int iStart)
{
	int iRoot, iLCh, iRCh, iSCh;
	CTentacle *lpTemp;
	iRoot = iStart;
	iLCh = iRoot * 2 + 1;
	iRCh = iLCh + 1;
	while (iLCh < m_nHeapSize)
	{
		if (iRCh < m_nHeapSize && m_lplpTenHeap[iRCh]->m_ftAccLen < m_lplpTenHeap[iLCh]->m_ftAccLen)
		{
			iSCh = iRCh;
		}
		else
		{
			iSCh = iLCh;
		}
		if (m_lplpTenHeap[iRoot]->m_ftAccLen < m_lplpTenHeap[iRoot]->m_ftAccLen)
		{
			break;
		}
		else
		{
			lpTemp = m_lplpTenHeap[iSCh];
			m_lplpTenHeap[iSCh] = m_lplpTenHeap[iRoot];
			m_lplpTenHeap[iRoot] = lpTemp;
			m_lplpTenHeap[iRoot]->iHeapPos = iRoot;
			m_lplpTenHeap[iSCh]->iHeapPos = iSCh;
			iRoot = iSCh;
			iLCh = iSCh * 2 + 1;
			iRCh = iLCh + 1;
		}
	}
}

void CTracer::m_fnCheck(int iNode)
{
	int iChild, iRoot;
	iChild = iNode;
	CTentacle *lpTemp;
	while (iChild != 0)
	{
		iRoot = (iChild - 1) / 2;
		if (m_lplpTenHeap[iChild]->m_ftAccLen < m_lplpTenHeap[iRoot]->m_ftAccLen)
		{
			lpTemp = m_lplpTenHeap[iChild];
			m_lplpTenHeap[iChild] = m_lplpTenHeap[iRoot];
			m_lplpTenHeap[iRoot] = lpTemp;
			m_lplpTenHeap[iChild]->iHeapPos = iChild;
			m_lplpTenHeap[iRoot]->iHeapPos = iRoot;
			iChild = iRoot;
		}
		else
		{
			break;
		}
	}
}

void CTracer::m_fnRemoveTentacle(int iHeapPos)
{
	CTentacle *lpTempTentacle;
	if (iHeapPos < m_nHeapSize && m_lplpTenHeap[iHeapPos]->bActive)
	{
		--m_nHeapSize;
		lpTempTentacle = m_lplpTenHeap[iHeapPos];
		lpTempTentacle->bActive = false;
		lpTempTentacle->m_lstJoints.clear();
		m_lplpTenHeap[iHeapPos] = m_lplpTenHeap[m_nHeapSize];
		m_lplpTenHeap[m_nHeapSize] = lpTempTentacle;
		m_lplpTenHeap[iHeapPos]->iHeapPos = iHeapPos;
		m_lplpTenHeap[m_nHeapSize]->iHeapPos = m_nHeapSize;
		m_fnSift(iHeapPos);
	}
}

void CTracer::m_fnCountPort(CMeshBase::CPolyhedron& mshSurface)
{
	CMeshBase::CPolyhedron::Vertex_iterator vi;
	CMeshBase::CPolyhedron::Halfedge_handle hh0, hh;
	CMeshBase::FT ftAngle, ftAngleDif;
	CMeshBase::CPolyhedron::Facet_iterator fi;
	int iDir;
	CPort newPort;
	newPort.bSatisfied = false;
	newPort.lpInducedPath = NULL;
	for (vi = mshSurface.vertices_begin(); vi != mshSurface.vertices_end(); ++vi)
	{
		if ((vi->iBordId == -1 && vi->nNumPorts != 0) || (vi->iBordId != -1 && vi->nNumPorts > 2))
		{
			if (vi->iBordId == -1)
			{
				hh0 = vi->halfedge();
				hh = hh0;
			}
			else
			{
				hh0 = vi->halfedge();
				while(!hh0->is_border())
				{
					hh0 = hh0->next()->opposite();
				}
				hh = hh0->next()->opposite();
			}
			hh0 = vi->halfedge();
			ftAngle = hh0->facet()->ftChartDir;
			switch (hh0->iLocalInd)
			{
			case 0:
				ftAngle -= hh0->ftAngle;
				break;
			case 1:
				ftAngle += hh0->next()->ftAngle;
				break;
			}
			ftAngleDif = 0;
			if (!hh0->next()->opposite()->is_border())
			{
				ftAngleDif = hh0->next()->opposite()->facet()->ftChartDir - hh0->facet()->ftChartDir + hh0->next()->ftPolarAxisDif - hh0->next()->iSeamType * ftHalfPi;
				ftAngleDif = ftAngleDif - CMeshBase::FT(int(ftAngleDif / ftTwoPi + 15.5) - 15) * ftTwoPi;
			}
			ftAngle += ftAngleDif / 2;
			iDir = (12 - int(ftAngle / ftHalfPi + 10)) % 4 - 2;
			hh = hh0;
			do
			{
				if (!hh->is_border())
				{
					ftAngle = hh->facet()->ftChartDir + iDir * ftHalfPi + ftPi;
					switch (hh->iLocalInd)
					{
					case 1:
						ftAngle -= hh->prev()->ftAngle;
						break;
					case 2:
						ftAngle += hh->ftAngle;
						break;
					}
					ftAngleDif = 0;
					if (!hh->opposite()->is_border())
					{
						ftAngleDif = hh->opposite()->facet()->ftChartDir - hh->facet()->ftChartDir + hh->ftPolarAxisDif - hh->iSeamType * ftHalfPi;
						ftAngleDif = ftAngleDif - CMeshBase::FT(int(ftAngleDif / ftTwoPi + 15.5) - 15) * ftTwoPi;
					}
					ftAngle += ftAngleDif / 2;
					ftAngle = ftAngle - CMeshBase::FT(int(ftAngle / ftTwoPi + 15.5) - 15) * ftTwoPi;
					if (ftAngle < 0)
					{
						newPort.hhStart = hh;
						newPort.iDir = iDir;
						newPort.ftAngle = ftPi + ftAngle - ftAngleDif / 2;
						if (newPort.ftAngle < 0)
						{
							newPort.ftAngle = 0;
						}
						if (newPort.ftAngle > ftPi)
						{
							newPort.ftAngle = ftPi;
						}
						m_lstPorts.push_back(newPort);
						iDir = (iDir + 3) % 4 - 2;
					}
				}
				iDir = (iDir - hh->iSeamType + 6) % 4 - 2;
				hh = hh->opposite()->prev();
			} while (hh != hh0);
		}
	}
	for (fi = mshSurface.facets_begin(); fi != mshSurface.facets_end(); ++fi)
	{
		for (iDir = 0; iDir < 4; ++iDir)
		{
			fi->lpiTraceInd[iDir] = -1;
		}
	}
	m_nNumFacets = mshSurface.size_of_facets();
}


void CTracer::m_fnAccurateSearch()
{
	std::list<CPort>::iterator iterPort;
	unsigned int i, j;
	CTenJoint newTenJoint;
	CMeshBase::FT ftAngle, ftProportion;
	CTentacle **lplplpOccupiedFacets[4];
	CTentacle *lpCurrentTentacle, *lpTestTentacle;
	CTenJoint *lpCurrentJoint;
	CMeshBase::CPolyhedron::Halfedge_handle hh0;
	CInducedPath newInducedPath;
	m_nNumTentacles = int(m_lstPorts.size());
	m_lpTentacles = new CTentacle[m_nNumTentacles];
	m_lplpTenHeap = new CTentacle*[m_nNumTentacles];
	for (i = 0; i < 4; ++i)
	{
		lplplpOccupiedFacets[i] = new CTentacle*[m_nNumFacets];
	}
	std::list<CTenJoint>::iterator iterJoint;
	std::list<CTenJoint>::reverse_iterator r_iterJoint;
	CPathElement newPathElement;
	CMeshBase::Point_3 pt1, pt2;
	for (i = 0; i < 4; ++i)
	{
		for (j = 0; j < m_nNumFacets; ++j)
		{
			lplplpOccupiedFacets[i][j] = NULL;
		}
	}
	i = 0;
	for (iterPort = m_lstPorts.begin(); iterPort != m_lstPorts.end(); ++iterPort)//(i = 0; i < m_nNumPorts; ++i)
	{
		m_lpTentacles[i].bActive = true;
		m_lpTentacles[i].m_lpSource = &(*iterPort);
		m_lpTentacles[i].m_ftAccLen = 0;
		m_lpTentacles[i].iHeapPos = i;
		m_lplpTenHeap[i] = &(m_lpTentacles[i]);
		newTenJoint.iDir = iterPort->iDir;
		newTenJoint.ftProportion = 1;
		newTenJoint.hhEdge = iterPort->hhStart;
		lplplpOccupiedFacets[(iterPort->iDir + 4) % 4][iterPort->hhStart->facet()->iIndex] = &(m_lpTentacles[i]);
		m_lpTentacles[i].m_lstJoints.push_back(newTenJoint);
		++i;
	}
	newInducedPath.m_bSelected = false;
	m_nHeapSize = m_nNumTentacles;
	lpCurrentTentacle = m_lplpTenHeap[0];
	while (lpCurrentTentacle->bActive)
	{
		lpCurrentJoint = &(*(lpCurrentTentacle->m_lstJoints.rbegin()));
		hh0 = lpCurrentJoint->hhEdge;
		ftProportion = lpCurrentJoint->ftProportion;
		pt1 = CGAL::barycenter(hh0->vertex()->point(), ftProportion, hh0->prev()->vertex()->point(), 1 - ftProportion);
		ftAngle = hh0->facet()->ftChartDir + lpCurrentJoint->iDir * ftHalfPi;
		if (hh0->iLocalInd == 1)
		{
			ftAngle -= hh0->prev()->ftAngle;
		}
		if (hh0->iLocalInd == 2)
		{
			ftAngle += hh0->ftAngle;
		}
		ftAngle -= CMeshBase::FT(int(ftAngle / ftTwoPi + 15.5) - 15) * ftTwoPi;
		if (ftAngle < 0)
		{
			if (ftAngle < -ftHalfPi)
			{
				ftAngle = ftPi;
			}
			else
			{
				ftAngle = 0;
			}
		}
		if (ftAngle + hh0->prev()->ftAngle < ftPi)
		{
			newTenJoint.hhEdge = hh0->next()->opposite();
			newTenJoint.ftProportion = 1 - hh0->ftLen * (1 - ftProportion) * sin(ftAngle) / (hh0->next()->ftLen * sin(hh0->ftAngle - ftAngle));
		}
		else if (ftAngle > hh0->ftAngle)
		{
			newTenJoint.hhEdge = hh0->prev()->opposite();
			newTenJoint.ftProportion = hh0->ftLen * ftProportion * sin(ftAngle) / (hh0->prev()->ftLen * sin(ftAngle + hh0->prev()->ftAngle - ftPi));
		}
		else
		{
			if (sin(hh0->ftAngle - ftAngle) * sin(hh0->prev()->ftAngle) * ftProportion > sin(ftAngle + hh0->prev()->ftAngle - ftPi) * sin(hh0->ftAngle) * (1 - ftProportion))
			{
				newTenJoint.hhEdge = hh0->next()->opposite();
				newTenJoint.ftProportion = 1 - hh0->ftLen * (1 - ftProportion) * sin(ftAngle) / (hh0->next()->ftLen * sin(hh0->ftAngle - ftAngle));
			}
			else
			{
				newTenJoint.hhEdge = hh0->prev()->opposite();
				newTenJoint.ftProportion = hh0->ftLen * ftProportion * sin(ftAngle) / (hh0->prev()->ftLen * sin(ftAngle + hh0->prev()->ftAngle - ftPi));
			}
		}
		pt2 = CGAL::barycenter<CMeshBase::FT>(newTenJoint.hhEdge->vertex()->point(), newTenJoint.ftProportion, newTenJoint.hhEdge->prev()->vertex()->point(), 1 - newTenJoint.ftProportion);
		newTenJoint.iDir = (lpCurrentJoint->iDir + newTenJoint.hhEdge->iSeamType + 6) % 4 - 2;
		lpCurrentTentacle->m_lstJoints.push_back(newTenJoint);
		lpCurrentTentacle->m_ftAccLen += sqrt(CGAL::squared_distance(pt1, pt2));
		m_fnSift(0);
		if (newTenJoint.hhEdge->is_border())
		{
			newInducedPath.lplpTerminals[0] = lpCurrentTentacle->m_lpSource;
			newInducedPath.lplpTerminals[1] = NULL;
			newInducedPath.m_bClosed = false;
			for (iterJoint = lpCurrentTentacle->m_lstJoints.begin(); iterJoint != lpCurrentTentacle->m_lstJoints.end(); ++iterJoint)
			{
				newPathElement.iDir = 0;
				newPathElement.hhEdge = iterJoint->hhEdge;
				newPathElement.ftProportion = iterJoint->ftProportion;
				newInducedPath.m_lstData.push_back(newPathElement);
			}
			m_lstIP.push_back(newInducedPath);
			newInducedPath.m_lstData.clear();
			lpCurrentTentacle->m_lpSource->lpInducedPath = &(*m_lstIP.rbegin());
			m_fnRemoveTentacle(lpCurrentTentacle->iHeapPos);
		}
		else
		{
			lpTestTentacle = lplplpOccupiedFacets[(newTenJoint.iDir + 2) % 4][newTenJoint.hhEdge->facet()->iIndex];
			if (lpTestTentacle != NULL)
			{
				if (lpTestTentacle->bActive)
				{
					lpCurrentTentacle->m_lpSource->bSatisfied = true;
					lpTestTentacle->m_lpSource->bSatisfied = true;
					newInducedPath.lplpTerminals[0] = lpCurrentTentacle->m_lpSource;
					newInducedPath.lplpTerminals[1] = lpTestTentacle->m_lpSource;
					newInducedPath.m_bClosed = true;
					for (iterJoint = lpCurrentTentacle->m_lstJoints.begin(); iterJoint != lpCurrentTentacle->m_lstJoints.end(); ++iterJoint)
					{
						newPathElement.iDir = 0;
						newPathElement.hhEdge = iterJoint->hhEdge;
						newPathElement.ftProportion = iterJoint->ftProportion;
						newInducedPath.m_lstData.push_back(newPathElement);
					}
					--iterJoint;
					for (r_iterJoint = lpTestTentacle->m_lstJoints.rbegin(); r_iterJoint != lpTestTentacle->m_lstJoints.rend(); ++r_iterJoint)
					{
						if (r_iterJoint->hhEdge->facet() == iterJoint->hhEdge->facet())
						{
							break;
						}
					}
					while (r_iterJoint != lpTestTentacle->m_lstJoints.rend())
					{
						newPathElement.iDir = 1;
						newPathElement.hhEdge = r_iterJoint->hhEdge;
						newPathElement.ftProportion = r_iterJoint->ftProportion;
						newInducedPath.m_lstData.push_back(newPathElement);
						++r_iterJoint;
					}
					lpCurrentTentacle->m_lstJoints.clear();
					lpTestTentacle->m_lstJoints.clear();
					m_lstIP.push_back(newInducedPath);
					newInducedPath.m_lstData.clear();
					lpCurrentTentacle->m_lpSource->lpInducedPath = &(*m_lstIP.rbegin());
					lpTestTentacle->m_lpSource->lpInducedPath = &(*m_lstIP.rbegin());
					m_fnRemoveTentacle(lpCurrentTentacle->iHeapPos);
					m_fnRemoveTentacle(lpTestTentacle->iHeapPos);
				}
				else
				{
					m_fnRemoveTentacle(lpCurrentTentacle->iHeapPos);
				}
			}
			else
			{
				lpTestTentacle = lplplpOccupiedFacets[(newTenJoint.iDir + 4) % 4][newTenJoint.hhEdge->facet()->iIndex];
				if (lpTestTentacle != NULL)
				{
					m_fnRemoveTentacle(lpCurrentTentacle->iHeapPos);
				}
				else
				{
					lplplpOccupiedFacets[(newTenJoint.iDir + 4) % 4][newTenJoint.hhEdge->facet()->iIndex] = lpCurrentTentacle;
				}
			}
		}
		lpCurrentTentacle = m_lplpTenHeap[0];
	}
	for (i = 0; i < 4; ++i)
	{
		delete[]lplplpOccupiedFacets[3 - i];
	}
	m_nNumTentacles = 0;
	delete[]m_lplpTenHeap;
	delete[]m_lpTentacles;
}


void CTracer::m_fnFillCoord(CTenJoint *lpJoint, bool bOpposite)
{
    CMeshBase::CPolyhedron::Halfedge_handle hh;
	CMeshBase::FT lpftCoord[3];
	hh = lpJoint->hhEdge;
	if (!hh->is_border())
	{
		lpftCoord[hh->iLocalInd] = lpJoint->ftProportion;
		lpftCoord[(hh->iLocalInd + 1) % 3] = 0.0;
		lpftCoord[(hh->iLocalInd + 2) % 3] = 1.0 - lpJoint->ftProportion;
		hh->facet()->lpvtTrace[(lpJoint->iDir + 4) % 4] = CMeshBase::Vector_3(lpftCoord[0], lpftCoord[1], lpftCoord[2]);
		hh->facet()->lpiTraceInd[(lpJoint->iDir + 4) % 4] = hh->iLocalInd;
	}
    if (bOpposite && !hh->opposite()->is_border())
    {
        hh = hh->opposite();
        lpftCoord[hh->iLocalInd] = 1.0 - lpJoint->ftProportion;
        lpftCoord[(hh->iLocalInd + 1) % 3] = 0.0;
        lpftCoord[(hh->iLocalInd + 2) % 3] = lpJoint->ftProportion;
        hh->facet()->lpvtTrace[(lpJoint->iDir + hh->iSeamType + 6) % 4] = CMeshBase::Vector_3(lpftCoord[0], lpftCoord[1], lpftCoord[2]);
        hh->facet()->lpiTraceInd[(lpJoint->iDir + hh->iSeamType + 6) % 4] = hh->iLocalInd;
    }
}

void CTracer::m_fnTShapeSearch()
{
	std::list<CPort>::iterator iterPort;
	unsigned int i;
	CTenJoint newTenJoint, *lpCurrentJoint;
	CTentacle *lpCurrentTentacle;
	CMeshBase::FT ftAngle, ftProportion, ftSum;
	CMeshBase::CPolyhedron::Halfedge_handle hh0;
	m_nNumTentacles = int(m_lstPorts.size());
	m_lpTentacles = new CTentacle[m_nNumTentacles];
	m_lplpTenHeap = new CTentacle*[m_nNumTentacles];
	CMeshBase::Vector_3 vtIntersect, *lpvtCoef;
	CMeshBase::Point_3 pt1, pt2;
	i = 0;
	for (iterPort = m_lstPorts.begin(); iterPort != m_lstPorts.end(); ++iterPort)//(i = 0; i < m_nNumPorts; ++i)
	{
		m_lpTentacles[i].bActive = true;
		m_lpTentacles[i].m_lpSource = &(*iterPort);
		m_lpTentacles[i].m_ftAccLen = 0;
		m_lpTentacles[i].iHeapPos = i;
		m_lplpTenHeap[i] = &(m_lpTentacles[i]);
		newTenJoint.iDir = iterPort->iDir;
		newTenJoint.ftProportion = 1;
		newTenJoint.hhEdge = iterPort->hhStart;
		m_fnFillCoord(&newTenJoint, false);
		m_lpTentacles[i].m_lstJoints.push_back(newTenJoint);
		++i;
	}
	m_nHeapSize = m_nNumTentacles;
	lpCurrentTentacle = m_lplpTenHeap[0];
	while (m_nHeapSize > 0)
	{
		lpCurrentJoint = &(*(lpCurrentTentacle->m_lstJoints.rbegin()));
		if (lpCurrentJoint->hhEdge->facet()->lpiTraceInd[(lpCurrentJoint->iDir + 6) % 4] == -1)
        {
            hh0 = lpCurrentJoint->hhEdge;
            ftProportion = lpCurrentJoint->ftProportion;
            pt1 = CGAL::barycenter(hh0->vertex()->point(), ftProportion, hh0->prev()->vertex()->point(), 1 - ftProportion);
            ftAngle = hh0->facet()->ftChartDir + lpCurrentJoint->iDir * ftHalfPi;
            if (hh0->iLocalInd == 1)
            {
                ftAngle -= hh0->prev()->ftAngle;
            }
            if (hh0->iLocalInd == 2)
            {
                ftAngle += hh0->ftAngle;
            }
            ftAngle -= CMeshBase::FT(int(ftAngle / ftTwoPi + 15.5) - 15) * ftTwoPi;
            if (ftAngle < 0)
            {
                if (ftAngle < -ftHalfPi)
                {
                    ftAngle = ftPi;
                }
                else
                {
                    ftAngle = 0;
                }
            }
            if (ftAngle + hh0->prev()->ftAngle < ftPi)
            {
                newTenJoint.hhEdge = hh0->next()->opposite();
                newTenJoint.ftProportion = 1 - hh0->ftLen * (1 - ftProportion) * sin(ftAngle) / (hh0->next()->ftLen * sin(hh0->ftAngle - ftAngle));
            }
            else if (ftAngle > hh0->ftAngle)
            {
                newTenJoint.hhEdge = hh0->prev()->opposite();
                newTenJoint.ftProportion = hh0->ftLen * ftProportion * sin(ftAngle) / (hh0->prev()->ftLen * sin(ftAngle + hh0->prev()->ftAngle - ftPi));
            }
            else
            {
                if (sin(hh0->ftAngle - ftAngle) * sin(hh0->prev()->ftAngle) * ftProportion > sin(ftAngle + hh0->prev()->ftAngle - ftPi) * sin(hh0->ftAngle) * (1 - ftProportion))
                {
                    newTenJoint.hhEdge = hh0->next()->opposite();
                    newTenJoint.ftProportion = 1 - hh0->ftLen * (1 - ftProportion) * sin(ftAngle) / (hh0->next()->ftLen * sin(hh0->ftAngle - ftAngle));
                }
                else
                {
                    newTenJoint.hhEdge = hh0->prev()->opposite();
                    newTenJoint.ftProportion = hh0->ftLen * ftProportion * sin(ftAngle) / (hh0->prev()->ftLen * sin(ftAngle + hh0->prev()->ftAngle - ftPi));
                }
            }
			if (newTenJoint.ftProportion < 1.0 / 65536.0)
			{
				newTenJoint.ftProportion = 0.0;
			}
			if (newTenJoint.ftProportion > 65535.0 / 65536.0)
			{
				newTenJoint.ftProportion = 1.0;
			}
            pt2 = CGAL::barycenter<CMeshBase::FT>(newTenJoint.hhEdge->vertex()->point(), newTenJoint.ftProportion, newTenJoint.hhEdge->prev()->vertex()->point(), 1 - newTenJoint.ftProportion);
            newTenJoint.iDir = (lpCurrentJoint->iDir + newTenJoint.hhEdge->iSeamType + 6) % 4 - 2;
			if (newTenJoint.hhEdge->is_border())
			{
				m_fnFillCoord(&newTenJoint, true);
				m_fnRemoveTentacle(lpCurrentTentacle->iHeapPos);
			}
			else if (newTenJoint.hhEdge->facet()->lpiTraceInd[(newTenJoint.iDir + 6) % 4] != -1)
			{
				m_fnRemoveTentacle(lpCurrentTentacle->iHeapPos);
			}
			else if (newTenJoint.hhEdge->facet()->lpiTraceInd[(newTenJoint.iDir + 4) % 4] == -1)
			{
				m_fnFillCoord(&newTenJoint, true);
				if (lpCurrentJoint->hhEdge->facet()->lpiTraceInd[(lpCurrentJoint->iDir + 5) % 4] != -1 &&
					lpCurrentJoint->hhEdge->facet()->lpiTraceInd[(lpCurrentJoint->iDir + 6) % 4] != -1 &&
					lpCurrentJoint->hhEdge->facet()->lpiTraceInd[(lpCurrentJoint->iDir + 7) % 4] != -1)
				{
					lpvtCoef = lpCurrentJoint->hhEdge->facet()->lpvtTrace;
					vtIntersect =  CGAL::cross_product(CGAL::cross_product(lpvtCoef[0], lpvtCoef[2]), CGAL::cross_product(lpvtCoef[1], lpvtCoef[3]));
					ftSum = vtIntersect.x() + vtIntersect.y() + vtIntersect.z();
					vtIntersect = vtIntersect / ftSum;
					if (vtIntersect.x() > 0 && vtIntersect.y() > 0 && vtIntersect.z() > 0)
					{
						lpCurrentJoint->hhEdge->facet()->lpvtTrace[(lpCurrentJoint->iDir + 6) % 4] = vtIntersect;
						m_fnRemoveTentacle(lpCurrentTentacle->iHeapPos);
					}
					else
					{
						lpCurrentTentacle->m_lstJoints.clear();
						lpCurrentTentacle->m_lstJoints.push_back(newTenJoint);
						lpCurrentTentacle->m_ftAccLen += sqrt(CGAL::squared_distance(pt1, pt2));
						m_fnSift(0);
					}
				}
				else
				{
					lpCurrentTentacle->m_lstJoints.clear();
					lpCurrentTentacle->m_lstJoints.push_back(newTenJoint);
					lpCurrentTentacle->m_ftAccLen += sqrt(CGAL::squared_distance(pt1, pt2));
					m_fnSift(0);
				}
			}
			else
			{
				m_fnRemoveTentacle(lpCurrentTentacle->iHeapPos);
			}
        }
        else
        {
            m_fnRemoveTentacle(lpCurrentTentacle->iHeapPos);
        }
		lpCurrentTentacle = m_lplpTenHeap[0];
	}
	m_nNumTentacles = 0;
	delete[]m_lplpTenHeap;
	delete[]m_lpTentacles;
}


void CTracer::m_fnBiasSearch(CMeshBase::FT ftBiasAngle)
{
	std::list<CPort>::iterator iterPort;
	unsigned int i, j;
	CTenJoint newTenJoint;
	CMeshBase::FT ftAngle, ftProportion;
	CTentacle **lplplpOccupiedFacets[8];
	CTentacle *lpCurrentTentacle, *lpTestTentacle;
	CTenJoint *lpCurrentJoint;
	CMeshBase::CPolyhedron::Halfedge_handle hh0;
	CInducedPath newInducedPath;
	std::list<CTenJoint>::iterator iterJoint;
	std::list<CTenJoint>::reverse_iterator r_iterJoint;
	CPathElement newPathElement;
	CMeshBase::Point_3 pt1, pt2;
	bool bVacant;
	std::list<CInducedPath>::iterator iterIP;
	m_nNumTentacles = 0;
	for (iterPort = m_lstPorts.begin(); iterPort != m_lstPorts.end(); ++iterPort)
	{
		if (!iterPort->bSatisfied)
		{
			m_nNumTentacles += 2;
		}
	}
	m_lpTentacles = new CTentacle[m_nNumTentacles];
	m_lplpTenHeap = new CTentacle*[m_nNumTentacles];
	for (i = 0; i < 8; ++i)
	{
		lplplpOccupiedFacets[i] = new CTentacle*[m_nNumFacets];
	}
	for (i = 0; i < 8; ++i)
	{
		for (j = 0; j < m_nNumFacets; ++j)
		{
			lplplpOccupiedFacets[i][j] = NULL;
		}
	}
	i = 0;
	for (iterPort = m_lstPorts.begin(); iterPort != m_lstPorts.end(); ++iterPort)
	{
		if (!iterPort->bSatisfied)
		{
			m_lpTentacles[i].bActive = true;
			m_lpTentacles[i].m_lpSource = &(*iterPort);
			m_lpTentacles[i].m_ftAccLen = 0;
			m_lpTentacles[i].iHeapPos = i;
			m_lpTentacles[i].iBias = -1;
			m_lplpTenHeap[i] = &(m_lpTentacles[i]);
			newTenJoint.iDir = iterPort->iDir;
			newTenJoint.ftProportion = 0;
			newTenJoint.hhEdge = iterPort->hhStart;
			ftAngle = newTenJoint.hhEdge->facet()->ftChartDir + newTenJoint.iDir * ftHalfPi - ftBiasAngle;
			if (newTenJoint.hhEdge->iLocalInd == 1)
			{
				ftAngle -= newTenJoint.hhEdge->prev()->ftAngle;
			}
			if (newTenJoint.hhEdge->iLocalInd == 2)
			{
				ftAngle += newTenJoint.hhEdge->ftAngle;
			}
			ftAngle -= (int(ftAngle / ftTwoPi + 15.5) - 15) * ftTwoPi;
			while (ftAngle < 0)
			{
				newTenJoint.iDir = (newTenJoint.iDir + 6 - newTenJoint.hhEdge->iSeamType) % 4 - 2;
				newTenJoint.hhEdge = newTenJoint.hhEdge->opposite()->next();
				ftAngle = newTenJoint.hhEdge->facet()->ftChartDir + newTenJoint.iDir * ftHalfPi - ftBiasAngle;
				if (newTenJoint.hhEdge->iLocalInd == 1)
				{
					ftAngle -= newTenJoint.hhEdge->prev()->ftAngle;
				}
				if (newTenJoint.hhEdge->iLocalInd == 2)
				{
					ftAngle += newTenJoint.hhEdge->ftAngle;
				}
				ftAngle -= (int(ftAngle / ftTwoPi + 15.5) - 15) * ftTwoPi;
			}
			m_lpTentacles[i].m_lstJoints.push_back(newTenJoint);
			lplplpOccupiedFacets[(newTenJoint.iDir + 4) % 4 + (m_lpTentacles[i].iBias + 1) * 2][newTenJoint.hhEdge->facet()->iIndex] = &(m_lpTentacles[i]);
			++i;
			m_lpTentacles[i].bActive = true;
			m_lpTentacles[i].m_lpSource = &(*iterPort);
			m_lpTentacles[i].m_ftAccLen = 0;
			m_lpTentacles[i].iHeapPos = i;
			m_lpTentacles[i].iBias = 1;
			m_lplpTenHeap[i] = &(m_lpTentacles[i]);
			newTenJoint.iDir = iterPort->iDir;
			newTenJoint.ftProportion = 0;
			newTenJoint.hhEdge = iterPort->hhStart;
			ftAngle = newTenJoint.hhEdge->facet()->ftChartDir + newTenJoint.iDir * ftHalfPi + ftBiasAngle;
			if (newTenJoint.hhEdge->iLocalInd == 1)
			{
				ftAngle -= newTenJoint.hhEdge->prev()->ftAngle;
			}
			if (newTenJoint.hhEdge->iLocalInd == 2)
			{
				ftAngle += newTenJoint.hhEdge->ftAngle;
			}
			ftAngle -= (int(ftAngle / ftTwoPi + 15.5) - 15) * ftTwoPi;
			while (ftAngle > ftPi - newTenJoint.hhEdge->prev()->ftAngle)
			{
				newTenJoint.iDir = (newTenJoint.iDir + 6 - newTenJoint.hhEdge->prev()->iSeamType) % 4 - 2;
				newTenJoint.hhEdge = newTenJoint.hhEdge->prev()->opposite();
				ftAngle = newTenJoint.hhEdge->facet()->ftChartDir + newTenJoint.iDir * ftHalfPi + ftBiasAngle;
				if (newTenJoint.hhEdge->iLocalInd == 1)
				{
					ftAngle -= newTenJoint.hhEdge->prev()->ftAngle;
				}
				if (newTenJoint.hhEdge->iLocalInd == 2)
				{
					ftAngle += newTenJoint.hhEdge->ftAngle;
				}
				ftAngle -= (int(ftAngle / ftTwoPi + 15.5) - 15) * ftTwoPi;
			}
			m_lpTentacles[i].m_lstJoints.push_back(newTenJoint);
			lplplpOccupiedFacets[(newTenJoint.iDir + 4) % 4 + (m_lpTentacles[i].iBias + 1) * 2][newTenJoint.hhEdge->facet()->iIndex] = &(m_lpTentacles[i]);
			++i;
		}
	}
	newInducedPath.m_bSelected = false;
	m_nHeapSize = m_nNumTentacles;
	lpCurrentTentacle = m_lplpTenHeap[0];
	while (m_nHeapSize > 0 && lpCurrentTentacle->bActive)
	{
		lpCurrentJoint = &(*(lpCurrentTentacle->m_lstJoints.rbegin()));
		hh0 = lpCurrentJoint->hhEdge;
		ftProportion = lpCurrentJoint->ftProportion;
		pt1 = CGAL::barycenter(hh0->vertex()->point(), ftProportion, hh0->prev()->vertex()->point(), 1 - ftProportion);
		ftAngle = hh0->facet()->ftChartDir + lpCurrentJoint->iDir * ftHalfPi + (lpCurrentTentacle->iBias) * ftBiasAngle;
		if (hh0->iLocalInd == 1)
		{
			ftAngle -= hh0->prev()->ftAngle;
		}
		if (hh0->iLocalInd == 2)
		{
			ftAngle += hh0->ftAngle;
		}
		ftAngle -= CMeshBase::FT(int(ftAngle / ftTwoPi + 15.5) - 15) * ftTwoPi;
		if (ftAngle < 0)
		{
			if (ftAngle < -ftHalfPi)
			{
				ftAngle = ftPi;
			}
			else
			{
				ftAngle = 0;
			}
		}
		if (ftAngle + hh0->prev()->ftAngle < ftPi)
		{
			newTenJoint.hhEdge = hh0->next()->opposite();
			newTenJoint.ftProportion = 1 - hh0->ftLen * (1 - ftProportion) * sin(ftAngle) / (hh0->next()->ftLen * sin(hh0->ftAngle - ftAngle));
		}
		else if (ftAngle > hh0->ftAngle)
		{
			newTenJoint.hhEdge = hh0->prev()->opposite();
			newTenJoint.ftProportion = hh0->ftLen * ftProportion * sin(ftAngle) / (hh0->prev()->ftLen * sin(ftAngle + hh0->prev()->ftAngle - ftPi));
		}
		else
		{
			if (sin(hh0->ftAngle - ftAngle) * sin(hh0->prev()->ftAngle) * ftProportion > sin(ftAngle + hh0->prev()->ftAngle - ftPi) * sin(hh0->ftAngle) * (1 - ftProportion))
			{
				newTenJoint.hhEdge = hh0->next()->opposite();
				newTenJoint.ftProportion = 1 - hh0->ftLen * (1 - ftProportion) * sin(ftAngle) / (hh0->next()->ftLen * sin(hh0->ftAngle - ftAngle));
			}
			else
			{
				newTenJoint.hhEdge = hh0->prev()->opposite();
				newTenJoint.ftProportion = hh0->ftLen * ftProportion * sin(ftAngle) / (hh0->prev()->ftLen * sin(ftAngle + hh0->prev()->ftAngle - ftPi));
			}
		}
		pt2 = CGAL::barycenter<CMeshBase::FT>(newTenJoint.hhEdge->vertex()->point(), newTenJoint.ftProportion, newTenJoint.hhEdge->prev()->vertex()->point(), 1 - newTenJoint.ftProportion);
		newTenJoint.iDir = (lpCurrentJoint->iDir + newTenJoint.hhEdge->iSeamType + 6) % 4 - 2;
		lpCurrentTentacle->m_lstJoints.push_back(newTenJoint);
		lpCurrentTentacle->m_ftAccLen += sqrt(CGAL::squared_distance(pt1, pt2));
		m_fnSift(0);
		if (newTenJoint.hhEdge->is_border())
		{
			m_fnRemoveTentacle(lpCurrentTentacle->iHeapPos);
			m_fnRemoveTentacle(lpCurrentTentacle[-lpCurrentTentacle->iBias].iHeapPos);
		}
		else
		{
			bVacant = true;
			for (i = 0; i < 2; ++i)
			{
				lpTestTentacle = lplplpOccupiedFacets[(newTenJoint.iDir + 6) % 4 + i * 4][newTenJoint.hhEdge->facet()->iIndex];
				if (lpTestTentacle != NULL)
				{
					bVacant = false;
					if (lpTestTentacle->bActive && !lpTestTentacle->m_lpSource->bSatisfied && lpCurrentTentacle->m_lpSource != lpTestTentacle->m_lpSource)
					{
						lpCurrentTentacle->m_lpSource->bSatisfied = true;
						lpTestTentacle->m_lpSource->bSatisfied = true;
						if (lpCurrentTentacle->m_lpSource->lpInducedPath != NULL)
						{
							lpCurrentTentacle->m_lpSource->lpInducedPath->lplpTerminals[0] = NULL;
							lpCurrentTentacle->m_lpSource->lpInducedPath->lplpTerminals[1] = NULL;
						}
						if (lpTestTentacle->m_lpSource->lpInducedPath != NULL)
						{
							lpTestTentacle->m_lpSource->lpInducedPath->lplpTerminals[0] = NULL;
							lpTestTentacle->m_lpSource->lpInducedPath->lplpTerminals[1] = NULL;
						}
						newInducedPath.lplpTerminals[0] = lpCurrentTentacle->m_lpSource;
						newInducedPath.lplpTerminals[1] = lpTestTentacle->m_lpSource;
						newInducedPath.m_bClosed = true;
						for (iterJoint = lpCurrentTentacle->m_lstJoints.begin(); iterJoint != lpCurrentTentacle->m_lstJoints.end(); ++iterJoint)
						{
							newPathElement.iDir = 0;
							newPathElement.hhEdge = iterJoint->hhEdge;
							newPathElement.ftProportion = iterJoint->ftProportion;
							newInducedPath.m_lstData.push_back(newPathElement);
						}
						--iterJoint;
						for (r_iterJoint = lpTestTentacle->m_lstJoints.rbegin(); r_iterJoint != lpTestTentacle->m_lstJoints.rend(); ++r_iterJoint)
						{
							if (r_iterJoint->hhEdge->facet() == iterJoint->hhEdge->facet())
							{
								break;
							}
						}
						while (r_iterJoint != lpTestTentacle->m_lstJoints.rend())
						{
							newPathElement.iDir = 1;
							newPathElement.hhEdge = r_iterJoint->hhEdge;
							newPathElement.ftProportion = r_iterJoint->ftProportion;
							newInducedPath.m_lstData.push_back(newPathElement);
							++r_iterJoint;
						}
						m_lstIP.push_back(newInducedPath);
						newInducedPath.m_lstData.clear();
						lpCurrentTentacle->m_lpSource->lpInducedPath = &(*m_lstIP.rbegin());
						lpTestTentacle->m_lpSource->lpInducedPath = &(*m_lstIP.rbegin());
						m_fnRemoveTentacle(lpCurrentTentacle->iHeapPos);
						m_fnRemoveTentacle(lpCurrentTentacle[-lpCurrentTentacle->iBias].iHeapPos);
						m_fnRemoveTentacle(lpTestTentacle->iHeapPos);
						m_fnRemoveTentacle(lpTestTentacle[-lpTestTentacle->iBias].iHeapPos);
					}
					else if (lpCurrentTentacle->m_lpSource != lpTestTentacle->m_lpSource)
					{
						m_fnRemoveTentacle(lpCurrentTentacle->iHeapPos);
						m_fnRemoveTentacle(lpCurrentTentacle[-lpCurrentTentacle->iBias].iHeapPos);
					}
				}
			}
			if (bVacant)
			{
				lpTestTentacle = lplplpOccupiedFacets[(newTenJoint.iDir + 4) % 4 + (lpCurrentTentacle->iBias + 1) * 2][newTenJoint.hhEdge->facet()->iIndex];
				if (lpTestTentacle != NULL)
				{
					m_fnRemoveTentacle(lpCurrentTentacle->iHeapPos);
					m_fnRemoveTentacle(lpCurrentTentacle[lpCurrentTentacle->iBias].iHeapPos);
				}
				else
				{
					lplplpOccupiedFacets[(newTenJoint.iDir + 4) % 4 + (lpCurrentTentacle->iBias + 1) * 2][newTenJoint.hhEdge->facet()->iIndex] = lpCurrentTentacle;
				}
			}
		}
		lpCurrentTentacle = m_lplpTenHeap[0];
	}
	for (i = 0; i < 8; ++i)
	{
		delete[]lplplpOccupiedFacets[7 - i];
	}
	m_nNumTentacles = 0;
	delete[]m_lplpTenHeap;
	delete[]m_lpTentacles;
	for (iterIP = m_lstIP.begin(); iterIP != m_lstIP.end(); ++iterIP)
	{
		if (iterIP->lplpTerminals[0] == NULL && iterIP->lplpTerminals[1] == NULL)
		{
			iterIP->m_lstData.clear();
			iterIP = m_lstIP.erase(iterIP);
			--iterIP;
		}
	}
}

void CTracer::m_fnDeletePath()
{
	std::list<CInducedPath>::iterator iterIP;
	int i;
	if (!m_lstIP.empty())
	{
		for (iterIP = m_lstIP.begin(); iterIP != m_lstIP.end(); ++iterIP)
		{
			if (iterIP->m_bSelected)
			{
				for (i = 0; i < 2; ++i)
				{
					if (iterIP->lplpTerminals[i] != NULL)
					{
						iterIP->lplpTerminals[i]->bSatisfied = false;
						iterIP->lplpTerminals[i]->lpInducedPath = NULL;
						iterIP->lplpTerminals[i] = NULL;
					}
				}
				iterIP->m_lstData.clear();
				iterIP = m_lstIP.erase(iterIP);
				--iterIP;
			}
		}
	}
}


void CTracer::m_fnAdjustPath()
{
	std::list<CInducedPath>::iterator iterIP;
	std::list<CPathElement>::iterator iterCurrent, iterPE1, iterPE2;
	CMeshBase::CPolyhedron::Halfedge_handle hh1, hh2;
	CMeshBase::Point_3 ptIns1, ptIns2, ptImage;
	CMeshBase::Vector_3 vtAxis0, vtAxis1, vtDist;
	CMeshBase::FT ftPropertion1, ftPropertion2, ftLastAngle, ftL0, ftL1, ftH0, ftH1;
	if (!m_lstIP.empty())
	{
		for (iterIP = m_lstIP.begin(); iterIP != m_lstIP.end(); ++iterIP)
		{
			if (iterIP->m_bClosed && iterIP->m_lstData.size() > 1)
			{
				iterPE2 = iterIP->m_lstData.begin();
				iterPE1 = iterPE2;
				++iterPE2;
				iterPE1->ftDistStart = 0;
				iterPE1->ftAngleStart = 0;
				while (iterPE2 != iterIP->m_lstData.end())
				{
					if (iterPE1->iDir == 0)
					{
						hh1 = iterPE1->hhEdge;
						ftPropertion1 = iterPE1->ftProportion;
					}
					else
					{
						hh1 = iterPE1->hhEdge->opposite();
						ftPropertion1 = 1 - iterPE1->ftProportion;
					}
					if (iterPE2->iDir == 0)
					{
						hh2 = iterPE2->hhEdge->opposite();
						ftPropertion2 = 1 - iterPE2->ftProportion;
					}
					else
					{
						hh2 = iterPE2->hhEdge;
						ftPropertion2 = iterPE2->ftProportion;
					}
					ptIns1 = CGAL::barycenter(hh1->prev()->vertex()->point(), 1 - ftPropertion1, hh1->vertex()->point(), ftPropertion1);
					ptIns2 = CGAL::barycenter(hh2->prev()->vertex()->point(), 1 - ftPropertion2, hh2->vertex()->point(), ftPropertion2);
					vtAxis0 = -hh1->vtVector / hh1->ftLen;
					vtAxis1 = CGAL::cross_product(hh1->facet()->vtNorm, vtAxis0) / hh1->facet()->ftTwArea;
					ftLastAngle = iterPE1->ftAngleStart - hh1->ftMom / 2;
					if (ftLastAngle > ftPi)
					{
						ftLastAngle = ftPi;
					}
					if (ftLastAngle < 0)
					{
						ftLastAngle = 0;
					}
					ptImage = ptIns1 + iterPE1->ftDistStart * (vtAxis0 * cos(ftLastAngle) + vtAxis1 * sin(ftLastAngle));
					vtDist = ptImage - ptIns2;
					iterPE2->ftAngleStart = atan2(CGAL::cross_product(hh2->vtVector, vtDist) * hh2->facet()->vtNorm / hh2->facet()->ftTwArea, hh2->vtVector * vtDist) + hh2->ftMom / 2;
					if (iterPE2->ftAngleStart > ftPi)
					{
						iterPE2->ftAngleStart = ftPi;
					}
					if (iterPE2->ftAngleStart < 0)
					{
						if (iterPE2->ftAngleStart < -ftHalfPi)
						{
							iterPE2->ftAngleStart = ftPi;
						}
						else
						{
							iterPE2->ftAngleStart = 0;
						}
					}
					iterPE2->ftDistStart = sqrt(vtDist * vtDist);
					iterPE1 = iterPE2;
					++iterPE2;
				};
				iterPE2 = iterIP->m_lstData.end();
				--iterPE2;
				iterPE1 = iterPE2;
				--iterPE1;
				iterPE2->ftDistEnd = 0;
				iterPE2->ftAngleEnd = 0;
				while (iterPE1 != iterIP->m_lstData.begin())
				{
					if (iterPE1->iDir == 0)
					{
						hh1 = iterPE1->hhEdge;
						ftPropertion1 = iterPE1->ftProportion;
					}
					else
					{
						hh1 = iterPE1->hhEdge->opposite();
						ftPropertion1 = 1 - iterPE1->ftProportion;
					}
					if (iterPE2->iDir == 0)
					{
						hh2 = iterPE2->hhEdge->opposite();
						ftPropertion2 = 1 - iterPE2->ftProportion;
					}
					else
					{
						hh2 = iterPE2->hhEdge;
						ftPropertion2 = iterPE2->ftProportion;
					}
					ptIns1 = CGAL::barycenter(hh1->prev()->vertex()->point(), 1 - ftPropertion1, hh1->vertex()->point(), ftPropertion1);
					ptIns2 = CGAL::barycenter(hh2->prev()->vertex()->point(), 1 - ftPropertion2, hh2->vertex()->point(), ftPropertion2);
					vtAxis0 = -hh2->vtVector / hh2->ftLen;
					vtAxis1 = CGAL::cross_product(hh2->facet()->vtNorm, vtAxis0) / hh2->facet()->ftTwArea;
					ftLastAngle = iterPE2->ftAngleEnd - hh2->ftMom / 2;
					ptImage = ptIns2 + iterPE2->ftDistEnd * (vtAxis0 * cos(ftLastAngle) + vtAxis1 * sin(ftLastAngle));
					vtDist = ptImage - ptIns1;
					iterPE1->ftAngleEnd = atan2(CGAL::cross_product(hh1->vtVector, vtDist) * hh1->facet()->vtNorm / hh1->facet()->ftTwArea, hh1->vtVector * vtDist) + hh1->ftMom / 2;
					if (iterPE1->ftAngleEnd > ftPi)
					{
						iterPE1->ftAngleEnd = ftPi;
					}
					if (iterPE1->ftAngleEnd < 0)
					{
						if (iterPE1->ftAngleEnd < -ftHalfPi)
						{
							iterPE1->ftAngleEnd = ftPi;
						}
						else
						{
							iterPE1->ftAngleEnd = 0;
						}
					}
					iterPE1->ftDistEnd = sqrt(vtDist * vtDist);
					iterPE2 = iterPE1;
					--iterPE1;
				}
				iterPE2 = iterIP->m_lstData.end();
				--iterPE2;
				iterPE1 = iterIP->m_lstData.begin();
				++iterPE1;
				while (iterPE1 != iterPE2)
				{
					ftL0 = iterPE1->ftDistStart * cos(iterPE1->ftAngleStart);
					ftH0 = iterPE1->ftDistStart * sin(iterPE1->ftAngleStart);
					ftL1 = iterPE1->ftDistEnd * cos(iterPE1->ftAngleEnd);
					ftH1 = iterPE1->ftDistEnd * sin(iterPE1->ftAngleEnd);
					if (ftH0 + ftH1 != 0)
					{
						if (iterPE1->iDir == 0)
						{
							iterPE1->ftProportion -= (ftH1 * ftL0 - ftH0 * ftL1) / (ftH0 + ftH1) / iterPE1->hhEdge->ftLen;
						}
						else
						{
							iterPE1->ftProportion += (ftH1 * ftL0 - ftH0 * ftL1) / (ftH0 + ftH1) / iterPE1->hhEdge->ftLen;
						}
					}
					if (iterPE1->ftProportion < 0)
					{
						iterPE1->ftProportion = 0;
					}
					if (iterPE1->ftProportion > 1)
					{
						iterPE1->ftProportion = 1;
					}
					++iterPE1;
				}
			}
		}
	}
}

void CTracer::m_fnRenewPath()
{
	std::list<CInducedPath>::iterator iterIP;
	std::list<CPathElement>::iterator iterPE, iterPE0, iterPE1, iterPELast;
	CMeshBase::CPolyhedron::Halfedge_handle hhStart, hh;
	CMeshBase::CPolyhedron::Vertex_handle vhPole;
	CPathElement newElement;
	if (!m_lstIP.empty())
	{
		for (iterIP = m_lstIP.begin(); iterIP != m_lstIP.end(); ++iterIP)
		{
			if (iterIP->m_bClosed && iterIP->m_lstData.size() > 2)
			{
				iterPE = iterIP->m_lstData.begin();
				++iterPE;
				iterPELast = iterIP->m_lstData.end();
				--iterPELast;
				while (iterPE != iterPELast && iterPE != iterIP->m_lstData.end())
				{
					if ((iterPE->iDir == 0 && iterPE->ftProportion == 0) || (iterPE->iDir == 1 && iterPE->ftProportion == 1))
					{
						if (iterPE->iDir == 0)
						{
							vhPole = iterPE->hhEdge->prev()->vertex();
						}
						else
						{
							vhPole = iterPE->hhEdge->vertex();
						}
						iterPE0 = iterPE;
						--iterPE0;
						while (((iterPE0->iDir == 0 && iterPE0->hhEdge->prev()->vertex() == vhPole) || (iterPE0->iDir == 1 && iterPE0->hhEdge->vertex() == vhPole)))
						{
							if (iterPE0 == iterIP->m_lstData.begin())
							{
								--iterPE0;
								break;
							}
							else
							{
								--iterPE0;
							}
						}
						++iterPE0;
						if (iterPE0->iDir == 0)
						{
							hhStart = iterPE0->hhEdge->opposite();
						}
						else
						{
							hhStart = iterPE0->hhEdge;
						}
						if (iterPE0 == iterIP->m_lstData.begin())
						{
							iterPE0->iDir = 1 - iterPE0->iDir;
							++iterPE0;
						}
						iterPE1 = iterPE;
						++iterPE1;
						while (((iterPE1->iDir == 0 && iterPE1->hhEdge->prev()->vertex() == vhPole) || (iterPE1->iDir == 1 && iterPE1->hhEdge->vertex() == vhPole)))
						{
							if (iterPE1 == iterPELast)
							{
								++iterPE1;
								break;
							}
							else
							{
								++iterPE1;
							}
						}
						iterPE = iterPE1;
						--iterPE;
						if (iterPE->iDir == 0)
						{
							hh = iterPE->hhEdge->prev();
						}
						else
						{
							hh = iterPE->hhEdge->opposite()->prev();
						}
						if (iterPE1 == iterIP->m_lstData.end())
						{
							iterPE->iDir = 1 - iterPE->iDir;
							--iterPE1;
						}
						iterPE = iterPE1;
						iterIP->m_lstData.erase(iterPE0, iterPE1);
						newElement.iDir = 0;
						newElement.ftProportion = 1;
						while (hh != hhStart)
						{
							newElement.hhEdge = hh;
							iterPE = iterIP->m_lstData.insert(iterPE, newElement);
							hh = hh->opposite()->prev();
						}
						iterPE = iterPE1;
						iterPELast = iterIP->m_lstData.end();
						--iterPELast;
					}
					else if ((iterPE->iDir == 0 && iterPE->ftProportion == 1) || (iterPE->iDir == 1 && iterPE->ftProportion == 0))
					{
						if (iterPE->iDir == 0)
						{
							vhPole = iterPE->hhEdge->vertex();
						}
						else
						{
							vhPole = iterPE->hhEdge->prev()->vertex();
						}
						iterPE0 = iterPE;
						--iterPE0;
						while (((iterPE0->iDir == 0 && iterPE0->hhEdge->vertex() == vhPole) || (iterPE0->iDir == 1 && iterPE0->hhEdge->prev()->vertex() == vhPole)))
						{
							if (iterPE0 == iterIP->m_lstData.begin())
							{
								--iterPE0;
								break;
							}
							else
							{
								--iterPE0;
							}
						}
						++iterPE0;
						if (iterPE0->iDir == 0)
						{
							hhStart = iterPE0->hhEdge->opposite();
						}
						else
						{
							hhStart = iterPE0->hhEdge;
						}
						if (iterPE0 == iterIP->m_lstData.begin())
						{
							iterPE0->iDir = 1 - iterPE0->iDir;
							++iterPE0;
						}
						iterPE1 = iterPE;
						++iterPE1;
						while (((iterPE1->iDir == 0 && iterPE1->hhEdge->vertex() == vhPole) || (iterPE1->iDir == 1 && iterPE1->hhEdge->prev()->vertex() == vhPole)))
						{
							if (iterPE1 == iterPELast)
							{
								++iterPE1;
								break;
							}
							else
							{
								++iterPE1;
							}
						}
						iterPE = iterPE1;
						--iterPE;
						if (iterPE->iDir == 0)
						{
							hh = iterPE->hhEdge->next();
						}
						else
						{
							hh = iterPE->hhEdge->opposite()->next();
						}
						if (iterPE1 == iterIP->m_lstData.end())
						{
							iterPE->iDir = 1 - iterPE->iDir;
							--iterPE1;
						}
						iterPE = iterPE1;
						iterIP->m_lstData.erase(iterPE0, iterPE1);
						newElement.iDir = 0;
						newElement.ftProportion = 0;
						while (hh != hhStart)
						{
							newElement.hhEdge = hh;
							iterPE = iterIP->m_lstData.insert(iterPE, newElement);
							hh = hh->opposite()->next();
						}
						iterPE = iterPE1;
						iterPELast = iterIP->m_lstData.end();
						--iterPELast;
					}
					else
					{
						++iterPE;
					}
				}
			}
		}
	}
}

void CTracer::m_fnMapToMesh()
{
	std::list<CInducedPath>::iterator iterPath;
	std::list<CPathElement>::iterator iterPE0, iterPE1;
	CMeshBase::CPolyhedron::Halfedge_handle hh0, hh1;
	CMeshBase::FT ftProportion0, ftProportion1;
	for (iterPath = m_lstIP.begin(); iterPath != m_lstIP.end(); ++iterPath)
	{
		iterPE1 = iterPath->m_lstData.begin();
		iterPE0 = iterPE1;
		++iterPE1;
		while (iterPE1 != iterPath->m_lstData.end())
		{
			if (iterPE0->iDir == 0)
			{
				hh0 = iterPE0->hhEdge;
				ftProportion0 = iterPE0->ftProportion;
			}
			else
			{
				hh0 = iterPE0->hhEdge->opposite();
				ftProportion0 = 1 - iterPE0->ftProportion;
			}
			if (iterPE1->iDir == 0)
			{
				hh1 = iterPE1->hhEdge->opposite();
				ftProportion1 = 1 - iterPE1->ftProportion;
			}
			else
			{
				hh1 = iterPE1->hhEdge;
				ftProportion1 = iterPE1->ftProportion;
			}
			if (hh0->next() == hh1)
			{
				if (ftProportion0 < 0.5 && ftProportion1 < 0.5)
				{
					hh0->bFreeEdge = true;
					hh0->opposite()->bFreeEdge = true;
				}
				if (ftProportion0 >= 0.5 && ftProportion1 >= 0.5)
				{
					hh1->bFreeEdge = true;
					hh1->opposite()->bFreeEdge = true;
				}
				if (ftProportion0 < 0.5 && ftProportion1 >= 0.5)
				{
					hh0->prev()->bFreeEdge = true;
					hh0->prev()->opposite()->bFreeEdge = true;
				}
			}
			else
			{
				if (ftProportion0 < 0.5 && ftProportion1 < 0.5)
				{
					hh1->bFreeEdge = true;
					hh1->opposite()->bFreeEdge = true;
				}
				if (ftProportion0 >= 0.5 && ftProportion1 >= 0.5)
				{
					hh0->bFreeEdge = true;
					hh0->opposite()->bFreeEdge = true;
				}
				if (ftProportion0 >= 0.5 && ftProportion1 < 0.5)
				{
					hh0->next()->bFreeEdge = true;
					hh0->next()->opposite()->bFreeEdge = true;
				}
			}
			iterPE0 = iterPE1;
			++iterPE1;
		}
	}
}

#include "CTShapeTracer.h"

CTShapeTracer::CTShapeTracer(CMeshBase::CPolyhedron& mshSurface, std::vector<CTrackKnot>& vecExtraKnots) : m_mshSurface(mshSurface), m_vecExtraKnots(vecExtraKnots)
{
	;
}

void CTShapeTracer::m_fnSift(std::vector<CTrack> &vecHeap, int iPos)
{
	int iRoot, iLCh, iRCh, iSCh, nLen;
	nLen = vecHeap.size();
	CTShapeTracer::CTrack trTemp;
	iRoot = iPos;
	iLCh = iRoot * 2 + 1;
	iRCh = iLCh + 1;
	while (iLCh < nLen)
	{
		if (iRCh < nLen && vecHeap[iRCh].ftAccLen < vecHeap[iLCh].ftAccLen)
		{
			iSCh = iRCh;
		}
		else
		{
			iSCh = iLCh;
		}
		if (vecHeap[iRoot].ftAccLen < vecHeap[iRoot].ftAccLen)
		{
			break;
		}
		else
		{
			trTemp = vecHeap[iSCh];
			vecHeap[iSCh] = vecHeap[iRoot];
			vecHeap[iRoot] = trTemp;
			iRoot = iSCh;
			iLCh = iSCh * 2 + 1;
			iRCh = iLCh + 1;
		}
	}
}

void CTShapeTracer::m_fnCheck(std::vector<CTShapeTracer::CTrack> &vecHeap, int iPos)
{
	int iChild, iRoot;
	iChild = iPos;
	CTShapeTracer::CTrack trTemp;
	while (iChild != 0)
	{
		iRoot = (iChild - 1) / 2;
		if (vecHeap[iChild].ftAccLen < vecHeap[iRoot].ftAccLen)
		{
			trTemp = vecHeap[iChild];
			vecHeap[iChild] = vecHeap[iRoot];
			vecHeap[iRoot] = trTemp;
			iChild = iRoot;
		}
		else
		{
			break;
		}
	}
}

int CTShapeTracer::m_fnOnPath(CMeshBase::CPolyhedron::Halfedge_handle hh0, int iDir)
{
	int iPathID;
	CMeshBase::CPolyhedron::Halfedge_handle hh;
	int iLocalInd;
	iPathID = -1;
	hh = hh0;
	do
	{
		iLocalInd = hh->facet()->lpiTraceLocalInd[(iDir + 4) % 4];
		if (iLocalInd == hh->iLocalInd + 3)
		{
			iPathID = hh->facet()->lpiPathId[(iDir + 4) % 4];
		}
		iDir = (iDir - hh->next()->iSeamType + 4) % 4;
		hh = hh->next()->opposite();
	} while (hh != hh0);
	return iPathID;
}


void CTShapeTracer::m_fnFillGridCoef(CMeshBase::CPolyhedron::Halfedge_handle hh0, int iDir, int iPathID)
{
	if (!hh0->is_border())
	{
		hh0->facet()->lpiTraceInd[(iDir + 4) % 4] = hh0->vertex()->iIndex;
		hh0->facet()->lpiTraceLocalInd[(iDir + 4) % 4] = hh0->iLocalInd + 3;
		hh0->facet()->lpiPathId[(iDir + 4) % 4] = iPathID;
		hh0->facet()->lpvtTrace[(iDir + 4) % 4].m_lptCoord[hh0->iLocalInd] = 1.0;
		hh0->facet()->lpvtTrace[(iDir + 4) % 4].m_lptCoord[hh0->next()->iLocalInd] = 0.0;
		hh0->facet()->lpvtTrace[(iDir + 4) % 4].m_lptCoord[hh0->prev()->iLocalInd] = 0.0;
	}
}

int CTShapeTracer::m_fnFillIndOut(CMeshBase::CPolyhedron::Halfedge_handle hh0, int iDir, CMeshBase::FT ftProportion, int iPathID)
{
	int iType, i, lpiLocalInd[4], iTerm;
	/*iType:
		0 = continue to extend from a vertex
		1 = terminate at a vertex
		2 = terminate on a parallel edge 
		3 = terminate in a facet
		4 = continue to extend from an edge
		5 = terminate on a cross edge
		6 = terminate on a boder edge
		7 = delete path
	*/
	CMeshBase::CPolyhedron::Halfedge_handle hhPath;
	CMeshBase::CPolyhedron::Vector_3 vtNew, vtIntersection;
	CTShapeTracer::CTrackKnot newKnot;
	iType = -1;
	for (i = 0; i < 4; ++i)
	{
		lpiLocalInd[i] = hh0->facet()->lpiTraceLocalInd[i];
	}
	if (ftProportion == 1.0)
	{
		lpiLocalInd[(iDir + 4) % 4] = hh0->iLocalInd + 3;
	}
	else
	{
		lpiLocalInd[(iDir + 4) % 4] = hh0->iLocalInd;
	}
	if (lpiLocalInd[(iDir + 4) % 4] >= 3 && lpiLocalInd[(iDir + 4) % 4] < 6) //path reach a vertex
	{
		iType = 0;
		if (lpiLocalInd[(iDir + 6) % 4] >= 3 && lpiLocalInd[(iDir + 6) % 4] < 6) // path goes along an edge
		{
			switch ((lpiLocalInd[(iDir + 4) % 4] + 3 - lpiLocalInd[(iDir + 6) % 4]) % 3) // find edge
			{
			case 1:
				hhPath = hh0;
				if (hhPath->facet()->lpiTraceLocalInd[(iDir + 3) % 4] == hhPath->iLocalInd) //path cross with another perpendicular path on edge
				{
					hhPath->facet()->lpiTraceLocalInd[(iDir + 4) % 4] = hhPath->iLocalInd;
					hhPath->facet()->lpvtTrace[(iDir + 4) % 4] = hhPath->facet()->lpvtTrace[(iDir + 3) % 4];
					hhPath->facet()->lpiTraceInd[(iDir + 4) % 4] = hhPath->facet()->lpiTraceInd[(iDir + 3) % 4];
					hhPath->facet()->lpiPathId[(iDir + 4) % 4] = iPathID;
					m_vecPaths[iPathID].iTermID = hhPath->facet()->lpiPathId[(iDir + 3) % 4];
					hhPath->opposite()->facet()->lpiTraceLocalInd[(iDir + 4 - hhPath->iSeamType) % 4] = hhPath->opposite()->iLocalInd;
					hhPath->opposite()->facet()->lpvtTrace[(iDir + 4 - hhPath->iSeamType) % 4] = hhPath->opposite()->facet()->lpvtTrace[(iDir + 5 - hhPath->iSeamType) % 4];
					hhPath->opposite()->facet()->lpiTraceInd[(iDir + 4 - hhPath->iSeamType) % 4] = hhPath->opposite()->facet()->lpiTraceInd[(iDir + 5 - hhPath->iSeamType) % 4];
					hhPath->opposite()->facet()->lpiPathId[(iDir + 4 - hhPath->iSeamType) % 4] = iPathID;
					m_fnFillGridCoef(hhPath->opposite(), (iDir - hhPath->iSeamType + 6) % 4, iPathID);
					iType = 2;
				}
				else if (hhPath->vertex()->iBordId != -1 || hhPath->vertex()->nNumPorts != 0) // path cross with another perpendicular path at vertex
				{
					iType = 1;
				}
				else
				{
					iTerm = m_fnOnPath(hhPath, (iDir + 5) % 4);
					if (iTerm != -1)
					{
						m_vecPaths[iPathID].iTermID = iTerm;
						iType = 1;
					}
					else
					{
						iTerm = m_fnOnPath(hhPath, (iDir + 3) % 4);
						if (iTerm != -1)
						{
							m_vecPaths[iPathID].iTermID = iTerm;
							iType = 1;
						}
					}
				}
				if (iType == 0 || iType == 1)
				{
					m_fnFillGridCoef(hhPath, iDir, iPathID);
				}
				break;
			case 2:
				hhPath = hh0->next();
				if (hhPath->facet()->lpiTraceLocalInd[(iDir + 5) % 4] == hhPath->iLocalInd) //path cross with another perpendicular path on edge
				{
					hhPath->facet()->lpiTraceLocalInd[(iDir + 4) % 4] = hhPath->iLocalInd;
					hhPath->facet()->lpvtTrace[(iDir + 4) % 4] = hhPath->facet()->lpvtTrace[(iDir + 5) % 4];
					hhPath->facet()->lpiTraceInd[(iDir + 4) % 4] = hhPath->facet()->lpiTraceInd[(iDir + 5) % 4];
					hhPath->facet()->lpiPathId[(iDir + 4) % 4] = iPathID;
					m_vecPaths[iPathID].iTermID = hhPath->facet()->lpiPathId[(iDir + 5) % 4];
					hhPath->opposite()->facet()->lpiTraceLocalInd[(iDir + 4 - hhPath->iSeamType) % 4] = hhPath->opposite()->iLocalInd;
					hhPath->opposite()->facet()->lpvtTrace[(iDir + 4 - hhPath->iSeamType) % 4] = hhPath->opposite()->facet()->lpvtTrace[(iDir + 3 - hhPath->iSeamType) % 4];
					hhPath->opposite()->facet()->lpiTraceInd[(iDir + 4 - hhPath->iSeamType) % 4] = hhPath->opposite()->facet()->lpiTraceInd[(iDir + 3 - hhPath->iSeamType) % 4];
					hhPath->opposite()->facet()->lpiPathId[(iDir + 4 - hhPath->iSeamType) % 4] = iPathID;
					m_fnFillGridCoef(hhPath->opposite()->prev(), (iDir - hhPath->iSeamType + 6) % 4, iPathID);
					iType = 2;
				}
				else if (hhPath->prev()->vertex()->iBordId != -1 || hhPath->prev()->vertex()->nNumPorts != 0) // path cross with another perpendicular path at vertex
				{
					iType = 1;
				}
				else
				{
					iTerm = m_fnOnPath(hhPath->prev(), (iDir + 5) % 4);
					if (iTerm != -1)
					{
						m_vecPaths[iPathID].iTermID = iTerm;
						iType = 1;
					}
					else
					{
						iTerm = m_fnOnPath(hhPath->prev(), (iDir + 3) % 4);
						if (iTerm != -1)
						{
							m_vecPaths[iPathID].iTermID = iTerm;
							iType = 1;
						}
					}
				}
				if (iType == 0 || iType == 1)
				{
					m_fnFillGridCoef(hhPath->prev(), iDir, iPathID);
				}
				break;
			}
		}
		else //path split facet through a vertex
		{
			if ((lpiLocalInd[(iDir + 3) % 4] >= 0 && lpiLocalInd[(iDir + 3) % 4] != hh0->iLocalInd + 3) && 
				(lpiLocalInd[(iDir + 5) % 4] >= 0 && lpiLocalInd[(iDir + 5) % 4] != hh0->iLocalInd + 3) &&
				(lpiLocalInd[(iDir + 5) % 4] < 3 || lpiLocalInd[(iDir + 3) % 4] < 3)) //path cross with another perpendicular path in facet
			{
				vtNew.m_lptCoord[hh0->iLocalInd] = 1.0;
				vtNew.m_lptCoord[hh0->prev()->iLocalInd] = 0.0;
				vtNew.m_lptCoord[hh0->next()->iLocalInd] = 0.0;
				vtIntersection = CGAL::cross_product(CGAL::cross_product(hh0->facet()->lpvtTrace[(iDir + 6) % 4], vtNew), CGAL::cross_product(hh0->facet()->lpvtTrace[(iDir + 3) % 4], hh0->facet()->lpvtTrace[(iDir + 5) % 4]));
				vtIntersection = vtIntersection / (vtIntersection.x() + vtIntersection.y() + vtIntersection.z());
				if ((vtNew - vtIntersection) * (hh0->facet()->lpvtTrace[(iDir + 6) % 4] - vtIntersection) < 0 &&
					(hh0->facet()->lpvtTrace[(iDir + 5) % 4] - vtIntersection) * (hh0->facet()->lpvtTrace[(iDir + 3) % 4] - vtIntersection) < 0)
				{
					hh0->facet()->lpvtTrace[(iDir + 4) % 4] = vtIntersection;
					hh0->facet()->lpiTraceLocalInd[(iDir + 4) % 4] = 6;
					hh0->facet()->lpiTraceInd[(iDir + 4) % 4] = m_vecExtraKnots.size() + m_mshSurface.size_of_vertices();
					hh0->facet()->lpiPathId[(iDir + 4) % 4] = iPathID;
					m_vecPaths[iPathID].iTermID = hh0->facet()->lpiPathId[(iDir + 3) % 4];
					newKnot.ptPosition = CGAL::barycenter(hh0->facet()->halfedge()->vertex()->point(), vtIntersection.x(), hh0->facet()->halfedge()->next()->vertex()->point(), vtIntersection.y(), hh0->facet()->halfedge()->prev()->vertex()->point(), vtIntersection.z());
					m_vecExtraKnots.push_back(newKnot);
					iType = 3;
				}
				else 
				{
					m_fnFillGridCoef(hh0, (iDir + 4) % 4, iPathID);
					if (hh0->vertex()->iBordId != -1 || hh0->vertex()->nNumPorts != 0) // path cross with another perpendicular path at vertex
					{
						iType = 1;
					}
					else
					{
						iTerm = m_fnOnPath(hh0, (iDir + 5) % 4);
						if (iTerm != -1)
						{
							m_vecPaths[iPathID].iTermID = iTerm;
							iType = 1;
						}
						else
						{
							iTerm = m_fnOnPath(hh0, (iDir + 3) % 4);
							if (iTerm != -1)
							{
								m_vecPaths[iPathID].iTermID = iTerm;
								iType = 1;
							}
						}
					}
				}
			}
			else			
			{
				m_fnFillGridCoef(hh0, (iDir + 4) % 4, iPathID);
				if (hh0->vertex()->iBordId != -1 || hh0->vertex()->nNumPorts != 0) // path cross with another perpendicular path at vertex
				{
					iType = 1;
				}
				else
				{
					iTerm = m_fnOnPath(hh0, (iDir + 5) % 4);
					if (iTerm != -1)
					{
						m_vecPaths[iPathID].iTermID = iTerm;
						iType = 1;
					}
					else
					{
						iTerm = m_fnOnPath(hh0, (iDir + 3) % 4);
						if (iTerm != -1)
						{
							m_vecPaths[iPathID].iTermID = iTerm;
							iType = 1;
						}
					}
				}
			}
		}
	}
	else //path cross an edge
	{
		iType = 4;
		vtNew.m_lptCoord[hh0->iLocalInd] = ftProportion;
		vtNew.m_lptCoord[hh0->prev()->iLocalInd] = 1.0 - ftProportion;
		vtNew.m_lptCoord[hh0->next()->iLocalInd] = 0.0;
		if ((lpiLocalInd[(iDir + 3) % 4] != -1 && lpiLocalInd[(iDir + 5) % 4] != -1) && (lpiLocalInd[(iDir + 3) % 4] < 3 || lpiLocalInd[(iDir + 5) % 4] < 3) &&
			(lpiLocalInd[(iDir + 2) % 4] < 3 || (lpiLocalInd[(iDir + 2) % 4] != lpiLocalInd[(iDir + 3) % 4] && lpiLocalInd[(iDir + 2) % 4] != lpiLocalInd[(iDir + 5) % 4]))) // no coincided vertex
		{

			vtIntersection = CGAL::cross_product(CGAL::cross_product(hh0->facet()->lpvtTrace[(iDir + 6) % 4], vtNew), CGAL::cross_product(hh0->facet()->lpvtTrace[(iDir + 5) % 4], hh0->facet()->lpvtTrace[(iDir + 3) % 4]));
			vtIntersection = vtIntersection / (vtIntersection.x() + vtIntersection.y() + vtIntersection.z());
			if ((vtNew - vtIntersection) * (hh0->facet()->lpvtTrace[(iDir + 6) % 4] - vtIntersection) < 0 &&
				(hh0->facet()->lpvtTrace[(iDir + 5) % 4] - vtIntersection) * (hh0->facet()->lpvtTrace[(iDir + 3) % 4] - vtIntersection) < 0)//path cross with another perpendicular path on facet
			{
				hh0->facet()->lpvtTrace[(iDir + 4) % 4] = vtIntersection;
				hh0->facet()->lpiTraceLocalInd[(iDir + 4) % 4] = 6;
				hh0->facet()->lpiTraceInd[(iDir + 4) % 4] = m_vecExtraKnots.size() + m_mshSurface.size_of_vertices();
				hh0->facet()->lpiPathId[(iDir + 4) % 4] = iPathID;
				m_vecPaths[iPathID].iTermID = hh0->facet()->lpiPathId[(iDir + 3) % 4];
				newKnot.ptPosition = CGAL::barycenter(hh0->facet()->halfedge()->vertex()->point(), vtIntersection.x(), hh0->facet()->halfedge()->next()->vertex()->point(), vtIntersection.y(), hh0->facet()->halfedge()->prev()->vertex()->point(), vtIntersection.z());
				m_vecExtraKnots.push_back(newKnot);
				iType = 3;
			}
			else
			{
				hh0->facet()->lpvtTrace[(iDir + 4) % 4] = vtNew;
				hh0->facet()->lpiTraceLocalInd[(iDir + 4) % 4] = hh0->iLocalInd;
				hh0->facet()->lpiTraceInd[(iDir + 4) % 4] = m_vecExtraKnots.size() + m_mshSurface.size_of_vertices();
				hh0->facet()->lpiPathId[(iDir + 4) % 4] = iPathID;
				newKnot.ptPosition = CGAL::barycenter(hh0->facet()->halfedge()->vertex()->point(), vtNew.x(), hh0->facet()->halfedge()->next()->vertex()->point(), vtNew.y(), hh0->facet()->halfedge()->prev()->vertex()->point(), vtNew.z());
				m_vecExtraKnots.push_back(newKnot);
			}
		}
		else
		{
			hh0->facet()->lpvtTrace[(iDir + 4) % 4] = vtNew;
			hh0->facet()->lpiTraceLocalInd[(iDir + 4) % 4] = hh0->iLocalInd;
			hh0->facet()->lpiTraceInd[(iDir + 4) % 4] = m_vecExtraKnots.size() + m_mshSurface.size_of_vertices();
			hh0->facet()->lpiPathId[(iDir + 4) % 4] = iPathID;
			newKnot.ptPosition = CGAL::barycenter(hh0->facet()->halfedge()->vertex()->point(), vtNew.x(), hh0->facet()->halfedge()->next()->vertex()->point(), vtNew.y(), hh0->facet()->halfedge()->prev()->vertex()->point(), vtNew.z());
			m_vecExtraKnots.push_back(newKnot);
			if (hh0->opposite()->is_border())
			{
				iType = 6;
			}
			else if ((hh0->facet()->lpiTraceLocalInd[(iDir + 3) % 4] == hh0->prev()->iLocalInd + 3 && hh0->facet()->lpiTraceLocalInd[(iDir + 5) % 4] == hh0->iLocalInd + 3) ||
				(hh0->opposite()->facet()->lpiTraceLocalInd[(iDir - hh0->iSeamType + 5) % 4] == hh0->opposite()->prev()->iLocalInd + 3 && hh0->opposite()->facet()->lpiTraceLocalInd[(iDir - hh0->iSeamType + 5) % 4] == hh0->opposite()->iLocalInd + 3))
			{
				m_vecPaths[iPathID].iTermID = hh0->facet()->lpiPathId[(iDir + 3) % 4];
				iType = 5;
			}
		}
	}
	return iType;
}

void CTShapeTracer::m_fnDeletePath(int iPathID)
{
	CMeshBase::CPolyhedron::Facet_iterator fi;
	std::queue<int> quPathID;
	std::vector<CTShapeTracer::CTrack>::iterator iterPath;
	quPathID.push(iPathID);
	int i, iDeleting;
	while (!quPathID.empty())
	{
		iDeleting = quPathID.front();
		quPathID.pop();
		if (iDeleting >= 0 && iDeleting < m_vecPaths.size())
		{
			for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
			{
				for (i = 0; i < 4; ++i)
				{
					if (fi->lpiPathId[i] == iDeleting)
					{
						fi->lpiPathId[i] = -1;
						fi->lpiTraceInd[i] = -1;
						fi->lpiTraceLocalInd[i] = -1;
						fi->lpvtTrace[i] = CMeshBase::Vector_3(-1, -1, -1);
					}
				}
			}
			m_vecPaths[iDeleting].iStatus = 0;
			m_vecPaths[iDeleting].iTermID = -1;
			for (iterPath = m_vecPaths.begin(); iterPath != m_vecPaths.end(); ++iterPath)
			{
				if (iterPath->iTermID == iDeleting)
				{
					quPathID.push(iterPath->iPathID);
				}
			}
		}
	}
}

void CTShapeTracer::m_fnRefillPort(std::vector<CTrack>& vecHeap)
{
	std::vector<CTrack>::iterator iterPath;
	for (iterPath = m_vecPaths.begin(); iterPath != m_vecPaths.end(); ++iterPath)
	{
		if (iterPath->iStatus == 0 && iterPath->hhEdge->facet()->lpiTraceLocalInd[(iterPath->iDir + 2) % 4] == -1)
		{
			m_fnFillGridCoef(iterPath->hhEdge, (iterPath->iDir + 2) % 4, iterPath->iPathID);
			vecHeap.push_back((*iterPath));
		}
	}
}

void CTShapeTracer::m_fnFindExtend(CMeshBase::CPolyhedron::Halfedge_handle& hh0, int& iDir)
{
	CMeshBase::CPolyhedron::Halfedge_handle hh, hhMin;
	CMeshBase::FT ftAngle, ftAngleSum, ftMinAngleSum;
	int iDirCur, iDirMin;
	hhMin = hh0;
	ftMinAngleSum = ftTwoPi;
	iDirMin = 0;
	iDirCur = iDir;
	hh = hh0;
	do
	{
		ftAngle = hh->facet()->ftChartDir + iDirCur * ftHalfPi;
		switch (hh->iLocalInd)
		{
		case 0:
			ftAngle -= hh->ftAngle;
			break;
		case 1:
			ftAngle += hh->next()->ftAngle;
			break;
		}
		ftAngle -= CMeshBase::FT(int(ftAngle / ftTwoPi + 15.5) - 15) * ftTwoPi;
		ftAngleSum = fabs(ftAngle * 2 - ftPi + hh->ftAngle);
		if (ftAngleSum < ftMinAngleSum)
		{
			iDirMin = iDirCur;
			ftMinAngleSum = ftAngleSum;
			hhMin = hh;
		}
		iDirCur = (iDirCur - hh->iSeamType + 6) % 4 - 2;
		hh = hh->opposite()->prev();
	} while (hh != hh0);
	ftAngle = hhMin->facet()->ftChartDir + iDirMin * ftHalfPi;
	switch (hhMin->iLocalInd)
	{
	case 0:
		ftAngle -= hhMin->ftAngle;
		break;
	case 1:
		ftAngle += hhMin->next()->ftAngle;
		break;
	}
	ftAngle -= CMeshBase::FT(int(ftAngle / ftTwoPi + 15.5) - 15) * ftTwoPi;
	if (ftAngle < 0)
	{
		hh0 = hhMin->next()->opposite();
		iDir = (iDirMin - hhMin->next()->iSeamType + 6) % 4 - 2;
	}
	else
	{
		hh0 = hhMin;
		iDir = iDirMin;
	}
}

void CTShapeTracer::m_fnCountPort(std::vector<CTShapeTracer::CTrack> &vecHeap)
{
	CMeshBase::CPolyhedron::Vertex_iterator vi;
	CMeshBase::CPolyhedron::Halfedge_handle hh0, hh;
	int iInd, iDir;
	CMeshBase::FT ftAngle, ftAngleDif;
	CTShapeTracer::CTrack newTrack;
	newTrack.ftAccLen = 0;
	newTrack.iTermID = -1;
	newTrack.iPathID = 0;
	iInd = 0;
	for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
	{
		vi->iIndex = iInd;
		++iInd;
	}
	for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
	{
		if (vi->iBordId == -1 && vi->nNumPorts != 0)
		{
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
			iDir = (12 - int((ftAngle + ftAngleDif / 2) / ftHalfPi + 10)) % 4 - 2;
			hh = hh0;
			do
			{
				ftAngle = hh->facet()->ftChartDir + iDir * ftHalfPi;
				switch (hh->iLocalInd)
				{
				case 0:
					ftAngle -= hh->ftAngle;
					break;
				case 1:
					ftAngle += hh->next()->ftAngle;
					break;
				}
				ftAngleDif = 0;
				if (!hh->opposite()->is_border())
				{
					ftAngleDif = hh->opposite()->facet()->ftChartDir - hh->facet()->ftChartDir + hh->ftPolarAxisDif - hh->iSeamType * ftHalfPi;
					ftAngleDif = ftAngleDif - CMeshBase::FT(int(ftAngleDif / ftTwoPi + 15.5) - 15) * ftTwoPi;
				}
				ftAngle = ftAngle - CMeshBase::FT(int(ftAngle / ftTwoPi + 15.5) - 15) * ftTwoPi;
				if (ftAngle + ftAngleDif / 2 < ftPi - hh->ftAngle)
				{
					if (ftAngle < 0)
					{
						newTrack.hhEdge = hh->next()->opposite();
						newTrack.iDir = (iDir + 6 - hh->next()->iSeamType) % 4 - 2;
						newTrack.ftProportion = 1.0;
						hh->next()->opposite()->facet()->lpiTraceLocalInd[(newTrack.iDir + 6) % 4] = hh->next()->opposite()->iLocalInd + 3;
						hh->next()->opposite()->facet()->lpiPathId[(newTrack.iDir + 6) % 4] = newTrack.iPathID;
						hh->next()->opposite()->facet()->lpvtTrace[(newTrack.iDir + 6) % 4].m_lptCoord[hh->next()->opposite()->iLocalInd] = 1.0;
						hh->next()->opposite()->facet()->lpvtTrace[(newTrack.iDir + 6) % 4].m_lptCoord[hh->next()->opposite()->prev()->iLocalInd] = 0.0;
						hh->next()->opposite()->facet()->lpvtTrace[(newTrack.iDir + 6) % 4].m_lptCoord[hh->next()->opposite()->next()->iLocalInd] = 0.0;
						hh->next()->opposite()->facet()->lpiTraceInd[(newTrack.iDir + 6) % 4] = hh->vertex()->iIndex;
					}
					else
					{
						newTrack.hhEdge = hh;
						newTrack.iDir = iDir;
						newTrack.ftProportion = 1.0;
						hh->facet()->lpiTraceLocalInd[(iDir + 6) % 4] = hh->iLocalInd + 3;
						hh->facet()->lpiPathId[(iDir + 6) % 4] = newTrack.iPathID;
						hh->facet()->lpvtTrace[(iDir + 6) % 4].m_lptCoord[hh->iLocalInd] = 1.0;
						hh->facet()->lpvtTrace[(iDir + 6) % 4].m_lptCoord[hh->prev()->iLocalInd] = 0.0;
						hh->facet()->lpvtTrace[(iDir + 6) % 4].m_lptCoord[hh->next()->iLocalInd] = 0.0;
						hh->facet()->lpiTraceInd[(iDir + 6) % 4] = hh->vertex()->iIndex;
					}
					newTrack.iSingID = vi->iIndex;
					newTrack.iStatus = 0;
					vecHeap.push_back(newTrack);
					m_vecPaths.push_back(newTrack);
					++newTrack.iPathID;
					iDir = (iDir + 3) % 4 - 2;
				}
				iDir = (iDir - hh->iSeamType + 6) % 4 - 2;
				hh = hh->opposite()->prev();
			} while (hh != hh0);
		}
	}
}

void CTShapeTracer::m_fnTrace(std::vector<CTrack>& vecHeap)
{
	CTShapeTracer::CTrack *lptrCur;
	CMeshBase::Point_3 ptIn, ptOut;
	CMeshBase::CPolyhedron::Halfedge_handle hhIn, hhOut;
	CMeshBase::FT ftAngle, ftProportionIn, ftProportionOut, ftDivideAngle, ftPosRound, ftNegRound;
	CMeshBase::Vector_3 vtDivideVector;
	int iType;
	int iDirIn, iDirOut;
	while (!vecHeap.empty())
	{
		lptrCur = &(vecHeap[0]);
		if (!lptrCur->hhEdge->is_border() && lptrCur->hhEdge->facet()->bSelected)
		{
			iDirIn = 0;
		}
		if (lptrCur->hhEdge->facet()->lpiTraceLocalInd[(lptrCur->iDir + 4) % 4] == -1)
		{
			hhIn = lptrCur->hhEdge;
			ftProportionIn = lptrCur->ftProportion;
			iDirIn = lptrCur->iDir;
			ptIn = CGAL::barycenter(hhIn->vertex()->point(), ftProportionIn, hhIn->prev()->vertex()->point(), 1 - ftProportionIn);
			vtDivideVector = hhIn->next()->vertex()->point() - ptIn;
			ftDivideAngle = atan2(sqrt(CGAL::squared_length(CGAL::cross_product(hhIn->vtVector, vtDivideVector))), hhIn->vtVector * vtDivideVector);
			ftAngle = hhIn->facet()->ftChartDir + iDirIn * ftHalfPi;
			if (hhIn->iLocalInd == 1)
			{
				ftAngle -= hhIn->prev()->ftAngle;
			}
			if (hhIn->iLocalInd == 2)
			{
				ftAngle += hhIn->ftAngle;
			}
			ftAngle -= (CMeshBase::FT(int(ftAngle / ftTwoPi + 15.25) - 15) * ftTwoPi);
			if (ftAngle + hhIn->prev()->ftAngle < ftPi)
			{
				hhOut = hhIn->next();
				if (ftAngle < 0.0)
				{
					ftProportionOut = 0.0;
				}
				else
				{
					ftProportionOut = hhIn->ftLen * (1 - ftProportionIn) * sin(ftAngle) / (hhIn->next()->ftLen * sin(hhIn->ftAngle - ftAngle));
				}
				ftPosRound = ftDivideAngle - ftAngle;
				ftNegRound = ftAngle;
			}
			else if (ftAngle > hhIn->ftAngle)
			{
				hhOut = hhIn->prev();
				if (ftAngle > ftPi)
				{
					ftProportionOut = 1.0;
				}
				else
				{
					ftProportionOut = 1.0 - hhIn->ftLen * ftProportionIn * sin(ftAngle) / (hhIn->prev()->ftLen * sin(ftAngle + hhIn->prev()->ftAngle - ftPi));
				}
				ftPosRound = ftPi - ftAngle;
				ftNegRound = ftAngle - ftDivideAngle;
			}
			else
			{
				if (sin(hhIn->ftAngle - ftAngle) * sin(hhIn->prev()->ftAngle) * ftProportionIn > sin(ftAngle + hhIn->prev()->ftAngle - ftPi) * sin(hhIn->ftAngle) * (1 - ftProportionIn))
				{
					hhOut = hhIn->next();
					ftProportionOut = hhIn->ftLen * (1 - ftProportionIn) * sin(ftAngle) / (hhIn->next()->ftLen * sin(hhIn->ftAngle - ftAngle));
					ftPosRound = ftDivideAngle - ftAngle;
					ftNegRound = ftAngle;
				}
				else
				{
					hhOut = hhIn->prev();
					ftProportionOut = 1.0 - hhIn->ftLen * ftProportionIn * sin(ftAngle) / (hhIn->prev()->ftLen * sin(ftAngle + hhIn->prev()->ftAngle - ftPi));
					ftPosRound = ftPi - ftAngle;
					ftNegRound = ftAngle - ftDivideAngle;
				}
			}
			if (ftNegRound < ftPosRound && ftNegRound < 1.0 / 8)
			{
				ftProportionOut = 1.0;
				hhOut = hhOut->prev();
			}
			else if (ftPosRound < ftNegRound && ftPosRound < 1.0 / 8)
			{
				ftProportionOut = 1.0;
			}
			iType = m_fnFillIndOut(hhOut, iDirIn, ftProportionOut, lptrCur->iPathID);
			/*iType:
			0 = continue to extend from a vertex
			1 = terminate at a vertex
			2 = terminate on a parallel edge
			3 = terminate in a facet
			4 = continue to extend from an edge
			5 = terminate on a cross edge
			6 = terminate on a boder edge
			7 = delete path
			*/
			if (iType == 0)
			{
				if (hhOut->vertex()->iBordId == -1 && hhOut->vertex()->nNumPorts == 0 && m_fnOnPath(hhOut, (iDirIn + 2) % 4) != -1)
				{
					m_fnDeletePath(lptrCur->iPathID);
					iType = 7;
				}
				else
				{
					hhIn = hhOut;
					iDirOut = iDirIn;
					m_fnFindExtend(hhIn, iDirOut);
					if (hhIn->facet()->lpiTraceLocalInd[(iDirOut + 2) % 4] != -1)
					{
						m_fnDeletePath(lptrCur->iPathID);
						iType = 7;
					}
					else
					{
						m_fnFillGridCoef(hhIn, (iDirOut + 2) % 4, lptrCur->iPathID);
					}
				}
			}
			else if (iType == 4 || iType == 5)
			{
				if (hhOut->opposite()->facet()->lpiTraceLocalInd[(iDirIn - hhOut->iSeamType + 6) % 4] != -1)
				{
					m_fnDeletePath(lptrCur->iPathID);
					iType = 7;
				}
				else
				{
					iDirOut = (iDirIn + 4 - hhOut->iSeamType) % 4;
					hhOut->opposite()->facet()->lpiTraceLocalInd[(iDirOut + 6) % 4] = hhOut->opposite()->iLocalInd;
					hhOut->opposite()->facet()->lpiTraceInd[(iDirOut + 6) % 4] = hhOut->facet()->lpiTraceInd[(iDirIn + 4) % 4];
					hhOut->opposite()->facet()->lpiPathId[(iDirOut + 6) % 4] = lptrCur->iPathID;
					hhOut->opposite()->facet()->lpvtTrace[(iDirOut + 6) % 4].m_lptCoord[hhOut->opposite()->iLocalInd] = 1.0 - ftProportionOut;
					hhOut->opposite()->facet()->lpvtTrace[(iDirOut + 6) % 4].m_lptCoord[hhOut->opposite()->prev()->iLocalInd] = ftProportionOut;
					hhOut->opposite()->facet()->lpvtTrace[(iDirOut + 6) % 4].m_lptCoord[hhOut->opposite()->next()->iLocalInd] = 0.0;
				}
			}
			else if (iType == 7)
			{
				m_fnDeletePath(lptrCur->iPathID);
			}
			if (iType == 0)
			{
				ptOut = hhOut->vertex()->point();
				lptrCur->ftProportion = 1.0;
				lptrCur->hhEdge = hhIn;
				lptrCur->iDir = iDirOut;
				lptrCur->ftAccLen += sqrt(CGAL::squared_distance(ptIn, ptOut));
				m_fnSift(vecHeap, 0);
			}
			else if (iType == 4)
			{
				ptOut = m_vecExtraKnots[hhOut->facet()->lpiTraceInd[(iDirIn + 4) % 4] - m_mshSurface.size_of_vertices()].ptPosition;
				lptrCur->ftProportion = 1.0 - ftProportionOut;
				lptrCur->hhEdge = hhOut->opposite();
				lptrCur->iDir = (iDirIn - hhOut->iSeamType + 6) % 4 - 2;
				lptrCur->ftAccLen += sqrt(CGAL::squared_distance(ptIn, ptOut));
				m_fnSift(vecHeap, 0);
			}
			else
			{
				if (iType != 7)
				{
					m_vecPaths[lptrCur->iPathID].iStatus = 1;
				}
				vecHeap[0] = *(vecHeap.rbegin());
				vecHeap.pop_back();
				m_fnSift(vecHeap, 0);
			}
		}
		else
		{
			m_fnDeletePath(lptrCur->iPathID);
			vecHeap[0] = *(vecHeap.rbegin());
			vecHeap.pop_back();
			m_fnSift(vecHeap, 0);
		}
	}
}

void CTShapeTracer::m_fnTidyKnots()
{
	CMeshBase::CPolyhedron::Facet_iterator fi;
	std::vector<CTrackKnot>::iterator iterKnot;
	int i = 0;
	for (iterKnot = m_vecExtraKnots.begin(); iterKnot != m_vecExtraKnots.end(); ++iterKnot)
	{
		iterKnot->iMark = -1;
	}
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
	{
		for (i = 0; i < 4; ++i)
		{
			if (fi->lpiTraceInd[i] != -1 && fi->lpiTraceInd[i] >= m_mshSurface.size_of_vertices())
			{
				m_vecExtraKnots[fi->lpiTraceInd[i] - m_mshSurface.size_of_vertices()].iMark = 0;
			}
		}
	}
	i = m_mshSurface.size_of_vertices();
	for (iterKnot = m_vecExtraKnots.begin(); iterKnot != m_vecExtraKnots.end(); ++iterKnot)
	{
		if (iterKnot->iMark != -1)
		{
			iterKnot->iMark = i;
			++i;
		}
	}
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
	{
		for (i = 0; i < 4; ++i)
		{
			if (fi->lpiTraceInd[i] != -1 && fi->lpiTraceInd[i] >= m_mshSurface.size_of_vertices())
			{
				fi->lpiTraceInd[i] = m_vecExtraKnots[fi->lpiTraceInd[i] - m_mshSurface.size_of_vertices()].iMark;
			}
		}
	}
	i = 0;
	for (iterKnot = m_vecExtraKnots.begin(); iterKnot != m_vecExtraKnots.end(); ++iterKnot)
	{
		if (iterKnot->iMark != -1)
		{
			m_vecExtraKnots[i] = (*iterKnot);
			++i;
		}
	}
	while (m_vecExtraKnots.size() > i)
	{
		m_vecExtraKnots.pop_back();
	}
}

void CTShapeTracer::m_fnPartMark(int* lpiMark, int nSize, int* lpiSep, int &nParts)
{
	int i, lpiBuffer[4];
	for (i = 0; i < 4; ++i)
	{
		lpiBuffer[i] = -1;
	}
	nParts = 0;
	for (i = 0; i < nSize; ++i)
	{
		if (lpiMark[i] != -1)
		{
			lpiBuffer[nParts] = i;
			++nParts;
		}
	}
	if (nParts < 4)
	{
		for (i = 0; i < 4; ++i)
		{
			lpiSep[i] = lpiBuffer[i];
		}
	}
	else
	{
		if ((lpiMark[lpiBuffer[0]] + 4 - lpiMark[lpiBuffer[1]]) % 2 == 0)
		{
			for (i = 0; i < 4; ++i)
			{
				lpiSep[i] = lpiBuffer[i];
			}
		}
		else
		{
			for (i = 0; i < 4; ++i)
			{
				lpiSep[i] = lpiBuffer[(i + 1) % 4];
			}
		}
	}
}

void CTShapeTracer::m_fnGenerateFinSurface(CMeshBase::CPolyhedron& mshFin)
{
	CMeshBase::CPolyhedron::Vertex_iterator vi;
	std::vector<CTShapeTracer::CTrackKnot>::iterator iterKnot;
	CMeshBase::CPolyhedron::Facet_iterator fi;
	CMeshBase::CPolyhedron::Halfedge_handle hh0, hh;
	unsigned int lpiIndex[8];
	CMeshBase::FT lpftProportion[8], ftTemp;
	int iCenter, i, nSize, iRepeat, iTemp, iPviot, nShape, iStart, j, lpiEnum[8], lpiMark[8], lpiSep[4], nParts, lpiDraw[4];
	bool bOccupy;
	for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
	{
		mshFin.append_vertex(vi->point().x(), vi->point().y(), vi->point().z());
	}
	for (iterKnot = m_vecExtraKnots.begin(); iterKnot != m_vecExtraKnots.end(); ++iterKnot)
	{
		mshFin.append_vertex(iterKnot->ptPosition.x(), iterKnot->ptPosition.y(), iterKnot->ptPosition.z());
	}
	mshFin.begin_facets();
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
	{
		if (fi->bSelected)
		{
			i = 0;
		}
		for (i = 0; i < 8; ++i)
		{
			lpiEnum[i] = -1;
			lpiMark[i] = -1;
			lpftProportion[i] = 0.0;
		}
		for (i = 0; i < 2; ++i)
		{
			lpiDraw[i] = 0;
			lpiDraw[i + 2] = 0;
			if (fi->lpiTraceLocalInd[i] == -1 || fi->lpiTraceLocalInd[i + 2] == -1 || fi->lpiTraceLocalInd[i] == 6 || fi->lpiTraceLocalInd[i + 2] == 6)
			{
				lpiDraw[i] = 1;
				lpiDraw[i + 2] = 1;
			}
			else if((fi->lpiTraceLocalInd[i] >= 3 && fi->lpiTraceLocalInd[i] < 6) && (fi->lpiTraceLocalInd[i + 2] >= 3 && fi->lpiTraceLocalInd[i + 2] < 6))
			{
				lpiDraw[i] = 0;
				lpiDraw[i + 2] = 0;
			}
			else
			{
				lpiDraw[i] = 1;
				lpiDraw[i + 2] = 1;
			}
		}
		nSize = 0;
		hh0 = fi->halfedge();
		hh = hh0;
		do
		{
			iRepeat = 0;
			bOccupy = false;
			for (i = 0; i < 4; ++i)
			{
				if (fi->lpiTraceLocalInd[i] == hh->iLocalInd && fi->lpiTraceInd[i] != -1)
				{
					lpiEnum[nSize] = fi->lpiTraceInd[i];
					if (lpiDraw[i] == 1)
					{
						lpiMark[nSize] = i;
					}
					lpftProportion[nSize] = fi->lpvtTrace[i].m_lptCoord[fi->lpiTraceLocalInd[i]];
					if (iRepeat > 0)
					{
						j = 0;
						while (j < iRepeat)
						{
							if (lpftProportion[nSize - j - 1] > lpftProportion[nSize - j])
							{
								iTemp = lpiEnum[nSize - j - 1];
								lpiEnum[nSize - j - 1] = lpiEnum[nSize - j];
								lpiEnum[nSize - j] = iTemp;
								iTemp = lpiMark[nSize - j - 1];
								lpiMark[nSize - j - 1] = lpiMark[nSize - j];
								lpiMark[nSize - j] = iTemp;
								ftTemp = lpftProportion[nSize - j - 1];
								lpftProportion[nSize - j - 1] = lpftProportion[nSize - j];
								lpftProportion[nSize - j] = ftTemp;
							}
							++j;
						}
					}
					++nSize;
					++iRepeat;
				}
				else if (fi->lpiTraceLocalInd[i] == hh->iLocalInd + 3 && !bOccupy)
				{
					lpiEnum[nSize] = fi->lpiTraceInd[i];
					if (lpiDraw[i] == 1)
					{
						lpiMark[nSize] = i;
					}
					lpftProportion[nSize] = 1.0;
					++nSize;
					++iRepeat;
					bOccupy = true;
				}
			}
			if (!bOccupy)
			{
				lpiEnum[nSize] = hh->vertex()->iIndex;
				lpftProportion[nSize] = 1.0;
				++nSize;
			}
			hh = hh->next();
		} while (hh != hh0);
		iCenter = -1;
		for (i = 0; i < 4; ++i)
		{
			if (fi->lpiTraceLocalInd[i] == 6)
			{
				iCenter = i;
			}
		}
		if (iCenter == -1)
		{
			m_fnPartMark(lpiMark, nSize, lpiSep, nParts);
			if (nParts < 2)
			{
				mshFin.append_facet(3, hh0->vertex()->iIndex, hh0->next()->vertex()->iIndex, hh0->prev()->vertex()->iIndex);
			}
			else if (nParts < 4)
			{
				for (i = 0; i < nParts; ++i)
				{
					nShape = 0;
					for (j = lpiSep[i]; j != lpiSep[(i + 1) % nParts]; j = (j + 1) % nSize)
					{
						lpiIndex[nShape] = lpiEnum[j];
						++nShape;
					}
					lpiIndex[nShape] = lpiEnum[lpiSep[(i + 1) % nParts]];
					++nShape;
					if (nShape > 2)
					{
						mshFin.append_facet(lpiIndex, nShape);
					}
				}
			}
			else
			{
				nShape = 0;
				for (j = lpiSep[0]; j != lpiSep[1]; j = (j + 1) % nSize)
				{
					lpiIndex[nShape] = lpiEnum[j];
					++nShape;
				}
				lpiIndex[nShape] = lpiEnum[lpiSep[1]];
				++nShape;
				if (nShape > 2)
				{
					mshFin.append_facet(lpiIndex, nShape);
				}
				nShape = 0;
				for (j = lpiSep[2]; j != lpiSep[3]; j = (j + 1) % nSize)
				{
					lpiIndex[nShape] = lpiEnum[j];
					++nShape;
				}
				lpiIndex[nShape] = lpiEnum[lpiSep[3]];
				++nShape;
				if (nShape > 2)
				{
					mshFin.append_facet(lpiIndex, nShape);
				}
				nShape = 0;
				for (j = lpiSep[3]; j != lpiSep[0]; j = (j + 1) % nSize)
				{
					lpiIndex[nShape] = lpiEnum[j];
					++nShape;
				}
				lpiIndex[nShape] = lpiEnum[lpiSep[0]];
				++nShape;
				for (j = lpiSep[1]; j != lpiSep[2]; j = (j + 1) % nSize)
				{
					lpiIndex[nShape] = lpiEnum[j];
					++nShape;
				}
				lpiIndex[nShape] = lpiEnum[lpiSep[2]];
				++nShape;
				if (nShape > 2)
				{
					mshFin.append_facet(lpiIndex, nShape);
				}
			}
		}
		else
		{
			for (i = 0; i < nSize; ++i)
			{
				if (lpiMark[i] != -1)
				{
					iPviot = i;
				}
			}
			iStart = iPviot;
			do
			{
				nShape = 0;
				lpiIndex[nShape] = fi->lpiTraceInd[iCenter];
				++nShape;
				for (i = 0; i < nSize; ++i)
				{
					lpiIndex[nShape] = lpiEnum[(iStart + i) % nSize];
					++nShape;
					if (i != 0 && lpiMark[(iStart + i) % nSize] != -1)
					{
						iStart = (iStart + i) % nSize;
						i = nSize;
					}
				}
				if (nShape > 2)
				{
					mshFin.append_facet(lpiIndex, nShape);
				}
			} while (iStart != iPviot);
		}
	}
	mshFin.end_facets(false);
}

void CTShapeTracer::m_fnMarkSeparate(CMeshBase::CPolyhedron& mshFin)
{
	CMeshBase::CPolyhedron::Vertex_iterator vi;
	CMeshBase::CPolyhedron::Edge_iterator ei;
	CMeshBase::CPolyhedron::Vertex_handle* lpVertexBuffer;
	CMeshBase::CPolyhedron::Facet_iterator fi;
	CMeshBase::CPolyhedron::Halfedge_handle hh0, hh, hh1;
	int iIndex, i, iCenter;
	bool bSatisfy;
	lpVertexBuffer = new CMeshBase::CPolyhedron::Vertex_handle[mshFin.size_of_vertices()];
	iIndex = 0;
	for (vi = mshFin.vertices_begin(); vi != mshFin.vertices_end(); ++vi)
	{
		vi->iIndex = iIndex;
		lpVertexBuffer[iIndex] = &(*vi);
		++iIndex;
	}
	for (ei = mshFin.edges_begin(); ei != mshFin.edges_end(); ++ei)
	{
		ei->bFreeEdge = false;
		ei->opposite()->bFreeEdge = false;
	}
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
	{
		iCenter = -1;
		for (i = 0; i < 4; ++i)
		{
			if (fi->lpiTraceLocalInd[i] == 6)
			{
				iCenter = i;
			}
		}
		if (iCenter == -1)
		{
			for (i = 0; i < 2; ++i)
			{
				if (fi->lpiTraceInd[i] != -1 && fi->lpiTraceInd[i + 2] != -1)
				{
					bSatisfy = false;
					if (fi->lpiTraceInd[i] < 0 || fi->lpiTraceInd[i] >= mshFin.size_of_vertices())
					{
						i = i;
					}
					hh0 = lpVertexBuffer[fi->lpiTraceInd[i]]->halfedge();
					hh = hh0;
					do
					{
						if (hh->prev()->vertex()->iIndex == fi->lpiTraceInd[i + 2])
						{
							hh->bFreeEdge = true;
							hh->opposite()->bFreeEdge = true;
							bSatisfy = true;
						}
						hh = hh->next()->opposite();
					} while (hh != hh0);
					if (!bSatisfy)
					{
						if (fi->lpiTraceInd[i] < 0 || fi->lpiTraceInd[i] >= mshFin.size_of_vertices())
						{
							i = i;
						}
						hh0 = lpVertexBuffer[fi->lpiTraceInd[i]]->halfedge();
						hh = hh0;
						do
						{
							if (hh->prev()->vertex()->iIndex >= m_mshSurface.size_of_vertices())
							{
								hh1 = hh->prev();
								while (hh1 != hh->opposite())
								{
									if (hh1->prev()->vertex()->iIndex == fi->lpiTraceInd[i + 2])
									{
										hh->bFreeEdge = true;
										hh->opposite()->bFreeEdge = true;
										hh1->bFreeEdge = true;
										hh1->opposite()->bFreeEdge = true;
									}
									hh1 = hh1->opposite()->prev();
								}
							}
							hh = hh->next()->opposite();
						} while (hh != hh0);
					}
				}
			}
		}
		else
		{
			for (i = 1; i < 4; ++i)
			{
				if (fi->lpiTraceInd[(iCenter + i) % 4] != -1)
				{
					hh0 = lpVertexBuffer[fi->lpiTraceInd[iCenter]]->halfedge();
					hh = hh0;
					do
					{
						if (hh->prev()->vertex()->iIndex == fi->lpiTraceInd[(iCenter + i) % 4])
						{
							hh->bFreeEdge = true;
							hh->opposite()->bFreeEdge = true;
						}
						hh = hh->next()->opposite();
					} while (hh != hh0);
				}
			}
		}
	}
	delete[]lpVertexBuffer;
}

int CTShapeTracer::m_fnMarkSubFacet(CMeshBase::CPolyhedron& mshFin)
{
	CMeshBase::CPolyhedron::Facet_iterator fi;
	CMeshBase::CPolyhedron::Halfedge_handle hh0, hh;
	CMeshBase::CPolyhedron::Facet_handle fhSeed;
	std::queue<CMeshBase::CPolyhedron::Facet_handle> qufh;
	int iIndex;
	for (fi = mshFin.facets_begin(); fi != mshFin.facets_end(); ++fi)
	{
		fi->iTempIndex = -1;
	}
	iIndex = 0;
	for (fi = mshFin.facets_begin(); fi != mshFin.facets_end(); ++fi)
	{
		if (fi->iTempIndex == -1)
		{
			fi->iTempIndex = iIndex;
			qufh.push(&(*fi));
			while (!qufh.empty())
			{
				fhSeed = qufh.front();
				qufh.pop();
				hh0 = fhSeed->halfedge();
				hh = hh0;
				do
				{
					if (!(hh->opposite()->is_border() || hh->bFreeEdge))
					{
						if (hh->opposite()->facet()->iTempIndex == -1)
						{
							hh->opposite()->facet()->iTempIndex = fhSeed->iTempIndex;
							qufh.push(hh->opposite()->facet());
						}
					}
					hh = hh->next();
				} while (hh != hh0);
			}
			++iIndex;
		}
	}
	return iIndex;
}

void CTShapeTracer::m_fnExportSubFacet(CMeshBase::CPolyhedron& mshFin, std::ofstream& ofs)
{
	CMeshBase::CPolyhedron::Vertex_iterator vi;
	CMeshBase::CPolyhedron::Facet_iterator fi;
	CMeshBase::CPolyhedron::Halfedge_handle hh0, hh;
	CMeshBase::CPolyhedron::Edge_iterator ei;
	int i, nParts;
	nParts = 0;
	for (fi = mshFin.facets_begin(); fi != mshFin.facets_end(); ++fi)
	{
		if (fi->iTempIndex >= nParts)
		{
			nParts = fi->iTempIndex + 1;
		}
	}
	for (vi = mshFin.vertices_begin(); vi != mshFin.vertices_end(); ++vi)
	{
		ofs << "v " << vi->point().x() << ' ' << vi->point().y() << ' ' << vi->point().z() << " \n";
	}
	for (i = 0; i < nParts; ++i)
	{
		ofs << "g part" << i << " \n";
		for (fi = mshFin.facets_begin(); fi != mshFin.facets_end(); ++fi)
		{
			if (fi->iTempIndex == i)
			{
				ofs << "f ";
				hh0 = fi->halfedge();
				hh = hh0;
				do
				{
					ofs << hh->vertex()->iIndex + 1 << ' ';
					hh = hh->next();
				} while (hh != hh0);
				ofs << '\n';
			}
		}
	}
	ofs << "g trace\n";
	for (ei = m_mshSurface.edges_begin(); ei != m_mshSurface.edges_end(); ++ei)
	{
		if (ei->bFreeEdge)
		{
			ofs << "l " << ei->prev()->vertex()->iIndex + 1 << ' ' << ei->vertex()->iIndex + 1 << " \n";
		}
	}
}

void CTShapeTracer::m_fnDelaunay(CMeshBase::CPolyhedron& mshTri)
{
	CMeshBase::Kernel::FT ftCross1, ftCross2, ftCross3, ftCross4, ftDot1, ftDot2, ftDot3, ftDot4;
	CMeshBase::Vector_3 vtNorm1, vtNorm2, vtNorm3, vtNorm4;
	bool bComplete, bEdgeFlip, bMeshSplit;
	CMeshBase::CPolyhedron::Edge_iterator ei;
	bComplete = false;
	bEdgeFlip = false;
	bMeshSplit = false;
	while (!bComplete)
	{
		bComplete = true;
		for (ei = mshTri.edges_begin(); ei != mshTri.edges_end(); ++ei)
		{
			if (!(ei->is_border_edge()))
			{
				ei->vtVector = ei->vertex()->point() - ei->prev()->vertex()->point();
				ei->opposite()->vtVector = -ei->vtVector;
				ei->next()->vtVector = ei->next()->vertex()->point() - ei->vertex()->point();
				ei->prev()->vtVector = ei->prev()->vertex()->point() - ei->next()->vertex()->point();
				ei->opposite()->next()->vtVector = ei->opposite()->next()->vertex()->point() - ei->opposite()->vertex()->point();
				ei->opposite()->prev()->vtVector = ei->opposite()->prev()->vertex()->point() - ei->opposite()->next()->vertex()->point();
				ftDot1 = -(ei->next()->vtVector * ei->prev()->vtVector);
				ftDot2 = -(ei->prev()->vtVector * ei->opposite()->next()->vtVector);
				ftDot3 = -(ei->opposite()->next()->vtVector * ei->opposite()->prev()->vtVector);
				ftDot4 = -(ei->opposite()->prev()->vtVector * ei->next()->vtVector);
				vtNorm1 = CGAL::cross_product(ei->next()->vtVector, ei->prev()->vtVector);
				vtNorm2 = CGAL::cross_product(ei->prev()->vtVector, ei->opposite()->next()->vtVector);
				vtNorm3 = CGAL::cross_product(ei->opposite()->next()->vtVector, ei->opposite()->prev()->vtVector);
				vtNorm4 = CGAL::cross_product(ei->opposite()->prev()->vtVector, ei->next()->vtVector);
				ftCross1 = sqrt(vtNorm1 * vtNorm1);
				ftCross2 = sqrt(vtNorm2 * vtNorm2);
				ftCross3 = sqrt(vtNorm3 * vtNorm3);
				ftCross4 = sqrt(vtNorm4 * vtNorm4);
				if (atan2(ftCross1, ftDot1) + atan2(ftCross3, ftDot3) > atan2(ftCross2, ftDot2) + atan2(ftCross4, ftDot4))
				{
					mshTri.flip_edge(&(*ei));
					bComplete = false;
				}
			}
		}
	}
}

void CTShapeTracer::m_fnExportTriangle(CMeshBase::CPolyhedron& mshFin, std::ofstream& ofs, int iPart)
{
	CMeshBase::CPolyhedron::Vertex_iterator vi;
	CMeshBase::CPolyhedron::Facet_iterator fi;
	CMeshBase::CPolyhedron::Halfedge_handle hh0, hh;
	int iVertexInd, nShape, i;
	CMeshBase::CPolyhedron::Vertex_handle lpvhPolygon[8];
	CMeshBase::CPolyhedron mshTri;
	for (vi = mshFin.vertices_begin(); vi != mshFin.vertices_end(); ++vi)
	{
		vi->iTempIndex = -1;
	}
	iVertexInd = 0;
	for (fi = mshFin.facets_begin(); fi != mshFin.facets_end(); ++fi)
	{
		if (fi->iTempIndex == iPart)
		{
			hh0 = fi->halfedge();
			hh = hh0;
			do
			{
				if (hh->vertex()->iTempIndex == -1)
				{
					hh->vertex()->iTempIndex = iVertexInd;
					mshTri.append_vertex(hh->vertex()->point().x(), hh->vertex()->point().y(), hh->vertex()->point().z());
					++iVertexInd;
				}
				hh = hh->next();
			} while (hh != hh0);
		}
	}
	mshTri.begin_facets();
	for (fi = mshFin.facets_begin(); fi != mshFin.facets_end(); ++fi)
	{
		if (fi->iTempIndex == iPart)
		{
			nShape = 0;
			hh0 = fi->halfedge();
			hh = hh0;
			do
			{
				lpvhPolygon[nShape] = hh->vertex();
				++nShape;
				hh = hh->next();
			} while (hh != hh0);
			for (i = 0; i < nShape - 2; ++i)
			{
				mshTri.append_facet(3, lpvhPolygon[0]->iTempIndex, lpvhPolygon[i + 1]->iTempIndex, lpvhPolygon[i + 2]->iTempIndex);
			}
		}
	}
	mshTri.end_facets();
	m_fnDelaunay(mshTri);
	ofs << mshTri;
	mshTri.clear();
}

void CTShapeTracer::m_fnExportOriginalMesh(std::ofstream &ofs)
{
	CMeshBase::CPolyhedron::Vertex_iterator vi;
	CMeshBase::CPolyhedron::Facet_iterator fi;
	CMeshBase::CPolyhedron::Halfedge_handle hh0, hh;
	CMeshBase::CPolyhedron::Edge_iterator ei;
	for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
	{
		ofs << "v " << vi->point().x() << ' ' << vi->point().y() << ' ' << vi->point().z() << " \n";
	}
	ofs << "vt 0 0 \n";
	ofs << "vt 0.965925826 0.2588190451 \n";
	ofs << "vt 0.2588190451 0.965925826 \n";
	ofs << "g facet\n";
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
	{
		ofs << "f ";
		hh0 = fi->halfedge();
		hh = hh0;
		do
		{
			ofs << hh->vertex()->iIndex + 1 << '/' << hh->iLocalInd + 1 << ' ';
			hh = hh->next();
		} while (hh != hh0);
		ofs << "\n";
	}
	ofs << "g boundary\n";
	for (ei = m_mshSurface.edges_begin(); ei != m_mshSurface.edges_end(); ++ei)
	{
		if (ei->is_border_edge())
		{
			ofs << "l " << ei->prev()->vertex()->iIndex + 1 << ' ' << ei->vertex()->iIndex + 1 << " \n";
		}
	}
	ofs << "g inner_edges\n";
	for (ei = m_mshSurface.edges_begin(); ei != m_mshSurface.edges_end(); ++ei)
	{
		if (!ei->is_border_edge())
		{
			ofs << "l " << ei->prev()->vertex()->iIndex + 1 << ' ' << ei->vertex()->iIndex + 1 << " \n";
		}
	}
}

void CTShapeTracer::m_fnExportExtraKnots(std::ofstream& ofs, CMeshBase::CPolyhedron& mshFin)
{
	int i;
	CMeshBase::CPolyhedron::Vertex_iterator vi;
	for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
	{
		if (vi->iBordId == -1 && vi->nNumPorts != 0)
		{
			ofs << vi->point().x() << '\t' << vi->point().y() << '\t' << vi->point().z() << '\t' << vi->nNumPorts << "\t\n";
		}
	}
	i = 0;
	for (vi = mshFin.vertices_begin(); vi != mshFin.vertices_end(); ++vi)
	{
		vi->iIndex = i;
		if (vi->iIndex >= m_mshSurface.size_of_vertices())
		{
			ofs << vi->point().x() << '\t' << vi->point().y() << '\t' << vi->point().z() << '\t' << 4 << "\t\n";
		}
		++i;
	}
}

void CTShapeTracer::m_fnExportPaths(std::ofstream& ofs, CMeshBase::CPolyhedron& mshFin)
{
	CMeshBase::CPolyhedron::Vertex_iterator vi, ovi;
	CMeshBase::CPolyhedron::Halfedge_handle hh0, hh;
	CMeshBase::CPolyhedron::Edge_iterator ei;
	bool bOnPath;
	int iInd;
	for (vi = mshFin.vertices_begin(); vi != mshFin.vertices_end(); ++vi)
	{
		vi->iTempIndex = -1;
	}
	iInd = 0;
	for (vi = mshFin.vertices_begin(); vi != mshFin.vertices_end(); ++vi)
	{
		if (vi->iTempIndex == -1)
		{
			bOnPath = false;
			hh0 = vi->halfedge();
			hh = hh0;
			do
			{
				if (hh->bFreeEdge)
				{
					bOnPath = true;
				}
				hh = hh->next()->opposite();
			} while (hh != hh0);
			if (bOnPath)
			{
				ofs << "v " << vi->point().x() << ' ' << vi->point().y() << ' ' << vi->point().z() << " \n";
				vi->iTempIndex = iInd;
				++iInd;
			}
		}
	}
	ofs << "g path\n";
	for (ei = mshFin.edges_begin(); ei != mshFin.edges_end(); ++ei)
	{
		if (ei->bFreeEdge)
		{
			ofs << "l " << ei->prev()->vertex()->iTempIndex + 1 << ' ' << ei->vertex()->iTempIndex + 1 << " \n";
		}
	}
	ofs << "g port\n";
	vi = mshFin.vertices_begin();
	for (ovi = m_mshSurface.vertices_begin(); ovi != m_mshSurface.vertices_end(); ++ovi)
	{
		if (ovi->nNumPorts != 0)
		{
			hh0 = vi->halfedge();
			hh = hh0;
			do
			{
				if (hh->bFreeEdge)
				{
					ofs << "l " << hh->prev()->vertex()->iTempIndex + 1 << ' ' << hh->vertex()->iTempIndex + 1 << " \n";
				}
				hh = hh->next()->opposite();
			} while (hh != hh0);
		}
		++vi;
	}
}
#include "CMeshSplitter.h"
CMeshSplitter::CMeshSplitter(CMeshBase::CPolyhedron &mshSurface) : m_mshSurface(mshSurface)
{
	;
}


void CMeshSplitter::m_fnExtendFacet(CMeshBase::CPolyhedron::Facet_handle fh)
{
	std::queue<CMeshBase::CPolyhedron::Facet_handle> qufh;
	CMeshBase::CPolyhedron::Facet_handle fhSeed;
	CMeshBase::CPolyhedron::Halfedge_handle hh0, hh;
	qufh.push(fh);
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
				if (hh->opposite()->facet()->iTempIndex != fhSeed->iTempIndex)
				{
					hh->opposite()->facet()->iTempIndex = fhSeed->iTempIndex;
					qufh.push(hh->opposite()->facet());
				}
			}
			hh = hh->next();
		} while (hh != hh0);
	}
}

void CMeshSplitter::m_fnActiveBoundaries()
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

void CMeshSplitter::m_fnPartVertices()
{
	CMeshBase::CPolyhedron::Vertex_iterator vi;
	std::vector<CMeshBase::CPolyhedron::Vertex_handle> vecVH;
	CMeshBase::CPolyhedron::Vertex_handle vhSeed;
	CMeshBase::CPolyhedron::Halfedge_handle hh0, hh;
	CMeshBase::FT ftDist;
	int iTempIndex;
	iTempIndex = -1;
	for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
	{
		vi->iHeapPos = -1;
		vi->hhCome = NULL;
		if (vi->iBordId != -1)
		{
			vi->ftDist = 0;
			vi->iTempIndex = vi->iBordId;
			vi->bDistConf = true;
			if (vi->iTempIndex >= iTempIndex)
			{
				iTempIndex = vi->iTempIndex + 1;
			}
		}
		else
		{
			vi->bDistConf = false;
		}
	}
	if (iTempIndex == -1)
	{
		iTempIndex = 0;
	}
	for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
	{
		if (vi->iBordId == -1)
		{
			if (vi->nNumPorts != 0)
			{
				vi->iTempIndex = iTempIndex;
				vi->ftDist = 0;
				vi->iHeapPos = vecVH.size();
				vecVH.push_back(&(*vi));
				++iTempIndex;
			}
			else
			{
				vi->iTempIndex = -1;
			}
		}
	}
	while (!vecVH.empty())
	{
		vhSeed = vecVH[0];
		if (vecVH.size() > 1)
		{
			vecVH[0] = *(vecVH.rbegin());
		}
		vecVH.pop_back();
		if (vecVH.size() > 0)
		{
			vecVH[0]->iHeapPos = 0;
			if (vecVH.size() > 1)
			{
				m_fnSift(vecVH, 0);
			}
		}
		vhSeed->iHeapPos = -1;
		vhSeed->bDistConf = true;
		hh0 = vhSeed->halfedge()->opposite();
		hh = hh0;
		do
		{
			ftDist = vhSeed->ftDist + hh->ftLen;
			if (!hh->vertex()->bDistConf)
			{
				if (hh->vertex()->iTempIndex == -1)
				{
					hh->vertex()->ftDist = ftDist;
					hh->vertex()->hhCome = hh;
					hh->vertex()->iTempIndex = vhSeed->iTempIndex;
					hh->vertex()->iHeapPos = vecVH.size();
					vecVH.push_back(hh->vertex());
					m_fnCheck(vecVH, hh->vertex()->iHeapPos);
				}
				else if (ftDist < hh->vertex()->ftDist)
				{
					hh->vertex()->ftDist = ftDist;
					hh->vertex()->hhCome = hh;
					hh->vertex()->iTempIndex = vhSeed->iTempIndex;
					m_fnCheck(vecVH, hh->vertex()->iHeapPos);
				}
			}
			hh = hh->opposite()->next();
		} while (hh != hh0);
	}
}

void CMeshSplitter::m_fnSimplifyConnection()
{
	CMeshBase::CPolyhedron::Edge_iterator ei;
	CMeshBase::CPolyhedron::Halfedge_handle *lphhLinkage;
	CMeshBase::FT ftMinDist, ftDist;
	CMeshBase::CPolyhedron::Halfedge_handle hh0, hh;
	int nNum, nRemain, iWrite, iRead, iNewLinkage;
	bool bLinkable;
	nNum = 0;
	for (ei = m_mshSurface.edges_begin(); ei != m_mshSurface.edges_end(); ++ei)
	{
		if (ei->prev()->vertex()->iTempIndex != ei->vertex()->iTempIndex)
		{
			ei->iTempIndex = 1;
			ei->opposite()->iTempIndex = 1;
			++nNum;
		}
		else
		{
			ei->iTempIndex = 0;
			ei->opposite()->iTempIndex = 0;
		}
	}
	lphhLinkage = new CMeshBase::CPolyhedron::Halfedge_handle[nNum];
	iWrite = 0;
	for (ei = m_mshSurface.edges_begin(); ei != m_mshSurface.edges_end(); ++ei)
	{
		if (ei->prev()->vertex()->iTempIndex != ei->vertex()->iTempIndex)
		{
			lphhLinkage[iWrite] = &(*ei);
			++iWrite;
		}
	}
	nRemain = nNum;
	while (nRemain > 0)
	{
		iNewLinkage = -1;
		ftMinDist = DBL_MAX;
		for (iRead = 0; iRead < nRemain; ++iRead)
		{
			ftDist = lphhLinkage[iRead]->prev()->vertex()->ftDist + lphhLinkage[iRead]->vertex()->ftDist + lphhLinkage[iRead]->ftLen;
			if (ftDist < ftMinDist)
			{
				ftMinDist = ftDist;
				iNewLinkage = iRead;
			}
		}
		hh0 = lphhLinkage[iNewLinkage];
		bLinkable = true;
		hh = hh0;
		while (hh != NULL)
		{
			if (hh->bFreeEdge)
			{
				bLinkable = false;
			}
			hh = hh->opposite()->vertex()->hhCome;
		}
		hh = hh0->vertex()->hhCome;
		while (hh != NULL)
		{
			if (hh->bFreeEdge)
			{
				bLinkable = false;
			}
			hh = hh->opposite()->vertex()->hhCome;
		}
		if (bLinkable)
		{
			hh = hh0;
			while (hh != NULL)
			{
				hh->bFreeEdge = true;
				hh->opposite()->bFreeEdge = true;
				hh = hh->opposite()->vertex()->hhCome;
			}
			hh = hh0->vertex()->hhCome;
			while (hh != NULL)
			{
				hh->bFreeEdge = true;
				hh->opposite()->bFreeEdge = true;
				hh = hh->opposite()->vertex()->hhCome;
			}
			hh = hh0;
			while (hh != NULL)
			{
				hh->iTempIndex = 0;
				hh->opposite()->iTempIndex = 0;
				if (hh->opposite()->is_border())
				{
					hh = NULL;
				}
				else
				{
					if (hh->opposite()->next()->vertex()->iTempIndex == hh0->prev()->vertex()->iTempIndex)
					{
						hh = hh->opposite()->prev();
					}
					else if (hh->opposite()->next()->vertex()->iTempIndex == hh0->vertex()->iTempIndex)
					{
						hh = hh->opposite()->next();
					}
					else
					{
						hh = NULL;
					}
				}
				if (hh == hh0)
				{
					hh = NULL;
				}
			}
			hh = hh0;
			while (hh != NULL)
			{
				if (hh->is_border())
				{
					hh = NULL;
				}
				else
				{
					hh->iTempIndex = 0;
					hh->opposite()->iTempIndex = 0;
					if (hh->next()->vertex()->iTempIndex == hh0->prev()->vertex()->iTempIndex)
					{
						hh = hh->next()->opposite();
					}
					else if (hh->next()->vertex()->iTempIndex == hh0->vertex()->iTempIndex)
					{
						hh = hh->prev()->opposite();
					}
					else
					{
						hh = NULL;
					}
				}
				if (hh == hh0)
				{
					hh = NULL;
				}
			}
		}
		else
		{
			hh0->iTempIndex = 0;
			hh0->opposite()->iTempIndex = 0;
		}
		iWrite = 0;
		for (iRead = 0; iRead < nRemain; ++iRead)
		{
			if (lphhLinkage[iRead]->iTempIndex == 1)
			{
				if (iRead > iWrite)
				{
					lphhLinkage[iWrite] = lphhLinkage[iRead];
				}
				++iWrite;
			}
		}
		nRemain = iWrite;
	}
	delete[]lphhLinkage;
}

void CMeshSplitter::m_fnMarkArea()
{
	CMeshBase::CPolyhedron::Facet_iterator fi;
	int iIndex;
	iIndex = 0;
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
	{
		fi->iTempIndex = -1;
	}
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
	{
		if (fi->iTempIndex == -1)
		{
			fi->iTempIndex = iIndex;
			m_fnExtendFacet(&(*fi));
			++iIndex;
		}
	}
}

void CMeshSplitter::m_fnMergeArea()
{
	CMeshBase::CPolyhedron::Edge_iterator ei;
	CMeshBase::CPolyhedron::Halfedge_handle *lphhLinkage;
	CMeshBase::FT ftMaxDist, ftDist;
	CMeshBase::CPolyhedron::Halfedge_handle hh0, hh;
	int nNum, nRemain, iWrite, iRead, iRedudantLinkage;
	nNum = 0;
	for (ei = m_mshSurface.edges_begin(); ei != m_mshSurface.edges_end(); ++ei)
	{
		ei->iTempIndex = 0;
		ei->opposite()->iTempIndex = 0;
		if (ei->bFreeEdge && ei->facet()->iTempIndex != ei->opposite()->facet()->iTempIndex && ei->prev()->vertex()->iTempIndex != ei->vertex()->iTempIndex)
		{
			ei->iTempIndex = 1;
			ei->opposite()->iTempIndex = 1;
			++nNum;
		}
		else
		{
			ei->iTempIndex = 0;
			ei->opposite()->iTempIndex = 0;
		}
	}
	lphhLinkage = new CMeshBase::CPolyhedron::Halfedge_handle[nNum];
	iWrite = 0;
	for (ei = m_mshSurface.edges_begin(); ei != m_mshSurface.edges_end(); ++ei)
	{
		if (ei->bFreeEdge && ei->facet()->iTempIndex != ei->opposite()->facet()->iTempIndex && ei->prev()->vertex()->iTempIndex != ei->vertex()->iTempIndex)
		{
			lphhLinkage[iWrite] = &(*ei);
			++iWrite;
		}
	}
	nRemain = nNum;
	while (nRemain > 0)
	{
		iRedudantLinkage = -1;
		ftMaxDist = 0;
		for (iRead = 0; iRead < nRemain; ++iRead)
		{
			ftDist = lphhLinkage[iRead]->prev()->vertex()->ftDist + lphhLinkage[iRead]->vertex()->ftDist + lphhLinkage[iRead]->ftLen;
			if (ftDist > ftMaxDist)
			{
				ftMaxDist = ftDist;
				iRedudantLinkage = iRead;
			}
		}
		hh0 = lphhLinkage[iRedudantLinkage];
		hh = hh0;
		while (hh != NULL)
		{
			hh->bFreeEdge = false;
			hh->opposite()->bFreeEdge = false;
			hh = hh->opposite()->vertex()->hhCome;
		}
		hh = hh0->vertex()->hhCome;
		while (hh != NULL)
		{
			hh->bFreeEdge = false;
			hh->opposite()->bFreeEdge = false;
			hh = hh->opposite()->vertex()->hhCome;
		}
		if (hh0->facet()->iTempIndex < hh0->opposite()->facet()->iTempIndex)
		{
			hh0->opposite()->facet()->iTempIndex = hh0->facet()->iTempIndex;
			m_fnExtendFacet(hh0->opposite()->facet());
		}
		else
		{
			hh0->facet()->iTempIndex = hh0->opposite()->facet()->iTempIndex;
			m_fnExtendFacet(hh0->facet());
		}
		hh0->iTempIndex = 0;
		hh0->opposite()->iTempIndex = 0;
		iWrite = 0;
		for (iRead = 0; iRead < nRemain; ++iRead)
		{
			if (lphhLinkage[iRead]->iTempIndex == 1)
			{
				if (lphhLinkage[iRead]->facet()->iTempIndex != lphhLinkage[iRead]->opposite()->facet()->iTempIndex)
				{
					if (iRead > iWrite)
					{
						lphhLinkage[iWrite] = lphhLinkage[iRead];
					}
					++iWrite;
				}
				else
				{
					lphhLinkage[iRead]->iTempIndex = 0;
					lphhLinkage[iRead]->opposite()->iTempIndex = 0;
				}
			}
		}
		nRemain = iWrite;
	}
	delete[]lphhLinkage;
}

void CMeshSplitter::m_fnSplitAlongPara()
{
	CMeshBase::CPolyhedron::Edge_iterator ei;
	for (ei = m_mshSurface.edges_begin(); ei != m_mshSurface.edges_end(); ++ei)
	{
		if (ei->is_border_edge())
		{
			ei->bFreeEdge = true;
			ei->opposite()->bFreeEdge = true;
		}
		else
		{
			if
			(
				ei->lpftGlobalPara != ei->opposite()->prev()->lpftGlobalPara ||
				ei->opposite()->lpftGlobalPara != ei->prev()->lpftGlobalPara
			)
			{
				ei->bFreeEdge = true;
				ei->opposite()->bFreeEdge = true;
			}
		}
	}
}
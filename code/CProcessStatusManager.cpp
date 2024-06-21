#include "CProcessStatusManager.h"


CProcessStatusManager::CProcessStatusManager()
{
	m_psProcessStatus = PS_EMPTY;
}

CProcessStatusManager::~CProcessStatusManager()
{
	m_fnStatusRetreat(PS_EMPTY);
}

void CProcessStatusManager::m_fnStatusRetreat(PROCESS_STATUS psNewProcessStatus)
{
	while (m_psProcessStatus > psNewProcessStatus)
	{
		switch (m_psProcessStatus)
		{
		case PS_NEW_COARSE:
			m_mshRemeshedSurface.clear();
			break;
		case PS_NEW_SCATTER:
			m_vecFacetLinkages.clear();
			break;
		case PS_NEW_LINKAGE:
			m_fnClearLinkage();
			break;
		case PS_NEW_POINTS:
			m_vecSampleBlocks.clear();
			m_vecSamplePoints.clear();
			break;
		case PS_CREVASSE_GLOBAL:
			m_fnClearGlobalPara();
			break;
		case PS_DIRECTION_UNIFORM:
			m_fnClearCrevasses();
			m_fnScatterDirection();
			break;
		case PS_DIRECTION_DECOMPOSE:
			m_fnClearFreeEdge();
			break;
		case PS_TRACK_PATH:
			m_Tracer.m_lstIP.clear();
			m_Tracer.m_lstPorts.clear();
			break;
		case PS_GEOMETRY_BLANK:
			m_mshSurface.clear();
			break;
		};
		--m_psProcessStatus;
	}
}

void CProcessStatusManager::m_fnClearLinkage()
{
	std::vector<CSamplePoint>::iterator iterSamplePoint;
	for (iterSamplePoint = m_vecSamplePoints.begin(); iterSamplePoint != m_vecSamplePoints.end(); ++iterSamplePoint)
	{
		iterSamplePoint->m_vecAdjacent.clear();
	}
}

void CProcessStatusManager::m_fnClearDensity()
{
	CMeshBase::CPolyhedron::Vertex_iterator vi;
	for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
	{
		vi->ftDensity = 0;
	}
}

void CProcessStatusManager::m_fnClearCrevasse()
{
	CMeshBase::CPolyhedron::Edge_iterator ei;
	for (ei = m_mshSurface.edges_begin(); ei != m_mshSurface.edges_end(); ++ei)
	{
		ei->iCrevasseId = -1;
		ei->opposite()->iCrevasseId = -1;
	}
	m_vecCrevasses.clear();
}
void CProcessStatusManager::m_fnScatterDirection()
{
	CMeshBase::CPolyhedron::Facet_iterator fi;
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
	{
		fi->bUnif = false;
	}
}

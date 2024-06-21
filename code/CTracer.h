#include <list>
#include "CParaConvertor.h"
//included by "CTShapeTracer.h"
#ifndef CTRACER_H
#define CTRACER_H
class CTracer
{
public:
	struct CInducedPath;
	struct CPort
	{
		CMeshBase::CPolyhedron::Halfedge_handle hhStart;
		CInducedPath *lpInducedPath;
		CMeshBase::FT ftAngle;
		int iDir;
		bool bSatisfied;
	};
	std::list<CPort> m_lstPorts;
	unsigned m_nNumFacets;
	struct CTenJoint
	{
		int iDir;
		CMeshBase::CPolyhedron::Halfedge_handle hhEdge;
		CMeshBase::FT ftProportion, ftDistant;
	};
	struct CTentacle
	{
		std::list<CTenJoint> m_lstJoints;
		CMeshBase::FT m_ftAccLen;
		CPort *m_lpSource;
		unsigned int iHeapPos;
		int iBias;
		bool bActive;
	};
	CTentacle *m_lpTentacles;
	int m_nNumTentacles;
	CTentacle **m_lplpTenHeap;
	int m_nHeapSize;
	struct CPathElement
	{
		int iDir;
		CMeshBase::CPolyhedron::Halfedge_handle hhEdge;
		CMeshBase::FT ftProportion;
		CMeshBase::FT ftDistStart, ftDistEnd, ftAngleStart, ftAngleEnd;
	};
	struct CInducedPath
	{
		CPort* lplpTerminals[2];
		bool m_bClosed;
		bool m_bSelected;
		std::list<CPathElement> m_lstData;
	};
	std::list<CInducedPath> m_lstIP;
	void m_fnSift(int iStart);
	void m_fnCheck(int iNode);
	void m_fnRemoveTentacle(int iHeapPos);
	void m_fnCountPort(CMeshBase::CPolyhedron &mshSurface);
	void m_fnAccurateSearch();
	void m_fnFillCoord(CTenJoint *lpJoint, bool bOpposite);
	void m_fnTShapeSearch();
	void m_fnBiasSearch(CMeshBase::FT ftBiasAngle);
	void m_fnDeletePath();
	void m_fnAdjustPath();
	void m_fnRenewPath();
	void m_fnMapToMesh();
};
#endif

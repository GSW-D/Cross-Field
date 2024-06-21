#include "CTracer.h"
//included by CPolyhedronModifier.h
#ifndef CTSHAPETRACER_H
#define CTSHAPETRACER_H
#include <queue>
class CTShapeTracer
{
public:
	struct CTrackKnot
	{
		CMeshBase::Point_3 ptPosition;
		int iMark;
	};
	struct CTrack
	{
		CMeshBase::FT ftAccLen;
		CMeshBase::CPolyhedron::Halfedge_handle hhEdge;
		int iDir, iSingID, iStatus, iPathID, iTermID;
		CMeshBase::FT ftProportion;
	};
	CTShapeTracer(CMeshBase::CPolyhedron& mshSurface, std::vector<CTrackKnot> &vecExtraKnots);
	CMeshBase::CPolyhedron &m_mshSurface;
	std::vector<CTrackKnot> &m_vecExtraKnots;
	std::vector<CTrack> m_vecPaths;
	void m_fnCountPort(std::vector<CTrack> &vecHeap);
	void m_fnRefillPort(std::vector<CTrack>& vecHeap);
	void m_fnTrace(std::vector<CTrack> &vecHeap);
	void m_fnTidyKnots();
	void m_fnGenerateFinSurface(CMeshBase::CPolyhedron& mshFin);
	void m_fnMarkSeparate(CMeshBase::CPolyhedron& mshFin);
	int m_fnMarkSubFacet(CMeshBase::CPolyhedron& mshFin);
	void m_fnExportSubFacet(CMeshBase::CPolyhedron& mshFin, std::ofstream &ofs);
	void m_fnExportTriangle(CMeshBase::CPolyhedron& mshFin, std::ofstream& ofs, int iPart);
	void m_fnExportOriginalMesh(std::ofstream &ofs);
	void m_fnExportExtraKnots(std::ofstream& ofs, CMeshBase::CPolyhedron& mshFin);
	void m_fnExportPaths(std::ofstream& ofs, CMeshBase::CPolyhedron& mshFin);
private:
	void m_fnSift(std::vector<CTrack> &vecHeap, int iPos);
	void m_fnCheck(std::vector<CTrack> &vecHeap, int iPos);
	int m_fnOnPath(CMeshBase::CPolyhedron::Halfedge_handle hh0, int iDir);
	void m_fnFillGridCoef(CMeshBase::CPolyhedron::Halfedge_handle hh0, int iDir, int iPathID);
	int m_fnFillIndOut(CMeshBase::CPolyhedron::Halfedge_handle hh0, int iDir, CMeshBase::FT ftProportion, int iPathId);
	void m_fnFindExtend(CMeshBase::CPolyhedron::Halfedge_handle &hh0, int &iDir);
	void m_fnDeletePath(int iPathID);
	void m_fnPartMark(int* lpiMark, int nSize, int* lpiSep, int &nParts);
	void m_fnDelaunay(CMeshBase::CPolyhedron &mshTri);
};
#endif
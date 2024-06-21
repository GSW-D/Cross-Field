#include "CDirectionAdjustor.h"
//included by CParemeterEquationFiller
#ifndef CMESHSPLITTER_H
#define CMESHSPLITTER_H
class CMeshSplitter: public CDistanceHeap
{
public:
	CMeshSplitter(CMeshBase::CPolyhedron &mshSurface);
private:
	void m_fnExtendFacet(CMeshBase::CPolyhedron::Facet_handle fh);
public:
	void m_fnActiveBoundaries();
	void m_fnPartVertices();
	void m_fnSimplifyConnection();
	void m_fnMarkArea();
	void m_fnMergeArea();
	void m_fnSplitAlongPara();
private:
	CMeshBase::CPolyhedron &m_mshSurface;

};
#endif

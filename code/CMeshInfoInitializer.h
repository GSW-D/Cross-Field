#include "CCrevasseProcessor.hpp"
//included by "CMeshInfoRepresentor.h"
#ifndef CMESHINFOINITIALIZER_H
#define CMESHINFOINITIALIZER_H
class CMeshInfoInitializer: public CCrevasseProcessor
{
public:
	CMeshBase::CPolyhedron m_mshSurface, m_mshFin;
	CMeshBase::CSamplerPolyhedron m_mshRemeshedSurface;
	bool m_bProtectGlobalPara;
	int m_nNumExpectedFacets;
	CMeshBase::FT m_ftMeanRadius, m_ftNewMeanRadius, m_ftGeometricDimension;
	CMeshBase::FT m_ftCurlFieldDimension, m_ftMaxCrossGrad, m_ftMaxElectricField, m_ftMaxCurlField;
	CTracer m_Tracer;
	CMeshBase::CPolyhedron::Halfedge_handle m_hhFirstInteger;
	CMeshBase::FT *m_lpftGlobalPara;
	int *m_lpiGlobalIndex;
	void m_fnCountNewMeshSize();
	void m_fnMarkNewMeshBorder();
	void m_fnMarkNewMeshCrevasse();
	void m_fnMarkBorder();
	void m_fnFillGeometricInfo();
	void m_fnDetectSeamType();
	void m_fnCorrectSeams();
	void m_fnMarkSingularity();
	void m_fnNormalizeMesh();
protected:
	void m_fnInitVector();
	void m_fnInitDot();
	void m_fnInitFacetNormal();
	void m_fnInitLength();
	void m_fnInitAngle();
	void m_fnInitHELocalInd();
	void m_fnInitAngleDefect();
	void m_fnInitLaplaceCoef();
	void m_fnInitPolarAxisDif();
	void m_fnInscribedCircle();
	void m_fnFacetPrincipleDirection(CMeshBase::CPolyhedron::Facet_handle fh);
	void m_fnInitMainDirection();
	void m_fnInitOtherDefaultInfo();
};
#endif

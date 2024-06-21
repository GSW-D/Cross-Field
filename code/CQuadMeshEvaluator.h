#include "CMeshOptimizer.h"
//included by "CProcessStatusManager.h" 
#ifndef CQUADMESHEVALUATOR_H
#define CQUADMESHEVALUATOR_H
class CQuadMeshEvaluator
{
public:
	CQuadMeshEvaluator(CMeshBase::CSamplerPolyhedron &mshSurface);
	~CQuadMeshEvaluator();
	double m_lpdblAngleDistribution[18], *m_lpAspectRatioDistribution, *m_lpAreaDistribution;
	void m_fnJumpSampling();
	void m_fnCountAngle();
	void m_fnCountAspectRatio(double dblUnitWidth, int nGropus);
	void m_fnCountFacetArea(double dblUnitWidth, int nGroups);
	double m_fnCountAngularDistortion();
	CMeshBase::FT m_ftAreaDistortion;
protected:
	CMeshBase::CSamplerPolyhedron &m_mshSurface;
};
#endif

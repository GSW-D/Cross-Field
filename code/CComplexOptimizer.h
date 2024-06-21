#include "CDirectionGenerator.h"
//included by CDirectionAdjustor.h
#ifndef CCOMPLEXOPTIMIZER_H
#define CCOMPLEXOPTIMIZER_H
class CComplexOptimizer
{
public:
	CComplexOptimizer(CMeshBase::CPolyhedron &mshSurface);
	void m_fnGenerateFreeMoment();
	void m_fnFillComplexMatrix(CMeshBase::CCmplxSv &CMS);
	void m_fnInitVector(CMeshBase::CCmplxSv &CMS);
	void m_fnNormalizeVector(CMeshBase::CCmplxSv &CMS);
	void m_fnFillDirection(CMeshBase::CCmplxSv &CMS);
protected:
	CMeshBase::CPolyhedron &m_mshSurface;
};
#endif

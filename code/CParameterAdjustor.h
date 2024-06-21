#include "CParameterGenerator.h"
//included by "CSampler.h"
#ifndef CPARAMETERADJUSTOR_H
#define CPARAMETERADJUSTOR_H

class CParameterAdjustor
{
private:
	void m_fnAssignCrevassePara(CCrevasseProcessor::CCrevasse *lpCrevasse, CMeshBase::FT *lpftStart);
public:
    CParameterAdjustor(CMeshBase::CPolyhedron &mshSurface, std::vector<CMeshCenter::CCrevasse> &vecCrevasses);
    ~CParameterAdjustor();

	void m_fnIntegerizeCrevasseSpan();
	void m_fnAssignCrevasseIntPara();
	void m_fnSolveHarmonicPara();
	void m_fnExtractCrevasseOffset();

protected:
	CMeshBase::CPolyhedron& m_mshSurface;
	std::vector<CMeshCenter::CCrevasse>& m_vecCrevasses;
	CMeshBase::FT *m_lpftRelativeOffset;
private:
};

#endif // CPARAMETERADJUSTOR_H

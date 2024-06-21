#include "CMeshSplitter.h"
//included by "CParameterAdjustor.h"
#ifndef CPARAMETERGENERATOR_H
#define CPARAMETERGENERATOR_H

class CParameterGenerator//: public CParameterEquationFiller
{
    public:
        CParameterGenerator(CMeshBase::CPolyhedron &mshSurface, std::vector<CMeshCenter::CCrevasse>&vecCrevasses);
		void m_fnEstimateParaCurl(CMeshBase::FT *lpftCurl);
		void m_fnMarkCrevasse();
		void m_fnExtractCrevasses();
		void m_fnUniformizeChart();
		int m_fnInitGlobalPara(CMeshBase::FT*&lpftGlobalParameter, int *&lpiGlobalIndex);
		void m_fnSolveCrevasseParaRange();
		void m_fnAssignCrevassePara();
		void m_fnRescaleCrevassePara(int nFacets);
		void m_fnSolveHarmonicPara();
		void m_fnRescaleGlobal(int nFacets);
    protected:
		CMeshBase::CPolyhedron& m_mshSurface;
		std::vector<CMeshCenter::CCrevasse>& m_vecCrevasses;
    private:
};

#endif // CPARAMETERGENERATOR_H

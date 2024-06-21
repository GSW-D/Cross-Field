#include <queue>
#include "CMeshCenter.h"
//#included by "CDirectionAdjustor.h"
#ifndef CDIRECTIONGENERATOR_H
#define CDIRECTIONGENERATOR_H

class CDirectionGenerator
{
    public:
        CDirectionGenerator(CMeshBase::CPolyhedron &mshSurface):m_mshSurface(mshSurface){};
        void m_fnMarkSharpFeature(CMeshBase::FT ftSharpBound);
        void  m_fnInitFixedDirection();
        void m_fnExtendDirection();
		void m_fnRandomizeDirection();
        void m_fnMergeSort(CMeshBase::FT* lpftData, int* lpiOrder, int nLen);
        void m_fnMarkIsolatedCones();
        void m_fnGenerateMetric();
        void m_fnFixCones(CMeshBase::FT ftThreshold);
        int m_fnFixSmallestErrorCone();
		CMeshBase::FT m_fnMIQEnergy();
		void m_fnFillComplex(CMeshBase::CCmplxSv &CmplxSv);
		void m_fnSolveEigen(CMeshBase::CCmplxSv &CmplxSv);
		void m_fnFillOptDir(CMeshBase::CCmplxSv &CmplxSv);
		void m_fnFillDirFromPara(bool bFrameField);
        void m_fnGenerateDirFromDensity();
    protected:
        CMeshBase::CPolyhedron &m_mshSurface;
    private:
        void m_fnFacetPrincipleDirection(CMeshBase::CPolyhedron::Facet_handle fh);
};

#endif // CDIRECTIONGENERATOR_H

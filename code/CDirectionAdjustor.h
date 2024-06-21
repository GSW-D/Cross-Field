#include "CComplexOptimizer.h"
//included by "CParameterEquationFiller.h"
#ifndef CDIRECTIONADJUSTOR_H
#define CDIRECTIONADJUSTOR_H
class CDirectionAdjustor: public CDistanceHeap
{
    public:
        CDirectionAdjustor(CMeshBase::CPolyhedron &msh_Surface);
        void m_fnGenerateMoment();
		void m_fnInitDirectionMatrix(CMeshBase::CLDLTSv &LDLTS);
		void m_fnFillMoment(CMeshBase::CLDLTSv &LDLTS);
		void m_fnAdjustDirection(CMeshBase::CLDLTSv &LDLTS);
		CMeshBase::FT m_fnMaxCurlFlux(CMeshBase::CPolyhedron::Facet_handle &fh);
		void m_fnAdjustCurlFlux(CMeshBase::CPolyhedron::Facet_handle fh, int iSign);
		void m_fnFindPoles(CMeshBase::FT ftChargeRadius, CMeshBase::CPolyhedron::Vertex_handle &vhMin, CMeshBase::CPolyhedron::Vertex_handle &vhMax);
		void m_fnLinkPoles(CMeshBase::CPolyhedron::Vertex_handle &vhMin, CMeshBase::CPolyhedron::Vertex_handle &vhMax);
		void m_fnInitDensityMatrix(CMeshBase::CLDLTSv &LDLTS);
		void m_fnFillPoisson(CMeshBase::CLDLTSv &LDLTS);
		void m_fnAssignDensity(CMeshBase::CLDLTSv &LDLTS);
		void m_fnModifyDensity(CMeshBase::FT ftChargeRadius, double &dblBottom, double &dblTop);
		void m_fnMoveSingularity();
        void m_fnMoveSingularity(CMeshBase::CPolyhedron::Vertex_handle vhMin, CMeshBase::CPolyhedron::Vertex_handle vhMax);
		CMeshBase::FT m_fnSmoothEnergy(CMeshBase::FT ftChargeRadius);
		CMeshBase::FT m_fnEnergy(CMeshBase::FT ftChargeRadius);
		CMeshBase::FT m_fnCurlEnergy();
		void m_fnMarkCurlPath();
		bool m_fnTestCurlPath();
		void m_fnReleaseCurl();
		void m_fnGenerateNablaField();
		void m_fnGenerateForces();
		void m_fnGenerateParaForce();
//		CMeshBase::FT m_fnMaxForce(int iType);
		int m_fnAdjustSingularitiesTheta();
		int m_fnAdjustSingularitiesPhi();
		int m_fnAdjustSingularitiesPara();
		int m_fnConstrainedDecurl();
		int m_fnConstrainedOptimize();
		int m_fnCurlAdjust();
		int m_fnMergeAdjacentSingularities();
		int m_fnNumberOfSingularities();
		CMeshBase::FT m_fnDiscreteCurl();
    protected:
        CMeshBase::CPolyhedron &m_mshSurface;
};

#endif // CDIRECTIONADJUSTOR_H

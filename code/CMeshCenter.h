#include <fstream>
#include <wx/string.h>
#include "CMeshInfoRepresentor.h"
//included by "CDirectionGenerator.h"
#ifndef CMESHCENTER_H
#define CMESHCENTER_H

class CMeshCenter : public CDistanceHeap, public CMeshInfoRepresentor
{
public:
    CMeshCenter();
	~CMeshCenter();

    void m_fnClearCrevasses();
	void m_fnClearGlobalPara();
    void m_fnClearBarrier();
	void m_fnClearFixedDirection();
	void m_fnDelaunayFlip();
	void m_fnBackupField();
	void m_fnRestoreField();
	void m_fnGenerateMoment();
	void m_fnSolveParaScale();
	void m_fnAdjustParaScale();
	CMeshBase::FT m_fnMIQEnergy();
	void m_fnEstimateAngleError(CMeshBase::FT *lpftErr);
	CMeshBase::FT m_fnAlignEnergy();
	void m_fnEvaluatePara(CMeshBase::FT &ftGrids, CMeshBase::FT &ftAreaDistortion, CMeshBase::FT &ftShapeDistortion);
	void m_fnClearFreeEdge();
	void m_fnDetectLocalPath(CMeshBase::CPolyhedron::Vertex_handle vh0, CMeshBase::CPolyhedron::Vertex_handle vh1);
	void m_fnExtendPath(CMeshBase::CPolyhedron::Halfedge_handle hhStart);
	void m_fnClearSelectedPath();
	void m_fnFreeBoundary();
	bool m_fnIsAligned();
	void m_fnFillVertexCrossField();
	void m_fnFillVertexFrameField();
	void m_fnFillVertexElectricField();
	void m_fnFillVertexNablaTheta();
	void m_fnFillVertexCurlField();
	void m_fnDetectFieldMaximum();
	void m_fnFixDirection();
	void m_fnMultiFixDirection();
	int m_nNumBorders;
	CMeshBase::FT m_ftChargeRadius;
    std::vector<CMeshBase::CPolyhedron::Vertex_handle> m_vecvhDensityExtrema;
	bool m_bExtremaSelecting;
    int m_iDensityMinimum, m_iDensityMaximum;
	bool m_bAligned;
    std::vector<CCrevasse> m_vecCrevasses;
    std::vector<CSamplePoint> m_vecSamplePoints;
    std::vector<CSampleBlock> m_vecSampleBlocks;
    std::vector<CFacetLinkage> m_vecFacetLinkages;
	int m_nGlobalParameterSize;
	bool m_bIsotropy;
	CMeshBase::CPolyhedron::Halfedge_handle m_hhDirEdge;
	CMeshBase::CPolyhedron::Facet_handle m_fhFacetFix;
protected:
private:
};



#endif // CMESHCENTER_H

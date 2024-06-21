#include "CParameterAdjustor.h"
//included by "CMeshOptimizer.h"
#ifndef CSAMPLER_H
#define CSAMPLER_H

class CSampler: public CParaConvertor
{
    public:
        CSampler(CMeshBase::CPolyhedron &mshSurface, std::vector<CMeshCenter::CCrevasse> &vecCrevasses,CMeshBase::CSamplerPolyhedron &mshRemeshedSurface, std::vector<CSamplePoint> &vecSamplePoints, std::vector<CSampleBlock> &vecSampleBlocks, std::vector<CFacetLinkage> &vecFacetLinkages);
        ~CSampler();
        void m_fnInitSampleBlocks();
        void m_fnGenerateLinkage();
        void m_fnClearLinkage();
        void m_fnGenerateFacets();
        void m_fnCorrectFacetDirection();
        void m_fnExtractNewMesh();
    protected:
        CMeshBase::CPolyhedron &m_mshSurface;
        std::vector<CMeshCenter::CCrevasse> &m_vecCrevasses;
        CMeshBase::CSamplerPolyhedron &m_mshRemeshedSurface;
        std::vector<CSamplePoint> &m_vecSamplePoints;
        std::vector<CSampleBlock> &m_vecSampleBlocks;
        std::vector<CFacetLinkage> &m_vecFacetLinkages;
    private:
        void m_fnParaToCoef(CMeshBase::FT *lpftPara, CMeshBase::FT *lpftTriangle, CMeshBase::FT *lpftCoef);
        void m_fnFacetCoef(CMeshBase::CPolyhedron::Halfedge_handle hhBase, CMeshBase::FT *lpftPara, CMeshBase::FT *lpftCoef);
        void m_fnSearchGate(CMeshBase::FT *lpftCoef, int &iMax, int &iMin);
        void m_fnParaMap(CMeshBase::FT *lpftParaOriginal, CMeshBase::CPolyhedron::Halfedge_handle hhCrevasse, CMeshBase::FT *lpftParaImage);
        void m_fnLocateFromPara(CMeshBase::CPolyhedron::Halfedge_handle hhStart, CMeshBase::FT *lpftParaOriginal, CMeshBase::FT *lpftParaImage, CMeshBase::CPolyhedron::Halfedge_handle &hhPosition, CMeshBase::FT *lpftCoef);
        void m_fnEnumIntegerParameters(CMeshBase::FT *lpftTriangle, std::vector<CIntegerParameter> &vecParameters);
        void m_fnVertexVicinityLocate(CMeshBase::CPolyhedron::Halfedge_handle hhStart, CMeshBase::FT *lpftParaOriginal, CMeshBase::FT *lpftParaImage, CMeshBase::CPolyhedron::Halfedge_handle &hhPosition, CMeshBase::FT *lpftCoef);
        void m_fnEdgeVicinityLocate(CMeshBase::CPolyhedron::Halfedge_handle hhStart, CMeshBase::FT *lpftParaOriginal, CMeshBase::FT *lpftParaImage, CMeshBase::CPolyhedron::Halfedge_handle &hhPosition, CMeshBase::FT *lpftCoef);
        void m_fnExtractSamplePoints(CSampleBlock &currentSampleBlock);
        int m_fnSearchPointFromPara(CMeshBase::CPolyhedron::Halfedge_handle hhBase, CMeshBase::FT *lpftPara);
};


#endif // CSAMPLER_H

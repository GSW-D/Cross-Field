#include "CSampler.h"
//included by "QuadMeshEvaluator.h" 
#ifndef CMESHOPTIMIZER_H
#define CMESHOPTIMIZER_H

class CMeshOptimizer
{
    public:
        CMeshOptimizer(CMeshBase::CPolyhedron &mshSurface, CMeshBase::CSamplerPolyhedron &mshRemeshedSurface);
        void m_fnMarkFreedom();
        void m_fnFillInfo();
        void m_fnGenerateAdjustment();
        void m_fnAdjustVertices();
        void m_fnProjectToOriginal();
    protected:
        CMeshBase::CPolyhedron &m_mshSurface;
        CMeshBase::CSamplerPolyhedron &m_mshRemeshedSurface;
    private:
        void m_fnFacetPointToCoef(CMeshBase::Point_3 ptObject, CMeshBase::Point_3 ptVertex0, CMeshBase::Point_3 ptVertex1, CMeshBase::Point_3 ptVertex2, CMeshBase::FT *lpftCoef, CMeshBase::FT &ftDistance, int &iNearest);
        void m_fnEdgePointToCoef(CMeshBase::Point_3 ptObject, CMeshBase::Point_3 ptVertex0, CMeshBase::Point_3 ptVertex1, CMeshBase::FT *lpftCoef, CMeshBase::FT &ftDistance, int &iNearest);
        void m_fnFacetLocalProject(CMeshBase::CSamplerPolyhedron::Vertex_handle vh);
        void m_fnEdgeLocalProject(CMeshBase::CSamplerPolyhedron::Vertex_handle vh);
};

#endif // CMESHOPTIMIZER_H

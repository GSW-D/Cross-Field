#include "CPolyhedronModifier.h"
//included by "CCrevasseProcessor.hpp"
#ifndef CDISTANCEHEAP_H
#define CDISTANCEHEAP_H

class CDistanceHeap
{
    public:
        CDistanceHeap();
    protected:
        void m_fnSift(std::vector<CMeshBase::CPolyhedron::Vertex_handle> &vecVH, int iStart);
        void m_fnCheck(std::vector<CMeshBase::CPolyhedron::Vertex_handle> &vecVH, int iNode);
    private:
};

#endif // CDISTANCEHEAP_H

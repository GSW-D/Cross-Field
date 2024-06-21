#include "CDistanceHeap.h"

CDistanceHeap::CDistanceHeap()
{
    //ctor
}

void CDistanceHeap::m_fnSift(std::vector<CMeshBase::CPolyhedron::Vertex_handle> &vecVH, int iStart)
{
    unsigned int iRoot, iLCh, iRCh, iSCh;
    CMeshBase::CPolyhedron::Vertex_handle vhTemp;
    iRoot = iStart;
    iLCh = iRoot * 2 + 1;
    iRCh = iLCh + 1;
    while(iLCh < vecVH.size())
    {
        if(iRCh < vecVH.size() && vecVH[iRCh]->ftDist < vecVH[iLCh]->ftDist)
        {
            iSCh = iRCh;
        }
        else
        {
            iSCh = iLCh;
        }
        if(vecVH[iRoot]->ftDist < vecVH[iSCh]->ftDist)
        {
            break;
        }
        else
        {
            vhTemp = vecVH[iSCh];
            vecVH[iSCh] = vecVH[iRoot];
            vecVH[iRoot] = vhTemp;
            vecVH[iSCh]->iHeapPos = iSCh;
            vecVH[iRoot]->iHeapPos = iRoot;
            iRoot = iSCh,
            iLCh = iSCh * 2 + 1;
            iRCh = iLCh + 1;
        }
    }
}

void CDistanceHeap::m_fnCheck(std::vector<CMeshBase::CPolyhedron::Vertex_handle> &vecVH, int iNode)
{
    unsigned int iChild, iRoot;
    iChild = iNode;
    CMeshBase::CPolyhedron::Vertex_handle vh;
    while(iChild != 0)
    {
        iRoot = (iChild - 1) / 2;
        if(vecVH[iChild]->ftDist < vecVH[iRoot]->ftDist)
        {
            vh = vecVH[iChild];
            vecVH[iChild] = vecVH[iRoot];
            vecVH[iRoot] = vh;
            vecVH[iChild]->iHeapPos = iChild;
            vecVH[iRoot]->iHeapPos = iRoot;
            iChild = iRoot;
        }
        else
        {
            break;
        }
    }
}

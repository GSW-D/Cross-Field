#include "CMeshOptimizer.h"

CMeshOptimizer::CMeshOptimizer(CMeshBase::CPolyhedron &mshSurface, CMeshBase::CSamplerPolyhedron &mshRemeshedSurface): m_mshSurface(mshSurface), m_mshRemeshedSurface(mshRemeshedSurface)
{
    //ctor
}

void CMeshOptimizer::m_fnMarkFreedom()
{
    CMeshBase::CSamplerPolyhedron::Vertex_iterator vi;
    CMeshBase::CSamplerPolyhedron::Halfedge_iterator hi;
    CMeshBase::Vector_3 vtVector0, vtVector1, vtCross;
    CMeshBase::CPolyhedron::Halfedge_handle hh0;
	int iTraverse;
    for(vi = m_mshRemeshedSurface.vertices_begin(); vi != m_mshRemeshedSurface.vertices_end(); ++vi)
    {
		if (vi->iBorderId == -1)
		{
			vi->iFreedom = 2;
		}
		else
		{
			if (vi->degree() == 3)
			{
				vi->iFreedom = 1;
			}
			else
			{
				vi->iFreedom = 0;
			}
		}
    }
    for(vi = m_mshRemeshedSurface.vertices_begin(); vi != m_mshRemeshedSurface.vertices_end(); ++vi)
    {
        if(vi->iFreedom == 1)
        {
            hh0 = vi->hhBase;
            if(!hh0->opposite()->is_border())
            {
                if(hh0->next()->opposite()->is_border())
                {
                    vi->hhBase = hh0->next();
                }
                else if(hh0->prev()->opposite()->is_border())
                {
                    vi->hhBase = hh0->prev();
                }
                else
                {
					iTraverse = 0;
                    while(hh0->vertex()->iBordId == -1 && iTraverse < 3)
                    {
						++iTraverse;
                        hh0 = hh0->next();
                    }
					iTraverse = 0;
                    while(!hh0->opposite()->is_border() && iTraverse < 3)
                    {
						++iTraverse;
                        hh0 = hh0->next()->opposite();
                    }
                    vi->hhBase = hh0;
                }
            }
        }
    }
}

void CMeshOptimizer::m_fnFillInfo()
{
    CMeshBase::CSamplerPolyhedron::Edge_iterator ei;
    for(ei = m_mshRemeshedSurface.edges_begin(); ei != m_mshRemeshedSurface.edges_end(); ++ei)
    {
        ei->vtVector = ei->vertex()->point() - ei->prev()->vertex()->point();
        ei->ftSqLen = (ei->vtVector) * (ei->vtVector);
        ei->ftLen = sqrt(ei->ftSqLen);
        ei->opposite()->vtVector = -ei->vtVector;
        ei->opposite()->ftSqLen = ei->ftSqLen;
        ei->opposite()->ftLen = ei->ftLen;
    }
}

void CMeshOptimizer::m_fnGenerateAdjustment()
{
    CMeshBase::CSamplerPolyhedron::Vertex_iterator vi;
    CMeshBase::CSamplerPolyhedron::Halfedge_handle hh0, hh;
    CMeshBase::FT ftTwArea;
    CMeshBase::Vector_3 vtNorm;
    for(vi = m_mshRemeshedSurface.vertices_begin(); vi != m_mshRemeshedSurface.vertices_end(); ++vi)
    {
        if(vi->iFreedom != 0)
        {
            vi->vtAdjustment = CMeshBase::Vector_3(0, 0, 0);
            hh0 = vi->halfedge();
            hh = hh0;
            do
            {
                if(!hh->is_border())
                {
                    vtNorm = CGAL::cross_product(hh->vtVector, hh->next()->vtVector);
                    ftTwArea = sqrt(vtNorm * vtNorm);
                    vi->vtAdjustment = vi->vtAdjustment + (hh->next()->vtVector - hh->vtVector) * ((hh->ftSqLen + hh->next()->ftSqLen) / (32 * ftTwArea));
                }
                hh = hh->next()->opposite();
            }while(hh != hh0);
        }
    }
}

void CMeshOptimizer::m_fnAdjustVertices()
{
    CMeshBase::CSamplerPolyhedron::Vertex_iterator vi;
    for(vi = m_mshRemeshedSurface.vertices_begin(); vi != m_mshRemeshedSurface.vertices_end(); ++vi)
    {
        if(vi->iFreedom != 0)
        {
            vi->point() = vi->point() + vi->vtAdjustment;
        }
    }
}

void CMeshOptimizer::m_fnFacetPointToCoef(CMeshBase::Point_3 ptObject, CMeshBase::Point_3 ptVertex0, CMeshBase::Point_3 ptVertex1, CMeshBase::Point_3 ptVertex2, CMeshBase::FT *lpftCoef, CMeshBase::FT &ftDistance, int &iNearest)
{
    CMeshBase::Vector_3 lpvtVector[3], lpvtNorm[3], vtNorm, vtEdge;
    CMeshBase::FT ftCoefSum;
    int i, iMax, iMin, iMiddle;
    lpvtVector[0] = ptVertex0 - ptObject;
    lpvtVector[1] = ptVertex1 - ptObject;
    lpvtVector[2] = ptVertex2 - ptObject;
    for(i = 0; i < 3; ++i)
    {
        lpvtNorm[i] = CGAL::cross_product(lpvtVector[(i + 1) % 3], lpvtVector[(i + 2) % 3]);
    }
    vtNorm = lpvtNorm[0] + lpvtNorm[1] + lpvtNorm[2];
    for(i = 0; i < 3; ++i)
    {
        lpftCoef[i] = lpvtNorm[i] * vtNorm;
    }
    ftCoefSum = lpftCoef[0] + lpftCoef[1] +lpftCoef[2];
    for(i = 0; i < 3; ++i)
    {
        lpftCoef[i] /= ftCoefSum;
    }
    iMax = 0;
    iMin = 0;
    for(i = 1; i < 3; ++i)
    {
        if(lpftCoef[i] < lpftCoef[iMin])
        {
            iMin = i;
        }
        if(lpftCoef[i] > lpftCoef[iMax])
        {
            iMax = i;
        }
    }
    if(iMin == iMax)
    {
        iMin = 0;
        iMax = 2;
    }
    iMiddle = 3 - iMin - iMax;
    if(lpftCoef[iMin] < 0)
    {
        if(lpftCoef[iMiddle] < 0)
        {
            ftDistance = sqrt(lpvtVector[iMax] * lpvtVector[iMax]);
            iNearest = -(iMax + 1);
        }
        else
        {
            vtEdge = lpvtVector[iMax] - lpvtVector[iMiddle];
            ftDistance = sqrt((lpvtNorm[iMin] * lpvtNorm[iMin]) / (vtEdge * vtEdge));
            iNearest = iMin + 1;
        }
    }
    else
    {
        ftDistance = fabs(lpvtVector[iMax] * vtNorm / sqrt(vtNorm * vtNorm));
        iNearest = 0;
    }
    if(iNearest > 0)
    {
        lpftCoef[iMin] = 0;
        lpftCoef[iMiddle] = lpvtVector[iMax] * (lpvtVector[iMax] - lpvtVector[iMiddle]);//sqrt(fabs(lpvtVector[iMax] * lpvtVector[iMax] - ftDistance * ftDistance));
        lpftCoef[iMax] = lpvtVector[iMiddle] * (lpvtVector[iMiddle] - lpvtVector[iMax]);//sqrt(fabs(lpvtVector[iMiddle] * lpvtVector[iMiddle] - ftDistance * ftDistance));
        ftCoefSum = lpftCoef[iMiddle] + lpftCoef[iMax];
        lpftCoef[iMiddle] /= ftCoefSum;
        lpftCoef[iMax] /= ftCoefSum;
    }
    else if(iNearest < 0)
    {
        lpftCoef[iMax] = 1;
        lpftCoef[iMiddle] = 0;
        lpftCoef[iMin] = 0;
    }
}

void CMeshOptimizer::m_fnEdgePointToCoef(CMeshBase::Point_3 ptObject, CMeshBase::Point_3 ptVertex0, CMeshBase::Point_3 ptVertex1, CMeshBase::FT *lpftCoef, CMeshBase::FT &ftDistance, int &iNearest)
{
    CMeshBase::Vector_3 vtVector0, vtVector1, vtVector, vtCross;
    CMeshBase::FT ftSumDot;
    vtVector0 = ptObject - ptVertex0;
    vtVector1 = ptVertex1 - ptObject;
    vtVector = ptVertex1 - ptVertex0;
    lpftCoef[0] = vtVector1 * vtVector;
    lpftCoef[1] = vtVector0 * vtVector;
    ftSumDot = lpftCoef[0] + lpftCoef[1];
    lpftCoef[0] /= ftSumDot;
    lpftCoef[1] /= ftSumDot;
    if(lpftCoef[0] < 0)
    {
        ftDistance = sqrt(vtVector1 * vtVector1);
        iNearest = 2;
    }
    else if(lpftCoef[1] < 0)
    {
        ftDistance = sqrt(vtVector0 * vtVector0);
        iNearest = 1;
    }
    else
    {
        vtCross = CGAL::cross_product(vtVector0, vtVector1);
        ftDistance = sqrt((vtCross * vtCross) / (vtVector * vtVector));
        iNearest = 0;
    }
    if(iNearest == 1)
    {
        lpftCoef[0] = 1;
        lpftCoef[1] = 0;
    }
    else if(iNearest == 2)
    {
        lpftCoef[0] = 0;
        lpftCoef[1] = 1;
    }
}

void CMeshOptimizer::m_fnFacetLocalProject(CMeshBase::CSamplerPolyhedron::Vertex_handle vh)
{
    CMeshBase::FT lpftCoef[3], ftMinDistance, lpftCoefTemp[3], ftDistance;
    CMeshBase::CPolyhedron::Halfedge_handle hh, hhTest, hhGate;
    int iNearest;
    m_fnFacetPointToCoef(vh->point(), vh->hhBase->vertex()->point(), vh->hhBase->next()->vertex()->point(), vh->hhBase->prev()->vertex()->point(), lpftCoef, ftMinDistance, iNearest);
    hh = vh->hhBase;
    hhGate = NULL;
    do
    {
        hhTest = hh->opposite();
        if(!hhTest->is_border())
        {
            m_fnFacetPointToCoef(vh->point(), hhTest->vertex()->point(), hhTest->next()->vertex()->point(), hhTest->prev()->vertex()->point(), lpftCoefTemp, ftDistance, iNearest);
            if(ftDistance < ftMinDistance)
            {
                hhGate = hh;
                memcpy(lpftCoef, lpftCoefTemp, 3 * sizeof(CMeshBase::FT));
                ftMinDistance = ftDistance;
            }
        }
        hh = hh->next();
    }while(hh != vh->hhBase);
    while(hhGate != NULL)
    {
        vh->hhBase = hhGate->opposite();
        hhGate = NULL;
        hh = vh->hhBase->next();
        while(hh != vh->hhBase)
        {
            hhTest = hh->opposite();
            if(!hhTest->is_border())
            {
                m_fnFacetPointToCoef(vh->point(), hhTest->vertex()->point(), hhTest->next()->vertex()->point(), hhTest->prev()->vertex()->point(), lpftCoefTemp, ftDistance, iNearest);
                if(ftDistance < ftMinDistance)
                {
                    hhGate = hh;
                    memcpy(lpftCoef, lpftCoefTemp, 3 * sizeof(CMeshBase::FT));
                    ftMinDistance = ftDistance;
                }
            }
            hh = hh->next();
        }
    }
    vh->point() = CGAL::barycenter(vh->hhBase->vertex()->point(), lpftCoef[0], vh->hhBase->next()->vertex()->point(), lpftCoef[1], vh->hhBase->prev()->vertex()->point(), lpftCoef[2]);
}

void CMeshOptimizer::m_fnEdgeLocalProject(CMeshBase::CSamplerPolyhedron::Vertex_handle vh)
{
    CMeshBase::FT lpftCoef[2], ftMinDistance, lpftCoefTemp[2], ftDistance;
    CMeshBase::CPolyhedron::Halfedge_handle hhTest;
    int iNearest, iDirection;
    bool bChanged;
    m_fnEdgePointToCoef(vh->point(), vh->hhBase->prev()->vertex()->point(), vh->hhBase->vertex()->point(), lpftCoef, ftMinDistance, iNearest);
    iDirection = 0;
    hhTest = vh->hhBase->opposite()->prev()->opposite();
    m_fnEdgePointToCoef(vh->point(), hhTest->prev()->vertex()->point(), hhTest->vertex()->point(), lpftCoefTemp, ftDistance, iNearest);
    if(ftDistance < ftMinDistance)
    {
        vh->hhBase = hhTest;
        memcpy(lpftCoef, lpftCoefTemp, 2 * sizeof(CMeshBase::FT));
        ftMinDistance = ftDistance;
        iDirection = 1;
    }
    hhTest = vh->hhBase->opposite()->next()->opposite();
    m_fnEdgePointToCoef(vh->point(), hhTest->prev()->vertex()->point(), hhTest->vertex()->point(), lpftCoefTemp, ftDistance, iNearest);
    if(ftDistance < ftMinDistance)
    {
        vh->hhBase = hhTest;
        memcpy(lpftCoef, lpftCoefTemp, 2 * sizeof(CMeshBase::FT));
        ftMinDistance = ftDistance;
        iDirection = -1;
    }
    if(iDirection != 0)
    {
        do
        {
            bChanged = false;
            if(iDirection > 0)
            {
                hhTest = vh->hhBase->opposite()->prev()->opposite();
            }
            else
            {
                hhTest = vh->hhBase->opposite()->next()->opposite();
            }
            m_fnEdgePointToCoef(vh->point(), hhTest->prev()->vertex()->point(), hhTest->vertex()->point(), lpftCoefTemp, ftDistance, iNearest);
            if(ftDistance < ftMinDistance)
            {
                vh->hhBase = hhTest;
                memcpy(lpftCoef, lpftCoefTemp, 2 * sizeof(CMeshBase::FT));
                ftMinDistance = ftDistance;
                bChanged = true;
            }
        }while(bChanged);
    }
    vh->point() = CGAL::barycenter(vh->hhBase->prev()->vertex()->point(), lpftCoef[0], vh->hhBase->vertex()->point(), lpftCoef[1]);
}

void CMeshOptimizer::m_fnProjectToOriginal()
{
    CMeshBase::CSamplerPolyhedron::Vertex_iterator vi;
    for(vi = m_mshRemeshedSurface.vertices_begin(); vi != m_mshRemeshedSurface.vertices_end(); ++vi)
    {
        switch(vi->iFreedom)
        {
        case 1:
            m_fnEdgeLocalProject(&(*vi));
            break;
        case 2:
            m_fnFacetLocalProject(&(*vi));
            break;
        }
    }
}


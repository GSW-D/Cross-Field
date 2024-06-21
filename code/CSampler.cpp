#include "CSampler.h"

CSampler::CSampler(CMeshBase::CPolyhedron &mshSurface, std::vector<CMeshCenter::CCrevasse> &vecCrevasses, CMeshBase::CSamplerPolyhedron &mshRemeshedSurface, std::vector<CSamplePoint> &vecSamplePoints, std::vector<CSampleBlock> &vecSampleBlocks, std::vector<CFacetLinkage> &vecFacetLinkages):
    m_mshSurface(mshSurface),
    m_vecCrevasses(vecCrevasses),
    m_mshRemeshedSurface(mshRemeshedSurface),
    m_vecSamplePoints(vecSamplePoints),
    m_vecSampleBlocks(vecSampleBlocks),
    m_vecFacetLinkages(vecFacetLinkages)
{
    //ctor
}

CSampler::~CSampler()
{
    //dtor
}

void CSampler::m_fnParaToCoef(CMeshBase::FT *lpftPara, CMeshBase::FT *lpftTriangle, CMeshBase::FT *lpftCoef)
{
    int i, j, l;
    CMeshBase::FT lpftParaBuffer[24], lpftCrossProduct[3], ftTemp, lpftDot[4];
    for(i = 0; i < 4; ++i)
    {
        memcpy(lpftParaBuffer + i * 6, lpftTriangle, 6 * sizeof(CMeshBase::FT));
    }
    for(i = 0; i < 3; ++i)
    {
        memcpy(lpftParaBuffer + i * 8, lpftPara, 2 * sizeof(CMeshBase::FT));
    }
    for(i = 0; i < 4; ++i)
    {
        for(j = 0; j < 3; ++j)
        {
            lpftCrossProduct[j] = lpftParaBuffer[i * 6 + ((j + 1) % 3) * 2] * lpftParaBuffer[i * 6 + ((j + 2) % 3) * 2 + 1] - lpftParaBuffer[i * 6 + ((j + 2) % 3) * 2] * lpftParaBuffer[i * 6 + ((j + 1) % 3) * 2 + 1];
        }
        for(j = 0; j < 2; ++j)
        {
            for(l = j + 1; l < 3; ++l)
            {
                if(fabs(lpftCrossProduct[l]) < fabs(lpftCrossProduct[j]))
                {
                    ftTemp = lpftCrossProduct[j];
                    lpftCrossProduct[j] = lpftCrossProduct[l];
                    lpftCrossProduct[l] = ftTemp;
                }
            }
        }
        lpftDot[i] = lpftCrossProduct[0] + lpftCrossProduct[1] + lpftCrossProduct[2];
    }
    for(i = 0; i < 3; ++i)
    {
        lpftCoef[i] = lpftDot[i] / lpftDot[3];
    }
}

void CSampler::m_fnFacetCoef(CMeshBase::CPolyhedron::Halfedge_handle hhBase, CMeshBase::FT *lpftPara, CMeshBase::FT *lpftCoef)
{
    CMeshBase::FT lpftTriangle[6];
//    CMeshBase::CPolyhedron::Halfedge_handle hhCrevasse;
    memcpy(lpftTriangle, hhBase->lpftGlobalPara, 2 * sizeof(CMeshBase::FT));
    memcpy(lpftTriangle + 2, hhBase->next()->lpftGlobalPara, 2 * sizeof(CMeshBase::FT));
    memcpy(lpftTriangle + 4, hhBase->prev()->lpftGlobalPara, 2 * sizeof(CMeshBase::FT));
    m_fnParaToCoef(lpftPara, lpftTriangle, lpftCoef);
}

void CSampler::m_fnSearchGate(CMeshBase::FT *lpftCoef, int &iMax, int &iMin)//, CMeshBase::CPolyhedron::Halfedge_handle hhBase)
{
    int i;
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
}

void CSampler::m_fnParaMap(CMeshBase::FT *lpftParaOriginal, CMeshBase::CPolyhedron::Halfedge_handle hhCrevasse, CMeshBase::FT *lpftParaImage)
{
    CMeshCenter::CCrevasse *lpCrevasse;
    if(hhCrevasse == NULL || !hhCrevasse->bFreeEdge)
    {
		memcpy(lpftParaImage, lpftParaOriginal, 2 * sizeof(CMeshBase::FT)); 
    }
    else
    {
        lpCrevasse = &(m_vecCrevasses[hhCrevasse->iCrevasseId]);
		lpCrevasse->parameter_convert(lpftParaOriginal, lpftParaImage);
    }
}


void CSampler::m_fnVertexVicinityLocate(CMeshBase::CPolyhedron::Halfedge_handle hhStart, CMeshBase::FT *lpftParaOriginal, CMeshBase::FT *lpftParaImage, CMeshBase::CPolyhedron::Halfedge_handle &hhPosition, CMeshBase::FT *lpftCoef)
{
    CMeshBase::FT lpftCoefCurrent[3], lpftParaCurrent[2], ftCurrentMinArea, ftMinArea;
    CMeshBase::CPolyhedron::Halfedge_handle hh0, hh;
    int iMax, iMin;

    memcpy(lpftParaCurrent, lpftParaOriginal, 2 * sizeof(CMeshBase::FT));

    if(hhStart->vertex()->iBordId == -1)
    {
        {
            hh0 = hhStart;
        }
    }
    else
    {
        hh0 = hhStart->hhCrevasseNext->prev();
        while(!hh0->next()->opposite()->is_border())
        {
            m_fnParaMap(lpftParaCurrent, hh0->next(), lpftParaCurrent);
            hh0 = hh0->next()->opposite()->hhCrevasseNext->prev();
        }
    }

    m_fnFacetCoef(hh0, lpftParaCurrent, lpftCoefCurrent);
    m_fnSearchGate(lpftCoefCurrent, iMax, iMin);
    ftCurrentMinArea = lpftCoefCurrent[iMin] * hh0->facet()->ftTwArea;

    ftMinArea = ftCurrentMinArea;
    hhPosition = hh0;
    memcpy(lpftParaImage, lpftParaCurrent, 2 * sizeof(CMeshBase::FT));
    memcpy(lpftCoef, lpftCoefCurrent, 3 * sizeof(CMeshBase::FT));

    if(hh0->bFreeEdge)
    {
        m_fnParaMap(lpftParaCurrent, hh0, lpftParaCurrent);
    }
    hh = hh0->opposite()->prev();
    while(hh!= hh0 && !hh->is_border())
    {
        m_fnFacetCoef(hh, lpftParaCurrent, lpftCoefCurrent);
        m_fnSearchGate(lpftCoefCurrent, iMax, iMin);
        ftCurrentMinArea = lpftCoefCurrent[iMin] * hh->facet()->ftTwArea;
        if(ftMinArea < ftCurrentMinArea || (ftMinArea == ftCurrentMinArea && hh->facet()->iIndex < hhPosition->facet()->iIndex))
        {
            ftMinArea = ftCurrentMinArea;
            hhPosition = hh;
            memcpy(lpftParaImage, lpftParaCurrent, 2 * sizeof(CMeshBase::FT));
            memcpy(lpftCoef, lpftCoefCurrent, 3 * sizeof(CMeshBase::FT));
        }
        if(hh->bFreeEdge && !hh->opposite()->is_border())//(hh->bFreeEdge && !hh->opposite()->is_border() && !hh->opposite()->facet()->bBearingSingularity)
        {
            m_fnParaMap(lpftParaCurrent, hh, lpftParaCurrent);
        }
        hh = hh->opposite()->prev();
    }
}

void CSampler::m_fnEdgeVicinityLocate(CMeshBase::CPolyhedron::Halfedge_handle hhStart, CMeshBase::FT *lpftParaOriginal, CMeshBase::FT *lpftParaImage, CMeshBase::CPolyhedron::Halfedge_handle &hhPosition, CMeshBase::FT *lpftCoef)
{
    CMeshBase::FT lpftCoefCurrent[3], lpftParaCurrent[2], ftMinArea, ftCurrentMinArea;

    m_fnFacetCoef(hhStart, lpftParaOriginal, lpftCoefCurrent);
    ftCurrentMinArea = lpftCoefCurrent[1] * hhStart->facet()->ftTwArea;

    hhPosition = hhStart;
    memcpy(lpftParaImage, lpftParaOriginal, 2 * sizeof(CMeshBase::FT));
    memcpy(lpftCoef, lpftCoefCurrent, 3 * sizeof(CMeshBase::FT));
    ftMinArea = ftCurrentMinArea;

	if (!hhStart->opposite()->is_border())
	{
		if (hhStart->bFreeEdge)
		{
			m_fnParaMap(lpftParaOriginal, hhStart, lpftParaCurrent);
		}
		else
		{
			memcpy(lpftParaCurrent, lpftParaOriginal, 2 * sizeof(CMeshBase::FT));
		}
		m_fnFacetCoef(hhStart->opposite(), lpftParaCurrent, lpftCoefCurrent);
		ftCurrentMinArea = lpftCoefCurrent[1] * hhStart->opposite()->facet()->ftTwArea;
		if(ftMinArea < ftCurrentMinArea || (ftMinArea == ftCurrentMinArea && hhStart->opposite()->facet()->iIndex < hhStart->facet()->iIndex))
		{
			hhPosition = hhStart->opposite();
			memcpy(lpftParaImage, lpftParaCurrent, 2 * sizeof(CMeshBase::FT));
			memcpy(lpftCoef, lpftCoefCurrent, 3 * sizeof(CMeshBase::FT));
			ftMinArea = ftCurrentMinArea;
		}
	}
}

void CSampler::m_fnLocateFromPara(CMeshBase::CPolyhedron::Halfedge_handle hhStart, CMeshBase::FT *lpftParaOriginal, CMeshBase::FT *lpftParaImage, CMeshBase::CPolyhedron::Halfedge_handle &hhPosition, CMeshBase::FT *lpftCoef)
{
    int iMin, iMax;
    CMeshBase::CPolyhedron::Halfedge_handle hhGate, hh0;
    bool bFirstTime;
    hhGate = NULL;
    hh0 = NULL;
    memcpy(lpftParaImage, lpftParaOriginal, 2 *sizeof(CMeshBase::FT));
    hhPosition = hhStart;
    bFirstTime = true;
    do
    {
        if(hhPosition->is_border())
        {
            hhPosition = NULL;
            break;
        }
        m_fnFacetCoef(hhPosition, lpftParaImage, lpftCoef);
        m_fnSearchGate(lpftCoef, iMax, iMin);
        if(lpftCoef[iMin] < -0.00390625)
        {
            switch(iMin)
            {
            case 0:
                {
                    hhGate = hhPosition->prev();
                }
                break;
            case 1:
                {
                    hhGate = hhPosition;
                }
                break;
            case 2:
                {
                    hhGate = hhPosition->next();
                }
                break;
            }
            if(bFirstTime)
            {
                bFirstTime = false;
            }
            else if(iMin == 1)
            {
                hhPosition = NULL;
                break;
            }
            if(hhGate->bFreeEdge && !hhGate->opposite()->is_border())
            {
                m_fnParaMap(lpftParaImage, hhGate, lpftParaImage);
            }
            hhPosition = hhGate->opposite();
        }
    }while(lpftCoef[iMin] < -0.00390625);
    if(hhPosition != NULL)
    {
        if(lpftCoef[iMax] > 0.99609375)
        {
            switch(iMax)
            {
            case 0:
                hh0 = hhPosition;
                break;
            case 1:
                hh0 = hhPosition->next();
                break;
            case 2:
                hh0 = hhPosition->prev();
                break;
            }
            m_fnVertexVicinityLocate(hh0, lpftParaImage, lpftParaImage, hhPosition, lpftCoef);
        }
        else if(lpftCoef[iMin] < 0.00390625)
        {
            switch(iMin)
            {
            case 0:
                hh0 = hhPosition->prev();
                break;
            case 1:
                hh0 = hhPosition;
                break;
            case 2:
                hh0 = hhPosition->next();
                break;
            }
            if(!hh0->opposite()->is_border())
            {
                m_fnEdgeVicinityLocate(hh0, lpftParaImage, lpftParaImage, hhPosition, lpftCoef);
            }
        }
    }
}

void CSampler::m_fnEnumIntegerParameters(CMeshBase::FT *lpftTriangle, std::vector<CIntegerParameter> &vecParameters)
{
    CMeshBase::FT lpftMax[2], lpftMin[2], lpftCoef[3];
    CIntegerParameter IntegerPara;
    int i;
    CMeshBase::FT ftU, ftV;
    vecParameters.clear();
    for(i = 0; i < 2; ++i)
    {
        lpftMax[i] = lpftTriangle[i];
        lpftMin[i] = lpftTriangle[i];
        if(lpftTriangle[i + 2] > lpftMax[i])
        {
            lpftMax[i] = lpftTriangle[i + 2];
        }
        if(lpftTriangle[i + 2] < lpftMin[i])
        {
            lpftMin[i] = lpftTriangle[i + 2];
        }
        if(lpftTriangle[i + 4] > lpftMax[i])
        {
            lpftMax[i] = lpftTriangle[i + 4];
        }
        if(lpftTriangle[i + 4] < lpftMin[i])
        {
            lpftMin[i] = lpftTriangle[i + 4];
        }
    }
    for(i = 0; i < 2; ++i)
    {
		lpftMin[i] = floor(lpftMin[i]);//floor(lpftMin[i] * 2);
		lpftMax[i] = ceil(lpftMax[i]) + 0.5;//ceil(lpftMax[i] * 2) + 0.5;
    }
    ftU = lpftMin[0];
    while(ftU < lpftMax[0])
    {
        ftV = lpftMin[1];
        while(ftV < lpftMax[1])
        {
			IntegerPara.lpftPara[0] = ftU;//ftU / 2;
			IntegerPara.lpftPara[1] = ftV;//ftV / 2;
            m_fnParaToCoef(IntegerPara.lpftPara, lpftTriangle, lpftCoef);
            if(lpftCoef[0] > -0.00390625 && lpftCoef[1] > -0.00390625 && lpftCoef[2] > -0.00390625)
            {
                vecParameters.push_back(IntegerPara);
            }
            ftV += 1.0;
        }
        ftU += 1.0;
    }
}

void CSampler::m_fnExtractSamplePoints(CSampleBlock &currentSampleBlock)
{
    CMeshBase::FT lpftTriangle[6], lpftParaTemp[2], lpftCoef[3];
    CMeshBase::CPolyhedron::Halfedge_handle hh0, hhTemp;
    CSamplePoint newSamplePoint;
    std::vector<CIntegerParameter> vecParaSamples;
    std::vector<CIntegerParameter>::iterator iterParaVisitor;
    std::vector<CSamplePoint>::iterator iterSampleVisitor;
    int iIndex;
    iIndex = currentSampleBlock.iPointStartIndex;
    hh0 = currentSampleBlock.hhBase;
    memcpy(lpftTriangle, hh0->lpftGlobalPara, 2 * sizeof(CMeshBase::FT));
    memcpy(lpftTriangle + 2, hh0->next()->lpftGlobalPara, 2 * sizeof(CMeshBase::FT));
    memcpy(lpftTriangle + 4, hh0->prev()->lpftGlobalPara, 2 * sizeof(CMeshBase::FT));
    m_fnEnumIntegerParameters(lpftTriangle, vecParaSamples);
    if(!vecParaSamples.empty())
    {
        for(iterParaVisitor = vecParaSamples.begin(); iterParaVisitor != vecParaSamples.end(); ++iterParaVisitor)
        {
            m_fnLocateFromPara(hh0, iterParaVisitor->lpftPara, lpftParaTemp, hhTemp, lpftCoef);
            if(hhTemp != NULL && hhTemp->facet() == hh0->facet())
            {
                newSamplePoint.iBlockId = currentSampleBlock.hhBase->facet()->iIndex;
                newSamplePoint.iIndex = iIndex;
                memcpy(newSamplePoint.lpftPara, iterParaVisitor->lpftPara, 2 * sizeof(CMeshBase::FT));
                newSamplePoint.ptCoordinate = CGAL::barycenter(hhTemp->vertex()->point(), lpftCoef[0], hhTemp->next()->vertex()->point(), lpftCoef[1], hhTemp->prev()->vertex()->point(), lpftCoef[2]);
                m_vecSamplePoints.push_back(newSamplePoint);
                ++iIndex;
            }
        }
        vecParaSamples.clear();
    }
    currentSampleBlock.iPointEndIndex = iIndex;
}

void CSampler::m_fnInitSampleBlocks()
{
    int iIndex;
    CMeshBase::CPolyhedron::Facet_iterator fi;
    CSampleBlock newSampleBlock;
    iIndex = 0;
    for(fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
    {
        fi->iIndex = iIndex;
        ++iIndex;
    }
    m_mshRemeshedSurface.clear();
    m_vecFacetLinkages.clear();
    m_vecSamplePoints.clear();
    m_vecSampleBlocks.clear();
    iIndex = 0;
    for(fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
    {
        newSampleBlock.hhBase = fi->halfedge();
        newSampleBlock.iPointStartIndex = iIndex;
        m_fnExtractSamplePoints(newSampleBlock);
        m_vecSampleBlocks.push_back(newSampleBlock);
        iIndex = newSampleBlock.iPointEndIndex;
    }
}
int CSampler::m_fnSearchPointFromPara(CMeshBase::CPolyhedron::Halfedge_handle hhBase, CMeshBase::FT *lpftPara)
{
    int iBlockId, i;
    CMeshBase::FT *lpftCurrentPara;
    CSampleBlock *lpSampleBlock;
    iBlockId = hhBase->facet()->iIndex;
    lpSampleBlock = &(m_vecSampleBlocks[iBlockId]);
    if(lpSampleBlock->iPointStartIndex == lpSampleBlock->iPointEndIndex)
    {
        i = -1;
    }
    else
    {
        for(i = lpSampleBlock->iPointStartIndex; i != lpSampleBlock->iPointEndIndex; ++i)
        {
            lpftCurrentPara = m_vecSamplePoints[i].lpftPara;
            if(fabs(lpftCurrentPara[0] - lpftPara[0]) < 0.25 && fabs(lpftCurrentPara[1] - lpftPara[1]) < 0.25)
            {
                break;
            }
        }
    }
    if(i == lpSampleBlock->iPointEndIndex)
    {
        i = -1;
    }
    return i;
}



void CSampler::m_fnGenerateLinkage()
{
    std::vector<CSamplePoint>::iterator iterSamplePoint;
    std::vector<int> vecAdjacentBlockId;
    std::vector<int> vecAdjacentVertexId;
    CMeshBase::FT lpftAdjacentPara[2], lpftParaImage[2], lpftCoef[4];
//    CMeshBase::CPolyhedron::Facet_handle fhPosition;
    CMeshBase::CPolyhedron::Halfedge_handle hhPosition, hhBase;
    int i, j;
    for(iterSamplePoint = m_vecSamplePoints.begin(); iterSamplePoint != m_vecSamplePoints.end(); ++iterSamplePoint)
    {
        hhBase = m_vecSampleBlocks[iterSamplePoint->iBlockId].hhBase;
        for(i = 0; i < 4; ++i)
        {
			lpftAdjacentPara[0] = iterSamplePoint->lpftPara[0] + CMeshBase::FT(((i + 1) % 2) * ((((i + 2) % 4) / 2) * 2 - 1));// *0.5;
			lpftAdjacentPara[1] = iterSamplePoint->lpftPara[1] + CMeshBase::FT((i % 2) * ((((i + 2) % 4) / 2) * 2 - 1));// *0.5;
            m_fnLocateFromPara(hhBase, lpftAdjacentPara, lpftParaImage, hhPosition, lpftCoef);
            if(hhPosition != NULL)
            {
                j = m_fnSearchPointFromPara(hhPosition, lpftParaImage);
                if(j != -1)
                {
                    if(std::find(iterSamplePoint->m_vecAdjacent.begin(), iterSamplePoint->m_vecAdjacent.end(), j) == iterSamplePoint->m_vecAdjacent.end())
                    {
                        iterSamplePoint->m_vecAdjacent.push_back(j);
                    }
                    if(std::find(m_vecSamplePoints[j].m_vecAdjacent.begin(), m_vecSamplePoints[j].m_vecAdjacent.end(), iterSamplePoint->iIndex) == m_vecSamplePoints[j].m_vecAdjacent.end())
                    {
                        m_vecSamplePoints[j].m_vecAdjacent.push_back(iterSamplePoint->iIndex);
                    }
                }
            }
        }
    }
}

void CSampler::m_fnGenerateFacets()
{
    std::vector<CSamplePoint>::iterator iterSamplePoint;
    std::vector<int>::iterator iterAdjacent1, iterAdjacent2, iterAdjacent3, iterAdjacentEnd;
    CFacetLinkage newFacet;
    m_vecFacetLinkages.clear();
    for(iterSamplePoint = m_vecSamplePoints.begin(); iterSamplePoint != m_vecSamplePoints.end(); ++iterSamplePoint)
    {
		if (!iterSamplePoint->m_vecAdjacent.empty())
		{
			iterAdjacentEnd = iterSamplePoint->m_vecAdjacent.end();
			--iterAdjacentEnd;
			for(iterAdjacent1 = iterSamplePoint->m_vecAdjacent.begin(); iterAdjacent1 != iterAdjacentEnd; ++iterAdjacent1)
			{
				if(*iterAdjacent1 > iterSamplePoint->iIndex)
				{
					iterAdjacent2 = iterAdjacent1;
					++iterAdjacent2;
					while(iterAdjacent2 != iterSamplePoint->m_vecAdjacent.end())
					{
						if(*iterAdjacent2 > iterSamplePoint->iIndex)
						{
							for(iterAdjacent3 = m_vecSamplePoints[*iterAdjacent1].m_vecAdjacent.begin(); iterAdjacent3 != m_vecSamplePoints[*iterAdjacent1].m_vecAdjacent.end(); ++iterAdjacent3)
							{
								if((*iterAdjacent3 > iterSamplePoint->iIndex) && (std::find(m_vecSamplePoints[*iterAdjacent2].m_vecAdjacent.begin(), m_vecSamplePoints[*iterAdjacent2].m_vecAdjacent.end(), (*iterAdjacent3)) != m_vecSamplePoints[*iterAdjacent2].m_vecAdjacent.end()))
								{
									newFacet.lpiVertices[0] = iterSamplePoint->iIndex;
									newFacet.lpiVertices[1] = *iterAdjacent1;
									newFacet.lpiVertices[2] = *iterAdjacent3;
									newFacet.lpiVertices[3] = *iterAdjacent2;
									m_vecFacetLinkages.push_back(newFacet);
								}
							}
						}
						++iterAdjacent2;
					}
				}
			}
		}
    }
}


void CSampler::m_fnCorrectFacetDirection()
{
    CMeshBase::Vector_3 vtOriginalNorm, vtRemeshedNorm;
    CMeshBase::Point_3 lpptPoints[4];
    CMeshBase::FT ftProduct;
    int i, iTemp;
    std::vector<CFacetLinkage>::iterator iterFacetLinkage;
    for(iterFacetLinkage = m_vecFacetLinkages.begin(); iterFacetLinkage != m_vecFacetLinkages.end(); ++iterFacetLinkage)
    {
        vtOriginalNorm = CMeshBase::Vector_3(0, 0, 0);
        for(i = 0; i < 3; ++i)
        {
            vtOriginalNorm = vtOriginalNorm + m_vecSampleBlocks[m_vecSamplePoints[iterFacetLinkage->lpiVertices[i]].iBlockId].hhBase->facet()->vtNorm;
        }
        for(i = 0; i < 4; ++i)
        {
            lpptPoints[i] = m_vecSamplePoints[iterFacetLinkage->lpiVertices[i]].ptCoordinate;
        }
        vtRemeshedNorm = CGAL::cross_product((lpptPoints[1] - lpptPoints[0]), (lpptPoints[3] - lpptPoints[0])) + CGAL::cross_product((lpptPoints[3] - lpptPoints[2]), (lpptPoints[1] - lpptPoints[2]));
        ftProduct = vtOriginalNorm * vtRemeshedNorm;
        if(ftProduct < 0)
        {
            iTemp = iterFacetLinkage->lpiVertices[1];
            iterFacetLinkage->lpiVertices[1] = iterFacetLinkage->lpiVertices[3];
            iterFacetLinkage->lpiVertices[3] = iTemp;
        }
    }
}

void CSampler::m_fnExtractNewMesh()
{
    int nNumBoundaries;
    CMeshBase::CPolyhedron::Vertex_iterator vi;
    CMeshBase::CSamplerPolyhedron::Vertex_iterator svi;
    std::vector<CSamplePoint>::iterator iterSamplePoint;
    nNumBoundaries = 0;
    for(vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
    {
        if(vi->iBordId > nNumBoundaries)
        {
            nNumBoundaries = vi->iBordId;
        }
    }
    ++nNumBoundaries;
	CSamplePolyhedronModifier SamplePolyhedronModifier(m_vecSamplePoints, m_vecFacetLinkages, m_mshRemeshedSurface);
	SamplePolyhedronModifier.GenerateMesh();
    iterSamplePoint = m_vecSamplePoints.begin();
    svi = m_mshRemeshedSurface.vertices_begin();
    while(svi != m_mshRemeshedSurface.vertices_end() && iterSamplePoint != m_vecSamplePoints.end())
    {
        svi->hhBase = m_vecSampleBlocks[iterSamplePoint->iBlockId].hhBase;
        ++svi;
        ++iterSamplePoint;
    }
}

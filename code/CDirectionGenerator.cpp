#include "CDirectionGenerator.h"

void CDirectionGenerator::m_fnMarkSharpFeature(CMeshBase::FT ftSharpBound)
{
    CMeshBase::CPolyhedron::Edge_iterator ei;
    CMeshBase::CPolyhedron::Facet_handle fh0, fh1;
	CMeshBase::FT ftCosAngle;
	CMeshBase::FT ftCosSharpBound;
	ftCosSharpBound = -cos(ftSharpBound);
    for(ei = m_mshSurface.edges_begin(); ei != m_mshSurface.edges_end(); ++ei)
    {
		if (!ei->is_border_edge())
		{
			fh0 = ei->facet();
			fh1 = ei->opposite()->facet();
			ftCosAngle = ((fh0->vtNorm) * (fh1->vtNorm)) / ((fh0->ftTwArea) * (fh1->ftTwArea));
			if(ftCosAngle < ftCosSharpBound)
			{
				fh0->ftChartDir = atan2(ei->vtVector * fh0->vtAuxiliaryAxis, ei->vtVector * fh0->vtPrincipalAxis);
				fh1->ftChartDir = atan2(ei->vtVector * fh1->vtAuxiliaryAxis, ei->vtVector * fh1->vtPrincipalAxis);
				fh0->bDirFixed = true;
				fh1->bDirFixed = true;
				ei->bSharp = true;
				ei->opposite()->bSharp = true;
			}
			else
			{
				ei->bSharp = false;
				ei->opposite()->bSharp = false;
			}
		}
		else
		{
			ei->bSharp = false;
			ei->opposite()->bSharp = false;
		}
    }
}

void  CDirectionGenerator::m_fnInitFixedDirection()
{
    CMeshBase::CPolyhedron::Facet_iterator fi;
    CMeshBase::CPolyhedron::Halfedge_handle hh0, hh;
	std::queue<CMeshBase::CPolyhedron::Facet_handle> qufh;
    int iFacetType;
    for(fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
    {
        fi->iDistOrd = -1;
    }
    for(fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
    {
		if (!fi->bDirFixed)
		{
			hh0 = fi->halfedge();
			hh = hh0;
			iFacetType = 0;
			do
			{
				if (hh->opposite()->is_border())//(hh->bFreeEdge)
				{
					++iFacetType;
				}
				hh = hh->next();
			}while(hh != hh0);
			if(iFacetType == 1)
			{
				fi->iDistOrd = 0;
				do
				{
					if(hh->opposite()->is_border())
					{
						break;
					}
					hh = hh->next();
				}while(hh != hh0);
				switch(hh->iLocalInd)
				{
				case 0:
					fi->ftChartDir = 0;
					break;
				case 1:
					fi->ftChartDir = hh->prev()->ftAngle;
					break;
				case 2:
					fi->ftChartDir = -hh->ftAngle;
					break;
				}
				fi->bDirFixed = true;
				qufh.push(&*(fi));
			}
			else if(iFacetType == 2)
			{
				fi->iDistOrd = 0;
				do
				{
					if(!hh->opposite()->is_border())
					{
						break;
					}
					hh = hh->next();
				}while(hh != hh0);
				switch(hh->iLocalInd)
				{
				case 0:
					fi->ftChartDir = hh->ftAngle;
					break;
				case 1:
					fi->ftChartDir = -hh->next()->ftAngle;
					break;
				case 2:
					fi->ftChartDir = 0;
					break;
				}
				fi->ftChartDir += (hh->next()->ftAngle - (CMeshBase::FT(int(hh->next()->ftAngle / ftHalfPi + 0.5))) * ftHalfPi) * 0.5;
				fi->bDirFixed = true;
				qufh.push(&(*fi));
			}
			else if (iFacetType == 3)
			{
				fi->bDirFixed = true;
				fi->iDistOrd = 0;
			}
			else
			{
				fi->bDirFixed = false;
			}
		}
		else
		{
			fi->iDistOrd = 0;
		}
	}
		//while (!qufh.empty())
		//{
		//	fh0 = qufh.front();
		//	qufh.pop();
		//	hh0 = fh0->halfedge();
		//	hh = hh0;
		//	do
		//	{
		//		if (hh->ftLaplaceCoef <= 0.0 && !hh->opposite()->is_border())
		//		{
		//			fh = hh->opposite()->facet();
		//			if (!fh->bDirFixed)
		//			{
		//				fh->ftChartDir = fh0->ftChartDir - hh->ftPolarAxisDif;
		//				fh->ftChartDir = fh->ftChartDir - (int(fh->ftChartDir / ftTwoPi + 15.5) - 15) * ftTwoPi;
		//				fh->bDirFixed = true;
		//				qufh.push(fh);
		//			}
		//		}
		//		hh = hh->next();
		//	} while (hh != hh0);
		//}
}

void CDirectionGenerator::m_fnFacetPrincipleDirection(CMeshBase::CPolyhedron::Facet_handle fh)
{
    CMeshBase::CPolyhedron::Halfedge_handle hh0, hh, hh1, hh2;
    int i, j, k, n;
    CMeshBase::FT ftMinCurDif = 4.0, ftMaxErr = 4.0;
    CMeshBase::FT lpftWeightedCurvature[3];
    std::vector<CMeshBase::FT> lpvecftQuadraticItem[3], vecftZ;
    std::vector<CMeshBase::FT>::iterator iterftI, iterftJ, iterftSqX, iterftXY, iterftSqY, iterftZ;
    CMeshBase::Vector_3 vtSample, lpvtRow[3], lpvtCross[3], vtProjection, ftDot;
    CMeshBase::FT ftXCur, ftYCur, ftZCur, ftVolume, ftErr, ftSumSqErr, ftDirU, ftDirV, ftTestScale, ftTestErr;
    CMeshBase::FT lpftSymmetryMatrix[9], lpftQuadraticCoef[3], lpftTestCoef[3];
    CMeshBase::FT ftEquationCoef0, ftEquationCoef1, ftSquareRootedDiscriminant, ftEigenValue0, ftEigenValue1;
    hh0 = fh->halfedge();
    hh = hh0;
    n = 0;
    do
    {
        hh1 = hh->next()->opposite();
        hh2 = hh->opposite()->prev();
        while(hh1 != hh2)
        {
            vtSample = hh1->facet()->ptIncent - fh->ptIncent;
            ftXCur = vtSample * fh->vtPrincipalAxis / fh->ftRadius;
            ftYCur = vtSample * fh->vtAuxiliaryAxis / fh->ftRadius;
            ftZCur = vtSample * fh->vtNorm / fh->ftTwArea;
            lpvecftQuadraticItem[0].push_back((ftXCur * ftXCur));
            lpvecftQuadraticItem[1].push_back((ftXCur * ftYCur));
            lpvecftQuadraticItem[2].push_back((ftYCur * ftYCur));
            vecftZ.push_back(ftZCur);
            ++n;
            hh1 = hh1->next()->opposite();
        }
        hh = hh->next();
    }while(hh != hh0);
    for(i = 0; i < 3; ++i)
    {
        for(j = i; j < 3; ++j)
        {
            iterftI = lpvecftQuadraticItem[i].begin();
            iterftJ = lpvecftQuadraticItem[j].begin();
            lpftSymmetryMatrix[i * 3 + j] = 0;
            for(k = 0; k < n; ++k)
            {
                lpftSymmetryMatrix[i * 3 + j] += (*iterftI) * (*iterftJ);
                ++iterftI;
                ++iterftJ;
            }
        }
        iterftI = lpvecftQuadraticItem[i].begin();
        iterftJ = vecftZ.begin();
        lpftWeightedCurvature[i] = 0;
        for(k = 0; k < n; ++k)
        {
            lpftWeightedCurvature[i] += (*iterftI) * (*iterftJ);
            ++iterftI;
            ++iterftJ;
        }
    }
    lpftSymmetryMatrix[3] = lpftSymmetryMatrix[1];
    lpftSymmetryMatrix[6] = lpftSymmetryMatrix[2];
    lpftSymmetryMatrix[7] = lpftSymmetryMatrix[5];
    for(i = 0; i < 3; ++i)
    {
        lpvtRow[i] = CMeshBase::Vector_3(lpftSymmetryMatrix[i * 3], lpftSymmetryMatrix[i * 3 + 1], lpftSymmetryMatrix[i * 3 + 2]);
    }
    for(i = 0; i < 3; ++i)
    {
        j = (i + 1) % 3;
        k = (i + 2) % 3;
        lpvtCross[i] = CGAL::cross_product(lpvtRow[j], lpvtRow[k]);
    }
    ftVolume = lpvtCross[0] * lpvtRow[0];
    vtProjection = CMeshBase::Vector_3(lpftWeightedCurvature[0], lpftWeightedCurvature[1], lpftWeightedCurvature[2]) / ftVolume;
    for(i = 0; i < 3; ++i)
    {
        lpftQuadraticCoef[i] = lpvtCross[i] * vtProjection;
    }
    iterftSqX = lpvecftQuadraticItem[0].begin();
    iterftXY = lpvecftQuadraticItem[1].begin();
    iterftSqY = lpvecftQuadraticItem[2].begin();
    iterftZ = vecftZ.begin();
    ftSumSqErr = 0;
    for(k = 0; k < n; ++k)
    {
        ftErr = (*iterftZ) - (*iterftSqX) * lpftQuadraticCoef[0] - (*iterftXY) * lpftQuadraticCoef[1] - (*iterftSqY) * lpftQuadraticCoef[2];
        ftSumSqErr += ftErr * ftErr;
        ++iterftSqX;
        ++iterftXY;
        ++iterftSqY;
        ++iterftZ;
    }
    lpftQuadraticCoef[1] /= 2;
    ftEquationCoef0 = lpftQuadraticCoef[0] * lpftQuadraticCoef[2] - lpftQuadraticCoef[1] * lpftQuadraticCoef[1];
    ftEquationCoef1 = (lpftQuadraticCoef[0] + lpftQuadraticCoef[2]) / 2;
    ftSquareRootedDiscriminant = sqrt(ftEquationCoef1 * ftEquationCoef1 - ftEquationCoef0);
    ftEigenValue0 = ftEquationCoef1 + ftSquareRootedDiscriminant;
    ftEigenValue1 = ftEquationCoef1 - ftSquareRootedDiscriminant;
    if(fabs(ftEigenValue0) > fabs(ftEigenValue1) * ftMinCurDif)
    {
        if(lpftQuadraticCoef[0] > 0)
        {
            ftDirU = ftEigenValue0 - (lpftQuadraticCoef[0] - lpftQuadraticCoef[1]);
            ftDirV = -ftEigenValue0 - (lpftQuadraticCoef[1] - lpftQuadraticCoef[2]);
        }
        else
        {
            ftDirU = ftEigenValue0 - (lpftQuadraticCoef[0] + lpftQuadraticCoef[1]);
            ftDirV = ftEigenValue0 - (lpftQuadraticCoef[1] + lpftQuadraticCoef[2]);
        }
        ftTestScale = ftEigenValue1 / (ftDirU * ftDirU + ftDirV * ftDirV);
        lpftTestCoef[0] = ftDirU * ftDirU * ftTestScale;
        lpftTestCoef[1] = 2 * ftDirU * ftDirV * ftTestScale;
        lpftTestCoef[2] = ftDirV * ftDirV * ftTestScale;
        ftTestErr = 0;
        iterftSqX = lpvecftQuadraticItem[0].begin();
        iterftXY = lpvecftQuadraticItem[1].begin();
        iterftSqY = lpvecftQuadraticItem[2].begin();
        iterftZ = vecftZ.begin();
        for(k = 0; k < n; ++k)
        {
            ftErr = *iterftZ - (*iterftSqX) * (lpftTestCoef[0]) - (*iterftXY) * (lpftTestCoef[1]) - (*iterftSqY) * (lpftTestCoef[2]);
            ftTestErr += ftErr * ftErr;
            ++iterftSqX;
            ++iterftXY;
            ++iterftSqY;
            ++iterftZ;
        }
        if(ftSumSqErr * ftMaxErr < ftTestErr)
        {
            fh->ftChartDir = atan2(ftDirV, ftDirU);
            fh->iDistOrd = 0;
        }
    }
    else if(fabs(ftEigenValue1) > fabs(ftEigenValue0) * ftMinCurDif)
    {
        if(lpftQuadraticCoef[1] < 0)
        {
            ftDirU = ftEigenValue1- (lpftQuadraticCoef[0] - lpftQuadraticCoef[1]);
            ftDirV = -ftEigenValue1 - (lpftQuadraticCoef[1] - lpftQuadraticCoef[2]);
        }
        else
        {
            ftDirU = ftEigenValue1- (lpftQuadraticCoef[0] + lpftQuadraticCoef[1]);
            ftDirV = ftEigenValue1 - (lpftQuadraticCoef[1] + lpftQuadraticCoef[2]);
        }
        ftTestScale = ftEigenValue0 / (ftDirU * ftDirU + ftDirV * ftDirV);
        lpftTestCoef[0] = ftDirU * ftDirU * ftTestScale;
        lpftTestCoef[1] = 2 * ftDirU * ftDirV * ftTestScale;
        lpftTestCoef[2] = ftDirV * ftDirV * ftTestScale;
        ftTestErr = 0;
        iterftSqX = lpvecftQuadraticItem[0].begin();
        iterftXY = lpvecftQuadraticItem[1].begin();
        iterftSqY = lpvecftQuadraticItem[2].begin();
        iterftZ = vecftZ.begin();
        for(k = 0; k < n; ++k)
        {
            ftErr = *iterftZ - (*iterftSqX) * (lpftTestCoef[0]) - (*iterftXY) * (lpftTestCoef[1]) - (*iterftSqY) * (lpftTestCoef[2]);
            ftTestErr += ftErr * ftErr;
            ++iterftSqX;
            ++iterftXY;
            ++iterftSqY;
            ++iterftZ;
        }
        if(ftSumSqErr * ftMaxErr < ftTestErr)
        {
            fh->ftChartDir = atan2(ftDirV, ftDirU);
            fh->iDistOrd = 0;
        }
    }
    for(i = 0; i < 3; ++i)
    {
        lpvecftQuadraticItem[i].clear();
    }
    vecftZ.clear();
}

void CDirectionGenerator::m_fnExtendDirection()
{
    std::queue<CMeshBase::CPolyhedron::Facet_handle> qufhDirFacet;
    CMeshBase::CPolyhedron::Facet_iterator fi;
    CMeshBase::CPolyhedron::Facet_handle fh0, fh;
    CMeshBase::CPolyhedron::Halfedge_handle hh0, hh, hh1;
    CMeshBase::FT ftSeam, ftArg1, ftArg2;
    int iFacetType;
    hh1 = NULL;
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
	{
		if (fi->bDirFixed)
		{
			fi->iDistOrd = 0;
		}
		else
		{
			fi->iDistOrd = -1;
		}
	}
    for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
    {
        if(fi->iDistOrd == 0)
        {
            hh0 = fi->halfedge();
            hh = hh0;
            do
            {
                if(!hh->opposite()->is_border())
                {
                    fh = hh->opposite()->facet();
                    if(fh->iDistOrd == -1)
                    {
                        fh->iDistOrd = 1;
                        qufhDirFacet.push(fh);
                    }
                }
                hh = hh->next();
            }while(hh != hh0);
        }
    }
    while(!qufhDirFacet.empty())
    {
        fh0 = qufhDirFacet.front();
        qufhDirFacet.pop();
        iFacetType = 0;
        hh0 = fh0->halfedge();
        hh = hh0;
        do
        {
			if (!hh->opposite()->is_border())
			{ 
				fh = hh->opposite()->facet();
				if(fh->iDistOrd == fh0->iDistOrd - 1)
				{
					++iFacetType;
				}
			}
            hh = hh->next();
        }while(hh != hh0);
        if(iFacetType == 1)
        {
            hh0 = fh0->halfedge();
            hh = hh0;
            do
            {
				if (!hh->opposite()->is_border())
				{
					fh = hh->opposite()->facet();
					if(fh->iDistOrd == fh0->iDistOrd - 1)
					{
						hh1 = hh;
					}
				}
                hh = hh->next();
            }while(hh != hh0);
            fh = hh1->opposite()->facet();
            fh0->ftChartDir = fh->ftChartDir + hh1->ftPolarAxisDif;
            fh0->ftChartDir = fh0->ftChartDir - (CMeshBase::FT(int(fh0->ftChartDir / ftTwoPi + 15.5)) - 15) * ftTwoPi;
            fh = hh1->next()->opposite()->facet();
            if(fh->iDistOrd == -1)
            {
                fh->iDistOrd = fh0->iDistOrd +1;
                qufhDirFacet.push(fh);
            }
            fh = hh1->prev()->opposite()->facet();
            if(fh->iDistOrd == -1)
            {
                fh->iDistOrd = fh0->iDistOrd +1;
                qufhDirFacet.push(fh);
            }
        }
        else if(iFacetType == 2)
        {
            hh0 = fh0->halfedge();
            hh = hh0;
            do
            {
				if (!hh->opposite()->is_border())
				{
					fh = hh->opposite()->facet();
					if(fh->iDistOrd != fh0->iDistOrd - 1)
					{
						hh1 = hh;
					}
				}
                hh = hh->next();
            }while(hh != hh0);
            fh = hh1->next()->opposite()->facet();
            fh0->ftChartDir = fh->ftChartDir + hh1->next()->ftPolarAxisDif;
            fh = hh1->prev()->opposite()->facet();
            ftSeam = fh->ftChartDir + hh1->prev()->ftPolarAxisDif - fh0 ->ftChartDir;
            ftSeam = ftSeam - (int(ftSeam / ftHalfPi + 15.5) - 15) * ftHalfPi;
            fh0->ftChartDir += ftSeam / 2;
            fh0->ftChartDir = fh0->ftChartDir - (CMeshBase::FT(int(fh0->ftChartDir / ftTwoPi + 15.5)) - 15) * ftTwoPi;
            fh = hh1->opposite()->facet();
            if(fh->iDistOrd == -1)
            {
                fh->iDistOrd = fh0->iDistOrd + 1;
                qufhDirFacet.push(fh);
            }
        }
        else
        {
            hh1 = hh0;
            fh = hh1->opposite()->facet();
            fh0->ftChartDir = fh->ftChartDir +  hh1->ftPolarAxisDif;
            hh1 = hh1->next();
            fh = hh1->opposite()->facet();
            ftArg1 = fh->ftChartDir + hh1->ftPolarAxisDif - fh0->ftChartDir;
            ftArg1 = ftArg1 - (int(ftArg1 / ftHalfPi + 15.5) - 15) * ftHalfPi;
            hh1 = hh1->next();
            fh = hh1->opposite()->facet();
            ftArg2 = fh->ftChartDir + hh1->ftPolarAxisDif - fh0->ftChartDir;
            ftArg2 = ftArg2 - (int(ftArg2 / ftHalfPi + 15.5) - 15) * ftHalfPi;
            if(ftArg1 > ftArg2)
            {
                ftSeam = ftArg2;
                ftArg2 = ftArg1;
                ftArg1 = ftSeam;
            }
            if(ftArg2 - ftArg1 > ftHalfPi)
            {
                if(ftArg2 + ftArg1 > 0)
                {
                    ftArg2 -= ftHalfPi;
                }
                else
                {
                    ftArg1 += ftHalfPi;
                }
            }
            fh0->ftChartDir += (ftArg1 + ftArg2) / 3;
        }
    }
}

void CDirectionGenerator::m_fnRandomizeDirection()
{
	CMeshBase::CPolyhedron::Facet_iterator fi;
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
	{
		if (!fi->bDirFixed)
		{
			fi->ftChartDir = 0;// (rand() / double(RAND_MAX) * 2.0 - 1.0) * ftPi;
			//fi->ftChartDir = fi->ftChartDir - (CMeshBase::FT(int(fi->ftChartDir / ftTwoPi + 15.5)) - 15) * ftTwoPi;
		}
	}
}

void CDirectionGenerator::m_fnMergeSort(CMeshBase::FT* lpftData, int* lpiOrder, int nLen)
{
    int* lpiBuffer, *lplpiArray[2], iSource, iDest, nSubLen, iReader1, iReader2, iWriter, iEnd1, iEnd2, i;
    lpiBuffer = new int[nLen];
    for (i = 0; i < nLen; ++i)
    {
        lpiOrder[i] = i;
    }
    lplpiArray[0] = lpiOrder;
    lplpiArray[1] = lpiBuffer;
    iSource = 0;
    nSubLen = 1;
    while (nSubLen < nLen)
    {
        iDest = (iSource + 1) % 2;
        iReader1 = 0;
        iWriter = 0;
        while (iReader1 < nLen)
        {
            iEnd1 = iReader2 = iReader1 + nSubLen;
            if (iEnd1 > nLen)
            {
                iReader2 = iEnd1 = nLen;
            }
            iEnd2 = iReader2 + nSubLen;
            if (iEnd2 > nLen)
            {
                iEnd2 = nLen;
            }
            while (iWriter < iEnd2)
            {
                if (iReader1 == iEnd1 || (iReader2 < iEnd2 && lpftData[lplpiArray[iSource][iReader2]] < lpftData[lplpiArray[iSource][iReader1]]))
                {
                    lplpiArray[iDest][iWriter] = lplpiArray[iSource][iReader2];
                    ++iReader2;
                }
                else
                {
                    lplpiArray[iDest][iWriter] = lplpiArray[iSource][iReader1];
                    ++iReader1;
                }
                ++iWriter;
            }
            iReader1 = iEnd2;
        }
        nSubLen *= 2;
        iSource = iDest;
    }
    if (iSource == 1)
    {
        for (iWriter = 0; iWriter < nLen; ++iWriter)
        {
            lpiOrder[iWriter] = lpiBuffer[iWriter];
        }
    }
    delete[]lpiBuffer;
}

void CDirectionGenerator::m_fnMarkIsolatedCones()
{
    CMeshBase::CPolyhedron::Vertex_iterator vi;
    CMeshBase::CPolyhedron::Vertex_handle *lpvhBuffer;
    CMeshBase::CPolyhedron::FT* lpftGaussCurvature, ftSumArea;
    CMeshBase::CPolyhedron::Vertex_handle vh;
    CMeshBase::CPolyhedron::Halfedge_handle hh0, hh;
    int iInd, *lpiOrder, nLen;
    bool bAdjacent;
    iInd = 0;
    for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
    {
        vi->nNumPorts = 0;
        vi->iDegreeBias = 0;
        if (vi->iBordId == -1)
        {
            vi->iEquationIndex = 0;
            vi->iTempIndex = iInd;
            ++iInd;
        }
        else
        {
            vi->iEquationIndex = -1;
        }
    }
    nLen = iInd;
    lpvhBuffer = new CMeshBase::CPolyhedron::Vertex_handle[nLen];
    lpftGaussCurvature = new CMeshBase::FT[nLen];
    lpiOrder = new int[nLen];
    iInd = 0;
    for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
    {
        if (vi->iBordId == -1)
        {
            lpvhBuffer[iInd] = &(*vi);
            ftSumArea = 0;
            hh0 = vi->halfedge();
            hh = hh0;
            do
            {
                ftSumArea += (hh->ftSqLen * hh->ftLaplaceCoef) / 8.0;
                hh = hh->next()->opposite();
            } while (hh != hh0);
            lpftGaussCurvature[iInd] = -abs(vi->ftAngleDefect / ftSumArea);
            ++iInd;
        }
    }
    m_fnMergeSort(lpftGaussCurvature, lpiOrder, nLen);
    bAdjacent = false;
    for (iInd = 0; iInd < nLen ; ++iInd)
    {
        vh = lpvhBuffer[lpiOrder[iInd]];
        bAdjacent = false;
        hh0 = vh->halfedge()->opposite();
        hh = hh0;
        do
        {
            if (hh->vertex()->iEquationIndex == -1)
            {
                bAdjacent = true;
            }
            hh = hh->opposite()->next();
        } while (hh != hh0 && !bAdjacent);
        if (!bAdjacent)
        {
            vh->iEquationIndex = -1;
            if (vh->ftAngleDefect > 0)
            {
                vh->nNumPorts = 3;
            }
            else
            {
                vh->nNumPorts = 5;
            }
        }
    }
    delete[]lpiOrder;
    delete[]lpftGaussCurvature;
    delete[]lpvhBuffer;
}

void CDirectionGenerator::m_fnGenerateMetric()
{
    CMeshBase::CPolyhedron::Vertex_iterator vi;
    CMeshBase::CPolyhedron::Halfedge_handle hh0, hh;
    MatElem<CMeshBase::FT> meDiag, meTri;
    CMeshBase::CLDLTSv LDLTS;
    int iEquationIndex;
    iEquationIndex = 0;
    for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
    {
        if (vi->iEquationIndex != -1)
        {
            vi->iEquationIndex = iEquationIndex;
            ++iEquationIndex;
        }
    }
    LDLTS.init(iEquationIndex);
    for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
    {
        if (vi->iEquationIndex != -1)
        {
            meDiag.row = vi->iEquationIndex;
            meDiag.col = vi->iEquationIndex;
            meDiag.data = 0;
            hh0 = vi->halfedge();
            hh = hh0;
            do
            {
                if (!hh->is_border_edge())
                {
                    meDiag.data += hh->ftLaplaceCoef;
                    if (hh->prev()->vertex()->iEquationIndex != -1)
                    {
                        meTri.row = vi->iEquationIndex;
                        meTri.col = hh->prev()->vertex()->iEquationIndex;
                        meTri.data = -hh->ftLaplaceCoef;
                        LDLTS.add_data(meTri.row, meTri.col, meTri.data);
                    }
                }
                hh = hh->next()->opposite();
            } while (hh != hh0);
            LDLTS.add_data(meDiag.row, meDiag.col, meDiag.data);
            if (vi->nNumPorts == 0)
            {
                LDLTS.m_b[vi->iEquationIndex] = -vi->ftAngleDefect;
            }
            else
            {
                LDLTS.m_b[vi->iEquationIndex] = -vi->ftAngleDefect + (4 - vi->nNumPorts) * ftHalfPi;
            }
        }
    }
    LDLTS.decompose();
    LDLTS.solve();
    for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
    {
        if (vi->iEquationIndex == -1)
        {
            vi->ftDensity = 0;
        }
        else
        {
            vi->ftDensity = LDLTS.m_x[vi->iEquationIndex];
        }
    }
    LDLTS.clear();
}

void CDirectionGenerator::m_fnFixCones(CMeshBase::FT ftThreshold)
{
    CMeshBase::CPolyhedron::Vertex_iterator vi;
    CMeshBase::CPolyhedron::Halfedge_handle hh0, hh;
    CMeshBase::FT ftPoisson, ftRoundErr;
    int iCone;
    for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
    {
        if (vi->iEquationIndex == -1)
        {
            ftPoisson = 0;
            hh0 = vi->halfedge()->opposite();
            hh = hh0;
            do
            {
                ftPoisson += hh->vertex()->ftDensity * hh->ftLaplaceCoef;
                hh = hh->opposite()->next();
            } while (hh != hh0);
            iCone = int((-ftPoisson + vi->ftAngleDefect) / ftHalfPi + 15.5) - 15;
            ftRoundErr = fabs(vi->ftAngleDefect - ftPoisson - iCone * ftHalfPi);
            if (ftRoundErr < ftThreshold)
            {
                vi->iEquationIndex = 0;
                if (iCone == 0)
                {
                    vi->nNumPorts = 0;
                }
                else
                {
                    vi->nNumPorts = 4 - iCone;
                }
            }
        }
    }
}

int CDirectionGenerator::m_fnFixSmallestErrorCone()
{
    CMeshBase::CPolyhedron::Vertex_iterator vi;
    CMeshBase::CPolyhedron::Halfedge_handle hh0, hh;
    CMeshBase::CPolyhedron::Vertex_handle vhNextCone;
    CMeshBase::FT ftPoisson, ftRoundErr, ftMinRoundErr;
    int iCone, iMinErrCone, nUnfixed;
    nUnfixed = 0;
    iMinErrCone = 0;
    vhNextCone = NULL;
    ftMinRoundErr = ftHalfPi;
    for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
    {
        if (vi->iEquationIndex == -1)
        {
            ftPoisson = 0;
            hh0 = vi->halfedge()->opposite();
            hh = hh0;
            do
            {
                ftPoisson += hh->vertex()->ftDensity * hh->ftLaplaceCoef;
                hh = hh->opposite()->next();
            } while (hh != hh0);
            iCone = int((-ftPoisson + vi->ftAngleDefect) / ftHalfPi + 15.5) - 15;
            ftRoundErr = fabs(vi->ftAngleDefect - ftPoisson - iCone * ftHalfPi);
            if (ftRoundErr < ftMinRoundErr)
            {
                vhNextCone = &(*vi);
                ftMinRoundErr = ftRoundErr;
                iMinErrCone = iCone;
            }
            ++nUnfixed;
        }
    }
    if (vhNextCone != NULL)
    {

        vhNextCone->iEquationIndex = 0;
        if (iMinErrCone == 0)
        {
            vhNextCone->nNumPorts = 0;
        }
        else
        {
            vhNextCone->nNumPorts = 4 - iMinErrCone;
        }
        --nUnfixed;
    }
    return nUnfixed;
}


CMeshBase::FT CDirectionGenerator::m_fnMIQEnergy()
{
    CMeshBase::FT ftEnergy;
    CMeshBase::CPolyhedron::Edge_iterator ei;
    ftEnergy = 0;
    for(ei = m_mshSurface.edges_begin(); ei != m_mshSurface.edges_end(); ++ei)
    {
        if(!ei->is_border_edge())
        {
            ftEnergy += (ei->ftMom * ei->ftMom);
        }
    }
    return ftEnergy;
}

void CDirectionGenerator::m_fnFillComplex(CMeshBase::CCmplxSv &CmplxSv)
{
	CMeshBase::CPolyhedron::Facet_iterator fi;
	CMeshBase::CPolyhedron::Halfedge_handle hh0, hh;
	MatElem<std::complex<CMeshBase::FT>> meDiag, meTri;
	int i;
	i = 0;
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
	{
		fi->iIndex = i;
		++i;
	}
	CmplxSv.init(i);
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
	{
		meDiag.row = meDiag.col = fi->iIndex;
		meDiag.data = std::complex<CMeshBase::FT>(0.0);
		meTri.row = fi->iIndex;
		hh0 = fi->halfedge();
		hh = hh0;
		do
		{
			if (!hh->opposite()->is_border())
			{
				meTri.col = hh->opposite()->facet()->iIndex;
				meTri.data = -std::complex<CMeshBase::FT>(cos(hh->ftPolarAxisDif * 4), sin(hh->ftPolarAxisDif * 4)) / hh->ftLaplaceCoef;
				CmplxSv.add_data(meTri.row, meTri.col, meTri.data);
			}
			meDiag.data += 1.0 / hh->ftLaplaceCoef;
			hh = hh->next();
		} while (hh != hh0);
		CmplxSv.add_data(meDiag.row, meDiag.col, meDiag.data);
	}
}

void CDirectionGenerator::m_fnSolveEigen(CMeshBase::CCmplxSv &CmplxSv)
{
	CmplxSv.decompose();
	int i, j;
	CMeshBase::FT ftNorm;
	for (i = 0; i < CmplxSv.m_dim; ++i)
	{
		CmplxSv.m_x[i] = std::complex<CMeshBase::FT>(1.0);
	}
	for (i = 0; i < 256; ++i)
	{
		ftNorm = 0;
		for (j = 0; j < CmplxSv.m_dim; ++j)
		{
			ftNorm += (CmplxSv.m_x[j].real() * CmplxSv.m_x[j].real() + CmplxSv.m_x[j].imag() * CmplxSv.m_x[j].imag());
		}
		ftNorm = sqrt(ftNorm / CmplxSv.m_dim);
		for (j = 0; j < CmplxSv.m_dim; ++j)
		{
			CmplxSv.m_b[j] = CmplxSv.m_x[j] / ftNorm;
		}
		CmplxSv.solve();
	}
}

void CDirectionGenerator::m_fnFillOptDir(CMeshBase::CCmplxSv &CmplxSv)
{
	CMeshBase::CPolyhedron::Facet_iterator fi;
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
	{
		fi->ftChartDir = atan2(CmplxSv.m_x[fi->iIndex].imag(), CmplxSv.m_x[fi->iIndex].real()) / 4;
	}
}

void CDirectionGenerator::m_fnFillDirFromPara(bool bFrameField)
{
	CMeshBase::CPolyhedron::Facet_iterator fi;
	CMeshBase::CPolyhedron::Halfedge_handle hh0, hh;
	CMeshBase::Vector_3 lpvtGrad[2], vtAve;
	CMeshBase::FT ftSqLen0, ftSqLen1;
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
	{
		lpvtGrad[1] = lpvtGrad[0] = CMeshBase::Vector_3(0, 0, 0);
		hh0 = fi->halfedge();
		hh = hh0;
		do
		{
			lpvtGrad[0] = lpvtGrad[0] + hh->vtVector * hh->next()->lpftGlobalPara[0];
			lpvtGrad[1] = lpvtGrad[1] + hh->vtVector * hh->next()->lpftGlobalPara[1];
			hh = hh->next();
		} while (hh != hh0);
		if (bFrameField)
		{
			fi->ftChartDir = atan2(fi->vtAuxiliaryAxis * lpvtGrad[0], fi->vtPrincipalAxis * lpvtGrad[0]);
			if (CGAL::cross_product(lpvtGrad[0], lpvtGrad[1]) * fi->vtNorm < 0)
			{
				fi->ftFrameDir = atan2(-fi->vtAuxiliaryAxis * lpvtGrad[1], -fi->vtPrincipalAxis * lpvtGrad[1]);
			}
			else
			{
				fi->ftFrameDir = atan2(fi->vtAuxiliaryAxis * lpvtGrad[1], fi->vtPrincipalAxis * lpvtGrad[1]);
			}
		}
		else
		{
			ftSqLen0 = CGAL::squared_length(lpvtGrad[0]);
			ftSqLen1 = CGAL::squared_length(lpvtGrad[1]);
			if (ftSqLen1 < ftSqLen0)
			{
				fi->ftChartDir = atan2(fi->vtAuxiliaryAxis * lpvtGrad[0], fi->vtPrincipalAxis * lpvtGrad[0]);
				fi->ftFrameDir = atan2(fi->vtPrincipalAxis * lpvtGrad[0], -fi->vtAuxiliaryAxis * lpvtGrad[0]);
			}
			else
			{
				fi->ftChartDir = atan2(-fi->vtPrincipalAxis * lpvtGrad[1], fi->vtAuxiliaryAxis * lpvtGrad[1]);
				fi->ftFrameDir = atan2(fi->vtAuxiliaryAxis * lpvtGrad[1], fi->vtPrincipalAxis * lpvtGrad[1]);
			}
		}
	}
}

void CDirectionGenerator::m_fnGenerateDirFromDensity()
{
    CMeshBase::CPolyhedron::Facet_iterator fi;
    CMeshBase::CPolyhedron::Halfedge_handle hh0, hh;
    CMeshBase::CPolyhedron::Facet_handle fh;
    std::queue<CMeshBase::CPolyhedron::Facet_handle> qufh;
    for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
    {
        fi->ftChartDir = -ftTwoPi * 2;
    }
    --fi;
    fi->ftChartDir = 0.0;
    qufh.push(&(*fi));
    while (!qufh.empty())
    {
        fh = qufh.front();
        qufh.pop();
        hh0 = fh->halfedge();
        hh = hh0;
        do
        {
            if (!hh->opposite()->is_border())
            {
                if (hh->opposite()->facet()->ftChartDir < -ftTwoPi)
                {
                    hh->opposite()->facet()->ftChartDir = fh->ftChartDir - hh->ftPolarAxisDif + (hh->vertex()->ftDensity - hh->prev()->vertex()->ftDensity) * hh->ftLaplaceCoef;
                    hh->opposite()->facet()->ftChartDir = hh->opposite()->facet()->ftChartDir - CMeshBase::FT(int(hh->opposite()->facet()->ftChartDir / ftHalfPi + 15.5) - 15) * ftHalfPi;
                    qufh.push(hh->opposite()->facet());
                }
            }
            hh = hh->next();
        } while (hh != hh0);
    }
}
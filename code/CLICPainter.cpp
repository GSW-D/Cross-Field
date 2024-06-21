#include "CLICPainter.h"
CLICPainter::CLICPainter(CGLCanvas *lpGLCanvas)
{
	m_lpdpBuffer = NULL;
	m_iXOffset = m_iYOffset = m_iHeight = m_iWidth = 0;
	m_lpGLCanvas = lpGLCanvas;
}
CLICPainter::~CLICPainter()
{
	if (m_lpdpBuffer != NULL)
	{
		delete[]m_lpdpBuffer;
	}
}
void CLICPainter::m_fnGenerateViewCoordinates()
{
	CMeshBase::CPolyhedron::Vertex_iterator vi;
	for (vi = m_lpGLCanvas->m_lpMeshCenter->m_mshSurface.vertices_begin(); vi != m_lpGLCanvas->m_lpMeshCenter->m_mshSurface.vertices_end(); ++vi)
	{
		m_lpGLCanvas->m_fnViewCoordinate(vi->point().x(), vi->point().y(), vi->point().z(), vi->lpdblViewCoord[0], vi->lpdblViewCoord[1], vi->lpdblViewCoord[2]);
	}
}
void CLICPainter::m_fnTestSize(double dblRadiusScale)
{
	CMeshBase::CPolyhedron::Vertex_iterator vi;
	int iMinX, iMaxX, iMinY, iMaxY, iX, iY, i;
	double lpdblCoordinates[3];
	iMinX = INT_MAX;
	iMaxX = -INT_MAX;
	iMinY = INT_MAX;
	iMaxY = -INT_MAX;
	for (vi = m_lpGLCanvas->m_lpMeshCenter->m_mshSurface.vertices_begin(); vi != m_lpGLCanvas->m_lpMeshCenter->m_mshSurface.vertices_end(); ++vi)
	{
		iX = int(floor(vi->lpdblViewCoord[0]));
		iY = int(floor(vi->lpdblViewCoord[1]));
		if (iX < iMinX)
		{
			iMinX = iX;
		}
		if (iY < iMinY)
		{
			iMinY = iY;
		}
		iX = int(ceil(vi->lpdblViewCoord[0]));
		iY = int(ceil(vi->lpdblViewCoord[1]));
		if (iX > iMaxX)
		{
			iMaxX = iX;
		}
		if (iY > iMaxY)
		{
			iMaxY = iY;
		}
	}
	if (dblRadiusScale > 0.0)
	{
		for (vi = m_lpGLCanvas->m_lpMeshCenter->m_mshSurface.vertices_begin(); vi != m_lpGLCanvas->m_lpMeshCenter->m_mshSurface.vertices_end(); ++vi)
		{
			if (vi->iBordId == -1 && vi->nNumPorts != 0)
			{
				for (i = 0; i < 162; ++i)
				{
					m_lpGLCanvas->m_fnViewCoordinate(\
						vi->point().x() + m_lpGLCanvas->m_lpMeshCenter->m_ftMeanRadius * dblRadiusScale * m_lpGLCanvas->m_lpdblSphereCoordinates[i * 3], \
						vi->point().y() + m_lpGLCanvas->m_lpMeshCenter->m_ftMeanRadius * dblRadiusScale * m_lpGLCanvas->m_lpdblSphereCoordinates[i * 3 + 1], \
						vi->point().z() + m_lpGLCanvas->m_lpMeshCenter->m_ftMeanRadius * dblRadiusScale * m_lpGLCanvas->m_lpdblSphereCoordinates[i * 3 + 2], \
						lpdblCoordinates[0], lpdblCoordinates[1], lpdblCoordinates[2]\
					);
					iX = int(floor(lpdblCoordinates[0]));
					iY = int(floor(lpdblCoordinates[1]));
					if (iX < iMinX)
					{
						iMinX = iX;
					}
					if (iY < iMinY)
					{
						iMinY = iY;
					}
					iX = int(ceil(lpdblCoordinates[0]));
					iY = int(ceil(lpdblCoordinates[1]));
					if (iX > iMaxX)
					{
						iMaxX = iX;
					}
					if (iY > iMaxY)
					{
						iMaxY = iY;
					}
				}
			}
		}
	}
	m_iXOffset = iMinX;
	m_iYOffset = iMinY;
	m_iWidth = iMaxX - iMinX + 1;
	m_iHeight = iMaxY - iMinY + 1;
}
void CLICPainter::m_fnVideoTestSize(double dblRadiusScale, int iWidth, int iHeight)
{
//	CMeshBase::CPolyhedron::Vertex_iterator vi;
//	int iMinX, iMaxX, iMinY, iMaxY, iX, iY, i;
//	double lpdblCoordinates[3];
//	iMinX = INT_MAX;
//	iMaxX = -INT_MAX;
//	iMinY = INT_MAX;
//	iMaxY = -INT_MAX;
//	for (vi = m_lpGLCanvas->m_lpMeshCenter->m_mshSurface.vertices_begin(); vi != m_lpGLCanvas->m_lpMeshCenter->m_mshSurface.vertices_end(); ++vi)
//	{
//		iX = int(floor(vi->lpdblViewCoord[0]));
//		iY = int(floor(vi->lpdblViewCoord[1]));
//		if (iX < iMinX)
//		{
//			iMinX = iX;
//		}
//		if (iY < iMinY)
//		{
//			iMinY = iY;
//		}
//		iX = int(ceil(vi->lpdblViewCoord[0]));
//		iY = int(ceil(vi->lpdblViewCoord[1]));
//		if (iX > iMaxX)
//		{
//			iMaxX = iX;
//		}
//		if (iY > iMaxY)
//		{
//			iMaxY = iY;
//		}
//	}
//	if (dblRadiusScale > 0.0)
//	{
//		for (vi = m_lpGLCanvas->m_lpMeshCenter->m_mshSurface.vertices_begin(); vi != m_lpGLCanvas->m_lpMeshCenter->m_mshSurface.vertices_end(); ++vi)
//		{
//			if (vi->iBordId == -1 && vi->nNumPorts != 0)
//			{
//				for (i = 0; i < 162; ++i)
//				{
//					m_lpGLCanvas->m_fnViewCoordinate(\
//						vi->point().x() + m_lpGLCanvas->m_lpMeshCenter->m_ftMeanRadius * dblRadiusScale * m_lpGLCanvas->m_lpdblSphereCoordinates[i * 3], \
//						vi->point().y() + m_lpGLCanvas->m_lpMeshCenter->m_ftMeanRadius * dblRadiusScale * m_lpGLCanvas->m_lpdblSphereCoordinates[i * 3 + 1], \
//						vi->point().z() + m_lpGLCanvas->m_lpMeshCenter->m_ftMeanRadius * dblRadiusScale * m_lpGLCanvas->m_lpdblSphereCoordinates[i * 3 + 2], \
//						lpdblCoordinates[0], lpdblCoordinates[1], lpdblCoordinates[2]\
//					);
//					iX = int(floor(lpdblCoordinates[0]));
//					iY = int(floor(lpdblCoordinates[1]));
//					if (iX < iMinX)
//					{
//						iMinX = iX;
//					}
//					if (iY < iMinY)
//					{
//						iMinY = iY;
//					}
//					iX = int(ceil(lpdblCoordinates[0]));
//					iY = int(ceil(lpdblCoordinates[1]));
//					if (iX > iMaxX)
//					{
//						iMaxX = iX;
//					}
//					if (iY > iMaxY)
//					{
//						iMaxY = iY;
//					}
//				}
//			}
//		}
//	}
	m_iXOffset = m_lpGLCanvas->m_iHalfWidth - iWidth / 2;//(iMaxX + iMinX - iWidth) / 2;
	m_iYOffset = m_lpGLCanvas->m_iHalfHeight - iHeight / 2;//(iMaxY + iMinY - iHeight) / 2;
	m_iWidth = iWidth;
	m_iHeight = iHeight;
}
void CLICPainter::m_fnInitBuffer()
{
	int nNumPixels, i;
	DirPixel *lpdpCur;
	if (m_lpdpBuffer != NULL)
	{
		delete[]m_lpdpBuffer;
	}
	nNumPixels = m_iWidth * m_iHeight;
	m_lpdpBuffer = new DirPixel[nNumPixels];
	lpdpCur = m_lpdpBuffer;
	for (i = 0; i < nNumPixels; ++i)
	{
		lpdpCur->fDepth = FLT_MAX;
		lpdpCur->hhLocation = NULL;
		lpdpCur->hhCrevasse = NULL;
		lpdpCur->lpiSeam[0] = lpdpCur->lpiSeam[1] = lpdpCur->lpiSeam[2] = lpdpCur->lpiSeam[3] = 2;
		lpdpCur->lpfDirection[0] = lpdpCur->lpfDirection[1] = 4.0f;
		lpdpCur->lpfColor[0] = lpdpCur->lpfColor[1] = lpdpCur->lpfColor[2] = 1.0f;
		++lpdpCur;
	}
}

void CLICPainter::m_fnSetLocations(bool bCrossField)
{
	CMeshBase::CPolyhedron::Facet_iterator fi;
	CMeshBase::CPolyhedron::Halfedge_handle hh0, hh;
	int iXStart, iXEnd, iYStart, iYEnd, iX, iY, i;
	double *lplpdblCoord[3], lpdblPixelCoord[2], lpdblAreaCoord[3], dblSumAreaCoord, lpdblViewCoord[6];
	float fDepth;
	CMeshBase::Vector_3 vtDirection;
	CMeshBase::Point_3 pntPosition;
	DirPixel *lpdpCurrent;
	for (fi = m_lpGLCanvas->m_lpMeshCenter->m_mshSurface.facets_begin(); fi != m_lpGLCanvas->m_lpMeshCenter->m_mshSurface.facets_end(); ++fi)
	{
		iXStart = INT_MAX;
		iXEnd = -INT_MAX;
		iYStart = INT_MAX;
		iYEnd = -INT_MAX;
		hh0 = fi->halfedge();
		hh = hh0;
		do
		{
			iX = int(floor(hh->vertex()->lpdblViewCoord[0]));
			iY = int(floor(hh->vertex()->lpdblViewCoord[1]));
			if (iX < iXStart)
			{
				iXStart = iX;
			}
			if (iY < iYStart)
			{
				iYStart = iY;
			}
			iX = int(ceil(hh->vertex()->lpdblViewCoord[0]));
			iY = int(ceil(hh->vertex()->lpdblViewCoord[1]));
			if (iX > iXEnd)
			{
				iXEnd = iX;
			}
			if (iY > iYEnd)
			{
				iYEnd = iY;
			}
			hh = hh->next();
		} while (hh != hh0);
		lplpdblCoord[0] = hh0->vertex()->lpdblViewCoord;
		lplpdblCoord[1] = hh0->next()->vertex()->lpdblViewCoord;
		lplpdblCoord[2] = hh0->prev()->vertex()->lpdblViewCoord;
		for (iX = iXStart; iX <= iXEnd; ++iX)
		{
			for (iY = iYStart; iY <= iYEnd; ++iY)
			{
				lpdblPixelCoord[0] = double(iX);
				lpdblPixelCoord[1] = double(iY);
				for (i = 0; i < 3; ++i)
				{
					lpdblAreaCoord[i] = m_fnPlanarArea(lpdblPixelCoord, lplpdblCoord[(i + 1) % 3], lplpdblCoord[(i + 2) % 3]);
				}
				dblSumAreaCoord = lpdblAreaCoord[0] + lpdblAreaCoord[1] + lpdblAreaCoord[2];
				if ((lpdblAreaCoord[0] >= 0 && lpdblAreaCoord[1] >= 0 && lpdblAreaCoord[2] >= 0 && dblSumAreaCoord > 0)||
					(lpdblAreaCoord[0] <= 0 && lpdblAreaCoord[1] <= 0 && lpdblAreaCoord[2] <= 0 && dblSumAreaCoord < 0))
				{
					lpdblAreaCoord[0] /= dblSumAreaCoord;
					lpdblAreaCoord[1] /= dblSumAreaCoord;
					lpdblAreaCoord[2] /= dblSumAreaCoord;
					fDepth = float(lplpdblCoord[0][2] * lpdblAreaCoord[0] + lplpdblCoord[1][2] * lpdblAreaCoord[1] + lplpdblCoord[2][2] * lpdblAreaCoord[2]);
					lpdpCurrent = m_lpdpBuffer + (iY - m_iYOffset) * m_iWidth + (iX - m_iXOffset);
					if (lpdpCurrent->hhLocation == NULL || lpdpCurrent->fDepth > fDepth)
					{
						lpdpCurrent->hhLocation = hh0;
						lpdpCurrent->fDepth = fDepth;
						pntPosition = CGAL::barycenter(hh0->vertex()->point(), lpdblAreaCoord[0], hh0->next()->vertex()->point(), lpdblAreaCoord[1], hh0->prev()->vertex()->point(), lpdblAreaCoord[2]);
						vtDirection = hh0->lpvtViewField[0] * lpdblAreaCoord[0] + hh0->next()->lpvtViewField[0] * lpdblAreaCoord[1] + hh0->prev()->lpvtViewField[0] * lpdblAreaCoord[2];
						m_lpGLCanvas->m_fnViewCoordinate(pntPosition.x(), pntPosition.y(), pntPosition.z(), lpdblViewCoord[0], lpdblViewCoord[1], lpdblViewCoord[2]);
						m_lpGLCanvas->m_fnViewCoordinate(pntPosition.x() + vtDirection.x(), pntPosition.y() + vtDirection.y(), pntPosition.z() + vtDirection.z(), lpdblViewCoord[3], lpdblViewCoord[4], lpdblViewCoord[5]);
						lpdpCurrent->lpfDirection[0] = atan2(lpdblViewCoord[4] - lpdblViewCoord[1], lpdblViewCoord[3] - lpdblViewCoord[0]);//float(m_lpGLCanvas->m_fnDirectionOnScreen(vtDirection.x(), vtDirection.y(), vtDirection.z()));
						if (bCrossField)
						{
							vtDirection = hh0->lpvtViewField[1] * lpdblAreaCoord[0] + hh0->next()->lpvtViewField[1] * lpdblAreaCoord[1] + hh0->prev()->lpvtViewField[1] * lpdblAreaCoord[2];
							m_lpGLCanvas->m_fnViewCoordinate(pntPosition.x(), pntPosition.y(), pntPosition.z(), lpdblViewCoord[0], lpdblViewCoord[1], lpdblViewCoord[2]);
							m_lpGLCanvas->m_fnViewCoordinate(pntPosition.x() + vtDirection.x(), pntPosition.y() + vtDirection.y(), pntPosition.z() + vtDirection.z(), lpdblViewCoord[3], lpdblViewCoord[4], lpdblViewCoord[5]);
							lpdpCurrent->lpfDirection[1] = atan2(lpdblViewCoord[4] - lpdblViewCoord[1], lpdblViewCoord[3] - lpdblViewCoord[0]);//float(m_lpGLCanvas->m_fnDirectionOnScreen(vtDirection.x(), vtDirection.y(), vtDirection.z()));
						}
						else
						{
							lpdpCurrent->fStrength = float(hh0->vertex()->ftFieldStrength * lpdblAreaCoord[0] + hh0->next()->vertex()->ftFieldStrength * lpdblAreaCoord[1] + hh0->prev()->vertex()->ftFieldStrength * lpdblAreaCoord[2]);
						}
					}
				}
			}
		}
	}
}

void CLICPainter::m_fnSetCrevasse(double dblWidth)
{
    CMeshBase::CPolyhedron::Edge_iterator ei;
	double lpdblDir[2], lpdblStart[2], lpdblEnd[2], dblLen, lpdblNorm[2], dblDist, lpdblProject[2];
	float fDepth;
	int i, j, iX, iY;
	bool bInside, bAroundStart, bAroundEnd;
	DirPixel *lpdpCurrent;
    for (ei = m_lpGLCanvas->m_lpMeshCenter->m_mshSurface.edges_begin(); ei != m_lpGLCanvas->m_lpMeshCenter->m_mshSurface.edges_end(); ++ei)
    {
        if (ei->bFreeEdge)
        {
            lpdblStart[0] = ei->prev()->vertex()->lpdblViewCoord[0];
            lpdblStart[1] = ei->prev()->vertex()->lpdblViewCoord[1];
            lpdblEnd[0] = ei->vertex()->lpdblViewCoord[0];
            lpdblEnd[1] = ei->vertex()->lpdblViewCoord[1];
            lpdblDir[0] = lpdblEnd[0] - lpdblStart[0];
            lpdblDir[1] = lpdblEnd[1] - lpdblStart[1];
            dblLen = sqrt(lpdblDir[0] * lpdblDir[0] + lpdblDir[1] * lpdblDir[1]);
            lpdblDir[0] /= dblLen;
            lpdblDir[1] /= dblLen;
            lpdblNorm[0] = -lpdblDir[1];
            lpdblNorm[1] = lpdblDir[0];
            for (i = -(ceil(dblWidth / 2.0) + 1); i < ceil(dblLen + dblWidth) + 1; ++i)
            {
                for (j = -(ceil(dblWidth / 2.0) + 1); j < ceil(dblWidth / 2.0) + 1; ++j)
                {
                    bAroundStart = false;
                    bAroundEnd = false;
                    bInside = false;
                    iX = round(lpdblStart[0] + lpdblDir[0] * i - lpdblDir[1] * j);
                    iY = round(lpdblStart[1] + lpdblDir[1] * i + lpdblDir[0] * j);
                    dblDist = (iX - lpdblStart[0]) * lpdblNorm[0] + (iY - lpdblStart[1]) * lpdblNorm[1];
                    lpdblProject[0] = (iX - lpdblStart[0]) * lpdblDir[0] + (iY - lpdblStart[1]) * lpdblDir[1];
                    lpdblProject[1] = (lpdblEnd[0] - iX) * lpdblDir[0] + (lpdblEnd[0] - iY) * lpdblDir[1];
                    if (dblDist <= dblWidth / 2 && dblDist >= -dblWidth / 2)
                    {
                        if (lpdblProject[0] < 0)
                        {
                            if (dblDist * dblDist + lpdblProject[0] * lpdblProject[0] <= dblWidth * dblWidth / 4)
                            {
                                bAroundStart = true;
                            }
                        }
                        else if(lpdblProject[1] < 0)
                        {
                            if (dblDist * dblDist + lpdblProject[1] * lpdblProject[1] <= dblWidth * dblWidth / 4)
                            {
                                bAroundEnd = true;
                            }
                        }
                        else
                        {
                            bInside = true;
                        }
                    }
                    if (bInside || bAroundStart || bAroundEnd)
                    {
                        lpdpCurrent = m_lpdpBuffer + (iY - m_iYOffset) * m_iWidth + (iX - m_iXOffset);
                        fDepth = float((ei->vertex()->lpdblViewCoord[2] * lpdblProject[0] + ei->prev()->vertex()->lpdblViewCoord[2] * lpdblProject[1]) / (lpdblProject[0] + lpdblProject[1]));
                        if (fDepth < lpdpCurrent->fDepth + 0.01f)
                        {
                            lpdpCurrent->hhCrevasse = &(*ei);
                            lpdpCurrent->fDepth = fDepth;
                        }
                    }
                }
            }
        }
    }
}


int CLICPainter::m_fnAdjacent(DirPixel *lpdp0, DirPixel *lpdp1)
{
	CMeshBase::CPolyhedron::Halfedge_handle hh0, hh1;
	int iSeamType, i;
	float fDif0, fDif1, fSumDif, fMinSumDif, fMinDif0, fMinDif1;
	if (lpdp0->hhLocation == NULL || lpdp1->hhLocation == NULL)
	{
		iSeamType = 2;
	}
	else
	{
		if (lpdp0->hhLocation->facet() == lpdp1->hhLocation->facet())
		{
			iSeamType = 0;
		}
		else
		{
			iSeamType = 2;
			hh0 = lpdp0->hhLocation;
			do
			{
				hh1 = lpdp1->hhLocation;
				do
				{
					if (hh1->vertex() == hh0->vertex())
					{
						iSeamType = 0;
					}
					hh1 = hh1->next();
				} while (hh1 != lpdp1->hhLocation);
				hh0 = hh0->next();
			} while (hh0 != lpdp0->hhLocation);
		}
	}
	if (iSeamType == 0)
	{
		fMinSumDif = FLT_MAX;
		for (i = -2; i < 2; ++i)
		{
			fDif0 = lpdp1->lpfDirection[(i + 2) % 2] + ftPi * (1 - ((i + 3) % 4) / 2) - lpdp0->lpfDirection[0];
			fDif0 = fDif0 - round(fDif0 / ftTwoPi) * ftTwoPi;
			fDif1 = lpdp1->lpfDirection[1 - (i + 2) % 2] + ftPi * (1 - (i + 2) / 2) - lpdp0->lpfDirection[1];
			fDif1 = fDif1 - round(fDif1 / ftTwoPi) * ftTwoPi;
			fSumDif = fDif0 * fDif0 + fDif1 * fDif1;
			if (fSumDif < fMinSumDif)
			{
				iSeamType = i;
				fMinDif0 = fDif0;
				fMinDif1 = fDif1;
				fMinSumDif = fSumDif;
			}
		}
		if (fMinDif0 > ftHalfPi || fMinDif0 < -ftHalfPi || fMinDif1 > ftHalfPi || fMinDif1 < -ftHalfPi)
		{
			iSeamType = 2;
		}
	}
	return iSeamType;
}

void CLICPainter::m_fnDetectSeam()
{
	int i, j;
	DirPixel *lpdpCur, *lpdpNext;
	for (i = 0; i < m_iHeight; ++i)
	{
		for (j = 0; j < m_iWidth - 1; ++j)
		{
			lpdpCur = m_lpdpBuffer + i * m_iWidth + j;
			lpdpNext = lpdpCur + 1;
			if (lpdpCur->hhLocation != NULL && lpdpNext->hhLocation != NULL)
			{
				lpdpCur->lpiSeam[0] = m_fnAdjacent(lpdpCur, lpdpNext);
				if (lpdpCur->lpiSeam[0] != 2)
				{
					lpdpNext->lpiSeam[2] = (2 - lpdpCur->lpiSeam[0]) % 4 - 2;
				}
				else
				{
					lpdpNext->lpiSeam[2] = 2;
				}
			}
		}
	}
	for (i = 0; i < m_iHeight - 1; ++i)
	{
		for (j = 0; j < m_iWidth; ++j)
		{
			lpdpCur = m_lpdpBuffer + i * m_iWidth + j;
			lpdpNext = lpdpCur + m_iWidth;
			if (lpdpCur->hhLocation != NULL && lpdpNext->hhLocation != NULL)
			{
				lpdpCur->lpiSeam[1] = m_fnAdjacent(lpdpCur, lpdpNext);
				if (lpdpCur->lpiSeam[1] != 2)
				{
					lpdpNext->lpiSeam[3] = (2 - lpdpCur->lpiSeam[1]) % 4 - 2;
				}
				else
				{
					lpdpNext->lpiSeam[3] = 2;
				}
			}
		}
	}
}

void CLICPainter::m_fnAssignOriginalColor(wxBitmap &bmpOriginal)
{
	int i, j, iRow, iCol, iOriginalWidth, iOriginalHeight;
	wxMemoryDC memdc;
	DirPixel *lpdpCur;
	iOriginalWidth = bmpOriginal.GetWidth();
	iOriginalHeight = bmpOriginal.GetHeight();
	lpdpCur = m_lpdpBuffer;
	memdc.SelectObject(bmpOriginal);
	wxColour clrCurrent;
	lpdpCur = m_lpdpBuffer;
	for (i = 0; i < m_iHeight; ++i)
	{
		for (j = 0; j < m_iWidth; ++j)
		{
			iRow = i % iOriginalHeight;
			iCol = j % iOriginalWidth;
			memdc.GetPixel(iCol, iRow, &clrCurrent);
			lpdpCur->lpfColor[0] = clrCurrent.Red() / 255.0f;
			lpdpCur->lpfColor[1] = clrCurrent.Green() / 255.0f;
			lpdpCur->lpfColor[2] = clrCurrent.Blue() / 255.0f;
			++lpdpCur;
		}
	}

}

void CLICPainter::m_fnCrossTrace(float fRadius)
{
	int i, j, k, iInlet, iDir;
	DirPixel *lpdpCur, *lpdpTr;
	float fAccLen, fInletPos, fIncrement, fAngle, fMaxTan, fMinTan, fRelativeAngle, fSinRelativeAngle, fCosRelativeAngle, fSumWeight, fDist, fWeight, fAngleDif;
	lpdpCur = m_lpdpBuffer;
	for (i = 0; i < m_iHeight; ++i)
	{
		for (j = 0; j < m_iWidth; ++j)
		{
			if (lpdpCur->hhLocation != NULL)
			{
				lpdpCur->lpfColor[3] = lpdpCur->lpfColor[4] = lpdpCur->lpfColor[5] = 0.0f;
				fSumWeight = 0.0f;
				for (k = -2; k < 2; ++k)
				{
					fAccLen = 0;
					iDir = k;
					lpdpTr = lpdpCur;
					if (lpdpTr->lpfDirection[(iDir + 2) % 2] != 4.0f)
					{
						fAngle = lpdpTr->lpfDirection[(iDir + 2) % 2] + ftPi * (1 - (iDir + 2)/2);
						fAngle = fAngle - round(fAngle / ftTwoPi) * ftTwoPi;
					}
					else
					{
						fAngle = 4.0f;
					}
					if ((fAngle > -3.5 && fAngle <= -ftPi * 3 / 4) || (fAngle >= ftPi * 3 / 4 && fAngle < 3.5))
					{
						iInlet = 2;
						fInletPos = sin(fAngle) / cos(fAngle);
						fIncrement = -0.5 / cos(fAngle);
					}
					else if (fAngle > -ftPi * 3 / 4 && fAngle < -ftPi / 4)
					{
						iInlet = 3;
						fInletPos = -cos(fAngle) / sin(fAngle);
						fIncrement = -0.5 / sin(fAngle);
					}
					else if (fAngle >= -ftPi / 4 && fAngle <= ftPi / 4)
					{
						iInlet = 0;
						fInletPos = tan(fAngle);
						fIncrement = 0.5 / cos(fAngle);
					}
					else if (fAngle > ftPi / 4 && fAngle <= ftPi * 3 / 4)
					{
						iInlet = 1;
						fInletPos = -cos(fAngle) / sin(fAngle);
						fIncrement = 0.5 / sin(fAngle);
					}
					else
					{
						fIncrement = 0.0;
						iInlet = -1;
					}
					if (iInlet != -1 && lpdpCur->lpiSeam[iInlet] == 2)
					{
						iInlet = -1;
					}
					if (iInlet != -1)
					{
						fDist = 0;
						fWeight = fIncrement;
					}
					else
					{
						fWeight = 1.0;
					}
					lpdpCur->lpfColor[3] += lpdpCur->lpfColor[0] * fWeight;
					lpdpCur->lpfColor[4] += lpdpCur->lpfColor[1] * fWeight;
					lpdpCur->lpfColor[5] += lpdpCur->lpfColor[2] * fWeight;
					fSumWeight += fWeight;
					fAccLen += fIncrement;
					while (iInlet != -1 && fAccLen < fRadius)
					{
						iDir = (iDir - lpdpTr->lpiSeam[iInlet] + 6) % 4 - 2;
						switch (iInlet)
						{
						case 0:
							++lpdpTr;
							break;
						case 1:
							lpdpTr += m_iWidth;
							break;
						case 2:
							--lpdpTr;
							break;
						case 3:
							lpdpTr -= m_iWidth;
							break;
						}
						fAngleDif = fAngle;
						if (lpdpTr->lpfDirection[(iDir + 2) % 2] != 4.0f)
						{
							fAngle = lpdpTr->lpfDirection[(iDir + 2) % 2] + ftPi * (1 - (iDir + 2) / 2);
							fAngle = fAngle - round(fAngle / ftTwoPi) * ftTwoPi;
							fAngleDif -= fAngle;
							fAngleDif = fAngleDif - round(fAngleDif / ftTwoPi) * ftTwoPi;
							if (fAngleDif > ftHalfPi || fAngleDif < -ftHalfPi)
							{
								fAngle = 4.0f;
							}
						}
						else
						{
							fAngle = 4.0f;
						}
						if (fAngle < -3.5 || fAngle > 3.5)
						{
							iInlet = -1;
						}
						else
						{
							fRelativeAngle = fAngle - ftHalfPi *iInlet;
							fSinRelativeAngle = sin(fRelativeAngle);
							fCosRelativeAngle = cos(fRelativeAngle);
							if (fCosRelativeAngle < 0)
							{
								fCosRelativeAngle = 0;
								if (fSinRelativeAngle < 0)
								{
									fSinRelativeAngle = -1;
								}
								else
								{
									fSinRelativeAngle = 1;
								}
							}
							fMaxTan = (1 - fInletPos) / 2.0f;
							fMinTan = -(fInletPos + 1) / 2.0f;
							if (fSinRelativeAngle < 0 && fCosRelativeAngle * fMinTan > fSinRelativeAngle)
							{
								iInlet = (iInlet + 3) % 4;
								fIncrement = -(fInletPos + 1) * 0.5f / fSinRelativeAngle;
								fInletPos = -(fInletPos + 1) * fCosRelativeAngle / fSinRelativeAngle - 1;
							}
							else if (fSinRelativeAngle > 0 && fCosRelativeAngle * fMaxTan < fSinRelativeAngle)
							{
								iInlet = (iInlet + 1) % 4;
								fIncrement = (1 - fInletPos) * 0.5f / fSinRelativeAngle;
								fInletPos = 1 - (1 - fInletPos) * fCosRelativeAngle / fSinRelativeAngle;
							}
							else
							{
								fIncrement = 1.0f / fCosRelativeAngle;
								fInletPos += 2 * fSinRelativeAngle / fCosRelativeAngle;
							}
							if (lpdpTr->lpiSeam[iInlet] == 2)
							{
								iInlet = -1;
							}
						}
						if (fAccLen + fIncrement > fRadius)
						{
							fDist = (fAccLen + fRadius) / 2.0f;
							fWeight = fRadius - fAccLen;
						}
						else
						{
							fDist = fAccLen + fIncrement / 2.0f;
							fWeight = fIncrement;
						}
						lpdpCur->lpfColor[3] += lpdpTr->lpfColor[0] * fWeight;
						lpdpCur->lpfColor[4] += lpdpTr->lpfColor[1] * fWeight;
						lpdpCur->lpfColor[5] += lpdpTr->lpfColor[2] * fWeight;
						fSumWeight += fWeight;
						fAccLen += fIncrement;
					}
				}
				lpdpCur->lpfColor[3] /= fSumWeight;
				lpdpCur->lpfColor[4] /= fSumWeight;
				lpdpCur->lpfColor[5] /= fSumWeight;
			}
			++lpdpCur;
		}
	}
}


void CLICPainter::m_fnLinearTrace(float fRadius)
{
	int i, j, k, iInlet;
	DirPixel *lpdpCur, *lpdpTr;
	float fAccLen, fInletPos, fIncrement, fAngle, fMaxTan, fMinTan, fRelativeAngle, fSinRelativeAngle, fCosRelativeAngle, fSumWeight, fDist, fWeight, fAngleDif;
	lpdpCur = m_lpdpBuffer;
	for (i = 0; i < m_iHeight; ++i)
	{
		for (j = 0; j < m_iWidth; ++j)
		{
			if (lpdpCur->hhLocation != NULL)
			{
				lpdpCur->lpfColor[3] = lpdpCur->lpfColor[4] = lpdpCur->lpfColor[5] = 0.0f;
				fSumWeight = 0.0f;
				for (k = 0; k < 2; ++k)
				{
					fAccLen = 0;
					lpdpTr = lpdpCur;
					if (lpdpTr->lpfDirection[0] != 4.0f)
					{
						fAngle = lpdpTr->lpfDirection[0] + ftPi * k;
						fAngle = fAngle - round(fAngle / ftTwoPi) * ftTwoPi;
					}
					else
					{
						fAngle = 4.0f;
					}
					if ((fAngle > -3.5 && fAngle <= -ftPi * 3 / 4) || (fAngle >= ftPi * 3 / 4 && fAngle < 3.5))
					{
						iInlet = 2;
						fInletPos = sin(fAngle) / cos(fAngle);
						fIncrement = -0.5 / cos(fAngle);
					}
					else if (fAngle > -ftPi * 3 / 4 && fAngle < -ftPi / 4)
					{
						iInlet = 3;
						fInletPos = -cos(fAngle) / sin(fAngle);
						fIncrement = -0.5 / sin(fAngle);
					}
					else if (fAngle >= -ftPi / 4 && fAngle <= ftPi / 4)
					{
						iInlet = 0;
						fInletPos = tan(fAngle);
						fIncrement = 0.5 / cos(fAngle);
					}
					else if (fAngle > ftPi / 4 && fAngle <= ftPi * 3 / 4)
					{
						iInlet = 1;
						fInletPos = -cos(fAngle) / sin(fAngle);
						fIncrement = 0.5 / sin(fAngle);
					}
					else
					{
						fIncrement = 0.0;
						iInlet = -1;
					}
					if (iInlet != -1)
					{
						fDist = 0;
						fWeight = fIncrement;
					}
					else
					{
						fWeight = 1.0;
					}
					lpdpCur->lpfColor[3] += lpdpCur->lpfColor[0] * fWeight;
					lpdpCur->lpfColor[4] += lpdpCur->lpfColor[1] * fWeight;
					lpdpCur->lpfColor[5] += lpdpCur->lpfColor[2] * fWeight;
					fSumWeight += fWeight;
					fAccLen += fIncrement;
					while (iInlet != -1 && fAccLen < fRadius)
					{
						switch (iInlet)
						{
						case 0:
							++lpdpTr;
							break;
						case 1:
							lpdpTr += m_iWidth;
							break;
						case 2:
							--lpdpTr;
							break;
						case 3:
							lpdpTr -= m_iWidth;
							break;
						}
						fAngleDif = fAngle;
						if (lpdpTr->lpfDirection[0] != 4.0f)
						{
							fAngle = lpdpTr->lpfDirection[0] + ftPi * k;
							fAngle = fAngle - round(fAngle / ftTwoPi) * ftTwoPi;
							fAngleDif -= fAngle;
							fAngleDif = fAngleDif - round(fAngleDif / ftTwoPi) * ftTwoPi;
							if (fAngleDif > ftHalfPi || fAngleDif < -ftHalfPi)
							{
								fAngle = 4.0f;
							}
						}
						else
						{
							fAngle = 4.0f;
						}
						if (fAngle < -3.5f || fAngle > 3.5f)
						{
							iInlet = -1;
						}
						else
						{
							fRelativeAngle = fAngle - ftHalfPi *iInlet;
							fSinRelativeAngle = sin(fRelativeAngle);
							fCosRelativeAngle = cos(fRelativeAngle);
							if (fCosRelativeAngle < 0)
							{
								fCosRelativeAngle = 0;
								if (fSinRelativeAngle < 0)
								{
									fSinRelativeAngle = -1;
								}
								else
								{
									fSinRelativeAngle = 1;
								}
							}
							fMaxTan = (1 - fInletPos) / 2.0f;
							fMinTan = -(fInletPos + 1) / 2.0f;
							if (fSinRelativeAngle < 0 && fCosRelativeAngle * fMinTan > fSinRelativeAngle)
							{
								iInlet = (iInlet + 3) % 4;
								fIncrement = -(fInletPos + 1) * 0.5f / fSinRelativeAngle;
								fInletPos = -(fInletPos + 1) * fCosRelativeAngle / fSinRelativeAngle - 1;
							}
							else if (fSinRelativeAngle > 0 && fCosRelativeAngle * fMaxTan < fSinRelativeAngle)
							{
								iInlet = (iInlet + 1) % 4;
								fIncrement = (1 - fInletPos) * 0.5f / fSinRelativeAngle;
								fInletPos = 1 - (1 - fInletPos) * fCosRelativeAngle / fSinRelativeAngle;
							}
							else
							{
								fIncrement = 1.0f / fCosRelativeAngle;
								fInletPos += 2 * fSinRelativeAngle / fCosRelativeAngle;
							}
						}
						if (fAccLen + fIncrement > fRadius)
						{
							fDist = (fAccLen + fRadius) / 2.0f;
							fWeight = fRadius - fAccLen;
						}
						else
						{
							fDist = fAccLen + fIncrement / 2.0f;
							fWeight = fIncrement;
						}
						lpdpCur->lpfColor[3] += lpdpTr->lpfColor[0] * fWeight;
						lpdpCur->lpfColor[4] += lpdpTr->lpfColor[1] * fWeight;
						lpdpCur->lpfColor[5] += lpdpTr->lpfColor[2] * fWeight;
						fSumWeight += fWeight;
						fAccLen += fIncrement;
					}
				}
				lpdpCur->lpfColor[3] /= fSumWeight;
				lpdpCur->lpfColor[4] /= fSumWeight;
				lpdpCur->lpfColor[5] /= fSumWeight;
			}
			++lpdpCur;
		}
	}
}


void CLICPainter::m_fnEnhanceBitmap(bool bShowStrength, float fMaxStrength)
{
	DirPixel *lpdpCur;
	int i, j, k, n;
	float lpfAve[3], lpfSigma[3], fDif;
	double lpdblColor[3];
	lpfAve[0] = lpfAve[1] = lpfAve[2] = 0.0f;
	n = 0;
	lpdpCur = m_lpdpBuffer;
	for (i = 0; i < m_iHeight; ++i)
	{
		for (j = 0; j < m_iWidth; ++j)
		{
			if (lpdpCur->hhLocation != NULL)
			{
				for (k = 0; k < 3; ++k)
				{
					lpfAve[k] += lpdpCur->lpfColor[k + 3];
				}
				++n;
			}
			++lpdpCur;
		}
	}
	for (i = 0; i < 3; ++i)
	{
		lpfAve[i] /= float(n);
	}
	lpfSigma[0] = lpfSigma[1] = lpfSigma[2] = 0.0f;
	lpdpCur = m_lpdpBuffer;
	for (i = 0; i < m_iHeight; ++i)
	{
		for (j = 0; j < m_iWidth; ++j)
		{
			if (lpdpCur->hhLocation != NULL)
			{
				for (k = 0; k < 3; ++k)
				{
					fDif = lpdpCur->lpfColor[k + 3] - lpfAve[k];
					lpfSigma[k] += fDif * fDif;
				}
			}
			++lpdpCur;
		}
	}
	for (i = 0; i < 3; ++i)
	{
		lpfSigma[i] = sqrt(lpfSigma[i] / float(n));
	}
	lpdpCur = m_lpdpBuffer;
	if (bShowStrength)
	{
		for (i = 0; i < m_iHeight; ++i)
		{
			for (j = 0; j < m_iWidth; ++j)
			{
				if (lpdpCur->hhLocation != NULL)
				{
					m_lpGLCanvas->m_lpMeshCenter->m_fnColorMap(double(lpdpCur->fStrength), fMaxStrength, 0, lpdblColor);
					for (k = 0; k < 3; ++k)
					{
						lpdpCur->lpfColor[k + 3] = 1.0f - ((lpdpCur->lpfColor[k + 3] - lpfAve[k]) / (lpfSigma[k] * 6.0f) + 0.5f) * float(1.0f - lpdblColor[k]);
						if (lpdpCur->lpfColor[k + 3] > 1.0f)
						{
							lpdpCur->lpfColor[k + 3] = 1.0f;
						}
						if (lpdpCur->lpfColor[k + 3] < 0.0f)
						{
							lpdpCur->lpfColor[k + 3] = 0.0f;
						}
					}
				}
				++lpdpCur;
			}
		}
	}
	else
	{
		for (i = 0; i < m_iHeight; ++i)
		{
			for (j = 0; j < m_iWidth; ++j)
			{
				if (lpdpCur->hhLocation != NULL)
				{
					for (k = 0; k < 3; ++k)
					{
						lpdpCur->lpfColor[k + 3] = (lpdpCur->lpfColor[k + 3] - lpfAve[k]) / (lpfSigma[k] * 6.0f) + 0.5f;
						if (lpdpCur->lpfColor[k + 3] > 1.0f)
						{
							lpdpCur->lpfColor[k + 3] = 1.0f;
						}
						if (lpdpCur->lpfColor[k + 3] < 0.0f)
						{
							lpdpCur->lpfColor[k + 3] = 0.0f;
						}
					}
				}
				++lpdpCur;
			}
		}
	}
}

void CLICPainter::m_fnWeakFieldRender(float fMaxStrength)
{
	DirPixel *lpdpCur;
	int i, j, k, n;
	float lpfAve[3], lpfSigma[3], fDif;
	double lpdblColor[3];
	lpfAve[0] = lpfAve[1] = lpfAve[2] = 0.0f;
	n = 0;
	lpdpCur = m_lpdpBuffer;
	for (i = 0; i < m_iHeight; ++i)
	{
		for (j = 0; j < m_iWidth; ++j)
		{
			if (lpdpCur->hhLocation != NULL)
			{
				for (k = 0; k < 3; ++k)
				{
					lpfAve[k] += lpdpCur->lpfColor[k + 3];
				}
				++n;
			}
			++lpdpCur;
		}
	}
	for (i = 0; i < 3; ++i)
	{
		lpfAve[i] /= float(n);
	}
	lpfSigma[0] = lpfSigma[1] = lpfSigma[2] = 0.0f;
	lpdpCur = m_lpdpBuffer;
	for (i = 0; i < m_iHeight; ++i)
	{
		for (j = 0; j < m_iWidth; ++j)
		{
			if (lpdpCur->hhLocation != NULL)
			{
				for (k = 0; k < 3; ++k)
				{
					fDif = lpdpCur->lpfColor[k + 3] - lpfAve[k];
					lpfSigma[k] += fDif * fDif;
				}
			}
			++lpdpCur;
		}
	}
	for (i = 0; i < 3; ++i)
	{
		lpfSigma[i] = sqrt(lpfSigma[i] / float(n));
	}
	lpdpCur = m_lpdpBuffer;
	for (i = 0; i < m_iHeight; ++i)
	{
		for (j = 0; j < m_iWidth; ++j)
		{
			if (lpdpCur->hhLocation != NULL)
			{
				m_lpGLCanvas->m_lpMeshCenter->m_fnColorMap(double(lpdpCur->fStrength), fMaxStrength, 0, lpdblColor);
				for (k = 0; k < 3; ++k)
				{
					lpdpCur->lpfColor[k + 3] = 1.0f - ((lpdpCur->lpfColor[k + 3] - lpfAve[k]) * (lpdpCur->fStrength / fMaxStrength) / (lpfSigma[k] * 6.0f) + 0.5f) * float(1.0f - lpdblColor[k]);
					if (lpdpCur->lpfColor[k + 3] > 1.0f)
					{
						lpdpCur->lpfColor[k + 3] = 1.0f;
					}
					if (lpdpCur->lpfColor[k + 3] < 0.0f)
					{
						lpdpCur->lpfColor[k + 3] = 0.0f;
					}
				}
			}
			++lpdpCur;
		}
	}
}

void CLICPainter::m_fnNoFieldRender(float fMaxStrength)
{
	DirPixel *lpdpCur;
	int i, j, k;
	double lpdblColor[3];
	lpdpCur = m_lpdpBuffer;
    for (i = 0; i < m_iHeight; ++i)
    {
        for (j = 0; j < m_iWidth; ++j)
        {
            if (lpdpCur->hhLocation != NULL)
            {
                m_lpGLCanvas->m_lpMeshCenter->m_fnColorMap(double(lpdpCur->fStrength), fMaxStrength, 0, lpdblColor);
                for (k = 0; k < 3; ++k)
                {
                    lpdpCur->lpfColor[k + 3] = 1.0f - 1.0f * float(1.0f - lpdblColor[k]);
                    if (lpdpCur->lpfColor[k + 3] > 1.0f)
                    {
                        lpdpCur->lpfColor[k + 3] = 1.0f;
                    }
                    if (lpdpCur->lpfColor[k + 3] < 0.0f)
                    {
                        lpdpCur->lpfColor[k + 3] = 0.0f;
                    }
                }
            }
            ++lpdpCur;
        }
    }
}

void CLICPainter::m_fnDrawDiscrete()
{
	CMeshBase::CPolyhedron::Facet_iterator fi;
	CMeshBase::CPolyhedron::Edge_iterator ei;
	CMeshBase::CPolyhedron::Halfedge_handle hh;
	CMeshBase::Point_3 ptChartEnd;
	DirPixel *lpdpCur;
	wxColour wxclrPixel;
	int i, j, lpiTerminals[4], iStart, iEnd, iAuxStart, iAuxEnd, iChange, iX, iY, iXStart, iXEnd, iYStart, iYEnd;
	double lpdblViewCoord[3];
	wxColour lpwxclrChart[4] = {wxColour(255, 0, 0), wxColour(0, 255, 0), wxColour(0, 0, 255), wxColour(0, 0, 0)};
	wxBitmap wxImgBmp(m_iWidth, m_iHeight);
	wxMemoryDC wxImgDC(wxImgBmp);
	wxPen wxpnEdge(wxColour(127, 127, 127), 1), lpwxpnChart[4] = {wxPen(lpwxclrChart[0], 1), wxPen(lpwxclrChart[1], 1), wxPen(lpwxclrChart[2], 1), wxPen(lpwxclrChart[3], 1)};
	wxImgDC.SetBackground(*wxWHITE_BRUSH);
	wxImgDC.Clear();
	wxImgDC.SetPen(wxpnEdge);
	for (ei = m_lpGLCanvas->m_lpMeshCenter->m_mshSurface.edges_begin(); ei != m_lpGLCanvas->m_lpMeshCenter->m_mshSurface.edges_end(); ++ei)
	{
		lpiTerminals[0] = int(ei->vertex()->lpdblViewCoord[0] + 0.5 - m_iXOffset);
		lpiTerminals[1] = int(ei->vertex()->lpdblViewCoord[1] + 0.5 - m_iYOffset);
		lpiTerminals[2] = int(ei->opposite()->vertex()->lpdblViewCoord[0] + 0.5) - m_iXOffset;
		lpiTerminals[3] = int(ei->opposite()->vertex()->lpdblViewCoord[1] + 0.5) - m_iYOffset;
		if (abs(lpiTerminals[2] - lpiTerminals[0]) < abs(lpiTerminals[3] - lpiTerminals[1]))
		{
			iChange = 1;
		}
		else
		{
			iChange = 0;
		}
		if (lpiTerminals[iChange] < lpiTerminals[iChange + 2])
		{
			iStart = lpiTerminals[iChange];
			iEnd = lpiTerminals[iChange + 2];
			iAuxStart = lpiTerminals[(iChange + 1) % 2];
			iAuxEnd = lpiTerminals[(iChange + 1) % 2 + 2];
		}
		else
		{
			iStart = lpiTerminals[iChange + 2];
			iEnd = lpiTerminals[iChange];
			iAuxStart = lpiTerminals[(iChange + 1) % 2 + 2];
			iAuxEnd = lpiTerminals[(iChange + 1) % 2];
		}
		iXStart = iYStart = iXEnd = iYEnd = -1;
		for (i = iStart; i <= iEnd; ++i)
		{
			if (iChange)
			{
				iY = i;
				iX = ((iEnd - i) * iAuxStart + (i - iStart) * iAuxEnd + (iEnd - iStart) / 2) / (iEnd - iStart);
			}
			else
			{
				iX = i;
				iY = ((iEnd - i) * iAuxStart + (i - iStart) * iAuxEnd + (iEnd - iStart) / 2) / (iEnd - iStart);
			}
			lpdpCur = m_lpdpBuffer + iY * m_iWidth + iX;
			if (ei->is_border_edge())
			{
				hh = &(*ei);
				if (hh->is_border())
				{
					hh = hh->opposite();
				}
				if (lpdpCur->hhLocation == NULL)
				{
					lpdpCur->hhLocation = hh;
				}
				if (lpdpCur->hhLocation->facet() == hh->facet())
				{
					if (iXStart == -1 && iYStart == -1)
					{
						iXStart = iX;
						iYStart = iY;
					}
					iXEnd = iX;
					iYEnd = iY;
				}
				else
				{
					if (!(iXStart == -1 || iYStart == -1 || iXEnd == -1 || iYEnd == -1))
					{
						wxImgDC.DrawLine(iXStart, iYStart, iXEnd, iYEnd);
					}
					iXStart = iYStart = iXEnd = iYEnd = -1;
				}
			}
			else
			{
				if (lpdpCur->hhLocation != NULL && (lpdpCur->hhLocation->facet() == ei->facet() || lpdpCur->hhLocation->facet() == ei->opposite()->facet()))
				{
					if (iXStart == -1 && iYStart == -1)
					{
						iXStart = iX;
						iYStart = iY;
					}
					iXEnd = iX;
					iYEnd = iY;
				}
				else
				{
					if (!(iXStart == -1 || iYStart == -1 || iXEnd == -1 || iYEnd == -1))
					{
						wxImgDC.DrawLine(iXStart, iYStart, iXEnd, iYEnd);
					}
					iXStart = iYStart = iXEnd = iYEnd = -1;
				}
			}
		}
		if (!(iXStart == -1 || iYStart == -1 || iXEnd == -1 || iYEnd == -1))
		{
			wxImgDC.DrawLine(iXStart, iYStart, iXEnd, iYEnd);
		}
	}
	for (fi = m_lpGLCanvas->m_lpMeshCenter->m_mshSurface.facets_begin(); fi != m_lpGLCanvas->m_lpMeshCenter->m_mshSurface.facets_end(); ++fi)
	{
		m_lpGLCanvas->m_fnViewCoordinate(fi->ptIncent.x(), fi->ptIncent.y(), fi->ptIncent.z(), lpdblViewCoord[0], lpdblViewCoord[1], lpdblViewCoord[2]);
		lpiTerminals[0] = int(lpdblViewCoord[0] + 0.5 - m_iXOffset);
		lpiTerminals[1] = int(lpdblViewCoord[1] + 0.5 - m_iYOffset);
		for (i = 0; i < 4; ++i)
		{
			wxImgDC.SetPen(lpwxpnChart[i]);
			ptChartEnd = fi->ptIncent + (fi->vtPrincipalAxis * cos(fi->ftChartDir + i * ftHalfPi) + fi->vtAuxiliaryAxis * sin(fi->ftChartDir + i * ftHalfPi));
			m_lpGLCanvas->m_fnViewCoordinate(ptChartEnd.x(), ptChartEnd.y(), ptChartEnd.z(), lpdblViewCoord[0], lpdblViewCoord[1], lpdblViewCoord[2]);
			lpiTerminals[2] = int(lpdblViewCoord[0] + 0.5 - m_iXOffset);
			lpiTerminals[3] = int(lpdblViewCoord[1] + 0.5 - m_iYOffset);
			if (abs(lpiTerminals[2] - lpiTerminals[0]) < abs(lpiTerminals[3] - lpiTerminals[1]))
			{
				iChange = 1;
			}
			else
			{
				iChange = 0;
			}
			if (lpiTerminals[iChange] < lpiTerminals[iChange + 2])
			{
				iStart = lpiTerminals[iChange];
				iEnd = lpiTerminals[iChange + 2];
				iAuxStart = lpiTerminals[(iChange + 1) % 2];
				iAuxEnd = lpiTerminals[(iChange + 1) % 2 + 2];
			}
			else
			{
				iStart = lpiTerminals[iChange + 2];
				iEnd = lpiTerminals[iChange];
				iAuxStart = lpiTerminals[(iChange + 1) % 2 + 2];
				iAuxEnd = lpiTerminals[(iChange + 1) % 2];
			}
			iXStart = iYStart = iXEnd = iYEnd = -1;
			for (j = iStart; j <= iEnd; ++j)
			{
				if (iChange)
				{
					iY = j;
					iX = ((iEnd - j) * iAuxStart + (j - iStart) * iAuxEnd + (iEnd - iStart) / 2) / (iEnd - iStart);
				}
				else
				{
					iX = j;
					iY = ((iEnd - j) * iAuxStart + (j - iStart) * iAuxEnd + (iEnd - iStart) / 2) / (iEnd - iStart);
				}
				lpdpCur = m_lpdpBuffer + iY * m_iWidth + iX;
				if (lpdpCur->hhLocation != NULL && lpdpCur->hhLocation->facet() == &(*fi))
				{
					if (iXStart == -1 && iYStart == -1)
					{
						iXStart = iX;
						iYStart = iY;
					}
					iXEnd = iX;
					iYEnd = iY;
				}
				else
				{
					if (!(iXStart == -1 || iYStart == -1 || iXEnd == -1 || iYEnd == -1))
					{
						wxImgDC.DrawLine(iXStart, iYStart, iXEnd, iYEnd);
					}
					iXStart = iYStart = iXEnd = iYEnd = -1;
				}
			}
			if (!(iXStart == -1 || iYStart == -1 || iXEnd == -1 || iYEnd == -1))
			{
				wxImgDC.DrawLine(iXStart, iYStart, iXEnd, iYEnd);
			}
		}
	}
	lpdpCur = m_lpdpBuffer;
	for (i = 0; i < m_iHeight; ++i)
	{
		for (j = 0; j < m_iWidth; ++j)
		{
			wxImgDC.GetPixel(j, i, &wxclrPixel);
			lpdpCur->lpfColor[3] = float(wxclrPixel.Red())/255.0f;
			lpdpCur->lpfColor[4] = float(wxclrPixel.Green())/255.0f;
			lpdpCur->lpfColor[5] = float(wxclrPixel.Blue())/255.0f;
			++lpdpCur;
		}
	}
}

void CLICPainter::m_fnAppendCrevasse(double dblWidth)
{
	DirPixel *lpdpCur;
	int i;
	lpdpCur = m_lpdpBuffer;
	for (i = 0; i < m_iHeight; ++i)
    {
        for (i = 0; i < m_iWidth; ++i)
        {
            if (lpdpCur->hhCrevasse != NULL)
            {
                lpdpCur->lpfColor[3] = 0.0;
                lpdpCur->lpfColor[4] = 1.0;
                lpdpCur->lpfColor[5] = 0.0;
            }
            ++lpdpCur;
        }
    }
}

void CLICPainter::m_fnAppendSingularities(double dblRadiusScale)
{
	CMeshBase::CPolyhedron::Vertex_iterator vi;
	double lpdblCoordinates[486], lpdblNormals[486], *lplpdblCoord[3], lpdblPixelCoord[2], lpdblAreaCoord[3], dblSumAreaCoord, *lplpdblNorm[3], lpdblCurrentNorm[3], dblNormLen;
	int iXStart, iXEnd, iYStart, iYEnd, iX, iY, i, j, k, iRow, iCol;
	float fDepth, fPoweredCos;
	DirPixel *lpdpCurrent;
	for (i = 0; i < 162; ++i)
	{
		m_lpGLCanvas->m_fnNormalMap(m_lpGLCanvas->m_lpdblSphereCoordinates + i * 3, lpdblNormals + i * 3);
	}
	for (vi = m_lpGLCanvas->m_lpMeshCenter->m_mshSurface.vertices_begin(); vi != m_lpGLCanvas->m_lpMeshCenter->m_mshSurface.vertices_end(); ++vi)
	{
		if (vi->iBordId == -1 && vi->nNumPorts != 0)
		{
			for (i = 0; i < 162; ++i)
			{
				m_lpGLCanvas->m_fnViewCoordinate(\
					vi->point().x() + m_lpGLCanvas->m_lpMeshCenter->m_ftMeanRadius * dblRadiusScale * m_lpGLCanvas->m_lpdblSphereCoordinates[i * 3], \
					vi->point().y() + m_lpGLCanvas->m_lpMeshCenter->m_ftMeanRadius * dblRadiusScale * m_lpGLCanvas->m_lpdblSphereCoordinates[i * 3 + 1], \
					vi->point().z() + m_lpGLCanvas->m_lpMeshCenter->m_ftMeanRadius * dblRadiusScale * m_lpGLCanvas->m_lpdblSphereCoordinates[i * 3 + 2], \
					lpdblCoordinates[i * 3], lpdblCoordinates[i * 3 + 1], lpdblCoordinates[i * 3 + 2]\
				);
			}
			for (i = 0; i < 320; ++i)
			{
				iXStart = INT_MAX;
				iXEnd = -INT_MAX;
				iYStart = INT_MAX;
				iYEnd = -INT_MAX;
				for (j = 0; j < 3; ++j)
				{
					lplpdblNorm[j] = lpdblNormals + m_lpGLCanvas->m_lpsSphereFacets[i * 3 + j] * 3;
					lplpdblCoord[j] = lpdblCoordinates + m_lpGLCanvas->m_lpsSphereFacets[i * 3 + j] * 3;
					iX = int(floor(lplpdblCoord[j][0]));
					iY = int(floor(lplpdblCoord[j][1]));
					if (iX < iXStart)
					{
						iXStart = iX;
					}
					if (iY < iYStart)
					{
						iYStart = iY;
					}
					iX = int(ceil(lplpdblCoord[j][0]));
					iY = int(ceil(lplpdblCoord[j][1]));
					if (iX > iXEnd)
					{
						iXEnd = iX;
					}
					if (iY > iYEnd)
					{
						iYEnd = iY;
					}
				}
				for (iX = iXStart; iX <= iXEnd; ++iX)
				{
					for (iY = iYStart; iY <= iYEnd; ++iY)
					{
						iRow = iY - m_iYOffset;
						iCol = iX - m_iXOffset;
						if (iRow >= 0 && iRow < m_iHeight && iCol >= 0 && iCol < m_iWidth)
						{
							lpdblPixelCoord[0] = double(iX);
							lpdblPixelCoord[1] = double(iY);
							for (k = 0; k < 3; ++k)
							{
								lpdblAreaCoord[k] = m_fnPlanarArea(lpdblPixelCoord, lplpdblCoord[(k + 1) % 3], lplpdblCoord[(k + 2) % 3]);
							}
							dblSumAreaCoord = lpdblAreaCoord[0] + lpdblAreaCoord[1] + lpdblAreaCoord[2];
							if ((lpdblAreaCoord[0] >= 0 && lpdblAreaCoord[1] >= 0 && lpdblAreaCoord[2] >= 0 && dblSumAreaCoord > 0) ||
								(lpdblAreaCoord[0] <= 0 && lpdblAreaCoord[1] <= 0 && lpdblAreaCoord[2] <= 0 && dblSumAreaCoord < 0))
							{
								lpdblAreaCoord[0] /= dblSumAreaCoord;
								lpdblAreaCoord[1] /= dblSumAreaCoord;
								lpdblAreaCoord[2] /= dblSumAreaCoord;
								lpdblCurrentNorm[0] = lplpdblNorm[0][0] * lpdblAreaCoord[0] + lplpdblNorm[1][0] * lpdblAreaCoord[1] + lplpdblNorm[2][0] * lpdblAreaCoord[2];
								lpdblCurrentNorm[1] = lplpdblNorm[0][1] * lpdblAreaCoord[0] + lplpdblNorm[1][1] * lpdblAreaCoord[1] + lplpdblNorm[2][1] * lpdblAreaCoord[2];
								lpdblCurrentNorm[2] = lplpdblNorm[0][2] * lpdblAreaCoord[0] + lplpdblNorm[1][2] * lpdblAreaCoord[1] + lplpdblNorm[2][2] * lpdblAreaCoord[2];
								dblNormLen = sqrt(lpdblCurrentNorm[0] * lpdblCurrentNorm[0] + lpdblCurrentNorm[1] * lpdblCurrentNorm[1] + lpdblCurrentNorm[2] * lpdblCurrentNorm[2]);
								lpdblCurrentNorm[0] /= dblNormLen;
								lpdblCurrentNorm[1] /= dblNormLen;
								lpdblCurrentNorm[2] /= dblNormLen;
								fDepth = float(lplpdblCoord[0][2] * lpdblAreaCoord[0] + lplpdblCoord[1][2] * lpdblAreaCoord[1] + lplpdblCoord[2][2] * lpdblAreaCoord[2]);
								lpdpCurrent = m_lpdpBuffer + (iY - m_iYOffset) * m_iWidth + (iX - m_iXOffset);
								if (lpdpCurrent->hhLocation == NULL || lpdpCurrent->fDepth > fDepth)
								{
									lpdpCurrent->fDepth = fDepth;
									fPoweredCos = pow(float(lpdblCurrentNorm[2]), 4.0f);
									if (vi->nNumPorts < 4)
									{
										lpdpCurrent->lpfColor[3] = 1.0f;
										lpdpCurrent->lpfColor[4] = fPoweredCos * 0.5f + 0.5f;
										lpdpCurrent->lpfColor[5] = fPoweredCos * 0.5f + 0.5f;
									}
									else
									{
										lpdpCurrent->lpfColor[3] = fPoweredCos * 0.5f + 0.5f;
										lpdpCurrent->lpfColor[4] = fPoweredCos * 0.5f + 0.5f;
										lpdpCurrent->lpfColor[5] = 1.0f;
									}
									if (lpdpCurrent->hhLocation == NULL)
									{
										lpdpCurrent->hhLocation = vi->halfedge();
									}
								}
							}
						}
					}
				}
			}
		}
	}
}

void CLICPainter::m_fnMarkForceToDraw()
{
	CMeshBase::CPolyhedron::Vertex_iterator vi;
	CMeshBase::CPolyhedron::Edge_iterator ei;
	CMeshBase::CPolyhedron::Halfedge_handle hh0, hh, hhPath;
	CMeshBase::FT ftMinDot, ftDot;
	for (ei = m_lpGLCanvas->m_lpMeshCenter->m_mshSurface.edges_begin(); ei != m_lpGLCanvas->m_lpMeshCenter->m_mshSurface.edges_end(); ++ei)
	{
		ei->ftForce = 0;
		ei->opposite()->ftForce = 0;
	}
	for (vi = m_lpGLCanvas->m_lpMeshCenter->m_mshSurface.vertices_begin(); vi != m_lpGLCanvas->m_lpMeshCenter->m_mshSurface.vertices_end(); ++vi)
	{
		vi->vtViewForce = CMeshBase::Vector_3(0, 0, 0);
		vi->iTempIndex = -1;
		if (vi->iDegreeBias != 0)
		{
			hh0 = vi->halfedge();
			hh = hh0;
			ftMinDot = 0;
			hhPath = NULL;
			do
			{
				if (hh->prev()->vertex()->iDegreeBias * vi->iDegreeBias < 0)
				{
					ftDot = hh->vtVector * vi->vtThetaForce / hh->ftLen;
					if (ftDot < ftMinDot && ftDot < 0.0)
					{
						ftMinDot = ftDot;
						hhPath = hh;
					}
				}
				hh = hh->next()->opposite();
			} while (hh != hh0);
			if (hhPath == NULL)
			{
				do
				{
					ftDot = hh->vtVector * vi->vtThetaForce / hh->ftLen;
					if (ftDot < ftMinDot)
					{
						ftMinDot = ftDot;
						hhPath = hh;
					}
					hh = hh->next()->opposite();
				} while (hh != hh0);
			}
			if (hhPath != NULL)
			{
				if (vi->iDegreeBias < 0)
				{
					hhPath->ftForce += ftMinDot;
					hhPath->opposite()->ftForce -= ftMinDot;
				}
				else
				{
					hhPath->ftForce -= ftMinDot;
					hhPath->opposite()->ftForce += ftMinDot;
				}
			}
		}
	}
	for (ei = m_lpGLCanvas->m_lpMeshCenter->m_mshSurface.edges_begin(); ei != m_lpGLCanvas->m_lpMeshCenter->m_mshSurface.edges_end(); ++ei)
	{
		if (ei->ftForce != 0.0)
		{
			if (fabs(ei->ftForce) > fabs(ei->ftForceBarrier))
			{
				if ((ei->prev()->vertex()->iDegreeBias * ei->vertex()->iDegreeBias) < 0)
				{
					if (ei->prev()->vertex()->iTempIndex == -1)
					{
						ei->prev()->vertex()->vtViewForce = ei->vtVector;
						ei->prev()->vertex()->iTempIndex = 0;
					}
					else
					{
						ei->prev()->vertex()->vtViewForce = ei->prev()->vertex()->vtThetaForce;
					}
					if (ei->vertex()->iTempIndex == -1)
					{
						ei->vertex()->vtViewForce = -ei->vtVector;
						ei->vertex()->iTempIndex = 0;
					}
					else
					{
						ei->vertex()->vtViewForce = ei->vertex()->vtThetaForce;
					}
				}
				else 
				{
					if (ei->prev()->vertex()->iDegreeBias != 0 && ei->prev()->vertex()->iTempIndex == -1)
					{
						ei->prev()->vertex()->vtViewForce = ei->prev()->vertex()->vtThetaForce;
					}
					if (ei->vertex()->iDegreeBias != 0 && ei->vertex()->iTempIndex == -1)
					{
						ei->vertex()->vtViewForce = ei->vertex()->vtThetaForce;
					}
				}
			}
		}
	}
}

void CLICPainter::m_fnDrawForce()
{
	int i, j;
	CMeshBase::Point_3 lpptArrow[2], lpptViewArrow[2], ptRatio, lpptAngles[2];
	CMeshBase::FT ftForceLen;
	CMeshBase::CPolyhedron::Vertex_iterator vi;
	DirPixel* lpdpCur;
	wxBitmap wxImgBmp(m_iWidth, m_iHeight);
	wxColour wxclrPixel;
	wxMemoryDC wxImgDC(wxImgBmp);
	wxPen wxpnPixel, wxpnArrow = wxPen(*wxBLACK, 3);
	wxBrush wxbrArrow = *wxBLACK_BRUSH;
	wxPoint lpwxpntTriangle[3];
	lpdpCur = m_lpdpBuffer;
	for (i = 0; i < m_iHeight; ++i)
	{
		for (j = 0; j < m_iWidth; ++j)
		{
			wxpnPixel.SetColour(unsigned char(lpdpCur->lpfColor[3] * 255), unsigned char(lpdpCur->lpfColor[4] * 255), unsigned char(lpdpCur->lpfColor[5] * 255));
			wxImgDC.SetPen(wxpnPixel);
			wxImgDC.DrawPoint(j, i);
			++lpdpCur;
		}
	}
	wxImgDC.SetPen(wxpnArrow);
	wxImgDC.SetBrush(wxbrArrow);
	for (vi = m_lpGLCanvas->m_lpMeshCenter->m_mshSurface.vertices_begin(); vi != m_lpGLCanvas->m_lpMeshCenter->m_mshSurface.vertices_end(); ++vi)
	{
		if (vi->iDegreeBias != 0 && CGAL::squared_length(vi->vtViewForce) > 1.0e-6)
		{
			ftForceLen = sqrt(CGAL::squared_length(vi->vtViewForce));
			lpptArrow[0] = vi->point() - vi->vtViewForce * (m_lpGLCanvas->m_lpMeshCenter->m_ftMeanRadius * 2 / ftForceLen);
			lpptArrow[1] = vi->point() + vi->vtViewForce * (m_lpGLCanvas->m_lpMeshCenter->m_ftMeanRadius * 2 / ftForceLen);
			m_lpGLCanvas->m_fnViewCoordinate(lpptArrow[0].x(), lpptArrow[0].y(), lpptArrow[0].z(), lpptViewArrow[0].m_tX, lpptViewArrow[0].m_tY, lpptViewArrow[0].m_tZ);
			m_lpGLCanvas->m_fnViewCoordinate(lpptArrow[1].x(), lpptArrow[1].y(), lpptArrow[1].z(), lpptViewArrow[1].m_tX, lpptViewArrow[1].m_tY, lpptViewArrow[1].m_tZ);
			ptRatio = CGAL::barycenter(lpptViewArrow[0], 0.25, lpptViewArrow[1], 0.75);
			lpptAngles[0] = ptRatio + CMeshBase::Vector_3(lpptViewArrow[0].y() - lpptViewArrow[1].y(), lpptViewArrow[1].x() - lpptViewArrow[0].x(), 0.0) * (1.0 / 8.0);
			lpptAngles[1] = ptRatio + CMeshBase::Vector_3(lpptViewArrow[1].y() - lpptViewArrow[0].y(), lpptViewArrow[0].x() - lpptViewArrow[1].x(), 0.0) * (1.0 / 8.0);
			lpwxpntTriangle[0] = wxPoint(int(lpptAngles[0].x() - m_iXOffset + 0.5), int(lpptAngles[0].y() - m_iYOffset + 0.5));
			lpwxpntTriangle[1] = wxPoint(int(lpptAngles[1].x() - m_iXOffset + 0.5), int(lpptAngles[1].y() - m_iYOffset + 0.5));
			lpwxpntTriangle[2] = wxPoint(int(lpptViewArrow[1].x() - m_iXOffset + 0.5), int(lpptViewArrow[1].y() - m_iYOffset + 0.5));
			wxImgDC.DrawLine(int(lpptViewArrow[0].x() - m_iXOffset + 0.5), int(lpptViewArrow[0].y() - m_iYOffset + 0.5), int(ptRatio.x() - m_iXOffset + 0.5), int(ptRatio.y() - m_iYOffset + 0.5));
			wxImgDC.DrawPolygon(3, lpwxpntTriangle);
		}
	}
	lpdpCur = m_lpdpBuffer;
	for (i = 0; i < m_iHeight; ++i)
	{
		for (j = 0; j < m_iWidth; ++j)
		{
			wxImgDC.GetPixel(j, i, &wxclrPixel);
			lpdpCur->lpfColor[3] = float(wxclrPixel.Red()) / 255.0f;
			lpdpCur->lpfColor[4] = float(wxclrPixel.Green()) / 255.0f;
			lpdpCur->lpfColor[5] = float(wxclrPixel.Blue()) / 255.0f;
			++lpdpCur;
		}
	}
}

void CLICPainter::m_fnWriteBitmap(wxBitmap &bmpLIC, int &iBmpX, int &iBmpY)
{
	int i, j;
	wxPen pnColoredPen;
	DirPixel *lpdpCur;
	iBmpX = m_iXOffset;
	iBmpY = m_iYOffset;
	bmpLIC.Create(m_iWidth, m_iHeight);
	wxMemoryDC memdc;
	memdc.SelectObject(bmpLIC);
	memdc.SetBackground(*wxWHITE_BRUSH);
	memdc.Clear();
	lpdpCur = m_lpdpBuffer;
	for (i = 0; i < m_iHeight; ++i)
	{
		for (j = 0; j < m_iWidth; ++j)
		{
			if (lpdpCur->hhLocation != NULL)
			{
				pnColoredPen.SetColour(char(lpdpCur->lpfColor[3] * 255), char(lpdpCur->lpfColor[4] * 255), char(lpdpCur->lpfColor[5] * 255));
				memdc.SetPen(pnColoredPen);
				memdc.DrawPoint(j, i);
			}
			++lpdpCur;
		}
	}
}

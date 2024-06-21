#include "CGLCanvas.h"

CGLCanvas::CGLCanvas(wxWindow *lpParent, int *lpAttribList) :wxGLCanvas(lpParent, wxID_ANY, lpAttribList), m_Timer(this)
{
    //ctor
	int i, j;
	PixelDataConvertor pdc;
	wxBitmap bmpSmallSphere;
	wxImage imgSmallSphere;
	m_lpGLContext = new CGLContext(this);
	std::ifstream ifs;
    m_dblNear = 0.125;
    m_dblDistance = 2;
    m_dblFar = 32;
    m_dblTanHalfFov = tan(56.1/360*3.14159265);
    m_dblRotateAngle = 0;
    m_lpdblRotateAxis[0] = 0;
    m_lpdblRotateAxis[1] = 0;
    m_lpdblRotateAxis[2] = 1;
    m_bFirstTime = true;
    m_lpMeshCenter = new CMeshInfoReserver;

	m_bShowEdge = true;
	m_bChartViewable = false;
	m_bFrameField = false;
    m_bSingularityViewable = false;
    //m_bElectricFieldViewable = false;
	m_iElectricFieldViewable = 0;
	m_bTrackViewable = false;
	m_iFacetJetColor = 0;

	m_bSelectPoint = false;
	m_bSelectTrack = false;
	m_bSelectCrevasse = false;
	m_bSelectFacet = false;
	m_bSelectEdge = false;

	m_iViewableLocalPara = -1;
    m_iViewableColor = 0;
    m_iViewableMeshInfo = -1;
    m_iViewableRay = 0;
    m_iViewableTriangle = 0;
    m_bProcessing = false;
	m_bMovieMode = false;
	m_bDiscrete = true;
	m_bRecordMode = false;
	m_bShowLIC = false;

	//wxImage imgSmallSphere(wxSize(44, 44));
	//ifs.open("SmallSphere.txt");
	//for (i = 0; i < 486; ++i)
	//{
	//	ifs >> m_lpdblSphereCoordinates[i];
	//	pdc.lpCoordinate[i] = m_lpdblSphereCoordinates[i];
	//}
	//for (i = 0; i < 960; ++i)
	//{
	//	ifs >> m_lpsSphereFacets[i];
	//	pdc.lpsIndex[i] = m_lpsSphereFacets[i];
	//}
	//ifs.close();
	//for (i = 0; i < 44; ++i)
	//{
	//	for (j = 0; j < 44; ++j)
	//	{
	//		imgSmallSphere.SetRGB(i, j, pdc.lpPixels[i * 44 * 3 + j * 3], pdc.lpPixels[i * 44 * 3 + j * 3 + 1], pdc.lpPixels[i * 44 * 3 + j * 3 + 2]);
	//	};
	//}
	//imgSmallSphere.SaveFile("SmallSphere.bmp", wxBITMAP_TYPE_BMP);

	bmpSmallSphere = wxBITMAP(SMALL_SPHERE_BMP);
	imgSmallSphere = bmpSmallSphere.ConvertToImage();
	for (i = 0; i < 44; ++i)
	{
		for (j = 0; j < 44; ++j)
		{
			pdc.lpPixels[i * 44 * 3 + j * 3] = imgSmallSphere.GetRed(i, j);
			pdc.lpPixels[i * 44 * 3 + j * 3 + 1] = imgSmallSphere.GetGreen(i, j);
			pdc.lpPixels[i * 44 * 3 + j * 3 + 2] = imgSmallSphere.GetBlue(i, j);
		}
	}
	for (i = 0; i < 486; ++i)
	{
		m_lpdblSphereCoordinates[i] = pdc.lpCoordinate[i];
	}
	for (i = 0; i < 960; ++i)
	{
		m_lpsSphereFacets[i] = pdc.lpsIndex[i];
	}

	m_bOrthodox = false;
	m_bmpOriginal = wxBITMAP(SHADER_BMP);
	m_iBmpX = 4096;
	m_iBmpY = 4096;
}

CGLCanvas::~CGLCanvas()
{
	delete m_lpMeshCenter;
	delete m_lpGLContext;
}

void CGLCanvas::m_fnRotate()
{
    CGeoPack::fnFillRotateMatrix(m_lpdblWorldMat, m_lpdblRotateAxis, m_dblRotateAngle);
}

void CGLCanvas::m_fnObserve()
{
    int i;
    memset(m_lpdblViewMat, 0, 16 * sizeof(double));
    for(i = 0; i < 4; ++i)
    {
        m_lpdblViewMat[i * 4 + i] = 1;
    }
    m_lpdblViewMat[14] = -m_dblDistance;
}

void CGLCanvas::m_fnCast()
{
    memset(m_lpdblProjectMat, 0, 16 * sizeof(double));
    m_lpdblProjectMat[0] = double(m_iHalfHeight) / double(m_iHalfWidth);
    m_lpdblProjectMat[5] = 1;
	m_lpdblProjectMat[10] = -m_dblFar * m_dblTanHalfFov / (m_dblFar - m_dblNear);
	m_lpdblProjectMat[11] = -m_dblTanHalfFov;//0;
	m_lpdblProjectMat[14] = -m_dblNear * m_dblFar * m_dblTanHalfFov / (m_dblFar - m_dblNear);
	m_lpdblProjectMat[15] = 0;
	if (m_bOrthodox)
	{
		m_lpdblProjectMat[10] = -m_dblDistance * m_dblTanHalfFov / (m_dblFar - m_dblNear);
		m_lpdblProjectMat[11] = 0;
		m_lpdblProjectMat[14] = -m_dblNear * m_dblDistance * m_dblTanHalfFov / (m_dblFar - m_dblNear);
		m_lpdblProjectMat[15] = m_dblDistance * m_dblTanHalfFov;
	}
}



void CGLCanvas::OnTimer(wxTimerEvent &event)
{
	if (!m_bProcessing)
	{
		m_Timer.Stop();
	}
	Refresh();
	//m_Monitor.m_Canvas.Refresh();
}


void CGLCanvas::OnMouseEvent(wxMouseEvent &event)
{
    int xMouse, yMouse;
	CMeshBase::CPolyhedron::Vertex_iterator vi;
	CMeshBase::CSamplerPolyhedron::Vertex_iterator svi;
	CMeshBase::CPolyhedron::Edge_iterator ei;
	CMeshBase::CPolyhedron::Facet_iterator fi;
	CMeshBase::CPolyhedron::Halfedge_handle hh0, hh;
	int iXDif, iYDif, iArea, iViewX0, iViewY0, iViewX1, iViewY1, i;
	double dblMinDepth, dblDepth, lpdblAreaCoord[3], lpdblVector[6], dblSumArea;
	int iMinSqDist, iSqDist;
	CMeshBase::CPolyhedron::Vertex_handle vhSelect;
	CMeshBase::CSamplerPolyhedron::Vertex_handle svhSelect;
	CMeshBase::CPolyhedron::Halfedge_handle hhSelect;
	CMeshBase::CPolyhedron::Facet_handle fhSelect;
	std::list<CTracer::CInducedPath>::iterator iterIP;
	std::list<CTracer::CPathElement>::iterator iterPE;
	CMeshBase::Point_3 ptPoint;
    if(event.LeftDown())
    {
        m_xMouse = event.GetX();
        m_yMouse = event.GetY();
        m_bLButtonDown = true;
		if (m_bSelectPoint)
		{
			dblMinDepth = 1.0;
			iMinSqDist = 25;
			if (m_iViewableMeshInfo == 3)
			{
				svhSelect = NULL;
				if (m_lpMeshCenter->m_mshRemeshedSurface.size_of_vertices() > 0)
				{
					for (svi = m_lpMeshCenter->m_mshRemeshedSurface.vertices_begin(); svi != m_lpMeshCenter->m_mshRemeshedSurface.vertices_end(); ++svi)
					{
						m_fnPositionOnScreen(svi->point().x(), svi->point().y(), svi->point().z(), svi->iViewX, svi->iViewY, svi->ftDepth);
						iXDif = svi->iViewX - m_xMouse;
						iYDif = svi->iViewY - m_yMouse;
						iSqDist = iXDif * iXDif + iYDif * iYDif;
						if (iSqDist < iMinSqDist)
						{
							if (svi->ftDepth < dblMinDepth)
							{
								dblMinDepth = svi->ftDepth;
								svhSelect = &(*svi);
							}
						}
						svi->iMark = 0;
					}
				}
				if (svhSelect != NULL)
				{
					svhSelect->iMark = 1;
				}
			}
			else
			{
				vhSelect = NULL;
				for (vi = m_lpMeshCenter->m_mshSurface.vertices_begin(); vi != m_lpMeshCenter->m_mshSurface.vertices_end(); ++vi)
				{
					m_fnPositionOnScreen(vi->point().x(), vi->point().y(), vi->point().z(), vi->iViewX, vi->iViewY, vi->ftDepth);
					iXDif = vi->iViewX - m_xMouse;
					iYDif = vi->iViewY - m_yMouse;
					iSqDist = iXDif * iXDif + iYDif * iYDif;
					if (iSqDist < iMinSqDist)
					{
						if (vi->ftDepth < dblMinDepth)
						{
							dblMinDepth = vi->ftDepth;
							vhSelect = &(*vi);
						}
					}
				}
				if (vhSelect != NULL)
				{
					for (vi = m_lpMeshCenter->m_mshSurface.vertices_begin(); vi != m_lpMeshCenter->m_mshSurface.vertices_end(); ++vi)
					{
						if (vi->iSelected == 1)
						{
							vi->iSelected = 0;
						}
						else if (vi->iSelected == 2)
						{
							vi->iSelected = 1;
						}
					}
					vhSelect->iSelected = 2;
				}
			}
		}
		else if (m_bSelectCrevasse)
		{
			for (vi = m_lpMeshCenter->m_mshSurface.vertices_begin(); vi != m_lpMeshCenter->m_mshSurface.vertices_end(); ++vi)
			{
				m_fnPositionOnScreen(vi->point().x(), vi->point().y(), vi->point().z(), vi->iViewX, vi->iViewY, vi->ftDepth);
			}
			dblMinDepth = 1.0;
			iMinSqDist = 25;
			hhSelect = NULL;
			for (ei = m_lpMeshCenter->m_mshSurface.edges_begin(); ei != m_lpMeshCenter->m_mshSurface.edges_end(); ++ei)
			{
				if (ei->bFreeEdge && !ei->is_border_edge() && !ei->bSelected)
				{
					iArea = m_xMouse * ei->prev()->vertex()->iViewY + ei->prev()->vertex()->iViewX * ei->vertex()->iViewY + ei->vertex()->iViewX * m_yMouse\
						- m_yMouse * ei->prev()->vertex()->iViewX - ei->prev()->vertex()->iViewY * ei->vertex()->iViewX - ei->vertex()->iViewY * m_xMouse;
					iXDif = ei->vertex()->iViewX - ei->prev()->vertex()->iViewX;
					iYDif = ei->vertex()->iViewY - ei->prev()->vertex()->iViewY;
					if (iArea * iArea < (iXDif * iXDif + iYDif * iYDif) * iMinSqDist && (iXDif * (m_xMouse - ei->vertex()->iViewX) + iYDif * (m_yMouse - ei->vertex()->iViewY)) < 0 && (iXDif * (m_xMouse - ei->prev()->vertex()->iViewX) + iYDif * (m_yMouse - ei->prev()->vertex()->iViewY)) > 0)
					{
						dblDepth = (ei->vertex()->ftDepth + ei->prev()->vertex()->ftDepth) / 2;
						if (dblDepth > 0.0 && dblDepth < dblMinDepth)
						{
							hhSelect = &(*ei);
							dblMinDepth = dblDepth;
						}
					}
				}
			}
			if (hhSelect != NULL)
			{
				hhSelect->bSelected = true;
				hhSelect->opposite()->bSelected = true;
				m_lpMeshCenter->m_fnExtendPath(hhSelect);
				m_lpMeshCenter->m_fnExtendPath(hhSelect->opposite());
			}
		}
		else if (m_bSelectTrack && m_bTrackViewable)
		{
			iMinSqDist = 25;
			if (!m_lpMeshCenter->m_Tracer.m_lstIP.empty())
			{
				for (iterIP = m_lpMeshCenter->m_Tracer.m_lstIP.begin(); iterIP != m_lpMeshCenter->m_Tracer.m_lstIP.end(); ++iterIP)
				{
					iterIP->m_bSelected = false;
					if (iterIP->m_lstData.size() > 1)
					{
						iterPE = iterIP->m_lstData.begin();
						ptPoint = CGAL::barycenter(iterPE->hhEdge->vertex()->point(), iterPE->ftProportion, iterPE->hhEdge->prev()->vertex()->point(), 1.0 - iterPE->ftProportion);
						m_fnPositionOnScreen(ptPoint.x(), ptPoint.y(), ptPoint.z(), iViewX1, iViewY1, dblDepth);
						iViewX0 = iViewX1;
						iViewY0 = iViewY1;
						++iterPE;
						while (iterPE != iterIP->m_lstData.end() && !iterIP->m_bSelected)
						{
							ptPoint = CGAL::barycenter(iterPE->hhEdge->vertex()->point(), iterPE->ftProportion, iterPE->hhEdge->prev()->vertex()->point(), 1.0 - iterPE->ftProportion);
							m_fnPositionOnScreen(ptPoint.x(), ptPoint.y(), ptPoint.z(), iViewX1, iViewY1, dblDepth);
							iArea = m_xMouse * iViewY0 + iViewX0 * iViewY1 + iViewX1 * m_yMouse - m_yMouse * iViewX0 - iViewY0 * iViewX1 - iViewY1 * m_xMouse;
							iXDif = iViewX1 - iViewX0;
							iYDif = iViewY1 - iViewY0;
							if (iArea * iArea < (iXDif * iXDif + iYDif * iYDif) * iMinSqDist && (iXDif * (m_xMouse - iViewX1) + iYDif * (m_yMouse - iViewY1)) < 0 && (iXDif * (m_xMouse - iViewX0) + iYDif * (m_yMouse - iViewY0)) > 0)
							{
								iterIP->m_bSelected = true;
							}
							iViewX0 = iViewX1;
							iViewY0 = iViewY1;
							++iterPE;
						}
					}
				}
			}
		}
		else if (m_bSelectFacet)
		{
			for (vi = m_lpMeshCenter->m_mshSurface.vertices_begin(); vi != m_lpMeshCenter->m_mshSurface.vertices_end(); ++vi)
			{
				m_fnPositionOnScreen(vi->point().x(), vi->point().y(), vi->point().z(), vi->iViewX, vi->iViewY, vi->ftDepth);
			}
			dblMinDepth = 1.0;
			fhSelect = NULL;
			for (fi = m_lpMeshCenter->m_mshSurface.facets_begin(); fi != m_lpMeshCenter->m_mshSurface.facets_end(); ++fi)
			{
				hh0 = fi->halfedge();
				hh = hh0;
				for (i = 0; i < 3; ++i)
				{
					lpdblVector[i * 2] = double(hh->vertex()->iViewX - m_xMouse);
					lpdblVector[i * 2 + 1] = double(hh->vertex()->iViewY - m_yMouse);
					hh = hh->next();
				}
				dblSumArea = 0;
				for (i = 0; i < 3; ++i)
				{
					lpdblAreaCoord[i] = lpdblVector[i * 2] * lpdblVector[(i + 1) % 3 * 2 + 1] - lpdblVector[i * 2 + 1] * lpdblVector[(i + 1) % 3 * 2];
					dblSumArea += lpdblAreaCoord[i];
				}
				for (i = 0; i < 3; ++i)
				{
					lpdblAreaCoord[i] /= dblSumArea;
				}
				if (lpdblAreaCoord[0] > 0 && lpdblAreaCoord[1] > 0 && lpdblAreaCoord[2] > 0)
				{
					dblDepth = hh0->prev()->vertex()->ftDepth * lpdblAreaCoord[0] + hh0->vertex()->ftDepth * lpdblAreaCoord[1] + hh0->next()->vertex()->ftDepth * lpdblAreaCoord[2];
					if (dblDepth > 0.0 && dblDepth < dblMinDepth)
					{
						dblMinDepth = dblDepth;
						fhSelect = &(*fi);
					}
				}
			}
			if (fhSelect != NULL)
			{
				fhSelect->bSelected = !fhSelect->bSelected;
				m_lpMeshCenter->m_fhFacetFix = fhSelect;
			}
		}
		else if (m_bSelectEdge)
		{
			for (vi = m_lpMeshCenter->m_mshSurface.vertices_begin(); vi != m_lpMeshCenter->m_mshSurface.vertices_end(); ++vi)
			{
				m_fnPositionOnScreen(vi->point().x(), vi->point().y(), vi->point().z(), vi->iViewX, vi->iViewY, vi->ftDepth);
			}
			dblMinDepth = 1.0;
			iMinSqDist = 1;
			hhSelect = NULL;
			for (ei = m_lpMeshCenter->m_mshSurface.edges_begin(); ei != m_lpMeshCenter->m_mshSurface.edges_end(); ++ei)
			{
				if (ei->vertex()->iViewX < m_iHalfWidth * 2 && ei->vertex()->iViewX > 0 && ei->vertex()->iViewY < m_iHalfHeight * 2 && ei->vertex()->iViewY > 0 &&\
					ei->prev()->vertex()->iViewX < m_iHalfWidth * 2 && ei->prev()->vertex()->iViewX > 0 && ei->prev()->vertex()->iViewY < m_iHalfHeight * 2 && ei->prev()->vertex()->iViewY > 0 &&\
					ei->vertex()->ftDepth > 0.0 && ei->vertex()->ftDepth < 1.0 && ei->prev()->vertex()->ftDepth > 0.0 && ei->prev()->vertex()->ftDepth < 1.0)
				{
					iArea = m_xMouse * ei->prev()->vertex()->iViewY + ei->prev()->vertex()->iViewX * ei->vertex()->iViewY + ei->vertex()->iViewX * m_yMouse\
						- m_yMouse * ei->prev()->vertex()->iViewX - ei->prev()->vertex()->iViewY * ei->vertex()->iViewX - ei->vertex()->iViewY * m_xMouse;
					iXDif = ei->vertex()->iViewX - ei->prev()->vertex()->iViewX;
					iYDif = ei->vertex()->iViewY - ei->prev()->vertex()->iViewY;
					if (iArea * iArea < (iXDif * iXDif + iYDif * iYDif) * iMinSqDist && (iXDif * (m_xMouse - ei->vertex()->iViewX) + iYDif * (m_yMouse - ei->vertex()->iViewY)) < 0 && (iXDif * (m_xMouse - ei->prev()->vertex()->iViewX) + iYDif * (m_yMouse - ei->prev()->vertex()->iViewY)) > 0)
					{
						dblDepth = (ei->vertex()->ftDepth + ei->prev()->vertex()->ftDepth) / 2;
						if (dblDepth < dblMinDepth)
						{
							hhSelect = &(*ei);
							dblMinDepth = dblDepth;
						}
					}
				}
			}
			if (hhSelect != NULL)
			{
				hhSelect->bSelected = !hhSelect->bSelected;
				hhSelect->opposite()->bSelected = hhSelect->bSelected;
				m_lpMeshCenter->m_hhDirEdge = hhSelect;
			}
		}
		Refresh();
    }
    if(event.LeftUp())
    {
        m_xMouse = event.GetX();
        m_yMouse = event.GetY();
        m_bLButtonDown = false;
    }
    if(event.Dragging() && m_bLButtonDown)
    {
        xMouse = event.GetX();
        yMouse = event.GetY();
        int xOffset1 = m_xMouse - m_iHalfWidth;
        int yOffset1 = m_iHalfHeight - m_yMouse;
        int xOffset2 = xMouse - m_iHalfWidth;
        int yOffset2 = m_iHalfHeight - yMouse;
        double lpdblVec1[3] = {double(xOffset1), double(yOffset1), sqrt(std::max<double>(0.0, double(m_iHalfHeight * m_iHalfHeight - xOffset1 * xOffset1 - yOffset1 * yOffset1)))};
        double lpdblVec2[3] = {double(xOffset2), double(yOffset2), sqrt(std::max<double>(0.0, double(m_iHalfHeight * m_iHalfHeight - xOffset2 * xOffset2 - yOffset2 * yOffset2)))};
        double lpdblCurrentAxis[3] = {0, 0, 1};
        double dblCurrentRotateAngle = 0.0;
        double lpdblNewAxis[3] = {0, 0, 1};
        double dblNewRotateAngle = 0.0;
        CGeoPack::fnMinRotation(lpdblCurrentAxis, dblCurrentRotateAngle, lpdblVec1, lpdblVec2);
        CGeoPack::fnMergeRotation(lpdblNewAxis, dblNewRotateAngle, m_lpdblRotateAxis, m_dblRotateAngle, lpdblCurrentAxis, dblCurrentRotateAngle);
        m_dblRotateAngle = dblNewRotateAngle;
        m_lpdblRotateAxis[0] = lpdblNewAxis[0];
        m_lpdblRotateAxis[1] = lpdblNewAxis[1];
        m_lpdblRotateAxis[2] = lpdblNewAxis[2];
        m_fnRotate();
        m_xMouse = xMouse;
        m_yMouse = yMouse;
        Refresh();
    }
    if(event.RightDown())
    {
		for (vi = m_lpMeshCenter->m_mshSurface.vertices_begin(); vi != m_lpMeshCenter->m_mshSurface.vertices_end(); ++vi)
		{
			vi->iSelected = 0;
		}
		for (ei = m_lpMeshCenter->m_mshSurface.edges_begin(); ei != m_lpMeshCenter->m_mshSurface.edges_end(); ++ei)
		{
			ei->bSelected = false;
			ei->opposite()->bSelected = false;
		}
		if (!m_lpMeshCenter->m_Tracer.m_lstIP.empty())
		{
			for (iterIP = m_lpMeshCenter->m_Tracer.m_lstIP.begin(); iterIP != m_lpMeshCenter->m_Tracer.m_lstIP.end(); ++iterIP)
			{
				iterIP->m_bSelected = false;
			}
		}
		for (fi = m_lpMeshCenter->m_mshSurface.facets_begin(); fi != m_lpMeshCenter->m_mshSurface.facets_end(); ++fi)
		{
			fi->bSelected = false;
		}
		if (m_iViewableMeshInfo == 3)
		{
			if (m_lpMeshCenter->m_mshRemeshedSurface.size_of_vertices() > 0)
			{
				for (svi = m_lpMeshCenter->m_mshRemeshedSurface.vertices_begin(); svi != m_lpMeshCenter->m_mshRemeshedSurface.vertices_end(); ++svi)
				{
					svi->iMark = 0;
				}
			}
		}
		m_lpMeshCenter->m_hhDirEdge = NULL;
		m_lpMeshCenter->m_fhFacetFix = NULL;
		Refresh();
    }
}


void CGLCanvas::OnKeyEvent(wxKeyEvent &event)
{
    wxKeyCode keyCode;
    keyCode = (wxKeyCode)event.GetKeyCode();
	CMeshBase::CPolyhedron::Vertex_handle vhMaximum, vhMinimum;
	CMeshBase::CPolyhedron::Vertex_iterator vi;
	CMeshBase::CPolyhedron::Facet_iterator fi;
	CMeshBase::CPolyhedron::Halfedge_handle hh0, hh;
    switch(keyCode)
	{

	case WXK_DELETE:
		{
			m_lpMeshCenter->m_fnClearSelectedPath();
			m_lpMeshCenter->m_Tracer.m_fnDeletePath();
			Refresh();
		}
		break;
	case WXK_HOME:
		{
		    if(m_lpMeshCenter->m_bExtremaSelecting && m_lpMeshCenter->m_vecvhDensityExtrema.size() - m_lpMeshCenter->m_iDensityMaximum > 1)
		    {
				++(m_lpMeshCenter->m_iDensityMaximum);
				vhMaximum = m_lpMeshCenter->m_vecvhDensityExtrema[m_lpMeshCenter->m_iDensityMaximum];
				vhMinimum = m_lpMeshCenter->m_vecvhDensityExtrema[m_lpMeshCenter->m_iDensityMinimum];
				m_lpMeshCenter->m_fnClearFreeEdge();
				m_lpMeshCenter->m_fnDetectLocalPath(vhMaximum, vhMinimum);
				Refresh();
		    };
		}
		break;
	case WXK_END:
		{
            if(m_lpMeshCenter->m_bExtremaSelecting && m_lpMeshCenter->m_iDensityMaximum - m_lpMeshCenter->m_iDensityMinimum > 1)
            {
                --(m_lpMeshCenter->m_iDensityMaximum);
				vhMaximum = m_lpMeshCenter->m_vecvhDensityExtrema[m_lpMeshCenter->m_iDensityMaximum];
				vhMinimum = m_lpMeshCenter->m_vecvhDensityExtrema[m_lpMeshCenter->m_iDensityMinimum];
				m_lpMeshCenter->m_fnClearFreeEdge();
				m_lpMeshCenter->m_fnDetectLocalPath(vhMaximum, vhMinimum);
				Refresh();
            }
		}
		break;
	case WXK_PAGEUP:
		{
            if(m_lpMeshCenter->m_bExtremaSelecting && m_lpMeshCenter->m_iDensityMaximum - m_lpMeshCenter->m_iDensityMinimum > 1)
            {
                ++(m_lpMeshCenter->m_iDensityMinimum);
				vhMaximum = m_lpMeshCenter->m_vecvhDensityExtrema[m_lpMeshCenter->m_iDensityMaximum];
				vhMinimum = m_lpMeshCenter->m_vecvhDensityExtrema[m_lpMeshCenter->m_iDensityMinimum];
				m_lpMeshCenter->m_fnClearFreeEdge();
				m_lpMeshCenter->m_fnDetectLocalPath(vhMaximum, vhMinimum);
				Refresh();
            }
		}
		break;
	case WXK_PAGEDOWN:
		{
            if(m_lpMeshCenter->m_bExtremaSelecting && m_lpMeshCenter->m_iDensityMinimum > 0)
            {
                --(m_lpMeshCenter->m_iDensityMinimum);
				vhMaximum = m_lpMeshCenter->m_vecvhDensityExtrema[m_lpMeshCenter->m_iDensityMaximum];
				vhMinimum = m_lpMeshCenter->m_vecvhDensityExtrema[m_lpMeshCenter->m_iDensityMinimum];
				m_lpMeshCenter->m_fnClearFreeEdge();
				m_lpMeshCenter->m_fnDetectLocalPath(vhMaximum, vhMinimum);
				Refresh();
            };
		}
		break;
    case WXK_UP:
        {
			m_dblDistance /= 1.2;
			m_dblFar /= 1.2;
			m_dblNear /= 1.2;
            m_fnObserve();
            m_fnCast();
            Refresh();
        }
        break;
    case WXK_DOWN:
        {
			m_dblDistance *= 1.2;
			m_dblFar *= 1.2;
			m_dblNear *= 1.2;
            m_fnObserve();
            m_fnCast();
            Refresh();
        }
        break;
    default:
		if (char(keyCode) == 'L' || char(keyCode) == 'l')
		{
			vhMaximum = NULL;
			vhMinimum = NULL;
			for (vi = m_lpMeshCenter->m_mshSurface.vertices_begin(); vi != m_lpMeshCenter->m_mshSurface.vertices_end(); ++vi)
			{
				if (vi->iSelected == 1)
				{
					vhMaximum = &(*vi);
				}
				else if (vi->iSelected == 2)
				{
					vhMinimum = &(*vi);
				}
			}
			if (vhMaximum != NULL && vhMinimum != NULL)
			{
				m_lpMeshCenter->m_fnDetectLocalPath(vhMaximum, vhMinimum);
				Refresh();
			}
		}
		else if (char(keyCode) == 'E' || char(keyCode) == 'e')
		{
			for (fi = m_lpMeshCenter->m_mshSurface.facets_begin(); fi != m_lpMeshCenter->m_mshSurface.facets_end(); ++fi)
			{
				fi->iTempIndex = -1;
			}
			for (fi = m_lpMeshCenter->m_mshSurface.facets_begin(); fi != m_lpMeshCenter->m_mshSurface.facets_end(); ++fi)
			{
				if (fi->bSelected)
				{
					hh0 = fi->halfedge();
					hh = hh0;
					do
					{
						if (!hh->bFreeEdge && !hh->opposite()->facet()->bSelected)
						{
							hh->opposite()->facet()->iTempIndex = 0;
						}
						hh = hh->next();
					} while (hh != hh0);
				}
			}
			for (fi = m_lpMeshCenter->m_mshSurface.facets_begin(); fi != m_lpMeshCenter->m_mshSurface.facets_end(); ++fi)
			{
				if (fi->iTempIndex == 0)
				{
					fi->bSelected = true;
					fi->iTempIndex = -1;
				}
			}
			Refresh();
		}
		else if (char(keyCode) == 'F' || char(keyCode) == 'f')
		{
			m_lpMeshCenter->m_fnFixDirection();
			Refresh();
		}
		else if (char(keyCode) == 'U' || char(keyCode) == 'u')
		{
			m_lpMeshCenter->m_fnMultiFixDirection();
			Refresh();
		}
		else if ((char(keyCode) == 'V' || char(keyCode) == 'v') && m_lpMeshCenter->m_fhFacetFix != NULL)
		{
			CDialogs::CDirectionConstraintDialog DCD(this);
			if (DCD.ShowModal() == wxID_OK)
			{
				CMeshBase::Vector_3 vtGuide(DCD.m_lpdblVector[0], DCD.m_lpdblVector[1], DCD.m_lpdblVector[2]);
				CMeshBase::FT ftMom;
				ftMom = atan2(m_lpMeshCenter->m_fhFacetFix->vtAuxiliaryAxis * vtGuide, m_lpMeshCenter->m_fhFacetFix->vtPrincipalAxis * vtGuide) - m_lpMeshCenter->m_fhFacetFix->ftChartDir;
				ftMom = ftMom - CMeshBase::FT(int(ftMom / ftHalfPi + 15.5) - 15) * ftHalfPi;
				m_lpMeshCenter->m_fhFacetFix->ftChartDir += ftMom;
				m_lpMeshCenter->m_fhFacetFix->bDirFixed = true;
				Refresh();
			}
		}
		else if ((char(keyCode) == 'R' || char(keyCode) == 'r') && m_lpMeshCenter->m_fhFacetFix != NULL)
		{
			CDialogs::CAngleBiasDialog ABD(this);
			if (ABD.ShowModal() == wxID_OK)
			{
				CMeshBase::FT ftMom;
				ftMom = ABD.m_dblAngleBias / 180.0 * ftPi;
				m_lpMeshCenter->m_fhFacetFix->ftChartDir += ftMom;
				m_lpMeshCenter->m_fhFacetFix->ftChartDir = m_lpMeshCenter->m_fhFacetFix->ftChartDir - CMeshBase::FT(int(m_lpMeshCenter->m_fhFacetFix->ftChartDir / ftTwoPi + 15.5) - 15) * ftTwoPi;
				Refresh();
			}
		}
        break;
    }
}


void CGLCanvas::OnSize(wxSizeEvent &event)
{
    wxRect rectSize;
    rectSize = event.GetSize();
    glViewport(0, 0, rectSize.width, rectSize.height);
    m_iHalfWidth = rectSize.width / 2;
    m_iHalfHeight = rectSize.height / 2;
    m_fnObserve();
    m_fnCast();
    Refresh();
}

void CGLCanvas::m_fnDrawSphere(double dblX, double dblY, double dblZ, double dblRadius, double dblR, double dblG, double dblB)
{
	int i;
	glBegin(GL_TRIANGLES);
	glColor3d(dblR, dblG, dblB);
	for (i = 0; i < 320; ++i)
	{
		glVertex3d(m_lpdblSphereCoordinates[m_lpsSphereFacets[i * 3] * 3] * dblRadius + dblX, m_lpdblSphereCoordinates[m_lpsSphereFacets[i * 3] * 3 + 1] * dblRadius + dblY, m_lpdblSphereCoordinates[m_lpsSphereFacets[i * 3] * 3 + 2] * dblRadius + dblZ);
		glVertex3d(m_lpdblSphereCoordinates[m_lpsSphereFacets[i * 3 + 1] * 3] * dblRadius + dblX, m_lpdblSphereCoordinates[m_lpsSphereFacets[i * 3 + 1] * 3 + 1] * dblRadius + dblY, m_lpdblSphereCoordinates[m_lpsSphereFacets[i * 3 + 1] * 3 + 2] * dblRadius + dblZ);
		glVertex3d(m_lpdblSphereCoordinates[m_lpsSphereFacets[i * 3 + 2] * 3] * dblRadius + dblX, m_lpdblSphereCoordinates[m_lpsSphereFacets[i * 3 + 2] * 3 + 1] * dblRadius + dblY, m_lpdblSphereCoordinates[m_lpsSphereFacets[i * 3 + 2] * 3 + 2] * dblRadius + dblZ);
	}
	glColor3d(1.0, 1.0, 1.0);
	glEnd();
}


void CGLCanvas::m_fnDrawHairyBall(double dblX, double dblY, double dblZ, double dblRadius, double dblR, double dblG, double dblB)
{

	double dblM = dblRadius * 0.618;
	double lpdblCoordinates[36] =
	{
		dblX, dblY + dblM, dblZ + dblRadius,
		dblX, dblY - dblM, dblZ + dblRadius,
		dblX, dblY - dblM, dblZ - dblRadius,
		dblX, dblY + dblM, dblZ - dblRadius,
		dblX + dblRadius, dblY, dblZ + dblM,
		dblX + dblRadius, dblY, dblZ - dblM,
		dblX - dblRadius, dblY, dblZ - dblM,
		dblX - dblRadius, dblY, dblZ + dblM,
		dblX + dblM, dblY + dblRadius, dblZ,
		dblX - dblM, dblY + dblRadius, dblZ,
		dblX - dblM, dblY - dblRadius, dblZ,
		dblX + dblM, dblY - dblRadius, dblZ
	};
	int i;
	glBegin(GL_LINES);
	glColor3d(dblR, dblG, dblB);
	for (i = 0; i < 12; ++i)
	{
		glVertex3d(dblX, dblY, dblZ);
		glVertex3d(lpdblCoordinates[i * 3], lpdblCoordinates[i * 3 + 1], lpdblCoordinates[i * 3 + 2]);
	}
	glEnd();
}

void CGLCanvas::m_fnDrawFacets()
{
    CMeshBase::CPolyhedron::Facet_iterator fi;
    CMeshBase::CPolyhedron::Halfedge_handle hh0, hh;
    if(m_lpMeshCenter->m_mshSurface.size_of_facets() > 0)
    {
        for(fi = m_lpMeshCenter->m_mshSurface.facets_begin(); fi != m_lpMeshCenter->m_mshSurface.facets_end(); ++fi)
		{
			if (m_iFacetJetColor == 1)
			{
				glBegin(GL_POLYGON);
				hh0 = fi->halfedge();
				hh = hh0;
				do
				{
					glColor3d(hh->vertex()->lpdblColor[0], hh->vertex()->lpdblColor[1], hh->vertex()->lpdblColor[2]);
					glVertex3d(hh->vertex()->point().x(), hh->vertex()->point().y(), hh->vertex()->point().z());
					hh = hh->next();
				} while (hh != hh0);
				glEnd();
			}
			else if (m_iFacetJetColor == 2)
			{
				glBegin(GL_POLYGON);
				hh0 = fi->halfedge();
				hh = hh0;
				do
				{
					glColor3d(fi->lpdblColor[0], fi->lpdblColor[1], fi->lpdblColor[2]);
					glVertex3d(hh->vertex()->point().x(), hh->vertex()->point().y(), hh->vertex()->point().z());
					hh = hh->next();
				} while (hh != hh0);
				glEnd();
			}
			else
			{
				if (fi->bSelected)
				{
					glColor3d(0.875, 0.875, 0.875);
				}
				else
				{
					glColor3d(1.0, 1.0, 1.0);
				}
				glBegin(GL_POLYGON);
				hh0 = fi->halfedge();
				hh = hh0;
				do
				{
					glVertex3d(hh->vertex()->point().x(), hh->vertex()->point().y(), hh->vertex()->point().z());
					hh = hh->next();
				} while (hh != hh0);
				glEnd();
			}
        }
    }
}


void CGLCanvas::m_fnDrawJetFacets()
{
	CMeshBase::CPolyhedron::Facet_iterator fi;
	CMeshBase::CPolyhedron::Halfedge_handle hh0, hh;
	int iPara;
	switch (m_iViewableColor)
	{
	case 1:
		iPara = 0;
		break;
	case 2:
		iPara = 1;
		break;
	default:
		iPara = -1;
	}
	if (m_lpMeshCenter->m_mshSurface.size_of_facets() > 0 && iPara != -1)
	{
		glBegin(GL_TRIANGLES);
		for (fi = m_lpMeshCenter->m_mshSurface.facets_begin(); fi != m_lpMeshCenter->m_mshSurface.facets_end(); ++fi)
		{
			hh0 = fi->halfedge();
			hh = hh0;
			do
			{
				glColor3d(hh->lpdblColor[iPara * 3 + 3], hh->lpdblColor[iPara * 3 + 4], hh->lpdblColor[iPara * 3 + 5]);
				glVertex3d(hh->vertex()->point().x(), hh->vertex()->point().y(), hh->vertex()->point().z());
				hh = hh->next();
			} while (hh != hh0);
		}
		glEnd();
	}
}

void CGLCanvas::m_fnDrawEdges()
{
    CMeshBase::CPolyhedron::Edge_iterator ei;
    if(m_lpMeshCenter->m_mshSurface.size_of_halfedges() > 0)
    {
        glBegin(GL_LINES);
        if(m_iViewableColor == 1 || m_iViewableColor == 2)
        {
            for(ei = m_lpMeshCenter->m_mshSurface.edges_begin(); ei != m_lpMeshCenter->m_mshSurface.edges_end(); ++ei)
            {
                if(ei->bFreeEdge)
                {
                    if(ei->is_border())
                    {
						glColor3d(ei->opposite()->prev()->lpdblColor[m_iViewableColor * 3], ei->opposite()->prev()->lpdblColor[m_iViewableColor * 3 + 1], ei->opposite()->prev()->lpdblColor[m_iViewableColor * 3 + 2]);
                    }
                    else if(ei->opposite()->is_border())
                    {
						glColor3d(ei->lpdblColor[m_iViewableColor * 3], ei->lpdblColor[m_iViewableColor * 3 + 1], ei->lpdblColor[m_iViewableColor * 3 + 2]);
                    }
                    else
                    {
                        glColor3d(0, 0, 0);
                    }
                }
                else
                {
					glColor3d(ei->lpdblColor[m_iViewableColor * 3], ei->lpdblColor[m_iViewableColor * 3 + 1], ei->lpdblColor[m_iViewableColor * 3 + 2]);
                }
                glVertex3d(ei->vertex()->point().x(), ei->vertex()->point().y(), ei->vertex()->point().z());
                if(ei->opposite()->bFreeEdge)
                {
                    if(ei->is_border())
                    {
						glColor3d(ei->opposite()->lpdblColor[m_iViewableColor * 3], ei->opposite()->lpdblColor[m_iViewableColor * 3 + 1], ei->opposite()->lpdblColor[m_iViewableColor * 3 + 2]);
                    }
                    else if(ei->opposite()->is_border())
                    {
						glColor3d(ei->prev()->lpdblColor[m_iViewableColor * 3], ei->prev()->lpdblColor[m_iViewableColor * 3 + 1], ei->prev()->lpdblColor[m_iViewableColor * 3 + 2]);
                    }
                    else
                    {
                        glColor3d(0, 0, 0);
                    }
                }
                else
                {
					glColor3d(ei->opposite()->lpdblColor[m_iViewableColor * 3], ei->opposite()->lpdblColor[m_iViewableColor * 3 + 1], ei->opposite()->lpdblColor[m_iViewableColor * 3 + 2]);
                }
                glVertex3d(ei->opposite()->vertex()->point().x(), ei->opposite()->vertex()->point().y(), ei->opposite()->vertex()->point().z());
            }
        }
		else
        {
            for(ei = m_lpMeshCenter->m_mshSurface.edges_begin(); ei != m_lpMeshCenter->m_mshSurface.edges_end(); ++ei)
            {
				if (ei->bFreeEdge)
                {
					//glColor3d(0.0, 0.0, 0.0);
					//glVertex3d(ei->vertex()->point().x(), ei->vertex()->point().y(), ei->vertex()->point().z());
					//glVertex3d(ei->opposite()->vertex()->point().x(), ei->opposite()->vertex()->point().y(), ei->opposite()->vertex()->point().z());
					if (ei->bSelected)
					{
						switch (ei->iSeamType)
						{
						case -2:
							glColor3d(1.0, 0.0, 1.0);
							glVertex3d(ei->vertex()->point().x(), ei->vertex()->point().y(), ei->vertex()->point().z());
							glVertex3d(ei->opposite()->vertex()->point().x(), ei->opposite()->vertex()->point().y(), ei->opposite()->vertex()->point().z());
							break;
						case -1:
							glColor3d(0.0, 1.0, 0.0);
							glVertex3d(ei->vertex()->point().x(), ei->vertex()->point().y(), ei->vertex()->point().z());
							glVertex3d(ei->opposite()->vertex()->point().x(), ei->opposite()->vertex()->point().y(), ei->opposite()->vertex()->point().z());
							break;
						case 1:
							glColor3d(0.0, 1.0, 0.0);
							glVertex3d(ei->vertex()->point().x(), ei->vertex()->point().y(), ei->vertex()->point().z());
							glVertex3d(ei->opposite()->vertex()->point().x(), ei->opposite()->vertex()->point().y(), ei->opposite()->vertex()->point().z());
							break;
						default:
							glColor3d(1.0, 0.0, 1.0);
							glVertex3d(ei->vertex()->point().x(), ei->vertex()->point().y(), ei->vertex()->point().z());
							glVertex3d(ei->opposite()->vertex()->point().x(), ei->opposite()->vertex()->point().y(), ei->opposite()->vertex()->point().z());
						}
					}
					else
					{
						switch (ei->opposite()->iSeamType)
						{
						case -2:
							glColor3d(0.0, 1.0, 0.0);
							glVertex3d(ei->vertex()->point().x(), ei->vertex()->point().y(), ei->vertex()->point().z());
							glVertex3d(ei->opposite()->vertex()->point().x(), ei->opposite()->vertex()->point().y(), ei->opposite()->vertex()->point().z());
							break;
						case -1:
							glColor3d(1.0, 0.0, 1.0);
							glVertex3d(ei->vertex()->point().x(), ei->vertex()->point().y(), ei->vertex()->point().z());
							glVertex3d(ei->opposite()->vertex()->point().x(), ei->opposite()->vertex()->point().y(), ei->opposite()->vertex()->point().z());
							break;
						case 1:
							glColor3d(1.0, 0.0, 1.0);
							glVertex3d(ei->vertex()->point().x(), ei->vertex()->point().y(), ei->vertex()->point().z());
							glVertex3d(ei->opposite()->vertex()->point().x(), ei->opposite()->vertex()->point().y(), ei->opposite()->vertex()->point().z());
							break;
						default:
							glColor3d(0.0, 1.0, 0.0);
							glVertex3d(ei->vertex()->point().x(), ei->vertex()->point().y(), ei->vertex()->point().z());
							glVertex3d(ei->opposite()->vertex()->point().x(), ei->opposite()->vertex()->point().y(), ei->opposite()->vertex()->point().z());
						}
					}
                }
                else
                {
					if (ei->bSelected)
					{
						glColor3d(1, 0, 0);
					}
					else
					{
						glColor3d(0.5, 0.5, 0.5);
					}
					glVertex3d(ei->vertex()->point().x(), ei->vertex()->point().y(), ei->vertex()->point().z());
					glVertex3d(ei->opposite()->vertex()->point().x(), ei->opposite()->vertex()->point().y(), ei->opposite()->vertex()->point().z());
				}
            }
        }
    glEnd();
    }


}

void CGLCanvas::m_fnDrawBoundary()
{
	CMeshBase::CPolyhedron::Edge_iterator ei;
	if (m_lpMeshCenter->m_mshSurface.size_of_edges() > 0)
	{
		glBegin(GL_LINES);
		glColor3d(0.0, 0.0, 0.0);
		for (ei = m_lpMeshCenter->m_mshSurface.edges_begin(); ei != m_lpMeshCenter->m_mshSurface.edges_end(); ++ei)
		{
			if (ei->is_border_edge())
			{
				glVertex3d(ei->prev()->vertex()->point().x(), ei->prev()->vertex()->point().y(), ei->prev()->vertex()->point().z());
				glVertex3d(ei->vertex()->point().x(), ei->vertex()->point().y(), ei->vertex()->point().z());
			}
		}
		glEnd();
	}
}

void CGLCanvas::m_fnDrawChart()
{
    CMeshBase::CPolyhedron::Facet_iterator fi;
    int i;
    double lpColorCompoR[4] = {1.0, 0, 0, 0.0};
    double lpColorCompoG[4] = {0.0, 1.0, 0.0, 0.0};
    double lpColorCompoB[4] = {0, 0, 1.0, 0.0};
    double dblCosTheta, dblSinTheta;
	CMeshBase::Vector_3 lpvtChart[4];
    if(m_lpMeshCenter->m_mshSurface.size_of_facets() > 0)
    {
        for(fi = m_lpMeshCenter->m_mshSurface.facets_begin(); fi != m_lpMeshCenter->m_mshSurface.facets_end(); ++fi)
        {
	        dblCosTheta = cos(fi->ftChartDir);
	        dblSinTheta = sin(fi->ftChartDir);
			lpvtChart[0] = fi->vtPrincipalAxis * dblCosTheta + fi->vtAuxiliaryAxis * dblSinTheta;
			lpvtChart[2] = -lpvtChart[0];
			if (m_bFrameField)
			{
				dblCosTheta = cos(fi->ftFrameDir);
				dblSinTheta = sin(fi->ftFrameDir);
				lpvtChart[1] = fi->vtPrincipalAxis * dblCosTheta + fi->vtAuxiliaryAxis * dblSinTheta;
				lpvtChart[3] = -lpvtChart[1];
			}
			else
			{
				lpvtChart[1] = -fi->vtPrincipalAxis * dblSinTheta + fi->vtAuxiliaryAxis * dblCosTheta;
				lpvtChart[3] = -lpvtChart[1];
			}
			if (fi->bDirFixed)
			{
				glLineWidth(3);
				glBegin(GL_LINES);
	            for(i = 0; i < 4; ++i)
	            {
	                glColor3d(lpColorCompoR[i], lpColorCompoG[i], lpColorCompoB[i]);
	                glVertex3d(fi->ptIncent.x(), fi->ptIncent.y(), fi->ptIncent.z());
	                glVertex3d
	                (
	                    (fi->ptIncent.x() + lpvtChart[i].x()),
	                    (fi->ptIncent.y() + lpvtChart[i].y()),
	                    (fi->ptIncent.z() + lpvtChart[i].z())
	                );
	            }
				glEnd();
			}
			else
			{
				glLineWidth(1);
				glBegin(GL_LINES);
				for(i = 0; i < 4; ++i)
				{
					//dblTheta = fi->ftChartDir + i * ftHalfPi;
					//dblCosTheta = cos(dblTheta);
					//dblSinTheta = sin(dblTheta);
					glColor3d(lpColorCompoR[i], lpColorCompoG[i], lpColorCompoB[i]);
					glVertex3d(fi->ptIncent.x(), fi->ptIncent.y(), fi->ptIncent.z());
					glVertex3d
					(
						(fi->ptIncent.x() + lpvtChart[i].x()),
						(fi->ptIncent.y() + lpvtChart[i].y()),
						(fi->ptIncent.z() + lpvtChart[i].z())
					);
				}
				glEnd();
			}
        }
		glLineWidth(1);
        glColor3d(lpColorCompoR[3], lpColorCompoG[3], lpColorCompoB[3]);
    }
}
//CMeshBase::CPolyhedron::Vertex_iterator vi;
	//CMeshBase::CPolyhedron::Halfedge_handle hh;
	//glBegin(GL_LINES);
	//for (vi = m_lpMeshCenter->m_mshSurface.vertices_begin(); vi != m_lpMeshCenter->m_mshSurface.vertices_end(); ++vi)
	//{
	//	hh = vi->halfedge();
	//	if (hh->is_border())
	//	{
	//		hh = hh->opposite()->prev();
	//	}
	//	for (i = 0; i < 4; ++i)
	//	{
	//		glColor3d(lpColorCompoR[i], lpColorCompoG[i], lpColorCompoB[i]);
	//		glVertex3d(vi->point().x(), vi->point().y(), vi->point().z());
	//		glVertex3d(
	//			vi->point().x() + hh->lpvtViewField[i % 2].x() * (1 - (i / 2) * 2) * hh->ftLen / 2,
	//			vi->point().y() + hh->lpvtViewField[i % 2].y() * (1 - (i / 2) * 2) * hh->ftLen / 2,
	//			vi->point().z() + hh->lpvtViewField[i % 2].z() * (1 - (i / 2) * 2) * hh->ftLen / 2);
	//	}
	//}
	//glEnd();

void CGLCanvas::m_fnDrawSingularities()
{
    CMeshBase::CPolyhedron::Vertex_iterator vi;
    if(m_lpMeshCenter->m_mshSurface.size_of_vertices() > 0)
    {
        for(vi = m_lpMeshCenter->m_mshSurface.vertices_begin(); vi != m_lpMeshCenter->m_mshSurface.vertices_end(); ++vi)
        {
			if (vi->nNumPorts != 0)
			{
				if (vi->iBordId == -1)
				{
					if (vi->iSelected == 0)
					{
						if(vi->nNumPorts > 4)
						{
							m_fnDrawSphere(vi->point().x(), vi->point().y(), vi->point().z(), m_lpMeshCenter->m_ftMeanRadius, 0.5, 0.5, 1.0);
						}
						if(vi->nNumPorts < 4)
						{
							m_fnDrawSphere(vi->point().x(), vi->point().y(), vi->point().z(), m_lpMeshCenter->m_ftMeanRadius, 1.0, 0.5, 0.5);
						}
					}
					else
					{
						if(vi->nNumPorts > 4)
						{
							m_fnDrawSphere(vi->point().x(), vi->point().y(), vi->point().z(), m_lpMeshCenter->m_ftMeanRadius, 0, 0, 1.0);
						}
						if(vi->nNumPorts < 4)
						{
							m_fnDrawSphere(vi->point().x(), vi->point().y(), vi->point().z(), m_lpMeshCenter->m_ftMeanRadius, 1.0, 0, 0);
						}
					}
				}
				else
				{
					if (vi->iSelected == 0)
					{
						if (vi->iDegreeBias < 0)
						{
							m_fnDrawSphere(vi->point().x(), vi->point().y(), vi->point().z(), m_lpMeshCenter->m_ftMeanRadius, 0.5, 0.5, 1.0);
						}
						if (vi->iDegreeBias > 0)
						{
							m_fnDrawSphere(vi->point().x(), vi->point().y(), vi->point().z(), m_lpMeshCenter->m_ftMeanRadius, 1.0, 0.5, 0.5);
						}
					}
					else
					{
						if (vi->iDegreeBias < 0)
						{
							m_fnDrawSphere(vi->point().x(), vi->point().y(), vi->point().z(), m_lpMeshCenter->m_ftMeanRadius, 0, 0, 1.0);
						}
						if (vi->iDegreeBias > 0)
						{
							m_fnDrawSphere(vi->point().x(), vi->point().y(), vi->point().z(), m_lpMeshCenter->m_ftMeanRadius, 1.0, 0, 0);
						}
					}
				}
			}
			else if (vi->iSelected != 0)
			{
				m_fnDrawSphere(vi->point().x(), vi->point().y(), vi->point().z(), m_lpMeshCenter->m_ftMeanRadius, 0, 0.5, 0);
			}
        }
    }
}

void CGLCanvas::m_fnDrawTrack()
{
	std::list<CTracer::CInducedPath>::iterator iterIP;
	std::list<CTracer::CPathElement>::iterator iterPE1, iterPE2;
	CMeshBase::Point_3 ptStart, ptEnd;
	if (!m_lpMeshCenter->m_Tracer.m_lstIP.empty())
	{
		glBegin(GL_LINES);
		for (iterIP = m_lpMeshCenter->m_Tracer.m_lstIP.begin(); iterIP != m_lpMeshCenter->m_Tracer.m_lstIP.end(); ++iterIP)
		{
			if (!iterIP->m_lstData.empty())
			{
				iterPE2 = iterIP->m_lstData.begin();
				iterPE1 = iterPE2;
				++iterPE2;
				while (iterPE2 != iterIP->m_lstData.end())
				{
					if (iterIP->m_bSelected)
					{
						glColor3d(1, 0, 0);
					}
					else
					{
						glColor3d(0, 0, 1);
					}
					ptStart = CGAL::barycenter(iterPE1->hhEdge->prev()->vertex()->point(), 1 - iterPE1->ftProportion, iterPE1->hhEdge->vertex()->point(), iterPE1->ftProportion);
					ptEnd = CGAL::barycenter(iterPE2->hhEdge->prev()->vertex()->point(), 1 - iterPE2->ftProportion, iterPE2->hhEdge->vertex()->point(), iterPE2->ftProportion);
					glVertex3d(ptStart.x(), ptStart.y(), ptStart.z());
					glVertex3d(ptEnd.x(), ptEnd.y(), ptEnd.z());
					iterPE1 = iterPE2;
					++iterPE2;
				}
			}
		}
		glEnd();
	}
}

void CGLCanvas::m_fnDrawElectricField(int iType)
{
    CMeshBase::CPolyhedron::Facet_iterator fi;
    double lpColorCompoR[4] = {1.0, 0.5, 0.5, 1.0};
    double lpColorCompoG[4] = {0.5, 1.0, 0.5, 1.0};
    double lpColorCompoB[4] = {0.5, 0.5, 1.0, 1.0};
	CMeshBase::FT ftCoef;
	CMeshBase::Vector_3 vtArrow;
    if(m_lpMeshCenter->m_mshSurface.size_of_facets() > 0)
    {
        glBegin(GL_LINES);
		switch (iType)
		{
		case 0:
			for(fi = m_lpMeshCenter->m_mshSurface.facets_begin(); fi != m_lpMeshCenter->m_mshSurface.facets_end(); ++fi)
			{
				ftCoef = m_lpMeshCenter->m_ftMeanRadius / (fi->ftRadius * m_lpMeshCenter->m_ftMaxElectricField);
				vtArrow = fi->vtPrincipalAxis * fi->lpftElectricFieldSum[0] + fi->vtAuxiliaryAxis * fi->lpftElectricFieldSum[1];
				glColor3d(lpColorCompoR[1], lpColorCompoG[1], lpColorCompoB[1]);
				glVertex3d(fi->ptIncent.x(), fi->ptIncent.y(), fi->ptIncent.z());
				glColor3d(lpColorCompoR[0], lpColorCompoG[0], lpColorCompoB[0]);
				glVertex3d
				(
					(fi->ptIncent.x() + vtArrow.x() * ftCoef),
					(fi->ptIncent.y() + vtArrow.y() * ftCoef),
					(fi->ptIncent.z() + vtArrow.z() * ftCoef)
				);
				glColor3d(lpColorCompoR[1], lpColorCompoG[1], lpColorCompoB[1]);
				glVertex3d(fi->ptIncent.x(), fi->ptIncent.y(), fi->ptIncent.z());
				glColor3d(lpColorCompoR[2], lpColorCompoG[2], lpColorCompoB[2]);
				glVertex3d
				(
					(fi->ptIncent.x() - vtArrow.x() * ftCoef),
					(fi->ptIncent.y() - vtArrow.y() * ftCoef),
					(fi->ptIncent.z() - vtArrow.z() * ftCoef)
				);
			}
			break;
		case 1:
			for(fi = m_lpMeshCenter->m_mshSurface.facets_begin(); fi != m_lpMeshCenter->m_mshSurface.facets_end(); ++fi)
			{
				ftCoef = m_lpMeshCenter->m_ftMeanRadius / (fi->ftRadius * m_lpMeshCenter->m_ftCurlFieldDimension);
				vtArrow = fi->vtPrincipalAxis * fi->lpftElectricFieldDif[0] + fi->vtAuxiliaryAxis * fi->lpftElectricFieldDif[1];
				glColor3d(lpColorCompoR[1], lpColorCompoG[1], lpColorCompoB[1]);
				glVertex3d(fi->ptIncent.x(), fi->ptIncent.y(), fi->ptIncent.z());
				glColor3d(lpColorCompoR[0], lpColorCompoG[0], lpColorCompoB[0]);
				glVertex3d
				(
					(fi->ptIncent.x() + vtArrow.x() * ftCoef),
					(fi->ptIncent.y() + vtArrow.y() * ftCoef),
					(fi->ptIncent.z() + vtArrow.z() * ftCoef)
				);
				glColor3d(lpColorCompoR[1], lpColorCompoG[1], lpColorCompoB[1]);
				glVertex3d(fi->ptIncent.x(), fi->ptIncent.y(), fi->ptIncent.z());
				glColor3d(lpColorCompoR[2], lpColorCompoG[2], lpColorCompoB[2]);
				glVertex3d
				(
					(fi->ptIncent.x() - vtArrow.x() * ftCoef),
					(fi->ptIncent.y() - vtArrow.y() * ftCoef),
					(fi->ptIncent.z() - vtArrow.z() * ftCoef)
				);
			}
			break;
		}
        glColor3d(lpColorCompoR[3], lpColorCompoG[3], lpColorCompoB[3]);
        glEnd();
    }
}

void CGLCanvas::m_fnDrawLocalPara()
{
    CMeshBase::CPolyhedron::Facet_iterator fi;
    double lpColorCompoR[4] = {1.0, 0.0, 0.0, 1.0};
    double lpColorCompoG[4] = {0.0, 0.5, 0.0, 1.0};
    double lpColorCompoB[4] = {0.0, 0.0, 1.0, 1.0};
    double dblCosTheta, dblSinTheta;
    if(m_lpMeshCenter->m_mshSurface.size_of_facets() > 0)
    {
        glBegin(GL_LINES);
        for(fi = m_lpMeshCenter->m_mshSurface.facets_begin(); fi != m_lpMeshCenter->m_mshSurface.facets_end(); ++fi)
        {
            dblCosTheta = fi->lpftLocalParaScl[m_iViewableLocalPara] * cos(fi->lpftLocalParaDir[m_iViewableLocalPara]);
            dblSinTheta = fi->lpftLocalParaScl[m_iViewableLocalPara] * sin(fi->lpftLocalParaDir[m_iViewableLocalPara]);
            glColor3d(lpColorCompoR[1], lpColorCompoG[1], lpColorCompoB[1]);
            glVertex3d(fi->ptIncent.x(), fi->ptIncent.y(), fi->ptIncent.z());
            glColor3d(lpColorCompoR[0], lpColorCompoG[0], lpColorCompoB[0]);
            glVertex3d
            (
                (fi->ptIncent.x() + fi->vtPrincipalAxis.x() * dblCosTheta + fi->vtAuxiliaryAxis.x() * dblSinTheta),
                (fi->ptIncent.y() + fi->vtPrincipalAxis.y() * dblCosTheta + fi->vtAuxiliaryAxis.y() * dblSinTheta),
                (fi->ptIncent.z() + fi->vtPrincipalAxis.z() * dblCosTheta + fi->vtAuxiliaryAxis.z() * dblSinTheta)
            );
            glColor3d(lpColorCompoR[1], lpColorCompoG[1], lpColorCompoB[1]);
            glVertex3d(fi->ptIncent.x(), fi->ptIncent.y(), fi->ptIncent.z());
            glColor3d(lpColorCompoR[2], lpColorCompoG[2], lpColorCompoB[2]);
            glVertex3d
            (
                (fi->ptIncent.x() - fi->vtPrincipalAxis.x() * dblCosTheta - fi->vtAuxiliaryAxis.x() * dblSinTheta),
                (fi->ptIncent.y() - fi->vtPrincipalAxis.y() * dblCosTheta - fi->vtAuxiliaryAxis.y() * dblSinTheta),
                (fi->ptIncent.z() - fi->vtPrincipalAxis.z() * dblCosTheta - fi->vtAuxiliaryAxis.z() * dblSinTheta)
            );
        }
        glColor3d(lpColorCompoR[3], lpColorCompoG[3], lpColorCompoB[3]);
        glEnd();
    }
}

void CGLCanvas::m_fnDrawGrids()
{
	CMeshBase::CPolyhedron::Facet_iterator fi;
	CMeshBase::CPolyhedron::Halfedge_handle hh0, hh;
	CMeshBase::FT ftMax, ftMin, ftStart, ftEnd, ftCur, ftCoef0, ftCoef1;
	int i, nNum, iPara;
	CMeshBase::CPolyhedron::Point_3 lpptTerminals[3];
	double lpdblColor[9];
	if (m_lpMeshCenter->m_mshSurface.size_of_facets() > 0)
	{
		for (iPara = 0; iPara < 2; ++iPara)
		{
			if ((m_iViewableColor & (1 << iPara)) != 0)
			{
				glBegin(GL_LINES);
				for (fi = m_lpMeshCenter->m_mshSurface.facets_begin(); fi != m_lpMeshCenter->m_mshSurface.facets_end(); ++fi)
				{
					hh0 = fi->halfedge();
					ftMax = hh0->lpftGlobalPara[iPara];
					ftMin = hh0->lpftGlobalPara[iPara];
					hh = hh0->next();
					while (hh != hh0)
					{
						if (hh->lpftGlobalPara[iPara] < ftMin)
						{
							ftMin = hh->lpftGlobalPara[iPara];
						}
						if (hh->lpftGlobalPara[iPara] > ftMax)
						{
							ftMax = hh->lpftGlobalPara[iPara];
						}
						hh = hh->next();
					}
					ftStart = floor(ftMin);//floor(ftMin * 2.0) / 2.0;
					ftEnd = ceil(ftMax);//ceil(ftMax * 2.0) / 2.0;
					if (ftStart < ftEnd)
					{
						for (ftCur = ftStart; ftCur < ftEnd; ftCur += 1.0)// (ftCur = ftStart; ftCur < ftEnd; ftCur += 0.5)
						{
							nNum = 0;
							hh = hh0;
							do
							{
								ftCoef0 = ftCur - hh->prev()->lpftGlobalPara[iPara];
								ftCoef1 = hh->lpftGlobalPara[iPara] - ftCur;
								if (((ftCoef0 <= 0 && ftCoef1 <= 0) || (ftCoef0 >= 0 && ftCoef1 >= 0)) && !(ftCoef0 == 0 && ftCoef1 == 0))
								{
									lpptTerminals[nNum] = CGAL::barycenter(hh->prev()->vertex()->point(), ftCoef1, hh->vertex()->point(), ftCoef0);
									for (i = 0; i < 3; ++i)
									{
										lpdblColor[nNum * 3 + i] = (hh->prev()->lpdblColor[iPara * 3 + 3 + i] * ftCoef1 + hh->lpdblColor[iPara * 3 + 3 + i] * ftCoef0) / (ftCoef0 + ftCoef1);
									}
									++nNum;
								}
								hh = hh->next();
							} while (hh != hh0);
							if (nNum > 1)
							{
								for (i = 1; i < nNum; ++i)
								{
									glColor3d(lpdblColor[(i - 1) * 3], lpdblColor[(i - 1) * 3 + 1], lpdblColor[(i - 1) * 3 + 2]);
									glVertex3d(lpptTerminals[i - 1].x(), lpptTerminals[i - 1].y(), lpptTerminals[i - 1].z());
									glColor3d(lpdblColor[i * 3], lpdblColor[i * 3 + 1], lpdblColor[i * 3 + 2]);
									glVertex3d(lpptTerminals[i].x(), lpptTerminals[i].y(), lpptTerminals[i].z());
								}
							}
						}
					}
				}
				glEnd();
			}
		}
	}
}

void CGLCanvas::m_fnDrawSamplePoints()
{
    std::vector<CSamplePoint>::iterator iterSamplePoint;
    if(m_lpMeshCenter->m_vecSamplePoints.size() > 0)
    {
        for(iterSamplePoint = m_lpMeshCenter->m_vecSamplePoints.begin(); iterSamplePoint != m_lpMeshCenter->m_vecSamplePoints.end(); ++iterSamplePoint)
        {
			m_fnDrawSphere(iterSamplePoint->ptCoordinate.x(), iterSamplePoint->ptCoordinate.y(), iterSamplePoint->ptCoordinate.z(), m_lpMeshCenter->m_ftMeanRadius, 0.25, 0.5, 0.25);
        }
    }
}

void CGLCanvas::m_fnDrawLinkages()
{
    std::vector<CSamplePoint>::iterator iterSamplePoint;
    CSamplePoint *lpAdjacentPoint;
    std::vector<int>::iterator iterAdjacent;
    if(m_lpMeshCenter->m_vecSamplePoints.size() > 0)
    {
        glBegin(GL_LINES);
        glColor3d(0.0, 0.0, 0.0);
        for(iterSamplePoint = m_lpMeshCenter->m_vecSamplePoints.begin(); iterSamplePoint != m_lpMeshCenter->m_vecSamplePoints.end(); ++iterSamplePoint)
        {
            if(iterSamplePoint->m_vecAdjacent.size() > 0)
            {
                for(iterAdjacent = iterSamplePoint->m_vecAdjacent.begin(); iterAdjacent != iterSamplePoint->m_vecAdjacent.end(); ++iterAdjacent)
                {
                    lpAdjacentPoint = &(m_lpMeshCenter->m_vecSamplePoints[*iterAdjacent]);
                    glVertex3d(iterSamplePoint->ptCoordinate.x(), iterSamplePoint->ptCoordinate.y(), iterSamplePoint->ptCoordinate.z());
                    glVertex3d(lpAdjacentPoint->ptCoordinate.x(), lpAdjacentPoint->ptCoordinate.y(), lpAdjacentPoint->ptCoordinate.z());
                }
            }
        }
        glEnd();
    }
}

void CGLCanvas::m_fnDrawFacetDirections()
{
    CSamplePoint* lplpSamplePoint[4];
    CMeshBase::Vector_3 vtNormal;
    CMeshBase::FT ftSqrtNormLen;
    CMeshBase::Point_3 ptBaryCenter, ptNormalTerminal;
    int i, j;
    std::vector<CFacetLinkage>::iterator iterFacet;
    if(m_lpMeshCenter->m_vecFacetLinkages.size() > 0)
    {
        for(iterFacet = m_lpMeshCenter->m_vecFacetLinkages.begin(); iterFacet != m_lpMeshCenter->m_vecFacetLinkages.end(); ++iterFacet)
        {
            for(i = 0; i < 4; ++i)
            {
                lplpSamplePoint[i] = &(m_lpMeshCenter->m_vecSamplePoints[iterFacet->lpiVertices[i]]);
            }
            ptBaryCenter = CGAL::barycenter<CMeshBase::FT>(lplpSamplePoint[0]->ptCoordinate, 1, lplpSamplePoint[1]->ptCoordinate, 1, lplpSamplePoint[2]->ptCoordinate, 1, lplpSamplePoint[3]->ptCoordinate, 1);
            vtNormal = CGAL::cross_product((lplpSamplePoint[1]->ptCoordinate - lplpSamplePoint[0]->ptCoordinate), (lplpSamplePoint[3]->ptCoordinate - lplpSamplePoint[0]->ptCoordinate)) + CGAL::cross_product((lplpSamplePoint[3]->ptCoordinate - lplpSamplePoint[2]->ptCoordinate), (lplpSamplePoint[1]->ptCoordinate - lplpSamplePoint[2]->ptCoordinate));
            ftSqrtNormLen = sqrt(sqrt(vtNormal * vtNormal));
            ptNormalTerminal = ptBaryCenter + vtNormal * 3.0 / ftSqrtNormLen;
            glBegin(GL_LINES);
            glColor3d(0, 0, 0);
            for(i = 0; i < 4; ++i)
            {
                j = (i + 1) % 4;
                glVertex3d(lplpSamplePoint[i]->ptCoordinate.x(), lplpSamplePoint[i]->ptCoordinate.y(), lplpSamplePoint[i]->ptCoordinate.z());
                glVertex3d(lplpSamplePoint[j]->ptCoordinate.x(), lplpSamplePoint[j]->ptCoordinate.y(), lplpSamplePoint[j]->ptCoordinate.z());
            }
            glColor3d(0.75, 0.875, 0.75);
            glVertex3d(ptBaryCenter.x(), ptBaryCenter.y(), ptBaryCenter.z());
            glVertex3d(ptNormalTerminal.x(), ptNormalTerminal.y(), ptNormalTerminal.z());
            glEnd();
        }
    }
}

void CGLCanvas::m_fnDrawNewMeshEdges()
{
    CMeshBase::CSamplerPolyhedron::Edge_iterator ei;
    if(m_lpMeshCenter->m_mshRemeshedSurface.size_of_halfedges() > 0)
    {
        glBegin(GL_LINES);
        for(ei = m_lpMeshCenter->m_mshRemeshedSurface.edges_begin(); ei != m_lpMeshCenter->m_mshRemeshedSurface.edges_end(); ++ei)
        {
			if (ei->bCrevasse && m_bTrackViewable)
			{
				glColor3d(1, 0, 0);
			}
			else
			{
				glColor3d(0, 0, 0);
			}
            glVertex3d(ei->vertex()->point().x(), ei->vertex()->point().y(), ei->vertex()->point().z());
            glVertex3d(ei->opposite()->vertex()->point().x(), ei->opposite()->vertex()->point().y(), ei->opposite()->vertex()->point().z());
        }
        glEnd();
    }
}

void CGLCanvas::m_fnDrawNewMeshFacets()
{
    CMeshBase::CSamplerPolyhedron::Facet_iterator fi;
    CMeshBase::CSamplerPolyhedron::Halfedge_handle hh;
    int i;
    if(m_lpMeshCenter->m_mshRemeshedSurface.size_of_facets() > 0)
    {
        glBegin(GL_QUADS);
        glColor3d(1, 1, 1);
        for(fi = m_lpMeshCenter->m_mshRemeshedSurface.facets_begin(); fi != m_lpMeshCenter->m_mshRemeshedSurface.facets_end(); ++fi)
        {
            hh = fi->halfedge();
            for(i = 0; i < 4; ++i)
            {
                glVertex3d(hh->vertex()->point().x(), hh->vertex()->point().y(), hh->vertex()->point().z());
                hh = hh->next();
            }
        }
        glEnd();
    }
}

void CGLCanvas::m_fnDrawNewMeshSingularities()
{
    CMeshBase::CSamplerPolyhedron::Vertex_iterator vi;
    if(m_lpMeshCenter->m_mshRemeshedSurface.size_of_vertices() > 0)
    {
        for(vi = m_lpMeshCenter->m_mshRemeshedSurface.vertices_begin(); vi != m_lpMeshCenter->m_mshRemeshedSurface.vertices_end(); ++vi)
        {
			if (vi->iBorderId == -1)
			{
				if(vi->degree() < 4)
				{
					m_fnDrawSphere(vi->point().x(), vi->point().y(), vi->point().z(), m_lpMeshCenter->m_ftNewMeanRadius * 2, 1.0, 0.5, 0.5);
				}
				else if(vi->degree() > 4)
				{
					m_fnDrawSphere(vi->point().x(), vi->point().y(), vi->point().z(), m_lpMeshCenter->m_ftNewMeanRadius * 2, 0.5, 0.5, 1.0);
				}
			}
			if (vi->iMark == 1)
			{
				m_fnDrawSphere(vi->point().x(), vi->point().y(), vi->point().z(), m_lpMeshCenter->m_ftNewMeanRadius * 2, 0.5, 1.0, 0.5);
			}
        }
    }
}

void CGLCanvas::OnPaint(wxPaintEvent&event)
{
	if(m_bFirstTime)
	{
		m_bFirstTime = false;
		m_fnRotate();
		m_fnObserve();
		m_fnCast();
		glEnable(GL_DEPTH_TEST);
		glHint(GL_LINE_SMOOTH_HINT, GL_DONT_CARE);
		glEnable(GL_POLYGON_OFFSET_FILL);
		glEnable(GL_POLYGON_OFFSET_POINT);
		glPolygonOffset(1.0, 1.0);
	}
	if (m_bShowLIC)
	{
		wxPaintDC dc(this);
		dc.DrawBitmap(m_bmpLIC, m_iBmpX, m_iBmpY);
	}
	else
	{
		glClearColor(1.0, 1.0, 1.0, 0.0);
		glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
		glLoadIdentity();
		glPushMatrix();
		glMultMatrixd(m_lpdblProjectMat);
		glMultMatrixd(m_lpdblViewMat);
		glMultMatrixd(m_lpdblWorldMat);
		if(m_iViewableMeshInfo == -1)
		{
			m_fnDrawFacets();
		}
		else if(m_iViewableMeshInfo == 3)
		{
			m_fnDrawNewMeshFacets();
		}
		if(m_iViewableMeshInfo != -1)
		{
			switch(m_iViewableMeshInfo)
			{
			case 0:
				m_fnDrawSamplePoints();
				break;
			case 1:
				m_fnDrawLinkages();
				break;
			case 2:
				m_fnDrawFacetDirections();
				break;
			case 3:
				m_fnDrawNewMeshEdges();
				if(m_bSingularityViewable)
				{
					m_fnDrawNewMeshSingularities();
				}
				break;
			}
		}
		else
		{
			if (m_bShowEdge)
			{
				m_fnDrawEdges();
			}
			else
			{
				m_fnDrawBoundary();
			}
			if ((m_iViewableColor == 1 || m_iViewableColor == 2 || m_iViewableColor == 3) && m_lpMeshCenter->m_lpftGlobalPara != NULL)
			{
				m_fnDrawGrids();
			}
			if(m_bChartViewable)
			{
				m_fnDrawChart();
			}
			if(m_bSingularityViewable)
			{
				m_fnDrawSingularities();
			}
			if (m_iElectricFieldViewable == 1 || m_iElectricFieldViewable == 3)
			{
				m_fnDrawElectricField(0);
			}
			if (m_iElectricFieldViewable == 2 || m_iElectricFieldViewable == 3)
			{
				m_fnDrawElectricField(1);
			}
			if (m_bTrackViewable)
			{
				m_fnDrawTrack();
			}
			if(m_iViewableLocalPara > -1 && m_iViewableLocalPara < 2)
			{
				m_fnDrawLocalPara();
			}
			if (m_iViewableColor == 4 || m_iViewableColor == 5)
			{
				//m_fnDrawWeight();
			}
		}
		if (m_bTrackViewable)
		{
			m_fnDrawTrack();
		}
		glPopMatrix();
		glFlush();
		SwapBuffers();
	}
}

void CGLCanvas::m_fnViewCoordinate(double dblX, double dblY, double dblZ, double &dblViewX, double &dblViewY, double &dblViewZ)
{
    int i, j;
    double lpdblCoordinateBuffer0[4], lpdblCoordinateBuffer1[4];
    lpdblCoordinateBuffer0[0] = dblX;
    lpdblCoordinateBuffer0[1] = dblY;
    lpdblCoordinateBuffer0[2] = dblZ;
    lpdblCoordinateBuffer0[3] = 1;
    for(i = 0; i < 4; ++i)
    {
        lpdblCoordinateBuffer1[i] = 0;
        for(j = 0; j < 4; ++j)
        {
            lpdblCoordinateBuffer1[i] += lpdblCoordinateBuffer0[j] * m_lpdblWorldMat[j * 4 + i];
        }
    }
    for(i = 0; i < 4; ++i)
    {
        lpdblCoordinateBuffer0[i] = 0;
        for(j = 0; j < 4; ++j)
        {
            lpdblCoordinateBuffer0[i] += lpdblCoordinateBuffer1[j] * m_lpdblViewMat[j * 4 + i];
        }
    }
    for(i = 0; i < 4; ++i)
    {
        lpdblCoordinateBuffer1[i] = 0;
        for(j = 0; j < 4; ++j)
        {
            lpdblCoordinateBuffer1[i] += lpdblCoordinateBuffer0[j] * m_lpdblProjectMat[j * 4 + i];
        }
    }
    dblViewX = lpdblCoordinateBuffer1[0] / lpdblCoordinateBuffer1[3] * double(m_iHalfWidth) + double(m_iHalfWidth);
    dblViewY = double(m_iHalfHeight) - lpdblCoordinateBuffer1[1] / lpdblCoordinateBuffer1[3] * double(m_iHalfHeight);
    dblViewZ = lpdblCoordinateBuffer1[2] / lpdblCoordinateBuffer1[3];
}

void CGLCanvas::m_fnPositionOnScreen(double dblX, double dblY, double dblZ, int &nX, int &nY, double &dblDepth)
{
	double dblViewX, dblViewY, dblViewZ;
	m_fnViewCoordinate(dblX, dblY, dblZ, dblViewX, dblViewY, dblViewZ);
	nX = int(round(dblViewX));
	nY = int(round(dblViewY));
	dblDepth = dblViewZ;
}

double CGLCanvas::m_fnDirectionOnScreen(double dblX, double dblY, double dblZ)
{
	int i, j;
	double lpdblCoordinateBuffer0[4], lpdblCoordinateBuffer1[4], dblTheta;
	lpdblCoordinateBuffer0[0] = dblX;
	lpdblCoordinateBuffer0[1] = dblY;
	lpdblCoordinateBuffer0[2] = dblZ;
	lpdblCoordinateBuffer0[3] = 0;
	for (i = 0; i < 4; ++i)
	{
		lpdblCoordinateBuffer1[i] = 0;
		for (j = 0; j < 4; ++j)
		{
			lpdblCoordinateBuffer1[i] += lpdblCoordinateBuffer0[j] * m_lpdblWorldMat[j * 4 + i];
		}
	}
	for (i = 0; i < 4; ++i)
	{
		lpdblCoordinateBuffer0[i] = 0;
		for (j = 0; j < 4; ++j)
		{
			lpdblCoordinateBuffer0[i] += lpdblCoordinateBuffer1[j] * m_lpdblViewMat[j * 4 + i];
		}
	}
	for (i = 0; i < 4; ++i)
	{
		lpdblCoordinateBuffer1[i] = 0;
		for (j = 0; j < 4; ++j)
		{
			lpdblCoordinateBuffer1[i] += lpdblCoordinateBuffer0[j] * m_lpdblProjectMat[j * 4 + i];
		}
	}
	if (lpdblCoordinateBuffer1[0] != 0 || lpdblCoordinateBuffer1[1] != 0)
	{
		dblTheta = atan2(-lpdblCoordinateBuffer1[1] * double(m_iHalfHeight), lpdblCoordinateBuffer1[0] * double(m_iHalfWidth));
	}
	else
	{
		dblTheta = 4.0;
	}
	return dblTheta;
}

void CGLCanvas::m_fnNormalMap(double *lpdblInput, double *lpdblOutput)
{
	int i, j;
	for (i = 0; i < 3; ++i)
	{
		lpdblOutput[i] = 0;
		for (j = 0; j < 3; ++j)
		{
			lpdblOutput[i] += lpdblInput[j] * m_lpdblWorldMat[j * 4 + i];
		}
	}
}

void CGLCanvas::m_fnDetectSelection(int xMouse, int yMouse, CMeshBase::CPolyhedron::Vertex_handle &vhSelecting)
{
    CMeshBase::CPolyhedron::Vertex_iterator vi;
    int xView, yView, xDist, yDist;
    double dblDepth, dblMinDepth;
    dblMinDepth = 1.0;
    vhSelecting = NULL;
    if(m_lpMeshCenter->m_mshSurface.size_of_vertices() > 0)
    {
        for(vi = m_lpMeshCenter->m_mshSurface.vertices_begin(); vi != m_lpMeshCenter->m_mshSurface.vertices_end(); ++vi)
        {
            m_fnPositionOnScreen(vi->point().x(), vi->point().y(), vi->point().z(), xView, yView, dblDepth);
            xDist = xView - xMouse;
            yDist = yView - yMouse;
            if(xDist * xDist + yDist * yDist < 25)
            {
                if(dblDepth < dblMinDepth)
                {
                    vhSelecting = &(*vi);
                    dblMinDepth = dblDepth;
                }
            }
        }
    }
}



BEGIN_EVENT_TABLE(CGLCanvas, wxGLCanvas)
    EVT_PAINT(CGLCanvas::OnPaint)
    EVT_MOUSE_EVENTS(CGLCanvas::OnMouseEvent)
    EVT_KEY_DOWN(CGLCanvas::OnKeyEvent)
    EVT_SIZE(CGLCanvas::OnSize)
    EVT_TIMER(wxID_ANY, CGLCanvas::OnTimer)
END_EVENT_TABLE()

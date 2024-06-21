#include <wx/dcclient.h>
#include <wx/timer.h>
#include <wx/wx.h>
#include <wx/scrolwin.h>
#include <wx/dcclient.h>
#include <wx/dcbuffer.h>
#include <wx/rawbmp.h>
#include "CGLContext.h"
#include "CGeoPack.hpp"
#include "CDialogs.h"
//included by "CLICPainter.h"
#ifndef CGLCANVAS_H
#define CGLCANVAS_H
union PixelDataConvertor
{
	unsigned char lpPixels[5808];
	struct
	{
		double lpCoordinate[486];
		short lpsIndex[960];
	};
};
class CGLCanvas : public wxGLCanvas
{
public:
    CGLCanvas(wxWindow *lpParent, int *lpAttribList);
	~CGLCanvas();
    CGLContext *m_lpGLContext;
    double m_lpdblWorldMat[16];
    double m_lpdblViewMat[16];
    double m_lpdblProjectMat[16];
    CMeshInfoReserver *m_lpMeshCenter;
    int m_xMouse, m_yMouse;
    bool m_bLButtonDown;
    void m_fnRotate();
    void m_fnObserve();
    void m_fnCast();
    int m_iHalfWidth;
    int m_iHalfHeight;
    double m_dblNear;
    double m_dblDistance;
    double m_dblFar;
    double m_dblRotateAngle;
    double m_lpdblRotateAxis[3];
    double m_dblTanHalfFov;
    bool m_bFirstTime;

	bool m_bOrthodox;
	bool m_bShowEdge;
    bool m_bChartViewable;
	bool m_bFrameField;
    bool m_bSingularityViewable;
    //bool m_bElectricFieldViewable;
	int m_iElectricFieldViewable;
	bool m_bTrackViewable;
	int m_iFacetJetColor;

	bool m_bSelectPoint;
	bool m_bSelectTrack;
	bool m_bSelectCrevasse;
	bool m_bSelectFacet;
	bool m_bSelectEdge;

	int m_iViewableLocalPara;
    int m_iViewableColor;
    int m_iViewableMeshInfo;
    int m_iViewableRay;
    int m_iViewableTriangle;


	bool m_bMovieMode, m_bRecordMode, m_bDiscrete;
	wxString m_wxstrMovieDir;
	wxString m_wxstrRecordFile;
	wxTimer m_Timer;
	bool m_bProcessing;
	void m_fnViewCoordinate(double dblX, double dblY, double dblZ, double &dblViewX, double &dblViewY, double &dblViewZ);
    void m_fnPositionOnScreen(double dblX, double dblY, double dblZ, int &nX, int &nY, double &dblDepth);
	double m_fnDirectionOnScreen(double dblX, double dblY, double dblZ);
	void m_fnNormalMap(double *lpdblInput, double *lpdblOutput);
	bool m_bShowLIC;
	wxBitmap m_bmpOriginal, m_bmpLIC;
	int m_iBmpX, m_iBmpY;
	double m_lpdblSphereCoordinates[486];
	short m_lpsSphereFacets[960];
protected:
private:
    void OnPaint(wxPaintEvent&event);
    void OnMouseEvent(wxMouseEvent &event);
    void OnKeyEvent(wxKeyEvent &event);
    void OnSize(wxSizeEvent &event);
    void OnTimer(wxTimerEvent &event);
	void m_fnDrawSphere(double dblX, double dblY, double dblZ, double dblRadius, double dblR, double dblG, double dblB);
	void m_fnDrawHairyBall(double dblX, double dblY, double dblZ, double dblRadius, double dblR, double dblG, double dblB);
    void m_fnDrawFacets();
	void m_fnDrawJetFacets();
    void m_fnDrawEdges();
	void m_fnDrawBoundary();
    void m_fnDrawChart();
	void m_fnDrawElectricField(int iType);
	void m_fnDrawTrack();
    void m_fnDrawLocalPara();
	void m_fnDrawGrids();
    void m_fnDrawSingularities();
    void m_fnDrawSamplePoints();
    void m_fnDrawLinkages();
    void m_fnDrawFacetDirections();
    void m_fnDrawNewMeshEdges();
    void m_fnDrawNewMeshFacets();
    void m_fnDrawNewMeshSingularities();
    void m_fnDetectSelection(int xMouse, int yMouse, CMeshBase::CPolyhedron::Vertex_handle &vhSelecting);

    DECLARE_EVENT_TABLE()
};

#endif // CGLCANVAS_H

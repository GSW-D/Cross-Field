#include <float.h>
#include "CGLCanvas.h"
//included by "CAlgorithmExecutor.h"
#ifndef CLICPAINTER_H
#define	CLICPAINTER_H
struct DirPixel
{
	float fDepth, lpfDirection[2];
	float lpfColor[6], fStrength;
	CMeshBase::CPolyhedron::Halfedge_handle hhLocation, hhCrevasse;
	short int lpiSeam[4];
};
class CLICPainter
{
public:
	CLICPainter(CGLCanvas *lpGLCanvas);
	~CLICPainter();
	DirPixel *m_lpdpBuffer;
	int m_iXOffset, m_iYOffset, m_iWidth, m_iHeight;
	void m_fnGenerateViewCoordinates();
	void m_fnTestSize(double dblRadiusScale);
	void m_fnVideoTestSize(double dblRadiusScale, int iWidth, int iHeight);
	void m_fnInitBuffer();
	void m_fnSetLocations(bool bCrossField);
	void m_fnSetCrevasse(double dblWidth);
	void m_fnDetectSeam();
	void m_fnAssignOriginalColor(wxBitmap &bmpOriginal);
	void m_fnCrossTrace(float fRadius);
	void m_fnLinearTrace(float fRadius);
	void m_fnEnhanceBitmap(bool bShowStrength, float fMaxStrength = 1.0);
	void m_fnWeakFieldRender(float fMaxStrength);
	void m_fnNoFieldRender(float fMaxStrength);
	void m_fnDrawDiscrete();
	void m_fnAppendCrevasse(double dblWidth);
	void m_fnAppendSingularities(double dblRadiusScale);
	void m_fnMarkForceToDraw();
	void m_fnDrawForce();
	void m_fnWriteBitmap(wxBitmap &bmpLIC, int &iBmpX, int &iBmpY);
protected:
	CGLCanvas *m_lpGLCanvas;
	double m_fnPlanarArea(double *lpdblPoint1, double *lpdblPoint2, double *lpdblPoint3)
	{
		return (
			lpdblPoint1[0] * lpdblPoint2[1] - lpdblPoint1[1] * lpdblPoint2[0] + \
			lpdblPoint2[0] * lpdblPoint3[1] - lpdblPoint2[1] * lpdblPoint3[0] + \
			lpdblPoint3[0] * lpdblPoint1[1] - lpdblPoint3[1] * lpdblPoint1[0]
			) / 2;
	}
	int m_fnAdjacent(DirPixel *lpdp0, DirPixel *lpdp1);
};
#endif

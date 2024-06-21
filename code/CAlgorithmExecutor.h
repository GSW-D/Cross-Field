#include <wx/thread.h>
#include <wx/filename.h>
#include <complex>
#include "CLICPainter.h"
//included by "CDialogs.h"
#ifndef CALGORITHMEXECUTOR_H
#define CALGORITHMEXECUTOR_H
enum PROCESS_TYPE
{
	PT_FLATTEN,
	PT_COMPLEX,
	PT_SINGULARITY,
	PT_GLOBAL,
	PT_COUNT_SINGULARITIES,
	PT_COARSE_DECURL,
	PT_FINE_DECURL,
	PT_CONSTRAINED_DECURL,
    PT_OPTIMIZE
};
class CAlgorithmExecutor: public wxThread
{
public:
    PROCESS_TYPE m_ptProcessType;
	int m_nIterateTimes, m_nTestTimes;
	//CMeshBase::FT m_ftChargeRadius;
	double m_dblStartLevel, m_dblDecrement;
	wxString m_wxstrInputField, m_wxstrOutputData;
	int m_nNumFacets, m_iAdjustType;

    CAlgorithmExecutor(CGLCanvas *lpGLCanvas, PROCESS_TYPE ptProcessType);
protected:
    CGLCanvas *m_lpGLCanvas;
    virtual ExitCode Entry();
	void m_fnFlattenCones();
	void m_fnComplexAdjust();
    void m_fnAdjustSingularities();
	void m_fnGlobalAdjust();
	void m_fnCountSingularities();
	void m_fnCoarseDecurl();
	void m_fnFineDecurl();
	void m_fnConstrainedDecurl();
    void m_fnOptimizeMesh();
};
#endif //CALGORITHMEXECUTOR_H

#include <wx/msgdlg.h>
#include <wx/filedlg.h>
#include <wx/textdlg.h>
#include "CAlgorithmExecutor.h"
//included by "CApp.h"
#ifndef CMAINFRAME_H
#define CMAINFRAME_H
class CMainFrame : public wxFrame
{
public:
    CMainFrame(wxFrame *lpParentFrame, const wxString &lpstrTitle);
    ~CMainFrame();
protected:
private:
    wxMenuBar *m_lpMenuBar;
	wxToolBar *m_lpToolBar;
    wxMenu *m_lpFileMenu;
    wxMenu *m_lpViewMenu;
    wxMenu *m_lpCoarseDirectionMenu;
    wxMenu *m_lpAccurateDirectionMenu;
	wxMenu *m_lpDirectionConstraintMenu;
	wxMenu *m_lpCurlEditMenu;
    wxMenu *m_lpParameterMenu;
    wxMenu *m_lpSampleMenu;
    CGLCanvas *m_lpGLCanvas;
	wxString m_wxstrMeshName;
	wxString m_wxstrMeshDir;

    void OnOpen(wxCommandEvent &event);
    void OnSave(wxCommandEvent &event);
    void OnSaveProcessStatus(wxCommandEvent &event);
    void OnOpenProcessStatus(wxCommandEvent &event);
	void OnExportMountain(wxCommandEvent &event);
	void OnExportDirectionField(wxCommandEvent &event);
	void OnExportParaDif(wxCommandEvent& event);
	void OnImportCones(wxCommandEvent& event);
	void OnExportSingularities(wxCommandEvent &event);
	void OnExportFixedDirection(wxCommandEvent &event);
	void OnExportHardFeature(wxCommandEvent &event);
	void OnExportFixedMark(wxCommandEvent &event);
	void OnImportFixedDirection(wxCommandEvent &event);
	void OnImportDirectionField(wxCommandEvent &event);
	void OnImportBitmap(wxCommandEvent &event);
	void OnExportLIC(wxCommandEvent &event);
	void OnExportTexturePara(wxCommandEvent &event);
	void OnExportAngleError(wxCommandEvent &event);
	void OnExportCrevasse(wxCommandEvent &event);
	void OnExportConfCoefColor(wxCommandEvent &event);
	void OnOpenMovieMode(wxCommandEvent &event);
	void OnOpenRecordMode(wxCommandEvent &event);
    void OnExportNewMesh(wxCommandEvent &event);
    void OnImportQuadMesh(wxCommandEvent &event);
	void OnExportQuadSingularities(wxCommandEvent &event);
	void OnExportQuadCylinders(wxCommandEvent &event);
	void OnExit(wxCommandEvent &event);

    void OnChangeViewPort(wxCommandEvent &event);
	void OnAppendRotation(wxCommandEvent &event);
	void OnOrthodox(wxCommandEvent &event);
	void OnShowEdge(wxCommandEvent &event);
	void OnShowLIC(wxCommandEvent &event);
    void OnChartViewable(wxCommandEvent &event);
	void OnFrameField(wxCommandEvent &event);
    void OnSingularityViewable(wxCommandEvent &event);
	void OnTrackViewable(wxCommandEvent &event);
    void OnShowElectricFieldSum(wxCommandEvent &event);
	void OnShowElectricFieldDif(wxCommandEvent &event);
	void OnShowGaussCurvature(wxCommandEvent &event);
	void OnShowConformalCoefficient(wxCommandEvent &event);
	void OnShowAngleError(wxCommandEvent &event);
	void OnShowParaCurl(wxCommandEvent& event);
    void OnDensityU(wxCommandEvent &event);
	void OnDensityV(wxCommandEvent &event);
	void OnLocalU(wxCommandEvent &event);
	void OnLocalV(wxCommandEvent &event);
	void OnGlobalU(wxCommandEvent &event);
	void OnGlobalV(wxCommandEvent &event);
    void OnShowSamplePoints(wxCommandEvent &event);
    void OnShowLinkages(wxCommandEvent &event);
    void OnShowFacetDirections(wxCommandEvent &event);
    void OnShowReconstructedMesh(wxCommandEvent &event);

	void OnNormalizeModal(wxCommandEvent &event);
    void OnInitilizeChart(wxCommandEvent &event);
	void OnEigenOptimal(wxCommandEvent &event);
	void OnFlatten(wxCommandEvent &event);
	void OnFixCones(wxCommandEvent &event);
    void OnDeletePath(wxCommandEvent &event);
	void OnComplexOptimize(wxCommandEvent &event);
	void OnStopIteration(wxCommandEvent &event);

    void OnDirectionAdj(wxCommandEvent &event);
    void OnAdjustSingularities(wxCommandEvent &event);
	void OnSetRadius(wxCommandEvent &event);
	void OnReportEnergy(wxCommandEvent &event);
	void OnReportSingularities(wxCommandEvent &event);
	void OnLinkPoles(wxCommandEvent &event);
    void OnMoveSingularity(wxCommandEvent &event);
    void OnClearBarrier(wxCommandEvent &event);
	void OnGlobalAdjust(wxCommandEvent &event);
	void OnCountSingularities(wxCommandEvent &event);


	void OnMarkSharpFeature(wxCommandEvent &event);
	void OnClearFix(wxCommandEvent &event);

	void OnDetectLoop(wxCommandEvent &event);
	void OnReleaseCurl(wxCommandEvent &event);
	void OnCoarseDecurl(wxCommandEvent &event);
	void OnFineDecurl(wxCommandEvent &event);
	void OnConstrainedDecurl(wxCommandEvent &event);
	void OnAdjustCurl(wxCommandEvent &event);

	void OnTraceDirection(wxCommandEvent &event);
	void OnBiasSearch(wxCommandEvent &event);
	void OnAdjustPath(wxCommandEvent &event);

	void OnSelectPoint(wxCommandEvent &event);
	void OnSelectTrack(wxCommandEvent &event);
	void OnSelectCrevasse(wxCommandEvent &event);
	void OnSelectFacet(wxCommandEvent &event);
	void OnSelectEdge(wxCommandEvent &event);

	void OnDecomposeMesh(wxCommandEvent &event);
    void OnUniformizeChart(wxCommandEvent &event);
	void OnSolveGlobalParameter(wxCommandEvent &event);
	void OnRescaleGlobalParameter(wxCommandEvent &event);
    void OnSeamIntegerize(wxCommandEvent &event);
	void OnEvaluatePara(wxCommandEvent &event);

    void OnExtractSamplePoints(wxCommandEvent &event);
    void OnGenerateLinkages(wxCommandEvent &event);
    void OnConstructFacets(wxCommandEvent &event);
    void OnUniformizeFacets(wxCommandEvent &event);
    void OnConstructNewMesh(wxCommandEvent &event);
    void OnOptimizeNewMesh(wxCommandEvent &event);
	void OnSimplifyNewMesh(wxCommandEvent &event);
	void OnEvaluateNewMesh(wxCommandEvent &event);

	enum
	{
		idMenuSaveProcessStatus = 1200,
		idMenuOpenProcessStatus,
		idMenuExportMountain,
		idMenuExportDirectionField,
		idMenuExportParaDif,
		idMenuImportCones,
		idMenuExportSingularities,
		idMenuExportFixedDirection,
		idMenuExportHardFeature,
		idMenuExportFixedMark,
		idMenuImportFixedDirection,
		idMenuImportDirectionField,
		idMenuImportBitmap,
		idMenuExportLIC,
		idMenuExportCrevasse,
		idMenuExportTexturePara,
		idMenuExportAngleError,
		idMenuExportConfCoefColor,
		idMenuOpenMovieMode,
		idMenuOpenRecordMode,
		idMenuExportNewMesh,
		idMenuImportQuadMesh,
		idMenuExportQuadSingularities,
		idMenuExportQuadCylinders,
		idMenuExit,

		idMenuChangeViewPort,
		idMenuAppendRotation,
		idOrthodox,
		idShowEdge,
		idShowLIC,
		idShowChart,
		idShowFrameField,
		idShowSingularity,
		idShowElectricFieldSum,
		idShowElectricFieldDif,
		idShowTrack,
		idSelectPoint,
		idSelectTrack,
		idSelectCrevasse,
		idSelectFacet,
		idSelectEdge,
		idShowGaussCurvature,
		idShowConformalCoefficient,
		idShowAngleError,
		idShowParaCurl,
		idShowDensityU,
		idShowDensityV,
		idShowLocalU,
		idShowLocalV,
		idShowGlobalU,
		idShowGlobalV,
		idMenuMarkedColorViewable,
		idShowSamplePoints,
		idShowLinkages,
		idShowFacetNormal,
		idShowNewMesh,

		idMenuNormalizeModal,
		idMenuInitilizeChart,
		idMenuEigenOptimal,
		idMenuFlatten,
		idMenuFixCones,
		idMenuDeletePath,
		idMenuOpenEnergyCalculator,
		idMenuComplexOptimize,
		idMenuStopIteration,

		idMenuDirectionAdj,
		idMenuAdjustSingularities,
		idMenuSetRadius,
		idMenuReportEnergy,
		idMenuReportSingularities,
		idMenuLinkPoles,
		idMenuMoveSingularity,
		idMenuClearBarrier,
		idMenuGlobalAdjust,
		idMenuCountSingularities,

		idMenuMarkSharpFeature,
		idMenuClearFix,

		idMenuDetectLoop,
		idMenuReleaseCurl,
		idMenuCoarseDecurl,
		idMenuFineDecurl,
		idMenuConstrainedDecurl,
		idMenuConstrainedOptimize,
		idMenuAdjustCurl,

		idMenuTraceDirection,
		idMenuBiasSearch,
		idMenuAdjustPath,
		idMenuDecomposeMesh,
		idMenuUniformizeChart,
		idMenuSolveGlobalParameter,
		idMenuRescaleGlobalParameter,
		idMenuSeamIntegerize,
		idMenuEvaluatePara,

		idMenuExtractSamplePoints,
		idMenuGenerateLinkages,
		idMenuConstructFacets,
		idMenuUniformizeFacets,
		idMenuConstructNewMesh,
		idMenuOptimizeNewMesh,
		idMenuSimplifyNewMesh,
		idMenuEvaluateNewMesh
	};
DECLARE_EVENT_TABLE()
};


#endif // CMAINFRAME_H

#include "CMainFrame.h"
CMainFrame::CMainFrame(wxFrame *lpParentFrame, const wxString &lpstrTitle):wxFrame(lpParentFrame, wxID_ANY, lpstrTitle, wxPoint(0, 0), wxSize(800, 650))
{
    //ctor
    m_lpMenuBar = new wxMenuBar();
	m_lpToolBar = new wxToolBar(this, wxID_ANY);
    m_lpFileMenu = new wxMenu(_T(""));
    m_lpFileMenu->Append(wxID_OPEN, _("&Open"), _(""));
    m_lpFileMenu->Append(wxID_SAVE, _("&Save"), _(""));
    m_lpFileMenu->Append(idMenuSaveProcessStatus, _("&SaveProcessStatus"), _(""));
    m_lpFileMenu->Append(idMenuOpenProcessStatus, _("&OpenProcessStatus"), _(""));
	m_lpFileMenu->Append(idMenuExportMountain, _("ExportMountain"), _(""));
    m_lpFileMenu->Append(idMenuExportDirectionField, _("&ExportDirectionField"), _(""));
	m_lpFileMenu->Append(idMenuExportParaDif, _("&ExportParaDif"), _(""));
	m_lpFileMenu->Append(idMenuImportCones, _("&ImportCones"), (""));
	m_lpFileMenu->Append(idMenuExportSingularities, _("Export Singularities"), _(""));
	m_lpFileMenu->Append(idMenuExportFixedDirection, _("ExportFixedDirection"), _(""));
	m_lpFileMenu->Append(idMenuExportHardFeature, _("ExportHardFeature"), _(""));
	m_lpFileMenu->Append(idMenuExportFixedMark, _("ExportFixedMark"), _(""));
	m_lpFileMenu->Append(idMenuImportFixedDirection, _("ImportFixedDirection"), _(""));
	m_lpFileMenu->Append(idMenuImportDirectionField, _("ImportDirectionField"), _(""));
	m_lpFileMenu->Append(idMenuImportBitmap, _("ImportBitmap"), _(""));
	m_lpFileMenu->Append(idMenuExportLIC, _("ExportLIC"), _(""));
	m_lpFileMenu->Append(idMenuExportCrevasse, _("Export Crevasse"), _(""));
	m_lpFileMenu->Append(idMenuExportTexturePara, _("&ExportTexturePara"), _(""));
	m_lpFileMenu->Append(idMenuExportAngleError, _("&ExportAngleError"), _(""));
	m_lpFileMenu->Append(idMenuExportConfCoefColor, _("Export Confermal Coefficient Color"), _(""));
	m_lpFileMenu->Append(idMenuOpenMovieMode, _("&OpenMovieMode"), _(""));
	m_lpFileMenu->Append(idMenuOpenRecordMode, _("&OpenRecordMode"), _(""));
    m_lpFileMenu->Append(idMenuExportNewMesh, _("&ExportNewMesh"), _(""));
    m_lpFileMenu->Append(idMenuImportQuadMesh, _("&ImportQuadMesh"), _(""));
	m_lpFileMenu->Append(idMenuExportQuadSingularities, _("ExportQuadSingularities"), _(""));
	m_lpFileMenu->Append(idMenuExportQuadCylinders, _("ExportQuadCylinders"), _(""));
    m_lpFileMenu->Append(idMenuExit, _("E&xit"), _(""));

    m_lpViewMenu = new wxMenu(_T(""));
    m_lpViewMenu->Append(idMenuChangeViewPort, _("ChangeViewPort"), _(""));
	m_lpViewMenu->Append(idMenuAppendRotation, _("AppendRotation"), _(""));

    m_lpCoarseDirectionMenu = new wxMenu(_T(""));
	m_lpCoarseDirectionMenu->Append(idMenuNormalizeModal, _("Normalize Modal"), _(""));
    m_lpCoarseDirectionMenu->Append(idMenuInitilizeChart, _("InitilizeChart"), _(""));
	m_lpCoarseDirectionMenu->Append(idMenuEigenOptimal, _("Eigen Optimal"), _(""));
	m_lpCoarseDirectionMenu->Append(idMenuFlatten, _("Flatten"), (""));
	m_lpCoarseDirectionMenu->Append(idMenuFixCones, _("FixCones"), (""));
    m_lpCoarseDirectionMenu->Append(idMenuDeletePath, _("&DeletePath"), _(""));
	m_lpCoarseDirectionMenu->Append(idMenuOpenEnergyCalculator, _("&OpenEnergyCalculator"), _(""));
	m_lpCoarseDirectionMenu->Append(idMenuComplexOptimize, _("&ComplexOptimize"), _(""));
	m_lpCoarseDirectionMenu->Append(idMenuStopIteration, _("Stop Iteration"), _(""));

    m_lpAccurateDirectionMenu = new wxMenu(_T(""));
    m_lpAccurateDirectionMenu->Append(idMenuDirectionAdj, _("Adjust&DirectionalDerivative"), _(""));
    m_lpAccurateDirectionMenu->Append(idMenuAdjustSingularities, _("Adjust Singularities"), _(""));
	m_lpAccurateDirectionMenu->Append(idMenuSetRadius, _("SetRadius"), _(""));
	m_lpAccurateDirectionMenu->Append(idMenuReportEnergy, _("Report Energy"),_(""));
	m_lpAccurateDirectionMenu->Append(idMenuReportSingularities, _("Report Singularities"), _(""));
	m_lpAccurateDirectionMenu->Append(idMenuLinkPoles, _("Link Poles"), _(""));
    m_lpAccurateDirectionMenu->Append(idMenuMoveSingularity, _("Move Singularity"), _(""));
    m_lpAccurateDirectionMenu->Append(idMenuClearBarrier, _("Clear Barriers"), _(""));
	m_lpAccurateDirectionMenu->Append(idMenuGlobalAdjust, _("Global Adjust"), _(""));
	m_lpAccurateDirectionMenu->Append(idMenuCountSingularities, _("CountSingularities"), _(""));

	m_lpDirectionConstraintMenu = new wxMenu(_T(""));
	m_lpDirectionConstraintMenu->Append(idMenuMarkSharpFeature, _("&Mark Sharp Feature"), _(""));
	m_lpDirectionConstraintMenu->Append(idMenuClearFix, _("&Clear Fix"), _(""));

	m_lpCurlEditMenu = new wxMenu(_T(""));
	m_lpCurlEditMenu->Append(idMenuDetectLoop, _("Detect Loop"), _(""));
	m_lpCurlEditMenu->Append(idMenuReleaseCurl, _("Release Curl"), _(""));
	m_lpCurlEditMenu->Append(idMenuCoarseDecurl, _("Coarse Decurl"), _(""));
	m_lpCurlEditMenu->Append(idMenuFineDecurl, _("Fine Decurl"), _(""));
	m_lpCurlEditMenu->Append(idMenuConstrainedDecurl, _("ConstrainedDecurl"), _(""));
	m_lpCurlEditMenu->Append(idMenuAdjustCurl, _("Adjust Curl"), _(""));

    m_lpParameterMenu = new wxMenu(_T(""));
	m_lpParameterMenu->Append(idMenuTraceDirection, _("&Trace"), _(""));
	m_lpParameterMenu->Append(idMenuBiasSearch, _("&BiasSearch"), _(""));
	m_lpParameterMenu->Append(idMenuAdjustPath, _("&AdjustPath"), _(""));
    m_lpParameterMenu->Append(idMenuDecomposeMesh, _("&Decompose"), _(""));
    m_lpParameterMenu->Append(idMenuUniformizeChart, _("&UniformizeChart"), _(""));
	m_lpParameterMenu->Append(idMenuSolveGlobalParameter, _("SolveGlobalParameter"), _(""));
	m_lpParameterMenu->Append(idMenuRescaleGlobalParameter, _("RescaleGlobalParameter"), _(""));
    m_lpParameterMenu->Append(idMenuSeamIntegerize, _("&SeamIntegerize"), _(""));
	m_lpParameterMenu->Append(idMenuEvaluatePara, _("&EvaluatePara"), _(""));

    m_lpSampleMenu = new wxMenu(_T(""));
    m_lpSampleMenu->Append(idMenuExtractSamplePoints, _("&ExtractSamplePoints"), _(""));
    m_lpSampleMenu->Append(idMenuGenerateLinkages, _("&GenerateLinkages"), _(""));
    m_lpSampleMenu->Append(idMenuConstructFacets, _("ConstructFacets"), _(""));
    m_lpSampleMenu->Append(idMenuUniformizeFacets, _("UniformizeFacets"), _(""));
    m_lpSampleMenu->Append(idMenuConstructNewMesh, _("ConstructNewMesh"), (""));
    m_lpSampleMenu->Append(idMenuOptimizeNewMesh, _("OptimizeNewMesh"), (""));
	m_lpSampleMenu->Append(idMenuSimplifyNewMesh, _("SimplifyNewMesh"), (""));
	m_lpSampleMenu->Append(idMenuEvaluateNewMesh, _("EvaluateNewMesh"), (""));

    m_lpMenuBar->Append(m_lpFileMenu, _("&File"));
    m_lpMenuBar->Append(m_lpViewMenu, _("&View"));
    m_lpMenuBar->Append(m_lpCoarseDirectionMenu, _("&CoarseDirection"));
    m_lpMenuBar->Append(m_lpAccurateDirectionMenu, _("&AccurateDirection"));
	m_lpMenuBar->Append(m_lpDirectionConstraintMenu, _("DirectionConstraint"));
	m_lpMenuBar->Append(m_lpCurlEditMenu, _("CurlEdit"));
    m_lpMenuBar->Append(m_lpParameterMenu, _("&Parameter"));
    m_lpMenuBar->Append(m_lpSampleMenu, _("Sample"));

	m_lpToolBar->AddTool(wxID_OPEN, wxEmptyString, wxBITMAP(OPEN_BMP), _("Open File"), wxITEM_NORMAL);
	m_lpToolBar->AddTool(wxID_SAVE, wxEmptyString, wxBITMAP(SAVE_BMP), _("Save File"), wxITEM_NORMAL);
	m_lpToolBar->AddTool(idOrthodox, wxEmptyString, wxBITMAP(ORTHODOX_BMP), _("Orthodox"), wxITEM_CHECK);
	m_lpToolBar->AddTool(idShowEdge, wxEmptyString, wxBITMAP(EDGE_BMP), _("Show Edge"), wxITEM_CHECK);
	m_lpToolBar->AddTool(idShowLIC, wxEmptyString, wxBITMAP(LIC_BMP), _("Show LIC"), wxITEM_CHECK);
	m_lpToolBar->AddTool(idShowChart, _("Show Chart"), wxBITMAP(CHART_BMP), _("Show Chart"), wxITEM_CHECK);
	m_lpToolBar->AddTool(idShowFrameField, _("Frame Field"), wxBITMAP(FRAME_BMP), _(""), wxITEM_CHECK);
	m_lpToolBar->AddTool(idShowSingularity, wxEmptyString, wxBITMAP(SINGULARITY_BMP), _("Show Singularity"), wxITEM_CHECK);
	m_lpToolBar->AddTool(idShowElectricFieldSum, wxEmptyString, wxBITMAP(ELEC_SUM_BMP), _("Show Electric Field Sum"), wxITEM_CHECK);
	m_lpToolBar->AddTool(idShowElectricFieldDif, wxEmptyString, wxBITMAP(ELEC_DIF_BMP), _("Show Electric Field Dif"), wxITEM_CHECK);
	m_lpToolBar->AddTool(idShowTrack, wxEmptyString, wxBITMAP(TRACK_BMP), _("Show Track"), wxITEM_CHECK);


	m_lpToolBar->AddTool(idSelectPoint, wxEmptyString, wxBITMAP(SELECT_PNT_BMP), _("Select Point"), wxITEM_CHECK);
	m_lpToolBar->AddTool(idSelectTrack, wxEmptyString, wxBITMAP(SELECT_TRACK_BMP), _("Select Track"), wxITEM_CHECK);
	m_lpToolBar->AddTool(idSelectCrevasse, wxEmptyString, wxBITMAP(SELECT_CREVASSE_BMP), _("Select Crevasse"), wxITEM_CHECK);
	m_lpToolBar->AddTool(idSelectFacet, wxEmptyString, wxBITMAP(SELECT_FACET_BMP), _("Select Facet"), wxITEM_CHECK);
	m_lpToolBar->AddTool(idSelectEdge, wxEmptyString, wxBITMAP(SELECT_EDGE_BMP), _("Select Edge"), wxITEM_CHECK);


	m_lpToolBar->AddTool(idShowGaussCurvature, wxEmptyString, wxBITMAP(GAUSS_CURVATURE_BMP), _("Show Gauss Curvature"), wxITEM_CHECK);
	m_lpToolBar->AddTool(idShowConformalCoefficient, wxEmptyString, wxBITMAP(CONF_COEF_BMP), _("Show Conformal Coefficient"), wxITEM_CHECK);
	m_lpToolBar->AddTool(idShowAngleError, wxEmptyString, wxBITMAP(ANGLE_ERR_BMP), _("Show Angle Error"), wxITEM_CHECK);
	m_lpToolBar->AddTool(idShowParaCurl, wxEmptyString, wxBITMAP(PARA_CURL_BMP), _("Show Para Curl"), wxITEM_CHECK);
	m_lpToolBar->AddTool(idShowDensityU, wxEmptyString, wxBITMAP(DENSITY_U_BMP), _("Show DensityU"), wxITEM_CHECK);
	m_lpToolBar->AddTool(idShowDensityV, wxEmptyString, wxBITMAP(DENSITY_V_BMP), _("Show DensityV"), wxITEM_CHECK);
	m_lpToolBar->AddTool(idShowLocalU, wxEmptyString, wxBITMAP(LOCAL_U_BMP), _("Show Local Parameter U"), wxITEM_CHECK);
	m_lpToolBar->AddTool(idShowLocalV, wxEmptyString, wxBITMAP(LOCAL_V_BMP), _("Show Local Parameter V"), wxITEM_CHECK);
	m_lpToolBar->AddTool(idShowGlobalU, wxEmptyString, wxBITMAP(GLOBAL_U_BMP), _("Show Gobal Parameter U"), wxITEM_CHECK);
	m_lpToolBar->AddTool(idShowGlobalV, wxEmptyString, wxBITMAP(GLOBAL_V_BMP), _("Show Gobal Parameter V"), wxITEM_CHECK);


	m_lpToolBar->AddTool(idShowSamplePoints, wxEmptyString, wxBITMAP(NEW_POINTS_BMP), _("Show Sample Points"), wxITEM_CHECK);
	m_lpToolBar->AddTool(idShowLinkages, wxEmptyString, wxBITMAP(NEW_LINKAGES_BMP), _("Show Linkages"), wxITEM_CHECK);
	m_lpToolBar->AddTool(idShowFacetNormal, wxEmptyString, wxBITMAP(NEW_NORMAL_BMP), _("Show Facet Normal"), wxITEM_CHECK);
	m_lpToolBar->AddTool(idShowNewMesh, wxEmptyString, wxBITMAP(NEW_MESH_BMP), _("Show New Mesh"), wxITEM_CHECK);

	m_lpToolBar->Realize();
    SetMenuBar(m_lpMenuBar);
	SetToolBar(m_lpToolBar);
	m_lpToolBar->ToggleTool(idShowEdge, true);
    int lpAttribList[] = { WX_GL_RGBA, WX_GL_DOUBLEBUFFER, WX_GL_STEREO, 0 };// {WX_GL_RGBA, WX_GL_DEPTH_SIZE, 1, 0};
    m_lpGLCanvas = new CGLCanvas(this, lpAttribList);
}

CMainFrame::~CMainFrame()
{
    //dtor
}

void CMainFrame::OnOpen(wxCommandEvent &event)
{
	std::ifstream ifs;
    wxFileDialog OpenFileDialog(this, _("OpenMeshFile"), "", "", "Mesh Files(*.off)|*.off|ObjFiles(*.obj)|*.obj|All Files(*.*)|*.*", wxFD_OPEN|wxFD_FILE_MUST_EXIST);
	wxMessageDialog MD(this,  _("Is the parameter from a frame field?"), _("Parametery Type"), wxYES_NO);
	char lpstrFileName[256], *lpstrFileType;
    if(OpenFileDialog.ShowModal() == wxID_OK)
    {
		strcpy(lpstrFileName, OpenFileDialog.GetPath().ToAscii());
		lpstrFileType = strrchr(lpstrFileName, '.');
		m_wxstrMeshName = OpenFileDialog.GetPath();
		m_wxstrMeshDir = OpenFileDialog.GetDirectory();
        m_lpGLCanvas->m_lpMeshCenter->m_fnStatusRetreat(PS_EMPTY);
		ifs.open(lpstrFileName);//(OpenFileDialog.GetPath().ToAscii());
		if (strcmp(lpstrFileType, ".off") == 0 || strcmp(lpstrFileType, ".OFF") == 0)
		{
			m_lpGLCanvas->m_lpMeshCenter->m_fnReadMesh(ifs);
			m_lpToolBar->EnableTool(wxID_SAVE, true);
			++(m_lpGLCanvas->m_lpMeshCenter->m_psProcessStatus);
			m_lpGLCanvas->m_lpMeshCenter->m_bProtectGlobalPara = false;
			m_lpGLCanvas->m_lpMeshCenter->m_fnFillGeometricInfo();
			m_lpGLCanvas->m_lpMeshCenter->m_fnDelaunayFlip();
			m_lpGLCanvas->m_lpMeshCenter->m_fnFillGeometricInfo();
			++(m_lpGLCanvas->m_lpMeshCenter->m_psProcessStatus);
		}
		if (strcmp(lpstrFileType, ".obj") == 0 || strcmp(lpstrFileType, ".OBJ") == 0)
		{
			if (m_lpGLCanvas->m_lpMeshCenter->m_fnImportParaObj(ifs) == 1)
			{
				m_lpToolBar->EnableTool(wxID_SAVE, true);
				++(m_lpGLCanvas->m_lpMeshCenter->m_psProcessStatus);
				m_lpGLCanvas->m_lpMeshCenter->m_bProtectGlobalPara = true;
				m_lpGLCanvas->m_lpMeshCenter->m_fnFillGeometricInfo();
				m_lpGLCanvas->m_lpMeshCenter->m_fnDelaunayFlip();
				m_lpGLCanvas->m_lpMeshCenter->m_fnFillGeometricInfo();
				m_lpGLCanvas->m_lpMeshCenter->m_bProtectGlobalPara = false;
				++(m_lpGLCanvas->m_lpMeshCenter->m_psProcessStatus);
				CDirectionGenerator DG(m_lpGLCanvas->m_lpMeshCenter->m_mshSurface);
				if (MD.ShowModal() == wxID_YES)
				{
					DG.m_fnFillDirFromPara(true);
				}
				else
				{
					DG.m_fnFillDirFromPara(false);
				}
				m_lpGLCanvas->m_lpMeshCenter->m_fnClearFreeEdge();
				m_lpGLCanvas->m_lpMeshCenter->m_fnDetectSeamType();
				m_lpGLCanvas->m_lpMeshCenter->m_fnGenerateMoment();
				m_lpGLCanvas->m_lpMeshCenter->m_fnMarkSingularity();
				++(m_lpGLCanvas->m_lpMeshCenter->m_psProcessStatus);
				++(m_lpGLCanvas->m_lpMeshCenter->m_psProcessStatus);
				++(m_lpGLCanvas->m_lpMeshCenter->m_psProcessStatus);
				CMeshSplitter MS(m_lpGLCanvas->m_lpMeshCenter->m_mshSurface);
				MS.m_fnSplitAlongPara();
				++(m_lpGLCanvas->m_lpMeshCenter->m_psProcessStatus);
				CParameterGenerator PG(m_lpGLCanvas->m_lpMeshCenter->m_mshSurface, m_lpGLCanvas->m_lpMeshCenter->m_vecCrevasses);
				PG.m_fnMarkCrevasse();
				PG.m_fnExtractCrevasses();
				++(m_lpGLCanvas->m_lpMeshCenter->m_psProcessStatus);
				m_lpGLCanvas->m_lpMeshCenter->m_fnRepresentGlobalPara();
			}
			else
			{
				++(m_lpGLCanvas->m_lpMeshCenter->m_psProcessStatus);
				m_lpGLCanvas->m_lpMeshCenter->m_fnFillGeometricInfo();
				m_lpGLCanvas->m_lpMeshCenter->m_fnDelaunayFlip();
				m_lpGLCanvas->m_lpMeshCenter->m_fnFillGeometricInfo();
				++(m_lpGLCanvas->m_lpMeshCenter->m_psProcessStatus);
			}
			m_wxstrMeshName.Replace(_(".obj"), _(".off"));
		}
		ifs.close();
        Refresh();
    }
}

void CMainFrame::OnSave(wxCommandEvent &event)
{
	std::ofstream ofs;
	char lpstrFileName[256], *lpstrFileType;
	wxFileDialog SaveFileDialog(this, _("SaveMeshFile"), "", "", "MeshFiles(*.off)|*.off|ObjFiles(*.obj)|*.obj", wxFD_SAVE | wxFD_OVERWRITE_PROMPT);
    if(m_lpGLCanvas->m_lpMeshCenter->m_psProcessStatus > PS_EMPTY)
    {
		if (SaveFileDialog.ShowModal() == wxID_OK)
		{
			strcpy(lpstrFileName, SaveFileDialog.GetPath().ToAscii());
			lpstrFileType = strrchr(lpstrFileName, '.');
			if (strcmp(lpstrFileType, ".off") == 0 || strcmp(lpstrFileType, ".OFF") == 0)
			{
				ofs.open(SaveFileDialog.GetPath().ToAscii());
				m_lpGLCanvas->m_lpMeshCenter->m_fnWriteMesh(ofs);
				ofs.close();
			}
			else if (strcmp(lpstrFileType, ".obj") == 0 || strcmp(lpstrFileType, ".OBJ") == 0)
			{
				ofs.open(SaveFileDialog.GetPath().ToAscii());
				m_lpGLCanvas->m_lpMeshCenter->m_fnExportTriangularPara(ofs);
				ofs.close();
			}
		}
    }
}

void CMainFrame::OnSaveProcessStatus(wxCommandEvent &event)
{
	std::ofstream ofs;
	wxFileDialog SaveFileDialog(this, _("SaveProcessStatus"), "", "", "ProcessStatusFiles(*.cfps)|*.cfps|All Files(*.*)|*.*", wxFD_SAVE | wxFD_OVERWRITE_PROMPT);
    if(SaveFileDialog.ShowModal() == wxID_OK)
    {
		ofs.open(SaveFileDialog.GetPath().ToAscii());
        m_lpGLCanvas->m_lpMeshCenter->m_fnWriteProcessStatus(ofs);
		ofs.close();
    }
}

void CMainFrame::OnOpenProcessStatus(wxCommandEvent &event)
{
	std::ifstream ifs;
    wxFileDialog OpenFileDialog(this, _("OpenProcessStatus"), "", "", "ProcessStatusFiles(*.cfps)|*.cfps|All Files(*.*)|*.*", wxFD_OPEN | wxFD_FILE_MUST_EXIST);
    if(OpenFileDialog.ShowModal() == wxID_OK)
    {
        m_lpGLCanvas->m_lpMeshCenter->m_fnStatusRetreat(PS_EMPTY);
		m_lpGLCanvas->m_lpMeshCenter->m_Tracer.m_lstIP.clear();
		ifs.open(OpenFileDialog.GetPath().ToAscii());
        m_lpGLCanvas->m_lpMeshCenter->m_fnReadProcessStatus(ifs);
		ifs.close();
		CDirectionGenerator DG(m_lpGLCanvas->m_lpMeshCenter->m_mshSurface);
		DG.m_fnInitFixedDirection();
		Refresh();
    }
}

void CMainFrame::OnExportMountain(wxCommandEvent &event)
{
	std::ofstream ofs;
	double dblBottom, dblTop, dblScale;
	CDirectionAdjustor DA(m_lpGLCanvas->m_lpMeshCenter->m_mshSurface);
	DA.m_fnModifyDensity(m_lpGLCanvas->m_lpMeshCenter->m_ftChargeRadius, dblBottom, dblTop);
	wxTextEntryDialog tedScale(this, "Scale:", _("Value Input"), _("1.0"), wxOK | wxCANCEL),
		tedMin(this, "MinColorValue:", _("Value Input"), wxString::FromDouble(dblBottom), wxOK | wxCANCEL),
		tedMax(this, "MaxColorValue:", _("Value Input"), wxString::FromDouble(dblTop), wxOK | wxCANCEL);
	wxFileDialog SaveFileDialog(this, _("Export Mountain"), "", "", "Wave Front Files(*.obj)|*.obj", wxFD_SAVE | wxFD_OVERWRITE_PROMPT);
	if (SaveFileDialog.ShowModal() == wxID_OK)
	{
		if (tedScale.ShowModal() == wxID_OK)
		{
			tedScale.GetValue().ToDouble(&dblScale);
			if (tedMin.ShowModal() == wxID_OK)
			{
				tedMin.GetValue().ToDouble(&dblBottom);
				if (tedMax.ShowModal() == wxID_OK)
				{
					tedMax.GetValue().ToDouble(&dblTop);
					ofs.open(SaveFileDialog.GetPath().ToAscii());
					m_lpGLCanvas->m_lpMeshCenter->m_fnExportMountain(ofs, dblBottom, dblTop, dblScale);
					ofs.close();
				}
			}
		}
	}
}


void CMainFrame::OnExportDirectionField(wxCommandEvent &event)
{
	std::ofstream ofs;
	wxString wxstrDefaultFile = m_wxstrMeshName;
	wxstrDefaultFile.Replace(".off", "_field.txt");

	wxFileDialog ExportFileDialog(this, _("ExportCrossField"), m_wxstrMeshDir, wxstrDefaultFile, "txt(*.txt)|*.txt", wxFD_SAVE | wxFD_OVERWRITE_PROMPT);
	if (ExportFileDialog.ShowModal() == wxID_OK)
	{
		ofs.open(ExportFileDialog.GetPath().ToAscii());
		m_lpGLCanvas->m_lpMeshCenter->m_fnExportCrossField(ofs);
		ofs.close();
	}
}
//{
//	std::ofstream ofs;
//	wxFileDialog ExportFileDialog(this, _("ExportCrossField"), "", "", "off(*.off)|*.off", wxFD_SAVE | wxFD_OVERWRITE_PROMPT);
//	if (ExportFileDialog.ShowModal() == wxID_OK)
//	{
//		ofs.open(ExportFileDialog.GetPath().ToAscii());
//		m_lpGLCanvas->m_lpMeshCenter->m_fnWriteMesh(ofs);
//		ofs.close();
//		wxString wxstrDirectionFieldFileName = ExportFileDialog.GetPath();
//		wxString wxstrSingularityFileName = ExportFileDialog.GetPath();
//		wxString wxstrElectricFieldFileName = ExportFileDialog.GetPath();
//		wxString wxstrCurlFieldFileName = ExportFileDialog.GetPath();
//		wxstrDirectionFieldFileName.Replace(_(".off"), _("_cross_field.txt"));
//		ofs.open(wxstrDirectionFieldFileName.ToAscii());
//		m_lpGLCanvas->m_lpMeshCenter->m_fnExportCrossField(ofs);
//		ofs.close();
//		wxstrSingularityFileName.Replace(_(".off"), _("_singularities.txt"));
//		ofs.open(wxstrSingularityFileName.ToAscii());
//		m_lpGLCanvas->m_lpMeshCenter->m_fnExportSingularities(ofs);
//		ofs.close();
//		wxstrElectricFieldFileName.Replace(_(".off"), _("_electric_field.txt"));
//		ofs.open(wxstrElectricFieldFileName.ToAscii());
//		m_lpGLCanvas->m_lpMeshCenter->m_fnExportElectricityField(ofs);
//		ofs.close();
//		wxstrCurlFieldFileName.Replace(_(".off"), _("_curl_field.txt"));
//		ofs.open(wxstrCurlFieldFileName.ToAscii());
//		m_lpGLCanvas->m_lpMeshCenter->m_fnExportCurlField(ofs);
//		ofs.close();
//	}
//}


void CMainFrame::OnExportParaDif(wxCommandEvent& event)
{
	std::ofstream ofs;
	wxFileDialog ExportFileDialog(this, _("Export ParaDif"), m_wxstrMeshDir, wxEmptyString, "txt(*.txt)|*.txt", wxFD_SAVE | wxFD_OVERWRITE_PROMPT);
	if (ExportFileDialog.ShowModal() == wxID_OK)
	{
		m_lpGLCanvas->m_lpMeshCenter->m_fnSolveParaScale();
		//m_lpGLCanvas->m_lpMeshCenter->m_fnAdjustParaScale();
		ofs.open(ExportFileDialog.GetPath().ToAscii());
		m_lpGLCanvas->m_lpMeshCenter->m_fnExportParaDif(ofs);
		ofs.close();
	}

}

void CMainFrame::OnImportCones(wxCommandEvent& event)
{
	std::ifstream ifs;
	wxFileDialog ImportFileDialog(this, _("ImportCones"), "", "", "TextFiles(*.txt)|*.txt", wxFD_OPEN | wxFD_FILE_MUST_EXIST);
	CDirectionAdjustor DA(m_lpGLCanvas->m_lpMeshCenter->m_mshSurface);
	CDirectionGenerator DG(m_lpGLCanvas->m_lpMeshCenter->m_mshSurface);
	CMeshBase::CLDLTSv LDLTS;
	if (ImportFileDialog.ShowModal() == wxID_OK)
	{
		ifs.open(ImportFileDialog.GetPath().ToAscii());
		m_lpGLCanvas->m_lpMeshCenter->m_fnImportIntCones(ifs);
		ifs.close();
		DA.m_fnInitDensityMatrix(LDLTS);
		LDLTS.decompose();
		DA.m_fnFillPoisson(LDLTS);
		LDLTS.solve();
		DA.m_fnAssignDensity(LDLTS);
		DG.m_fnGenerateDirFromDensity();
		m_lpGLCanvas->m_lpMeshCenter->m_fnDetectSeamType();
		m_lpGLCanvas->m_lpMeshCenter->m_fnMarkSingularity();
		Refresh();
	}
}

void CMainFrame::OnExportSingularities(wxCommandEvent &event)
{
	std::ofstream ofs;
	wxString wxstrDefaultFile = m_wxstrMeshName;
	wxstrDefaultFile.Replace(".off", "_singularities.txt");

	wxFileDialog ExportFileDialog(this, _("ExportSingularities"), m_wxstrMeshDir, wxstrDefaultFile, "txt(*.txt)|*.txt", wxFD_SAVE | wxFD_OVERWRITE_PROMPT);
	if (ExportFileDialog.ShowModal() == wxID_OK)
	{
		ofs.open(ExportFileDialog.GetPath().ToAscii());
		m_lpGLCanvas->m_lpMeshCenter->m_fnExportSingularities(ofs);
		ofs.close();
	}
}

void CMainFrame::OnExportFixedDirection(wxCommandEvent &event)
{
	std::ofstream ofs;
	wxFileDialog ExportFileDialog(this, _("ExportFixedDirection"), "", "", "Text Files(*.txt)|*.txt", wxFD_SAVE | wxFD_OVERWRITE_PROMPT);
	if (ExportFileDialog.ShowModal() == wxID_OK)
	{
		ofs.open(ExportFileDialog.GetPath().ToAscii());
		m_lpGLCanvas->m_lpMeshCenter->m_fnExportFixedDirection(ofs);
		ofs.close();
	}
}

void CMainFrame::OnExportHardFeature(wxCommandEvent &event)
{
	std::ofstream ofs;
	wxFileDialog ExportFileDialog(this, _("ExportHardFeature"), "", "", "Text Files(*.txt)|*.txt", wxFD_SAVE | wxFD_OVERWRITE_PROMPT);
	if (ExportFileDialog.ShowModal() == wxID_OK)
	{
		ofs.open(ExportFileDialog.GetPath().ToAscii());
		m_lpGLCanvas->m_lpMeshCenter->m_fnExportHardFeature(ofs);
		ofs.close();
	}
}

void CMainFrame::OnExportFixedMark(wxCommandEvent &event)
{
	std::ofstream ofs;
	wxFileDialog ExportFileDialog(this, _("ExportFixedMark"), "", "", "Text Files(*.txt)|*.txt", wxFD_SAVE | wxFD_OVERWRITE_PROMPT);
	if (ExportFileDialog.ShowModal() == wxID_OK)
	{
		ofs.open(ExportFileDialog.GetPath().ToAscii());
		m_lpGLCanvas->m_lpMeshCenter->m_fnExportFixedMark(ofs);
		ofs.close();
	}
}

void CMainFrame::OnImportFixedDirection(wxCommandEvent &event)
{
	std::ifstream ifs;
	wxFileDialog ImportFileDialog(this, _("ImportFixedDirection"), "", "", "TextFiles(*.txt)|*.txt", wxFD_OPEN | wxFD_FILE_MUST_EXIST);
	if (ImportFileDialog.ShowModal() == wxID_OK)
	{
		ifs.open(ImportFileDialog.GetPath().ToAscii());
		m_lpGLCanvas->m_lpMeshCenter->m_fnImportFixedDirection(ifs);
		ifs.close();
	}
	Refresh();
}


void CMainFrame::OnImportDirectionField(wxCommandEvent &event)
{
	std::ifstream ifs;
	wxFileDialog OpenFileDialog(this, _("ImportDirection"), "", "", "Text Files(*.txt)|*.txt", wxFD_OPEN | wxFD_FILE_MUST_EXIST);
	CDialogs::CFieldFileFormatChoiceDialog FFFCD(this);
	if (OpenFileDialog.ShowModal() == wxID_OK)
	{
		if (FFFCD.ShowModal() == wxID_OK)
		{
			ifs.open(OpenFileDialog.GetPath().ToAscii());
			m_lpGLCanvas->m_lpMeshCenter->m_fnImportDirectionField(ifs, FFFCD.bWithIndex, FFFCD.bFrameField);
			ifs.close();
			CDirectionAdjustor DA(m_lpGLCanvas->m_lpMeshCenter->m_mshSurface);
			m_lpGLCanvas->m_lpMeshCenter->m_fnDetectSeamType();
			m_lpGLCanvas->m_lpMeshCenter->m_fnMarkSingularity();
			m_lpGLCanvas->m_lpMeshCenter->m_fnGenerateMoment();
			DA.m_fnGenerateNablaField();
			Refresh();
		}
	}
}

void CMainFrame::OnImportBitmap(wxCommandEvent &event)
{
	wxFileDialog OpenFileDialog(this, _("ImportBitmap"), "", "", "Bmp Files (*.bmp)|*.bmp", wxFD_OPEN | wxFD_FILE_MUST_EXIST);
	if (OpenFileDialog.ShowModal() == wxID_OK)
	{
		m_lpGLCanvas->m_bmpOriginal = wxBitmap(OpenFileDialog.GetPath(), wxBITMAP_TYPE_BMP);
	}
	Refresh();
}

void CMainFrame::OnExportLIC(wxCommandEvent &event)
{
	wxFileDialog SaveFileDialog(this, _("Export LIC"), "", "", "Bmp Files (*.bmp)|*.bmp", wxFD_SAVE | wxFD_OVERWRITE_PROMPT);
	if (SaveFileDialog.ShowModal() == wxID_OK)
	{
		m_lpGLCanvas->m_bmpLIC.SaveFile(SaveFileDialog.GetPath(), wxBITMAP_TYPE_BMP);
	}
}

//void CMainFrame::OnExportElectricField(wxCommandEvent &event)
//{
//	std::ofstream ofs;
//	wxFileDialog ExportFileDialog(this, _("ExportElectricField"), "", "", "MFiles(*.m)|*.m", wxFD_SAVE | wxFD_OVERWRITE_PROMPT);
//	if (ExportFileDialog.ShowModal() == wxID_OK)
//	{
//		ofs.open(ExportFileDialog.GetPath().ToAscii());
//		m_lpGLCanvas->m_lpMeshCenter->m_fnExportMFile(ofs);
//		ofs.close();
//		wxString wxstrFieldFileName = ExportFileDialog.GetPath();
//		wxstrFieldFileName.Replace(_(".m"), _("_electric_field.txt"));
//		ofs.open(wxstrFieldFileName.ToAscii());
//		m_lpGLCanvas->m_lpMeshCenter->m_fnExportElectricityField(ofs);
//		ofs.close();
//	}
//}

void CMainFrame::OnExportCrevasse(wxCommandEvent &event)
{
	std::ofstream ofs;
	wxFileDialog ExportFileDialog(this, _("Export Texture Para"), "", "", "OBJFiles(*.obj)|*.obj", wxFD_SAVE | wxFD_OVERWRITE_PROMPT);
    if (ExportFileDialog.ShowModal() == wxID_OK)
    {
        ofs.open(ExportFileDialog.GetPath().ToAscii());
        m_lpGLCanvas->m_lpMeshCenter->m_fnExportCrevasse(ofs);
        ofs.close();
    }
}

void CMainFrame::OnExportTexturePara(wxCommandEvent &event)
{
	std::ofstream ofs;
	CDialogs::CTextureDialog TD(this);
	wxString wxstrDefaultFile;
	wxstrDefaultFile = m_wxstrMeshName;
	wxstrDefaultFile.Replace(".off", ".obj");
	wxFileDialog ExportFileDialog(this, _("Export Texture Para"), m_wxstrMeshDir + _("\\..\\field"), wxstrDefaultFile, "OBJFiles(*.obj)|*.obj", wxFD_SAVE | wxFD_OVERWRITE_PROMPT);
	if (TD.ShowModal() == wxID_OK)
	{
		if (ExportFileDialog.ShowModal() == wxID_OK)
		{
			ofs.open(ExportFileDialog.GetPath().ToAscii());
			if (TD.m_bSelectedAreaTexture)
			{
				if (TD.m_bIncludingOtherVertices)
				{
					m_lpGLCanvas->m_lpMeshCenter->m_fnExportPartPara(ofs, TD.m_bNormalizeTexture);
				}
				else
				{
					m_lpGLCanvas->m_lpMeshCenter->m_fnExportSelected(ofs, TD.m_bNormalizeTexture);
				}
			}
			else
			{
				m_lpGLCanvas->m_lpMeshCenter->m_fnExportTexturePara(ofs, TD.m_bNormalizeTexture);
			}
			ofs.close();
		}
	}
}

void CMainFrame::OnExportAngleError(wxCommandEvent &event)
{
	std::ofstream ofs;
	CMeshBase::FT lpftErr[9];
	wxFileDialog ExportFileDialog(this, _("Export Angle Error"), "", "", "OBJFiles(*.obj)|*.obj", wxFD_SAVE | wxFD_OVERWRITE_PROMPT);
	if (ExportFileDialog.ShowModal() == wxID_OK)
	{
		m_lpGLCanvas->m_lpMeshCenter->m_fnEstimateAngleError(lpftErr);
		CDialogs::CAngleErrDialog AED(this, lpftErr);
		if (AED.ShowModal() == wxID_OK)
		{
			ofs.open(ExportFileDialog.GetPath().ToAscii());
			m_lpGLCanvas->m_lpMeshCenter->m_fnExportAngleError(ofs, AED.m_iSelection, AED.m_lpdblErr[AED.m_iSelection * 3], AED.m_lpdblErr[AED.m_iSelection * 3 + 2]);
			ofs.close();
		}
	}
}

void CMainFrame::OnExportConfCoefColor(wxCommandEvent &event)
{
	std::ofstream ofs;
	wxFileDialog ExportFileDialog(this, _("Export Conformal Coefficient Color"), "", "", "OBJFiles(*.obj)|*.obj", wxFD_SAVE | wxFD_OVERWRITE_PROMPT);
	if (ExportFileDialog.ShowModal() == wxID_OK)
	{
		ofs.open(ExportFileDialog.GetPath().ToAscii());
		m_lpGLCanvas->m_lpMeshCenter->m_fnExportConfCoefColor(ofs);
		ofs.close();
	}
}

void CMainFrame::OnOpenMovieMode(wxCommandEvent &event)
{
	wxMessageDialog MD(this, _("Draw Discrete Bitmap?"), _("Warning"), wxYES_NO | wxCANCEL);
	wxDirDialog DirDlg(this, _("Choose Movie File Directory"), _(""), wxDD_DEFAULT_STYLE);
	int iMovieType;
	if (m_lpGLCanvas->m_bMovieMode)
	{
		m_lpGLCanvas->m_bMovieMode = false;
		m_lpGLCanvas->m_wxstrMovieDir.clear();
	}
	else
	{
		iMovieType = MD.ShowModal();
		if (iMovieType != wxCANCEL)
		{
			if (DirDlg.ShowModal() == wxID_OK)
			{
				m_lpGLCanvas->m_bMovieMode = true;
				m_lpGLCanvas->m_wxstrMovieDir = DirDlg.GetPath();
				if (iMovieType == wxID_YES)
				{
					m_lpGLCanvas->m_bDiscrete = true;
				}
				else
				{
					m_lpGLCanvas->m_bDiscrete = false;
				}
			}
		}
	}
}

void CMainFrame::OnOpenRecordMode(wxCommandEvent &event)
{
	if (m_lpGLCanvas->m_bMovieMode)
	{
		m_lpGLCanvas->m_bRecordMode = false;
		m_lpGLCanvas->m_wxstrRecordFile.clear();
	}
	else
	{
		wxFileDialog SaveFileDialog(this, _("Choose Record File"), "", "", "Text Files(*.txt)|*.txt", wxFD_SAVE | wxFD_OVERWRITE_PROMPT);
		if (SaveFileDialog.ShowModal() == wxID_OK)
		{
			m_lpGLCanvas->m_bRecordMode = true;
			m_lpGLCanvas->m_wxstrRecordFile = SaveFileDialog.GetPath();
		}
	}
}

void CMainFrame::OnExportNewMesh(wxCommandEvent &event)
{
	std::ofstream ofs;
	char lpstrFileName[256], *lpstrFileType;
	wxFileDialog SaveFileDialog(this, _("SaveMeshFile"), "", "", "MeshFiles(*.off)|*.off|ObjFiles(*.obj)|*.obj", wxFD_SAVE | wxFD_OVERWRITE_PROMPT);
    if(m_lpGLCanvas->m_lpMeshCenter->m_mshRemeshedSurface.size_of_facets() > 0 && SaveFileDialog.ShowModal() == wxID_OK)
    {
		strcpy(lpstrFileName, SaveFileDialog.GetPath().To8BitData());
		lpstrFileType = strrchr(lpstrFileName, '.');
		if (strcmp(lpstrFileType, ".off") == 0 || strcmp(lpstrFileType, ".OFF") == 0)
		{
			ofs.open(SaveFileDialog.GetPath().ToAscii());
			m_lpGLCanvas->m_lpMeshCenter->m_fnWriteNewMesh(ofs);
			ofs.close();
		}
		if (strcmp(lpstrFileType, ".obj") == 0 || strcmp(lpstrFileType, ".OBJ") == 0)
		{
			ofs.open(SaveFileDialog.GetPath().ToAscii());
			m_lpGLCanvas->m_lpMeshCenter->m_fnExportQuadObjFile(ofs);
			ofs.close();
		}
    }
}

void CMainFrame::OnImportQuadMesh(wxCommandEvent &event)
{
	std::ifstream ifs;
    wxFileDialog OpenFileDialog(this, _("OpenQuadMesh"), "", "", "MeshFiles(*.off)|*.off", wxFD_OPEN | wxFD_FILE_MUST_EXIST);
    if(OpenFileDialog.ShowModal() == wxID_OK)
    {
		ifs.open(OpenFileDialog.GetPath().ToAscii());
		m_lpGLCanvas->m_lpMeshCenter->m_fnReadNewMesh(ifs);
		ifs.close();
		m_lpGLCanvas->m_lpMeshCenter->m_fnCountNewMeshSize();
		m_lpGLCanvas->m_lpMeshCenter->m_mshRemeshedSurface.refresh_degree();
		m_lpGLCanvas->m_lpMeshCenter->m_fnMarkNewMeshBorder();
		m_lpGLCanvas->m_lpMeshCenter->m_fnMarkNewMeshCrevasse();
		Refresh();
    }
}

void CMainFrame::OnExportQuadSingularities(wxCommandEvent &event)
{
	std::ofstream ofs;
	wxFileDialog SaveFileDialog(this, _("Export Singularities"), _(""), _(""), _("TextFiles(*.txt)|*.txt"), wxFD_SAVE | wxFD_OVERWRITE_PROMPT);
	if (SaveFileDialog.ShowModal() == wxID_OK)
	{
		ofs.open(SaveFileDialog.GetPath().To8BitData());
		m_lpGLCanvas->m_lpMeshCenter->m_fnExportQuadSingularities(ofs);
		ofs.close();
	}
}

void CMainFrame::OnExportQuadCylinders(wxCommandEvent &event)
{
	std::ofstream ofs;
	wxFileDialog SaveFileDialog(this, _("Export Cylinder"), _(""), _(""), _("TextFiles(*.txt)|*.txt"), wxFD_SAVE | wxFD_OVERWRITE_PROMPT);
	if (SaveFileDialog.ShowModal() == wxID_OK)
	{
		ofs.open(SaveFileDialog.GetPath().To8BitData());
		m_lpGLCanvas->m_lpMeshCenter->m_fnExportQuadCylinder(ofs);
		ofs.close();
	}
}

void CMainFrame::OnExit(wxCommandEvent &event)
{
    Destroy();
}

void CMainFrame::OnChangeViewPort(wxCommandEvent &event)
{
    CDialogs::CViewPortDialog::ViewPortData vpData;
    double dblAxisLength;
    vpData.dblAxisX = m_lpGLCanvas->m_lpdblRotateAxis[0];
    vpData.dblAxisY = m_lpGLCanvas->m_lpdblRotateAxis[1];
    vpData.dblAxisZ = m_lpGLCanvas->m_lpdblRotateAxis[2];
    vpData.dblRotateAngle = m_lpGLCanvas->m_dblRotateAngle * 180 / ftPi;
    vpData.dblDistance = m_lpGLCanvas->m_dblDistance;
    vpData.dblFarPlane = m_lpGLCanvas->m_dblFar;
    vpData.dblNearPlane = m_lpGLCanvas->m_dblNear;
    CDialogs::CViewPortDialog ViewPortDialog(this, &vpData);
    if(ViewPortDialog.ShowModal() == wxID_OK)
    {
        dblAxisLength = sqrt(vpData.dblAxisX * vpData.dblAxisX + vpData.dblAxisY * vpData.dblAxisY + vpData.dblAxisZ * vpData.dblAxisZ);
        m_lpGLCanvas->m_lpdblRotateAxis[0] = vpData.dblAxisX / dblAxisLength;
        m_lpGLCanvas->m_lpdblRotateAxis[1] = vpData.dblAxisY / dblAxisLength;
        m_lpGLCanvas->m_lpdblRotateAxis[2] = vpData.dblAxisZ / dblAxisLength;
        m_lpGLCanvas->m_dblRotateAngle = vpData.dblRotateAngle * ftPi / 180;
        m_lpGLCanvas->m_dblDistance = vpData.dblDistance;
        m_lpGLCanvas->m_dblFar = vpData.dblFarPlane;
        m_lpGLCanvas->m_dblNear = vpData.dblNearPlane;
        m_lpGLCanvas->m_fnRotate();
        m_lpGLCanvas->m_fnObserve();
        m_lpGLCanvas->m_fnCast();
        Refresh();
    }
}

void CMainFrame::OnAppendRotation(wxCommandEvent &event)
{
    CDialogs::CViewPortDialog::ViewPortData vpData;
    double dblAxisLength, lpdblRotateAxis[3], dblRotateAngle, lpdblNewAxis[3], dblNewRotateAngle;
    vpData.dblAxisX = 0;
    vpData.dblAxisY = 0;
    vpData.dblAxisZ = 1;
    vpData.dblRotateAngle = 0;
    vpData.dblDistance = m_lpGLCanvas->m_dblDistance;
    vpData.dblFarPlane = m_lpGLCanvas->m_dblFar;
    vpData.dblNearPlane = m_lpGLCanvas->m_dblNear;
    CDialogs::CViewPortDialog ViewPortDialog(this, &vpData);
    if(ViewPortDialog.ShowModal() == wxID_OK)
    {
        dblAxisLength = sqrt(vpData.dblAxisX * vpData.dblAxisX + vpData.dblAxisY * vpData.dblAxisY + vpData.dblAxisZ * vpData.dblAxisZ);
        lpdblRotateAxis[0] = vpData.dblAxisX / dblAxisLength;
        lpdblRotateAxis[1] = vpData.dblAxisY / dblAxisLength;
        lpdblRotateAxis[2] = vpData.dblAxisZ / dblAxisLength;
        dblRotateAngle = vpData.dblRotateAngle * ftPi / 180;
		CGeoPack::fnMergeRotation(lpdblNewAxis, dblNewRotateAngle, m_lpGLCanvas->m_lpdblRotateAxis, m_lpGLCanvas->m_dblRotateAngle, lpdblRotateAxis, dblRotateAngle);
		m_lpGLCanvas->m_dblRotateAngle = dblNewRotateAngle;
		m_lpGLCanvas->m_lpdblRotateAxis[0] = lpdblNewAxis[0];
		m_lpGLCanvas->m_lpdblRotateAxis[1] = lpdblNewAxis[1];
		m_lpGLCanvas->m_lpdblRotateAxis[2] = lpdblNewAxis[2];
        m_lpGLCanvas->m_fnRotate();
        m_lpGLCanvas->m_fnObserve();
        m_lpGLCanvas->m_fnCast();
        Refresh();
    }
}

void CMainFrame::OnOrthodox(wxCommandEvent &event)
{
	if (m_lpGLCanvas->m_bOrthodox)
	{
		m_lpGLCanvas->m_bOrthodox = false;
	}
	else
	{
		m_lpGLCanvas->m_bOrthodox = true;
	}
	m_lpGLCanvas->m_fnCast();
	Refresh();
}

void CMainFrame::OnShowEdge(wxCommandEvent &event)
{
	if (m_lpGLCanvas->m_bShowEdge)
	{
		m_lpGLCanvas->m_bShowEdge = false;
	}
	else
	{
		m_lpGLCanvas->m_bShowEdge = true;
	}
	Refresh();
}

void CMainFrame::OnShowLIC(wxCommandEvent &event)
{
	wxMessageDialog MD(this, _("Is this a frame field?"), _("Field Type"), wxYES_NO);
	if (m_lpGLCanvas->m_bShowLIC)
	{
		m_lpGLCanvas->m_bShowLIC = false;
	}
	else
	{
		CDirectionAdjustor DA(m_lpGLCanvas->m_lpMeshCenter->m_mshSurface);
		m_lpGLCanvas->m_lpMeshCenter->m_fnDetectSeamType();
		DA.m_fnGenerateMoment();
		DA.m_fnGenerateNablaField();
		//m_lpGLCanvas->m_lpMeshCenter->m_fnFillVertexElectricField();
		//m_lpGLCanvas->m_lpMeshCenter->m_fnFillVertexNablaTheta();
		//m_lpGLCanvas->m_lpMeshCenter->m_fnFillVertexCurlField();
		m_lpGLCanvas->m_lpMeshCenter->m_fnDetectFieldMaximum();
		CDialogs::CLICDialog LD(this);
		LD.dblMaxCrossGrad = m_lpGLCanvas->m_lpMeshCenter->m_ftMaxCrossGrad;
		LD.dblMaxElectricField = m_lpGLCanvas->m_lpMeshCenter->m_ftMaxElectricField;
		LD.dblMaxCurlField = m_lpGLCanvas->m_lpMeshCenter->m_ftMaxCurlField;
		if (LD.ShowModal() == wxID_OK)
		{
			CLICPainter LICPainter(m_lpGLCanvas);
			if (LD.bDiscrete)
			{
				LICPainter.m_fnGenerateViewCoordinates();
				LICPainter.m_fnTestSize(LD.dblSingularityScale);
				LICPainter.m_fnInitBuffer();
				LICPainter.m_fnSetLocations(false);
				LICPainter.m_fnDrawDiscrete();
			}
			else
			{
				switch (LD.iFieldChoice)
				{
				case 0:
					if (MD.ShowModal() == wxID_YES)
					{
						m_lpGLCanvas->m_lpMeshCenter->m_fnFillVertexFrameField();
					}
					else
					{
						m_lpGLCanvas->m_lpMeshCenter->m_fnFillVertexCrossField();
					}
					break;
				case 1:
					m_lpGLCanvas->m_lpMeshCenter->m_fnFillVertexNablaTheta();
					break;
				case 2:
					m_lpGLCanvas->m_lpMeshCenter->m_fnFillVertexElectricField();
					break;
				case 3:
					m_lpGLCanvas->m_lpMeshCenter->m_fnFillVertexCurlField();
					break;
				}
				LICPainter.m_fnGenerateViewCoordinates();
				LICPainter.m_fnTestSize(LD.dblSingularityScale);
				LICPainter.m_fnInitBuffer();
				if (LD.iFieldChoice == 0)
				{
					LICPainter.m_fnSetLocations(true);
					LICPainter.m_fnDetectSeam();
				}
				else
				{
					LICPainter.m_fnSetLocations(false);
				}
				LICPainter.m_fnAssignOriginalColor(m_lpGLCanvas->m_bmpOriginal);
				switch (LD.iFieldChoice)
				{
				case 0:
					LICPainter.m_fnCrossTrace(LD.dblTraceRadius);
					LICPainter.m_fnEnhanceBitmap(false);
					break;
				case 1:
					switch (LD.iScaleOption)
					{
					case 0:
						LICPainter.m_fnLinearTrace(LD.dblTraceRadius);
						LICPainter.m_fnEnhanceBitmap(false);
						break;
					case 1:
						LICPainter.m_fnLinearTrace(LD.dblTraceRadius);
						LICPainter.m_fnEnhanceBitmap(true, LD.dblMaxCrossGrad);
						break;
					case 2:
						LICPainter.m_fnLinearTrace(LD.dblTraceRadius);
						LICPainter.m_fnWeakFieldRender(LD.dblMaxCrossGrad);
						break;
					}
					break;
				case 2:
					switch (LD.iScaleOption)
					{
					case 0:
						LICPainter.m_fnLinearTrace(LD.dblTraceRadius);
						LICPainter.m_fnEnhanceBitmap(false);
						break;
					case 1:
						LICPainter.m_fnLinearTrace(LD.dblTraceRadius);
						LICPainter.m_fnEnhanceBitmap(true, LD.dblMaxElectricField);
						break;
					case 2:
						LICPainter.m_fnLinearTrace(LD.dblTraceRadius);
						LICPainter.m_fnWeakFieldRender(LD.dblMaxElectricField);
						break;
					}
					break;
				case 3:
					switch (LD.iScaleOption)
					{
					case 0:
						LICPainter.m_fnLinearTrace(LD.dblTraceRadius);
						LICPainter.m_fnEnhanceBitmap(false);
						break;
					case 1:
						LICPainter.m_fnLinearTrace(LD.dblTraceRadius);
						LICPainter.m_fnEnhanceBitmap(true, LD.dblMaxCurlField);
						break;
					case 2:
						LICPainter.m_fnLinearTrace(LD.dblTraceRadius);
						LICPainter.m_fnWeakFieldRender(LD.dblMaxCurlField);
						break;
					}
					break;
				}
			}
			if (LD.bWithSingularities)
			{
				LICPainter.m_fnAppendSingularities(LD.dblSingularityScale);
			}
			if (LD.bDiscrete)
			{
				LICPainter.m_fnDrawForce();
			}
			LICPainter.m_fnWriteBitmap(m_lpGLCanvas->m_bmpLIC, m_lpGLCanvas->m_iBmpX, m_lpGLCanvas->m_iBmpY);
		}
		m_lpGLCanvas->m_bShowLIC = true;
	}
	Refresh();
}

void CMainFrame::OnChartViewable(wxCommandEvent &event)
{
	m_lpGLCanvas->m_bChartViewable = m_lpToolBar->GetToolState(idShowChart);
    Refresh();
}

void CMainFrame::OnFrameField(wxCommandEvent &event)
{
	m_lpGLCanvas->m_bFrameField = m_lpToolBar->GetToolState(idShowFrameField);
    Refresh();
}


void CMainFrame::OnSingularityViewable(wxCommandEvent &event)
{
    if(m_lpGLCanvas->m_bSingularityViewable)
    {
        m_lpGLCanvas->m_bSingularityViewable = false;
    }
    else
    {
        m_lpGLCanvas->m_bSingularityViewable = true;
    }
    Refresh();
}

void CMainFrame::OnShowElectricFieldSum(wxCommandEvent &event)
{
	switch (m_lpGLCanvas->m_iElectricFieldViewable)
	{
	case 0:
	case 2:
		m_lpGLCanvas->m_iElectricFieldViewable = m_lpGLCanvas->m_iElectricFieldViewable + 1;
		break;
	case 1:
	case 3:
		m_lpGLCanvas->m_iElectricFieldViewable = m_lpGLCanvas->m_iElectricFieldViewable - 1;
		break;
	}
    Refresh();
}


void CMainFrame::OnShowElectricFieldDif(wxCommandEvent &event)
{
	switch (m_lpGLCanvas->m_iElectricFieldViewable)
	{
	case 0:
	case 1:
		m_lpGLCanvas->m_iElectricFieldViewable = m_lpGLCanvas->m_iElectricFieldViewable + 2;
		break;
	case 2:
	case 3:
		m_lpGLCanvas->m_iElectricFieldViewable = m_lpGLCanvas->m_iElectricFieldViewable - 2;
		break;
	}
	Refresh();
}

void CMainFrame::OnTrackViewable(wxCommandEvent &event)
{
	if (m_lpGLCanvas->m_bTrackViewable)
	{
		m_lpGLCanvas->m_bTrackViewable = false;
	}
	else
	{
		m_lpGLCanvas->m_bTrackViewable = true;
	}
	Refresh();
}

void CMainFrame::OnShowGaussCurvature(wxCommandEvent &event)
{
	if (m_lpToolBar->GetToolState(idShowGaussCurvature))
	{
		m_lpToolBar->ToggleTool(idShowConformalCoefficient, false);
		m_lpToolBar->ToggleTool(idShowAngleError, false);
		m_lpGLCanvas->m_lpMeshCenter->m_fnRepresentGaussCurvature();
		m_lpGLCanvas->m_iFacetJetColor = 1;
	}
	else
	{
		m_lpGLCanvas->m_iFacetJetColor = 0;
	}
	Refresh();
}


void CMainFrame::OnShowConformalCoefficient(wxCommandEvent &event)
{
	if (m_lpToolBar->GetToolState(idShowConformalCoefficient))
	{
		m_lpToolBar->ToggleTool(idShowGaussCurvature, false);
		m_lpToolBar->ToggleTool(idShowAngleError, false);
		m_lpGLCanvas->m_lpMeshCenter->m_fnRepresentConformalCoefficient();
		m_lpGLCanvas->m_iFacetJetColor = 1;
	}
	else
	{
		m_lpGLCanvas->m_iFacetJetColor = 0;
	}
	Refresh();
}

void CMainFrame::OnShowAngleError(wxCommandEvent &event)
{
	CMeshBase::FT lpftErr[9];
	if (m_lpToolBar->GetToolState(idShowAngleError))
	{
		m_lpGLCanvas->m_lpMeshCenter->m_fnEstimateAngleError(lpftErr);
		CDialogs::CAngleErrDialog AED(this, lpftErr);
		if (AED.ShowModal() == wxID_OK)
		{
			m_lpToolBar->ToggleTool(idShowConformalCoefficient, false);
			m_lpToolBar->ToggleTool(idShowGaussCurvature, false);
			m_lpGLCanvas->m_lpMeshCenter->m_fnRepresentParaCurl(AED.m_iSelection, AED.m_lpdblErr);
			m_lpGLCanvas->m_iFacetJetColor = 2;
		}
		else
		{
			m_lpToolBar->ToggleTool(idShowAngleError, false);
		}
	}
	else
	{
		m_lpGLCanvas->m_iFacetJetColor = 0;
	}
	Refresh();
}

void CMainFrame::OnShowParaCurl(wxCommandEvent& event)
{
	CMeshBase::FT lpftCurl[9];
	CParameterGenerator PG(m_lpGLCanvas->m_lpMeshCenter->m_mshSurface, m_lpGLCanvas->m_lpMeshCenter->m_vecCrevasses);
	if (m_lpToolBar->GetToolState(idShowParaCurl))
	{
		PG.m_fnEstimateParaCurl(lpftCurl);
		CDialogs::CAngleErrDialog AED(this, lpftCurl);
		if (AED.ShowModal() == wxID_OK)
		{
			m_lpToolBar->ToggleTool(idShowConformalCoefficient, false);
			m_lpToolBar->ToggleTool(idShowGaussCurvature, false);
			m_lpGLCanvas->m_lpMeshCenter->m_fnRepresentParaCurl(AED.m_iSelection, AED.m_lpdblErr);
			m_lpGLCanvas->m_iFacetJetColor = 2;
		}
		else
		{
			m_lpToolBar->ToggleTool(idShowParaCurl, false);
		}
	}
	else
	{
		m_lpGLCanvas->m_iFacetJetColor = 0;
	}
	Refresh();
}

void CMainFrame::OnDensityU(wxCommandEvent &event)
{
	if (m_lpToolBar->GetToolState(idShowDensityU))
	{
		m_lpToolBar->ToggleTool(idShowDensityV, false);
		m_lpToolBar->ToggleTool(idShowGlobalU, false);
		m_lpToolBar->ToggleTool(idShowGlobalV, false);
		m_lpGLCanvas->m_iViewableColor = 4;
	}
	else
	{
		m_lpGLCanvas->m_iViewableColor = 0;
	}
	Refresh();
}


void CMainFrame::OnDensityV(wxCommandEvent &event)
{
	if (m_lpToolBar->GetToolState(idShowDensityV))
	{
		m_lpToolBar->ToggleTool(idShowDensityU, false);
		m_lpToolBar->ToggleTool(idShowGlobalU, false);
		m_lpToolBar->ToggleTool(idShowGlobalV, false);
		m_lpGLCanvas->m_iViewableColor = 5;
	}
	else
	{
		m_lpGLCanvas->m_iViewableColor = 0;
	}
	Refresh();
}

void CMainFrame::OnLocalU(wxCommandEvent &event)
{
	if (m_lpToolBar->GetToolState(idShowLocalU))
	{
		m_lpToolBar->ToggleTool(idShowLocalV, false);
		m_lpGLCanvas->m_iViewableLocalPara = 0;
	}
	else
	{
		m_lpGLCanvas->m_iViewableLocalPara = -1;
	}
	Refresh();
}

void CMainFrame::OnLocalV(wxCommandEvent &event)
{
	if (m_lpToolBar->GetToolState(idShowLocalV))
	{
		m_lpToolBar->ToggleTool(idShowLocalU, false);
		m_lpGLCanvas->m_iViewableLocalPara = 1;
	}
	else
	{
		m_lpGLCanvas->m_iViewableLocalPara = -1;
	}
	Refresh();
}


void CMainFrame::OnGlobalU(wxCommandEvent &event)
{
	if (m_lpToolBar->GetToolState(idShowGlobalU))
	{
		m_lpToolBar->ToggleTool(idShowDensityU, false);
		m_lpToolBar->ToggleTool(idShowDensityV, false);
		if ((m_lpGLCanvas->m_iViewableColor & 4) != 0)
		{
			m_lpGLCanvas->m_iViewableColor -= 4;
		}
		m_lpGLCanvas->m_iViewableColor += 1;
	}
	else
	{
		m_lpGLCanvas->m_iViewableColor -= 1;
	}
	Refresh();
}

void CMainFrame::OnGlobalV(wxCommandEvent &event)
{
	if (m_lpToolBar->GetToolState(idShowGlobalV))
	{
		m_lpToolBar->ToggleTool(idShowDensityU, false);
		m_lpToolBar->ToggleTool(idShowDensityV, false);
		if ((m_lpGLCanvas->m_iViewableColor & 4) != 0)
		{
			m_lpGLCanvas->m_iViewableColor -= 4;
		}
		m_lpGLCanvas->m_iViewableColor += 2;
	}
	else
	{
		m_lpGLCanvas->m_iViewableColor -= 2;
	}
	Refresh();
}


void CMainFrame::OnShowSamplePoints(wxCommandEvent &event)
{
	if (m_lpToolBar->GetToolState(idShowSamplePoints))
	{
		m_lpToolBar->ToggleTool(idShowLinkages, false);
		m_lpToolBar->ToggleTool(idShowFacetNormal, false);
		m_lpToolBar->ToggleTool(idShowNewMesh, false);
		m_lpGLCanvas->m_iViewableMeshInfo = 0;
	}
	else
	{
		m_lpGLCanvas->m_iViewableMeshInfo = -1;
	}
    Refresh();
}

void CMainFrame::OnShowLinkages(wxCommandEvent &event)
{
	if (m_lpToolBar->GetToolState(idShowLinkages))
	{
		m_lpToolBar->ToggleTool(idShowSamplePoints, false);
		m_lpToolBar->ToggleTool(idShowFacetNormal, false);
		m_lpToolBar->ToggleTool(idShowNewMesh, false);
		m_lpGLCanvas->m_iViewableMeshInfo = 1;
	}
	else
	{
		m_lpGLCanvas->m_iViewableMeshInfo = -1;
	}
	Refresh();
}


void CMainFrame::OnShowFacetDirections(wxCommandEvent &event)
{
	if (m_lpToolBar->GetToolState(idShowFacetNormal))
	{
		m_lpToolBar->ToggleTool(idShowSamplePoints, false);
		m_lpToolBar->ToggleTool(idShowLinkages, false);
		m_lpToolBar->ToggleTool(idShowNewMesh, false);
		m_lpGLCanvas->m_iViewableMeshInfo = 2;
	}
	else
	{
		m_lpGLCanvas->m_iViewableMeshInfo = -1;
	}
	Refresh();
}

void CMainFrame::OnShowReconstructedMesh(wxCommandEvent &event)
{
	if (m_lpToolBar->GetToolState(idShowNewMesh))
	{
		m_lpToolBar->ToggleTool(idShowSamplePoints, false);
		m_lpToolBar->ToggleTool(idShowLinkages, false);
		m_lpToolBar->ToggleTool(idShowFacetNormal, false);
		m_lpGLCanvas->m_iViewableMeshInfo = 3;
	}
	else
	{
		m_lpGLCanvas->m_iViewableMeshInfo = -1;
	}
    Refresh();
}

void CMainFrame::OnNormalizeModal(wxCommandEvent &event)
{
	wxMessageDialog MD(this, "Are you sure to normalize this Modal?", _("Command Confirmation"), wxOK | wxCANCEL);
	if (MD.ShowModal() == wxID_OK)
	{
		if (m_lpGLCanvas->m_lpMeshCenter->m_psProcessStatus >= PS_GEOMETRY_FILLED)
		{
			m_lpGLCanvas->m_lpMeshCenter->m_fnStatusRetreat(PS_GEOMETRY_FILLED);
			m_lpGLCanvas->m_lpMeshCenter->m_fnNormalizeMesh();
			m_lpGLCanvas->m_lpMeshCenter->m_fnFillGeometricInfo();
		}
		Refresh();
	}
}


void CMainFrame::OnInitilizeChart(wxCommandEvent &event)
{
	CDialogs::CChartInitDialog CID(this);
    if(m_lpGLCanvas->m_lpMeshCenter->m_psProcessStatus >= PS_GEOMETRY_FILLED)
    {
		if (CID.ShowModal() == wxID_OK)
		{
			m_lpGLCanvas->m_lpMeshCenter->m_fnStatusRetreat(PS_GEOMETRY_FILLED);
			CDirectionGenerator DG(m_lpGLCanvas->m_lpMeshCenter->m_mshSurface);
			if (CID.m_iChoice == 0)
			{
				//DG.m_fnMarkSharpFeature();
				DG.m_fnInitFixedDirection();
				DG.m_fnExtendDirection();
			}
			else
			{
				DG.m_fnRandomizeDirection();
			}
			m_lpGLCanvas->m_lpMeshCenter->m_fnClearFreeEdge();
			m_lpGLCanvas->m_lpMeshCenter->m_fnDetectSeamType();
			m_lpGLCanvas->m_lpMeshCenter->m_fnMarkSingularity();
			++(m_lpGLCanvas->m_lpMeshCenter->m_psProcessStatus);
			Refresh();
		}
    }
}

void CMainFrame::OnEigenOptimal(wxCommandEvent &event)
{
	clock_t iStart, iEnd;
	wxMessageDialog MD(this, _("Are you sure to fill optimal direction?"), _("Command Confirmation"), wxOK | wxCANCEL);
	if (MD.ShowModal() == wxID_OK)
	{
		CDirectionGenerator DG(m_lpGLCanvas->m_lpMeshCenter->m_mshSurface);
		CMeshBase::CCmplxSv CmplxSv;
		iStart = clock();
		DG.m_fnFillComplex(CmplxSv);
		DG.m_fnSolveEigen(CmplxSv);
		DG.m_fnFillOptDir(CmplxSv);
		iEnd = clock();
		m_lpGLCanvas->m_lpMeshCenter->m_fnDetectSeamType();
		m_lpGLCanvas->m_lpMeshCenter->m_fnMarkSingularity();
		Refresh();
		wxMessageBox(wxString::Format("Complete!\nTime consumption: %d milliseconds", int(iEnd - iStart)), _("Tip"), wxOK, m_lpGLCanvas);
	}
}

void CMainFrame::OnFlatten(wxCommandEvent& event)
{
	wxMessageDialog MD(this, _("Are you sure to apply flatten algorithm?"), _("Command Confirmation"), wxOK | wxCANCEL);
	if (MD.ShowModal() == wxID_OK)
	{
		CDirectionGenerator DG(m_lpGLCanvas->m_lpMeshCenter->m_mshSurface);
		DG.m_fnMarkIsolatedCones();
		DG.m_fnGenerateMetric();
		Refresh();
	}
}

void CMainFrame::OnFixCones(wxCommandEvent& event)
{
	wxMessageDialog MD(this, _("Are you sure to apply flatten algorithm?"), _("Command Confirmation"), wxOK | wxCANCEL);
	if (MD.ShowModal() == wxID_OK)
	{
		CAlgorithmExecutor* lpAE = new CAlgorithmExecutor(m_lpGLCanvas, PT_FLATTEN);
		m_lpGLCanvas->m_Timer.Start(4);
		lpAE->Run();
	}
}

void CMainFrame::OnDeletePath(wxCommandEvent &event)
{
	wxMessageDialog MD(this, _("Are you sure to clear all linkages?"), _("Command Confirmation"), wxOK | wxCANCEL);
	if (MD.ShowModal() == wxID_OK)
	{
		m_lpGLCanvas->m_lpMeshCenter->m_fnClearFreeEdge();
		Refresh();
	}
}

void CMainFrame::OnComplexOptimize(wxCommandEvent &event)
{
	if (m_lpGLCanvas->m_lpMeshCenter->m_psProcessStatus >= PS_DIRECTION_SCATTER)
	{
		CDialogs::CIterDlg ID(this);
		if (ID.ShowModal() == wxID_OK)
		{
			m_lpGLCanvas->m_lpMeshCenter->m_fnStatusRetreat(PS_DIRECTION_SCATTER);
			CAlgorithmExecutor *lpAlgorithmExecutor = new CAlgorithmExecutor(m_lpGLCanvas, PT_COMPLEX);
			lpAlgorithmExecutor->m_nIterateTimes = ID.m_nIterTimes;
			m_lpGLCanvas->m_Timer.Start(4);
			lpAlgorithmExecutor->Run();
		}
	}
}


void CMainFrame::OnStopIteration(wxCommandEvent &event)
{
	wxMessageDialog MD(this, _("Are you sure to stop iteration?"), _("Command Confirmation"), wxOK | wxCANCEL);
	if (MD.ShowModal() == wxID_OK)
	{
		m_lpGLCanvas->m_bProcessing = false;
		Refresh();
	}
}


void CMainFrame::OnDirectionAdj(wxCommandEvent &event)
{
	CDirectionAdjustor DA(m_lpGLCanvas->m_lpMeshCenter->m_mshSurface);
	CMeshBase::CLDLTSv LDLTS;
	wxMessageDialog MD(this, _("Are you sure to adjust direction?"), _("Command Confirmation"), wxOK | wxCANCEL);
    if(m_lpGLCanvas->m_lpMeshCenter->m_psProcessStatus >= PS_DIRECTION_SCATTER)
    {
		if (MD.ShowModal() == wxID_OK)
        {
			m_lpGLCanvas->m_lpMeshCenter->m_fnStatusRetreat(PS_DIRECTION_SCATTER);
			DA.m_fnGenerateMoment();
			DA.m_fnInitDirectionMatrix(LDLTS);
			DA.m_fnFillMoment(LDLTS);
			LDLTS.decompose();
			LDLTS.solve();
			DA.m_fnAdjustDirection(LDLTS);
			LDLTS.clear();
			DA.m_fnGenerateMoment();
			m_lpGLCanvas->m_lpMeshCenter->m_fnDetectSeamType();
			m_lpGLCanvas->m_lpMeshCenter->m_fnMarkSingularity();
			DA.m_fnInitDensityMatrix(LDLTS);
			DA.m_fnFillPoisson(LDLTS);
			LDLTS.decompose();
			LDLTS.solve();
			DA.m_fnAssignDensity(LDLTS);
			LDLTS.clear();
			DA.m_fnGenerateNablaField();
			m_lpGLCanvas->m_lpMeshCenter->m_fnRepresentElectricField();
			m_lpGLCanvas->m_lpMeshCenter->m_fnFillVertexElectricField();
			m_lpGLCanvas->m_lpMeshCenter->m_fnFillVertexCurlField();
			Refresh();
        }
    }
}


void CMainFrame::OnAdjustSingularities(wxCommandEvent &event)
{
	CDialogs::CAdjustTypeDialog ATD(this);
	if (m_lpGLCanvas->m_lpMeshCenter->m_psProcessStatus >= PS_DIRECTION_SCATTER)
	{
        if(ATD.ShowModal() == wxID_OK)
        {
			m_lpGLCanvas->m_lpMeshCenter->m_fnStatusRetreat(PS_DIRECTION_SCATTER);
            m_lpGLCanvas->m_Timer.Start(4);
            CAlgorithmExecutor *lpAlgorithmExecutor = new CAlgorithmExecutor(m_lpGLCanvas, PT_SINGULARITY);
			lpAlgorithmExecutor->m_iAdjustType = ATD.m_iLocalOptimizeType;
            lpAlgorithmExecutor->Run();
        }
    }
}

void CMainFrame::OnSetRadius(wxCommandEvent &event)
{
    double dblScale;
	wxTextEntryDialog TED(this, _("Radius"), _("Number Input"), _("-2.0"), wxOK | wxCANCEL);
	if (TED.ShowModal() == wxID_OK)
	{
		TED.GetValue().ToDouble(&dblScale);
		m_lpGLCanvas->m_lpMeshCenter->m_ftChargeRadius = m_lpGLCanvas->m_lpMeshCenter->m_ftGeometricDimension * exp(dblScale * log(10));
	}
}

void CMainFrame::OnReportEnergy(wxCommandEvent &event)
{
	CDialogs::CEnergyDialog ED(this);
	CMeshBase::FT ftEnergy;
	CDirectionAdjustor DA(m_lpGLCanvas->m_lpMeshCenter->m_mshSurface);
	if (ED.ShowModal() == wxID_OK)
	{
		switch (ED.m_iEnergyType)
		{
		case 0:
			ftEnergy = DA.m_fnEnergy(m_lpGLCanvas->m_lpMeshCenter->m_ftChargeRadius);
			break;
		case 1:
			ftEnergy = DA.m_fnSmoothEnergy(m_lpGLCanvas->m_lpMeshCenter->m_ftChargeRadius);
			break;
		case 2:
			ftEnergy = DA.m_fnCurlEnergy();
			break;
		case 3:
			DA.m_fnGenerateMoment();
			ftEnergy = DA.m_fnDiscreteCurl();
			break;
		case 4:
			ftEnergy = m_lpGLCanvas->m_lpMeshCenter->m_fnMIQEnergy();
			break;
		case 5:
			ftEnergy = m_lpGLCanvas->m_lpMeshCenter->m_fnAlignEnergy();
			break;
		}
		wxString wxstrEnergy = wxString::FromDouble(ftEnergy);
		wxTextEntryDialog TED(this, _("Energy"), _("Energy Report"), wxstrEnergy, wxOK);
		TED.ShowModal();
	}
}

void CMainFrame::OnReportSingularities(wxCommandEvent &event)
{
	int iPositive, iNegative;
	wxString wxstrMessage;
	m_lpGLCanvas->m_lpMeshCenter->m_fnReportSingularities(iPositive, iNegative);
	wxstrMessage = wxString::Format(wxString("Positive: %d\nNegative: %d"), iPositive, iNegative);
	wxMessageDialog MD(this, wxstrMessage, _("Number Report"), wxOK | wxCANCEL);
	MD.ShowModal();
}

void CMainFrame::OnLinkPoles(wxCommandEvent &event)
{
	wxMessageDialog MD(this, _("Are you sure to link poles?"), _("Command Confirmation"), wxOK | wxCANCEL);
	CMeshBase::CPolyhedron::Vertex_handle vhMin, vhMax;
	if (m_lpGLCanvas->m_lpMeshCenter->m_psProcessStatus >= PS_DIRECTION_SCATTER)
	{
		if (MD.ShowModal() == wxID_OK)
		{
			CDirectionAdjustor DA(m_lpGLCanvas->m_lpMeshCenter->m_mshSurface);
			DA.m_fnFindPoles(m_lpGLCanvas->m_lpMeshCenter->m_ftChargeRadius, vhMin, vhMax);
			vhMin->iSelected = 2;
			vhMax->iSelected = 1;
			DA.m_fnLinkPoles(vhMin, vhMax);
			Refresh();
		}
	}
}

void CMainFrame::OnMoveSingularity(wxCommandEvent &event)
{
	wxMessageDialog MD(this, _("Are you sure to move singularity?"), _("Command Confirmation"), wxOK | wxCANCEL);
    if(m_lpGLCanvas->m_lpMeshCenter->m_psProcessStatus >= PS_DIRECTION_SCATTER)
    {
		if (MD.ShowModal() == wxID_OK)
		{
			CDirectionAdjustor DA(m_lpGLCanvas->m_lpMeshCenter->m_mshSurface);
			m_lpGLCanvas->m_lpMeshCenter->m_fnStatusRetreat(PS_DIRECTION_SCATTER);
			DA.m_fnMoveSingularity();
			m_lpGLCanvas->m_lpMeshCenter->m_fnClearFreeEdge();
			m_lpGLCanvas->m_lpMeshCenter->m_fnGenerateMoment();
			m_lpGLCanvas->m_lpMeshCenter->m_fnMarkSingularity();
			Refresh();
		}
    }
}

void CMainFrame::OnClearBarrier(wxCommandEvent &event)
{
	wxMessageDialog MD(this, _("Are you sure to clear Barrier?"), _("Command Confirmation"), wxOK | wxCANCEL);
	if (MD.ShowModal() == wxID_OK)
	{
		m_lpGLCanvas->m_lpMeshCenter->m_fnClearBarrier();
	}
}

void CMainFrame::OnGlobalAdjust(wxCommandEvent &event)
{
	if (m_lpGLCanvas->m_lpMeshCenter->m_psProcessStatus >= PS_DIRECTION_SCATTER)
	{
		CDialogs::CGlobalAdjustDialog dlgGA(this);
		if (dlgGA.ShowModal() == wxID_OK)
		{
			m_lpGLCanvas->m_lpMeshCenter->m_fnStatusRetreat(PS_DIRECTION_SCATTER);
			CAlgorithmExecutor *lpAlgorithmExecutor = new CAlgorithmExecutor(m_lpGLCanvas, PT_GLOBAL);
			m_lpGLCanvas->m_lpMeshCenter->m_ftChargeRadius = m_lpGLCanvas->m_lpMeshCenter->m_ftGeometricDimension * exp(dlgGA.m_dblChargeRadius * log(10));
			lpAlgorithmExecutor->m_nIterateTimes = dlgGA.m_nTestTimes;
			lpAlgorithmExecutor->m_iAdjustType = dlgGA.m_iLocalOptimizeType;
			m_lpGLCanvas->m_Timer.Start(4);
			lpAlgorithmExecutor->Run();
		}
	}
}

void CMainFrame::OnCountSingularities(wxCommandEvent &event)
{
	CDialogs::CCountSingularityDialog CSD(this);
    wxFileDialog OpenFileDialog(this, _("Open Original CrossField"), "", "", "Cross Field(*.txt)|*.txt|All Files(*.*)|*.*", wxFD_OPEN | wxFD_FILE_MUST_EXIST);
	wxFileDialog SaveFileDialog(this, _("Export Data"), "", "", "Data Files(*.txt)|*.txt", wxFD_SAVE | wxFD_OVERWRITE_PROMPT);
	if (OpenFileDialog.ShowModal() == wxID_OK)
    {
        if (SaveFileDialog.ShowModal() == wxID_OK)
        {
            if (CSD.ShowModal() == wxID_OK)
            {
                CAlgorithmExecutor *lpAlgorithmExecutor = new CAlgorithmExecutor(m_lpGLCanvas, PT_COUNT_SINGULARITIES);
                lpAlgorithmExecutor->m_nTestTimes = (int)CSD.m_iTestTimes;
                lpAlgorithmExecutor->m_dblStartLevel = CSD.m_dblStartLevel;
                lpAlgorithmExecutor->m_dblDecrement = CSD.m_dblDecrement;
                lpAlgorithmExecutor->m_wxstrInputField = OpenFileDialog.GetPath();
                lpAlgorithmExecutor->m_wxstrOutputData = SaveFileDialog.GetPath();
                lpAlgorithmExecutor->m_iAdjustType = 1;
                m_lpGLCanvas->m_Timer.Start(4);
                lpAlgorithmExecutor->Run();
            }
        }
    }
}

void CMainFrame::OnMarkSharpFeature(wxCommandEvent &event)
{
	CDialogs::CAngleBiasDialog ABD(this);
	if (ABD.ShowModal() == wxID_OK)
	{
		CDirectionGenerator DG(m_lpGLCanvas->m_lpMeshCenter->m_mshSurface);
		DG.m_fnMarkSharpFeature(ABD.m_dblAngleBias * ftPi / 180);
		Refresh();
	}
}


void CMainFrame::OnClearFix(wxCommandEvent &event)
{
	wxMessageDialog MD(this, _("Are you sure to clear all the fix direction?"), _("Command Confirmation"), wxOK | wxCANCEL);
	if (MD.ShowModal() == wxID_OK)
	{
		m_lpGLCanvas->m_lpMeshCenter->m_fnClearFixedDirection();
	}
	Refresh();
}


void CMainFrame::OnDetectLoop(wxCommandEvent &event)
{
	CDirectionAdjustor DA(m_lpGLCanvas->m_lpMeshCenter->m_mshSurface);
	wxMessageDialog MD(this, _("Are you sure to detect loop?"), _("Command Confirmation"), wxOK | wxCANCEL);
	if (MD.ShowModal() == wxID_OK)
	{
		DA.m_fnMarkCurlPath();
		DA.m_fnTestCurlPath();
		Refresh();
	}
}

void CMainFrame::OnReleaseCurl(wxCommandEvent &event)
{
	CDirectionAdjustor DA(m_lpGLCanvas->m_lpMeshCenter->m_mshSurface);
	wxMessageDialog MD(this, _("Are you sure to release curl?"), _("Command Confirmation"), wxOK | wxCANCEL);
	if (MD.ShowModal() == wxID_OK)
	{
		DA.m_fnReleaseCurl();
		Refresh();
	}
}

void CMainFrame::OnCoarseDecurl(wxCommandEvent &event)
{
	if (m_lpGLCanvas->m_lpMeshCenter->m_psProcessStatus >= PS_DIRECTION_SCATTER)
	{
		CDialogs::CIterDlg ID(this);
		if (ID.ShowModal() == wxID_OK)
		{
			m_lpGLCanvas->m_lpMeshCenter->m_fnStatusRetreat(PS_DIRECTION_SCATTER);
			CAlgorithmExecutor *lpAlgorithmExecutor = new CAlgorithmExecutor(m_lpGLCanvas, PT_COARSE_DECURL);
			lpAlgorithmExecutor->m_nIterateTimes = ID.m_nIterTimes;
			m_lpGLCanvas->m_Timer.Start(4);
			lpAlgorithmExecutor->Run();
		}
	}
}

void CMainFrame::OnFineDecurl(wxCommandEvent &event)
{
	if (m_lpGLCanvas->m_lpMeshCenter->m_psProcessStatus >= PS_DIRECTION_SCATTER)
	{
		CDialogs::CIterDlg ID(this);
		if (ID.ShowModal() == wxID_OK)
		{
			m_lpGLCanvas->m_lpMeshCenter->m_fnStatusRetreat(PS_DIRECTION_SCATTER);
			CAlgorithmExecutor *lpAlgorithmExecutor = new CAlgorithmExecutor(m_lpGLCanvas, PT_FINE_DECURL);
			lpAlgorithmExecutor->m_nIterateTimes = ID.m_nIterTimes;
			m_lpGLCanvas->m_Timer.Start(4);
			lpAlgorithmExecutor->Run();
		}
	}
}


void CMainFrame::OnConstrainedDecurl(wxCommandEvent &event)
{
	if (m_lpGLCanvas->m_lpMeshCenter->m_psProcessStatus >= PS_DIRECTION_SCATTER)
	{
		CDialogs::CIterDlg ID(this);
		if (ID.ShowModal() == wxID_OK)
		{
			m_lpGLCanvas->m_lpMeshCenter->m_fnStatusRetreat(PS_DIRECTION_SCATTER);
			CAlgorithmExecutor *lpAlgorithmExecutor = new CAlgorithmExecutor(m_lpGLCanvas, PT_CONSTRAINED_DECURL);
			lpAlgorithmExecutor->m_nIterateTimes = ID.m_nIterTimes;
			m_lpGLCanvas->m_Timer.Start(4);
			lpAlgorithmExecutor->Run();
		}
	}
}

void CMainFrame::OnAdjustCurl(wxCommandEvent &event)
{
	CDirectionAdjustor DA(m_lpGLCanvas->m_lpMeshCenter->m_mshSurface);
	CMeshBase::CLDLTSv LDLTS0, LDLTS1;
	DA.m_fnInitDirectionMatrix(LDLTS0);
	LDLTS0.decompose();
	DA.m_fnInitDensityMatrix(LDLTS1);
	LDLTS1.decompose();
	DA.m_fnGenerateMoment();
	DA.m_fnFillMoment(LDLTS0);
	LDLTS0.solve();
	DA.m_fnAdjustDirection(LDLTS0);
	DA.m_fnGenerateMoment();
	DA.m_fnFillPoisson(LDLTS1);
	LDLTS1.solve();
	DA.m_fnAssignDensity(LDLTS1);
	DA.m_fnGenerateNablaField();
	DA.m_fnGenerateForces();
	DA.m_fnCurlAdjust();
	DA.m_fnGenerateMoment();
	DA.m_fnFillMoment(LDLTS0);
	LDLTS0.solve();
	DA.m_fnAdjustDirection(LDLTS0);
	DA.m_fnGenerateMoment();
	m_lpGLCanvas->m_lpMeshCenter->m_fnMarkSingularity();
	DA.m_fnFillPoisson(LDLTS1);
	LDLTS1.solve();
	DA.m_fnAssignDensity(LDLTS1);
	DA.m_fnGenerateNablaField();
	LDLTS1.clear();
	LDLTS0.clear();
	m_lpGLCanvas->m_lpMeshCenter->m_fnRepresentElectricField();
	Refresh();
}

void CMainFrame::OnTraceDirection(wxCommandEvent &event)
{
	wxMessageDialog MD(this, _("Are you sure to trace direction?"), _("Command Confirmation"), wxOK | wxCANCEL);
	if (m_lpGLCanvas->m_lpMeshCenter->m_psProcessStatus >= PS_DIRECTION_SCATTER)
	{
		if (MD.ShowModal() == wxID_OK)
		{
			m_lpGLCanvas->m_lpMeshCenter->m_fnStatusRetreat(PS_DIRECTION_SCATTER);
			m_lpGLCanvas->m_lpMeshCenter->m_Tracer.m_fnCountPort(m_lpGLCanvas->m_lpMeshCenter->m_mshSurface);
			m_lpGLCanvas->m_lpMeshCenter->m_Tracer.m_fnAccurateSearch();
			CDirectionAdjustor DA(m_lpGLCanvas->m_lpMeshCenter->m_mshSurface);
			DA.m_fnGenerateMoment();
			++(m_lpGLCanvas->m_lpMeshCenter->m_psProcessStatus);
			Refresh();
		}
	}
}

void CMainFrame::OnBiasSearch(wxCommandEvent &event)
{
	if (m_lpGLCanvas->m_lpMeshCenter->m_psProcessStatus >= PS_TRACK_PATH)
	{
		CDialogs::CAngleBiasDialog ABDlg(this);
		if (ABDlg.ShowModal() == wxID_OK)
		{
			m_lpGLCanvas->m_lpMeshCenter->m_fnStatusRetreat(PS_TRACK_PATH);
			m_lpGLCanvas->m_lpMeshCenter->m_Tracer.m_fnBiasSearch(ABDlg.m_dblAngleBias);
			Refresh();
		}
	}
}

void CMainFrame::OnAdjustPath(wxCommandEvent &event)
{
	long i;
	if (m_lpGLCanvas->m_lpMeshCenter->m_psProcessStatus >= PS_TRACK_PATH)
	{
		CDialogs::CIterDlg iterDlg(this, 10);
		if (iterDlg.ShowModal() == wxID_OK)
		{
			m_lpGLCanvas->m_lpMeshCenter->m_fnStatusRetreat(PS_TRACK_PATH);
			if (iterDlg.m_nIterTimes > 0)
			{
				for (i = 0; i < iterDlg.m_nIterTimes; ++i)
				{
					m_lpGLCanvas->m_lpMeshCenter->m_Tracer.m_fnAdjustPath();
					m_lpGLCanvas->m_lpMeshCenter->m_Tracer.m_fnRenewPath();
				}
			}
			Refresh();
		}
	}
}

void CMainFrame::OnSelectPoint(wxCommandEvent &event)
{
	m_lpGLCanvas->m_bSelectPoint = m_lpToolBar->GetToolState(idSelectPoint);
	if (m_lpGLCanvas->m_bSelectPoint)
	{
		m_lpToolBar->ToggleTool(idSelectTrack, false);
		m_lpToolBar->ToggleTool(idSelectCrevasse, false);
		m_lpToolBar->ToggleTool(idSelectFacet, false);
		m_lpToolBar->ToggleTool(idSelectEdge, false);
		m_lpGLCanvas->m_bSelectTrack = false;
		m_lpGLCanvas->m_bSelectCrevasse = false;
		m_lpGLCanvas->m_bSelectFacet = false;
		m_lpGLCanvas->m_bSelectEdge = false;
	}
	Refresh();
}


void CMainFrame::OnSelectTrack(wxCommandEvent &event)
{
	m_lpGLCanvas->m_bSelectTrack = m_lpToolBar->GetToolState(idSelectTrack);
	if (m_lpGLCanvas->m_bSelectTrack)
	{
		m_lpToolBar->ToggleTool(idSelectPoint, false);
		m_lpToolBar->ToggleTool(idSelectCrevasse, false);
		m_lpToolBar->ToggleTool(idSelectFacet, false);
		m_lpToolBar->ToggleTool(idSelectEdge, false);
		m_lpGLCanvas->m_bSelectPoint = false;
		m_lpGLCanvas->m_bSelectCrevasse = false;
		m_lpGLCanvas->m_bSelectFacet = false;
		m_lpGLCanvas->m_bSelectEdge = false;
	}
	Refresh();
}


void CMainFrame::OnSelectCrevasse(wxCommandEvent &event)
{
	m_lpGLCanvas->m_bSelectCrevasse = m_lpToolBar->GetToolState(idSelectCrevasse);
	if (m_lpGLCanvas->m_bSelectCrevasse)
	{
		m_lpToolBar->ToggleTool(idSelectPoint, false);
		m_lpToolBar->ToggleTool(idSelectTrack, false);
		m_lpToolBar->ToggleTool(idSelectFacet, false);
		m_lpToolBar->ToggleTool(idSelectEdge, false);
		m_lpGLCanvas->m_bSelectPoint = false;
		m_lpGLCanvas->m_bSelectTrack = false;
		m_lpGLCanvas->m_bSelectFacet = false;
		m_lpGLCanvas->m_bSelectEdge = false;
	}
	Refresh();
}

void CMainFrame::OnSelectFacet(wxCommandEvent &event)
{
	m_lpGLCanvas->m_bSelectFacet = m_lpToolBar->GetToolState(idSelectFacet);
	if (m_lpGLCanvas->m_bSelectFacet)
	{
		m_lpToolBar->ToggleTool(idSelectPoint, false);
		m_lpToolBar->ToggleTool(idSelectTrack, false);
		m_lpToolBar->ToggleTool(idSelectCrevasse, false);
		m_lpToolBar->ToggleTool(idSelectEdge, false);
		m_lpGLCanvas->m_bSelectPoint = false;
		m_lpGLCanvas->m_bSelectTrack = false;
		m_lpGLCanvas->m_bSelectCrevasse = false;
		m_lpGLCanvas->m_bSelectEdge = false;
	}
	Refresh();
}

void CMainFrame::OnSelectEdge(wxCommandEvent &event)
{
	m_lpGLCanvas->m_bSelectEdge = m_lpToolBar->GetToolState(idSelectEdge);
	if (m_lpGLCanvas->m_bSelectEdge)
	{
		m_lpToolBar->ToggleTool(idSelectPoint, false);
		m_lpToolBar->ToggleTool(idSelectTrack, false);
		m_lpToolBar->ToggleTool(idSelectCrevasse, false);
		m_lpToolBar->ToggleTool(idSelectFacet, false);
		m_lpGLCanvas->m_bSelectPoint = false;
		m_lpGLCanvas->m_bSelectTrack = false;
		m_lpGLCanvas->m_bSelectCrevasse = false;
		m_lpGLCanvas->m_bSelectFacet = false;
	}
	Refresh();
}

void CMainFrame::OnDecomposeMesh(wxCommandEvent &event)
{
	int iPositiveSingularities, iNegativeSingularities;
	wxMessageDialog MD(this, _("Are you sure to decompse mesh?"), _("Command Confirmation"), wxOK | wxCANCEL);
	if (m_lpGLCanvas->m_lpMeshCenter->m_psProcessStatus == PS_TRACK_PATH)
	{
		if (MD.ShowModal() == wxID_OK)
		{
			m_lpGLCanvas->m_lpMeshCenter->m_Tracer.m_fnMapToMesh();
			m_lpGLCanvas->m_lpMeshCenter->m_fnFreeBoundary();
			CParameterGenerator PG(m_lpGLCanvas->m_lpMeshCenter->m_mshSurface, m_lpGLCanvas->m_lpMeshCenter->m_vecCrevasses);
			++(m_lpGLCanvas->m_lpMeshCenter->m_psProcessStatus);
			Refresh();
		}
	}
	else if (m_lpGLCanvas->m_lpMeshCenter->m_psProcessStatus >= PS_DIRECTION_SCATTER)
	{
		if (MD.ShowModal() == wxID_OK)
		{
			m_lpGLCanvas->m_lpMeshCenter->m_fnStatusRetreat(PS_DIRECTION_SCATTER);
			CMeshSplitter MS(m_lpGLCanvas->m_lpMeshCenter->m_mshSurface);
			m_lpGLCanvas->m_lpMeshCenter->m_fnGenerateMoment();
			m_lpGLCanvas->m_lpMeshCenter->m_fnReportSingularities(iPositiveSingularities, iNegativeSingularities);
			if (iPositiveSingularities + iNegativeSingularities > 0)
			{
				MS.m_fnPartVertices();
				MS.m_fnSimplifyConnection();
				MS.m_fnMarkArea();
				MS.m_fnMergeArea();
			}
			MS.m_fnActiveBoundaries();
			m_lpGLCanvas->m_lpMeshCenter->m_psProcessStatus = PS_DIRECTION_DECOMPOSE;
			Refresh();
		}
	}
}

void CMainFrame::OnUniformizeChart(wxCommandEvent &event)
{
	wxMessageDialog MD(this, _("Are you sure to uniformize chart?"), _("Command Confirmation"), wxOK | wxCANCEL);
    if(m_lpGLCanvas->m_lpMeshCenter->m_psProcessStatus >= PS_DIRECTION_DECOMPOSE)
    {
		if (MD.ShowModal() == wxID_OK)
		{
			m_lpGLCanvas->m_lpMeshCenter->m_fnStatusRetreat(PS_DIRECTION_DECOMPOSE);
			CParameterGenerator PG(m_lpGLCanvas->m_lpMeshCenter->m_mshSurface, m_lpGLCanvas->m_lpMeshCenter->m_vecCrevasses);
			PG.m_fnMarkCrevasse();
			PG.m_fnExtractCrevasses();
			m_lpGLCanvas->m_lpMeshCenter->m_fnIsAligned();
			PG.m_fnUniformizeChart();
			++(m_lpGLCanvas->m_lpMeshCenter->m_psProcessStatus);
			Refresh();
		}
    }
}

void CMainFrame::OnSolveGlobalParameter(wxCommandEvent& event)
{
	CDialogs::CFacetControlDialog FCD(this);
	if (m_lpGLCanvas->m_lpMeshCenter->m_psProcessStatus >= PS_DIRECTION_UNIFORM)
	{
		if (FCD.ShowModal() == wxID_OK)
		{
			m_lpGLCanvas->m_lpMeshCenter->m_fnGenerateMoment();
			CParameterGenerator PG(m_lpGLCanvas->m_lpMeshCenter->m_mshSurface, m_lpGLCanvas->m_lpMeshCenter->m_vecCrevasses);
			PG.m_fnSolveCrevasseParaRange();
			m_lpGLCanvas->m_lpMeshCenter->m_nGlobalParameterSize = PG.m_fnInitGlobalPara(m_lpGLCanvas->m_lpMeshCenter->m_lpftGlobalPara, m_lpGLCanvas->m_lpMeshCenter->m_lpiGlobalIndex);
			PG.m_fnAssignCrevassePara();
			PG.m_fnRescaleCrevassePara(FCD.m_nNumFacets);
			PG.m_fnSolveHarmonicPara();
			m_lpGLCanvas->m_lpMeshCenter->m_fnRepresentLocalPara();
			m_lpGLCanvas->m_lpMeshCenter->m_fnRepresentGlobalPara();
			m_lpGLCanvas->m_lpMeshCenter->m_psProcessStatus = PS_CREVASSE_GLOBAL;
		}
	}
}

void CMainFrame::OnRescaleGlobalParameter(wxCommandEvent &event)
{
	long nFacets;
	wxTextEntryDialog TED(this, _("Number of Facets"), _("Rescale Global Parameter"));
	if (TED.ShowModal() == wxID_OK)
	{
		TED.GetValue().ToLong(&nFacets);
		CParameterGenerator PG(m_lpGLCanvas->m_lpMeshCenter->m_mshSurface, m_lpGLCanvas->m_lpMeshCenter->m_vecCrevasses);
		PG.m_fnRescaleGlobal(int(nFacets));
		m_lpGLCanvas->m_lpMeshCenter->m_fnRepresentLocalPara();
		m_lpGLCanvas->m_lpMeshCenter->m_fnRepresentGlobalPara();
		Refresh();
	}
}


void CMainFrame::OnSeamIntegerize(wxCommandEvent& event)
{
	CParameterAdjustor PA(m_lpGLCanvas->m_lpMeshCenter->m_mshSurface, m_lpGLCanvas->m_lpMeshCenter->m_vecCrevasses);
	wxMessageDialog MD(this, _("Are you sure to integrize global parameters?"), _("Command Confirmation"), wxOK | wxCANCEL);
	if (m_lpGLCanvas->m_lpMeshCenter->m_psProcessStatus >= PS_CREVASSE_GLOBAL)
	{
		if (MD.ShowModal() == wxID_OK)
		{
			PA.m_fnIntegerizeCrevasseSpan();
			PA.m_fnAssignCrevasseIntPara();
			PA.m_fnSolveHarmonicPara();
			PA.m_fnExtractCrevasseOffset();
			m_lpGLCanvas->m_lpMeshCenter->m_fnRepresentLocalPara();
			m_lpGLCanvas->m_lpMeshCenter->m_fnRepresentGlobalPara();
			++m_lpGLCanvas->m_lpMeshCenter->m_psProcessStatus;
			Refresh();
		}
	}

}

void CMainFrame::OnEvaluatePara(wxCommandEvent &event)
{
	std::ofstream ofs;
	wxString wxstrDefaultFile;
	wxstrDefaultFile = m_wxstrMeshName;
	wxstrDefaultFile.Replace(_(".off"), _("_evaluation.txt"));
	wxFileDialog SaveFileDialog(this, _("SaveEvaluation"), _(""), wxstrDefaultFile, "TxtFiles(*.txt)|*.txt", wxFD_SAVE | wxFD_OVERWRITE_PROMPT);
	double dblGrids, dblShapeDistortion, dblAreaDistortion;
	if (SaveFileDialog.ShowModal() == wxID_OK)
	{
		ofs.open(SaveFileDialog.GetPath().To8BitData());
		m_lpGLCanvas->m_lpMeshCenter->m_fnEvaluatePara(dblGrids, dblAreaDistortion, dblShapeDistortion);
		ofs << "Grids: " << dblGrids << '\n';
		ofs << "AreaDistortion: " << dblAreaDistortion << '\n';
		ofs << "ShapeDistortion: " << dblShapeDistortion << '\n';
		ofs.close();
	}
}

void CMainFrame::OnExtractSamplePoints(wxCommandEvent &event)
{
	wxMessageDialog MD(this, _("Are you sure to take samples according parameter?"), _("Command Confirmation"), wxOK | wxCANCEL);
    if(m_lpGLCanvas->m_lpMeshCenter->m_psProcessStatus >= PS_PARAMETER_INT)
    {
		if (MD.ShowModal() == wxID_OK)
		{
			m_lpGLCanvas->m_lpMeshCenter->m_fnStatusRetreat(PS_PARAMETER_INT);
			CSampler Sampler
			(
				m_lpGLCanvas->m_lpMeshCenter->m_mshSurface,
				m_lpGLCanvas->m_lpMeshCenter->m_vecCrevasses,
				m_lpGLCanvas->m_lpMeshCenter->m_mshRemeshedSurface,
				m_lpGLCanvas->m_lpMeshCenter->m_vecSamplePoints,
				m_lpGLCanvas->m_lpMeshCenter->m_vecSampleBlocks,
				m_lpGLCanvas->m_lpMeshCenter->m_vecFacetLinkages
			);
			Sampler.m_fnInitSampleBlocks();
			++(m_lpGLCanvas->m_lpMeshCenter->m_psProcessStatus);
			Refresh();
		}
    }
}

void CMainFrame::OnGenerateLinkages(wxCommandEvent &event)
{
	wxMessageDialog MD(this, _("Are you sure to generate linkages between samples?"), _("Command Confirmation"), wxOK | wxCANCEL);
    if(m_lpGLCanvas->m_lpMeshCenter->m_psProcessStatus >= PS_NEW_POINTS)
    {
		if (MD.ShowModal() == wxID_OK)
		{
			m_lpGLCanvas->m_lpMeshCenter->m_fnStatusRetreat(PS_NEW_POINTS);
			CSampler Sampler
				(
				m_lpGLCanvas->m_lpMeshCenter->m_mshSurface,
				m_lpGLCanvas->m_lpMeshCenter->m_vecCrevasses,
				m_lpGLCanvas->m_lpMeshCenter->m_mshRemeshedSurface,
				m_lpGLCanvas->m_lpMeshCenter->m_vecSamplePoints,
				m_lpGLCanvas->m_lpMeshCenter->m_vecSampleBlocks,
				m_lpGLCanvas->m_lpMeshCenter->m_vecFacetLinkages
				);
			Sampler.m_fnGenerateLinkage();
			++(m_lpGLCanvas->m_lpMeshCenter->m_psProcessStatus);
			Refresh();
		}
    }
}

void CMainFrame::OnConstructFacets(wxCommandEvent &event)
{
	wxMessageDialog MD(this, _("Are you sure to construct facets?"), _("Command Confirmation"), wxOK | wxCANCEL);
    if(m_lpGLCanvas->m_lpMeshCenter->m_psProcessStatus >= PS_NEW_LINKAGE)
    {
		m_lpGLCanvas->m_lpMeshCenter->m_fnStatusRetreat(PS_NEW_LINKAGE);
		if (MD.ShowModal() == wxID_OK)
		{
			CSampler Sampler
				(
				m_lpGLCanvas->m_lpMeshCenter->m_mshSurface,
				m_lpGLCanvas->m_lpMeshCenter->m_vecCrevasses,
				m_lpGLCanvas->m_lpMeshCenter->m_mshRemeshedSurface,
				m_lpGLCanvas->m_lpMeshCenter->m_vecSamplePoints,
				m_lpGLCanvas->m_lpMeshCenter->m_vecSampleBlocks,
				m_lpGLCanvas->m_lpMeshCenter->m_vecFacetLinkages
				);
			Sampler.m_fnGenerateFacets();
			++(m_lpGLCanvas->m_lpMeshCenter->m_psProcessStatus);
			Refresh();
		}
    }
}

void CMainFrame::OnUniformizeFacets(wxCommandEvent &event)
{
	wxMessageDialog MD(this, _("Are you sure to uniformize facets?"), _("Command Confirmation"), wxOK | wxCANCEL);
    if(m_lpGLCanvas->m_lpMeshCenter->m_psProcessStatus >= PS_NEW_SCATTER)
    {
		if (MD.ShowModal() == wxID_OK)
		{
			m_lpGLCanvas->m_lpMeshCenter->m_fnStatusRetreat(PS_NEW_SCATTER);
			CSampler Sampler
				(
				m_lpGLCanvas->m_lpMeshCenter->m_mshSurface,
				m_lpGLCanvas->m_lpMeshCenter->m_vecCrevasses,
				m_lpGLCanvas->m_lpMeshCenter->m_mshRemeshedSurface,
				m_lpGLCanvas->m_lpMeshCenter->m_vecSamplePoints,
				m_lpGLCanvas->m_lpMeshCenter->m_vecSampleBlocks,
				m_lpGLCanvas->m_lpMeshCenter->m_vecFacetLinkages
				);
			Sampler.m_fnCorrectFacetDirection();
			++(m_lpGLCanvas->m_lpMeshCenter->m_psProcessStatus);
			Refresh();
		}
    }
}

void CMainFrame::OnConstructNewMesh(wxCommandEvent &event)
{
	wxMessageDialog MD(this, _("Are you sure to construct new mesh?"), _("Command Confirmation"), wxOK | wxCANCEL);
    if(m_lpGLCanvas->m_lpMeshCenter->m_psProcessStatus >= PS_NEW_UNIFORM)
    {
		if (MD.ShowModal() == wxID_OK)
		{
			m_lpGLCanvas->m_lpMeshCenter->m_fnStatusRetreat(PS_NEW_UNIFORM);
			CSampler Sampler
				(
				m_lpGLCanvas->m_lpMeshCenter->m_mshSurface,
				m_lpGLCanvas->m_lpMeshCenter->m_vecCrevasses,
				m_lpGLCanvas->m_lpMeshCenter->m_mshRemeshedSurface,
				m_lpGLCanvas->m_lpMeshCenter->m_vecSamplePoints,
				m_lpGLCanvas->m_lpMeshCenter->m_vecSampleBlocks,
				m_lpGLCanvas->m_lpMeshCenter->m_vecFacetLinkages
				);
			Sampler.m_fnExtractNewMesh();
			m_lpGLCanvas->m_lpMeshCenter->m_mshRemeshedSurface.refresh_degree();
			m_lpGLCanvas->m_lpMeshCenter->m_fnCountNewMeshSize();
			m_lpGLCanvas->m_lpMeshCenter->m_fnMarkNewMeshBorder();
			m_lpGLCanvas->m_lpMeshCenter->m_fnMarkNewMeshCrevasse();
			CMeshOptimizer Optimizer(m_lpGLCanvas->m_lpMeshCenter->m_mshSurface, m_lpGLCanvas->m_lpMeshCenter->m_mshRemeshedSurface);
			Optimizer.m_fnMarkFreedom();
			++(m_lpGLCanvas->m_lpMeshCenter->m_psProcessStatus);
			Refresh();
		}
    }
}

void CMainFrame::OnOptimizeNewMesh(wxCommandEvent &event)
{
    CDialogs::CIterDlg dlgIter(this);
	if (m_lpGLCanvas->m_lpMeshCenter->m_psProcessStatus >= PS_NEW_COARSE)
	{
		if(dlgIter.ShowModal() == wxID_OK)
		{
			m_lpGLCanvas->m_lpMeshCenter->m_fnStatusRetreat(PS_NEW_COARSE);
			m_lpGLCanvas->m_Timer.Start(4);
			CAlgorithmExecutor *lpAlgorithmExecutor = new CAlgorithmExecutor(m_lpGLCanvas, PT_OPTIMIZE);
			lpAlgorithmExecutor->m_nIterateTimes = dlgIter.m_nIterTimes;
			lpAlgorithmExecutor->Run();
			Refresh();
		}
	}
}

void CMainFrame::OnSimplifyNewMesh(wxCommandEvent &event)
{
	CQuadMeshEvaluator QME(m_lpGLCanvas->m_lpMeshCenter->m_mshRemeshedSurface);
	QME.m_fnJumpSampling();
	Refresh();
}

void CMainFrame::OnEvaluateNewMesh(wxCommandEvent &event)
{
	CDialogs::CEstimateDialog ED(this);
	wxFileDialog FD(this, _("SaveMeshFile"), "", "", "TextFiles(*.txt)|*.txt", wxFD_SAVE | wxFD_OVERWRITE_PROMPT);
	std::ofstream ofs;
	CQuadMeshEvaluator QME(m_lpGLCanvas->m_lpMeshCenter->m_mshRemeshedSurface);
	long i;
	if (ED.ShowModal() == wxID_OK)
	{
		if (FD.ShowModal() == wxID_OK)
		{
			QME.m_fnCountAngle();
			QME.m_fnCountAspectRatio(ED.m_dblAspectRatioUnit, int(ED.m_nAspectRatioGroups));
			QME.m_fnCountFacetArea(ED.m_dblLogAreaUnit, int(ED.m_nAreaGroups));
			ofs.open(FD.GetPath().ToAscii());
			for (i = 0; i < 18; ++i)
			{
				ofs << QME.m_lpdblAngleDistribution[i] << '\t';
			}
			ofs << '\n';
			for (i = 0; i < ED.m_nAspectRatioGroups; ++i)
			{
				ofs << QME.m_lpAspectRatioDistribution[i] << '\t';
			}
			ofs << '\n';
			for (i = 0; i < ED.m_nAreaGroups; ++i)
			{
				ofs << QME.m_lpAreaDistribution[i] << '\t';
			}
			ofs << '\n';
			ofs << "Area Distortion: " << QME.m_ftAreaDistortion << '\n';
			ofs << "Shape Distortion: " << QME.m_fnCountAngularDistortion() << '\n';
			ofs << "Genus: " << (int(m_lpGLCanvas->m_lpMeshCenter->m_mshRemeshedSurface.size_of_vertices()) - int(m_lpGLCanvas->m_lpMeshCenter->m_mshRemeshedSurface.size_of_edges()) + int(m_lpGLCanvas->m_lpMeshCenter->m_mshRemeshedSurface.size_of_facets())) << '\n';
			ofs.close();
		}
	}
}

BEGIN_EVENT_TABLE(CMainFrame, wxFrame)
	EVT_MENU(wxID_OPEN, CMainFrame::OnOpen)
	EVT_TOOL(wxID_OPEN, CMainFrame::OnOpen)
	EVT_MENU(wxID_SAVE, CMainFrame::OnSave)
	EVT_TOOL(wxID_SAVE, CMainFrame::OnSave)
	EVT_MENU(idMenuSaveProcessStatus, CMainFrame::OnSaveProcessStatus)
	EVT_MENU(idMenuOpenProcessStatus, CMainFrame::OnOpenProcessStatus)
	EVT_MENU(idMenuExportMountain, CMainFrame::OnExportMountain)
	EVT_MENU(idMenuExportDirectionField, CMainFrame::OnExportDirectionField)
	EVT_MENU(idMenuExportParaDif, CMainFrame::OnExportParaDif)
	EVT_MENU(idMenuImportCones, CMainFrame::OnImportCones)
	EVT_MENU(idMenuExportSingularities, CMainFrame::OnExportSingularities)
	EVT_MENU(idMenuExportFixedDirection, CMainFrame::OnExportFixedDirection)
	EVT_MENU(idMenuExportHardFeature, CMainFrame::OnExportHardFeature)
	EVT_MENU(idMenuExportFixedMark, CMainFrame::OnExportFixedMark)
	EVT_MENU(idMenuImportFixedDirection, CMainFrame::OnImportFixedDirection)
	EVT_MENU(idMenuImportDirectionField, CMainFrame::OnImportDirectionField)
	EVT_MENU(idMenuImportBitmap, CMainFrame::OnImportBitmap)
	EVT_MENU(idMenuExportLIC, CMainFrame::OnExportLIC)
	EVT_MENU(idMenuExportCrevasse, CMainFrame::OnExportCrevasse)
	EVT_MENU(idMenuExportTexturePara, CMainFrame::OnExportTexturePara)
	EVT_MENU(idMenuExportAngleError, CMainFrame::OnExportAngleError)
	EVT_MENU(idMenuExportConfCoefColor, CMainFrame::OnExportConfCoefColor)
	EVT_MENU(idMenuOpenMovieMode, CMainFrame::OnOpenMovieMode)
	EVT_MENU(idMenuOpenRecordMode, CMainFrame::OnOpenRecordMode)
	EVT_MENU(idMenuExportNewMesh, CMainFrame::OnExportNewMesh)
	EVT_MENU(idMenuImportQuadMesh, CMainFrame::OnImportQuadMesh)
	EVT_MENU(idMenuExportQuadSingularities, CMainFrame::OnExportQuadSingularities)
	EVT_MENU(idMenuExportQuadCylinders, CMainFrame::OnExportQuadCylinders)
	EVT_MENU(idMenuExit, CMainFrame::OnExit)

	EVT_MENU(idMenuChangeViewPort, CMainFrame::OnChangeViewPort)
	EVT_MENU(idMenuAppendRotation, CMainFrame::OnAppendRotation)
	EVT_TOOL(idOrthodox, CMainFrame::OnOrthodox)
	EVT_TOOL(idShowEdge, CMainFrame::OnShowEdge)
	EVT_TOOL(idShowFrameField, CMainFrame::OnFrameField)
	EVT_TOOL(idShowLIC, CMainFrame::OnShowLIC)
	EVT_TOOL(idShowChart, CMainFrame::OnChartViewable)
	EVT_TOOL(idShowSingularity, CMainFrame::OnSingularityViewable)
	EVT_TOOL(idShowElectricFieldSum, CMainFrame::OnShowElectricFieldSum)
	EVT_TOOL(idShowElectricFieldDif, CMainFrame::OnShowElectricFieldDif)
	EVT_TOOL(idShowTrack, CMainFrame::OnTrackViewable)
	EVT_TOOL(idShowLocalU, CMainFrame::OnLocalU)
	EVT_TOOL(idShowLocalV, CMainFrame::OnLocalV)
	EVT_TOOL(idShowGaussCurvature, CMainFrame::OnShowGaussCurvature)
	EVT_TOOL(idShowConformalCoefficient, CMainFrame::OnShowConformalCoefficient)
	EVT_TOOL(idShowAngleError, CMainFrame::OnShowAngleError)
	EVT_TOOL(idShowParaCurl, CMainFrame::OnShowParaCurl)
	EVT_TOOL(idShowDensityU, CMainFrame::OnDensityU)
	EVT_TOOL(idShowDensityV, CMainFrame::OnDensityV)
	EVT_TOOL(idShowGlobalU, CMainFrame::OnGlobalU)
	EVT_TOOL(idShowGlobalV, CMainFrame::OnGlobalV)
	EVT_TOOL(idShowSamplePoints, CMainFrame::OnShowSamplePoints)
	EVT_TOOL(idShowLinkages, CMainFrame::OnShowLinkages)
	EVT_TOOL(idShowFacetNormal, CMainFrame::OnShowFacetDirections)
	EVT_TOOL(idShowNewMesh, CMainFrame::OnShowReconstructedMesh)

	EVT_MENU(idMenuNormalizeModal, CMainFrame::OnNormalizeModal)
	EVT_MENU(idMenuInitilizeChart, CMainFrame::OnInitilizeChart)
	EVT_MENU(idMenuEigenOptimal, CMainFrame::OnEigenOptimal)
	EVT_MENU(idMenuFlatten, CMainFrame::OnFlatten)
	EVT_MENU(idMenuFixCones, CMainFrame::OnFixCones)
	EVT_MENU(idMenuDeletePath, CMainFrame::OnDeletePath)
	EVT_MENU(idMenuComplexOptimize, CMainFrame::OnComplexOptimize)
	EVT_MENU(idMenuStopIteration, CMainFrame::OnStopIteration)

	EVT_MENU(idMenuDirectionAdj, CMainFrame::OnDirectionAdj)
	EVT_MENU(idMenuAdjustSingularities, CMainFrame::OnAdjustSingularities)
	EVT_MENU(idMenuSetRadius, CMainFrame::OnSetRadius)
	EVT_MENU(idMenuReportEnergy, CMainFrame::OnReportEnergy)
	EVT_MENU(idMenuReportSingularities, CMainFrame::OnReportSingularities)
	EVT_MENU(idMenuLinkPoles, CMainFrame::OnLinkPoles)
	EVT_MENU(idMenuMoveSingularity, CMainFrame::OnMoveSingularity)
	EVT_MENU(idMenuClearBarrier, CMainFrame::OnClearBarrier)
	EVT_MENU(idMenuGlobalAdjust, CMainFrame::OnGlobalAdjust)
	EVT_MENU(idMenuCountSingularities, CMainFrame::OnCountSingularities)

	EVT_MENU(idMenuMarkSharpFeature, CMainFrame::OnMarkSharpFeature)
	EVT_MENU(idMenuClearFix, CMainFrame::OnClearFix)

	EVT_MENU(idMenuDetectLoop, CMainFrame::OnDetectLoop)
	EVT_MENU(idMenuReleaseCurl, CMainFrame::OnReleaseCurl)
	EVT_MENU(idMenuCoarseDecurl, CMainFrame::OnCoarseDecurl)
	EVT_MENU(idMenuFineDecurl, CMainFrame::OnFineDecurl)
	EVT_MENU(idMenuConstrainedDecurl, CMainFrame::OnConstrainedDecurl)
	EVT_MENU(idMenuAdjustCurl, CMainFrame::OnAdjustCurl)

	EVT_MENU(idMenuTraceDirection, CMainFrame::OnTraceDirection)
	EVT_MENU(idMenuBiasSearch, CMainFrame::OnBiasSearch)
	EVT_MENU(idMenuAdjustPath, CMainFrame::OnAdjustPath)

	EVT_TOOL(idSelectPoint, CMainFrame::OnSelectPoint)
	EVT_TOOL(idSelectTrack, CMainFrame::OnSelectTrack)
	EVT_TOOL(idSelectCrevasse, CMainFrame::OnSelectCrevasse)
	EVT_TOOL(idSelectFacet, CMainFrame::OnSelectFacet)
	EVT_TOOL(idSelectEdge, CMainFrame::OnSelectEdge)

	EVT_MENU(idMenuDecomposeMesh, CMainFrame::OnDecomposeMesh)
    EVT_MENU(idMenuUniformizeChart, CMainFrame::OnUniformizeChart)
	EVT_MENU(idMenuSolveGlobalParameter, CMainFrame::OnSolveGlobalParameter)
	EVT_MENU(idMenuRescaleGlobalParameter, CMainFrame::OnRescaleGlobalParameter)
    EVT_MENU(idMenuSeamIntegerize, CMainFrame::OnSeamIntegerize)
	EVT_MENU(idMenuEvaluatePara, CMainFrame::OnEvaluatePara)

    EVT_MENU(idMenuExtractSamplePoints, CMainFrame::OnExtractSamplePoints)
    EVT_MENU(idMenuGenerateLinkages, CMainFrame::OnGenerateLinkages)
    EVT_MENU(idMenuConstructFacets, CMainFrame::OnConstructFacets)
    EVT_MENU(idMenuUniformizeFacets, CMainFrame::OnUniformizeFacets)
    EVT_MENU(idMenuConstructNewMesh, CMainFrame::OnConstructNewMesh)
    EVT_MENU(idMenuOptimizeNewMesh, CMainFrame::OnOptimizeNewMesh)
	EVT_MENU(idMenuSimplifyNewMesh, CMainFrame::OnSimplifyNewMesh)
	EVT_MENU(idMenuEvaluateNewMesh, CMainFrame::OnEvaluateNewMesh)
END_EVENT_TABLE()

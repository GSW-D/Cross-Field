#include <time.h>
#include "CAlgorithmExecutor.h"

CAlgorithmExecutor::CAlgorithmExecutor(CGLCanvas *lpGLCanvas, PROCESS_TYPE ptProcessType)
{
    m_lpGLCanvas = lpGLCanvas;
    m_ptProcessType = ptProcessType;
}

wxThread::ExitCode CAlgorithmExecutor::Entry()
{
	clock_t iStart, iEnd;
    m_lpGLCanvas->m_bProcessing = true;
	iStart = clock();
    switch(m_ptProcessType)
    {
	case PT_FLATTEN:
		m_fnFlattenCones();
		break;
	case PT_COMPLEX:
		m_fnComplexAdjust();
		break;
    case PT_SINGULARITY:
        m_fnAdjustSingularities();
        break;
	case PT_GLOBAL:
		m_fnGlobalAdjust();
		break;
	case PT_COUNT_SINGULARITIES:
		m_fnCountSingularities();
		break;
	case PT_COARSE_DECURL:
		m_fnCoarseDecurl();
		break;
	case PT_FINE_DECURL:
		m_fnFineDecurl();
		break;
	case PT_CONSTRAINED_DECURL:
		m_fnConstrainedDecurl();
		break;
    case PT_OPTIMIZE:
        m_fnOptimizeMesh();
        break;
    }
	iEnd = clock();
    m_lpGLCanvas->m_bProcessing = false;
	wxMessageBox(wxString::Format("Complete!\nTime consumption: %d milliseconds", int(iEnd - iStart)), _("Tip"), wxOK, m_lpGLCanvas);
    return (wxThread::ExitCode)0;
}

void CAlgorithmExecutor::m_fnFlattenCones()
{
	CDirectionGenerator DG(m_lpGLCanvas->m_lpMeshCenter->m_mshSurface);
	do
	{
		m_nIterateTimes = DG.m_fnFixSmallestErrorCone();
		DG.m_fnGenerateMetric();
	} while (m_nIterateTimes > 1 && m_lpGLCanvas->m_bProcessing);
	DG.m_fnGenerateDirFromDensity();
	m_lpGLCanvas->m_lpMeshCenter->m_fnDetectSeamType();
	m_lpGLCanvas->m_lpMeshCenter->m_fnMarkSingularity();
}

void CAlgorithmExecutor::m_fnComplexAdjust()
{
	int i;
	CDirectionAdjustor DA(m_lpGLCanvas->m_lpMeshCenter->m_mshSurface);
	CComplexOptimizer CO(m_lpGLCanvas->m_lpMeshCenter->m_mshSurface);
	CMeshBase::CCmplxSv CMS;
	CO.m_fnGenerateFreeMoment();
	CO.m_fnFillComplexMatrix(CMS);
	CO.m_fnInitVector(CMS);
	CMS.decompose();
	for (i = 0; i < m_nIterateTimes && m_lpGLCanvas->m_bProcessing; ++i)
	{
		CMS.solve();
		CO.m_fnNormalizeVector(CMS);
	}
	CO.m_fnFillDirection(CMS);
	m_lpGLCanvas->m_lpMeshCenter->m_fnDetectSeamType();
	DA.m_fnGenerateMoment();
	m_lpGLCanvas->m_lpMeshCenter->m_fnMarkSingularity();
}


void CAlgorithmExecutor::m_fnAdjustSingularities()
{
	CDirectionAdjustor DA(m_lpGLCanvas->m_lpMeshCenter->m_mshSurface);
	CMeshBase::CLDLTSv LDLTS0, LDLTS1;
	wxString wxstrFileName;
	std::ofstream ofs;
	std::vector<CMeshBase::FT> vecEnergy;
	std::vector<CMeshBase::FT>::iterator iterEnergy;
	int i, nAdjust;
    CLICPainter LICPainter(m_lpGLCanvas);
	i = 0;
	if (m_lpGLCanvas->m_bMovieMode)
	{
        LICPainter.m_fnVideoTestSize(2.0, 600, 600);
	}
	if (m_lpGLCanvas->m_bRecordMode)
    {
        ofs.open(m_lpGLCanvas->m_wxstrRecordFile.ToAscii());
    }
	DA.m_fnInitDirectionMatrix(LDLTS0);
	DA.m_fnInitDensityMatrix(LDLTS1);
	LDLTS0.decompose();
	LDLTS1.decompose();
	do
	{
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
		DA.m_fnGenerateForces();
		if (m_lpGLCanvas->m_bMovieMode)
		{
			if (m_lpGLCanvas->m_bDiscrete)
			{
				LICPainter.m_fnGenerateViewCoordinates();
				LICPainter.m_fnTestSize(1.0);
				LICPainter.m_fnInitBuffer();
				LICPainter.m_fnSetLocations(false);
				LICPainter.m_fnDrawDiscrete();
				LICPainter.m_fnAppendSingularities(1.0);
				LICPainter.m_fnMarkForceToDraw();
				LICPainter.m_fnDrawForce();
			}
			else
			{
				LICPainter.m_fnInitBuffer();
				LICPainter.m_fnGenerateViewCoordinates();
				m_lpGLCanvas->m_lpMeshCenter->m_fnFillVertexCrossField();
				LICPainter.m_fnSetLocations(true);
				LICPainter.m_fnDetectSeam();
				LICPainter.m_fnAssignOriginalColor(m_lpGLCanvas->m_bmpOriginal);
				LICPainter.m_fnCrossTrace(200);
				LICPainter.m_fnEnhanceBitmap(false);
				LICPainter.m_fnAppendSingularities(1.5);
			}
			LICPainter.m_fnWriteBitmap(m_lpGLCanvas->m_bmpLIC, m_lpGLCanvas->m_iBmpX, m_lpGLCanvas->m_iBmpY);
			wxstrFileName = wxString::Format("\\DirectionField%d%d%d.bmp", (i / 100), ((i / 10) % 10), (i % 10));
			m_lpGLCanvas->m_bmpLIC.SaveFile(m_lpGLCanvas->m_wxstrMovieDir + wxstrFileName, wxBITMAP_TYPE_BMP);
		}
		if (m_lpGLCanvas->m_bRecordMode)
        {
            ofs << i << '\t' << DA.m_fnNumberOfSingularities() << '\t' << DA.m_fnEnergy(m_lpGLCanvas->m_lpMeshCenter->m_ftChargeRadius) << "\t\n";
        }
		++i;
		switch (m_iAdjustType)
		{
		case 0:
			nAdjust = DA.m_fnAdjustSingularitiesTheta();
			break;
		case 1:
			nAdjust = DA.m_fnAdjustSingularitiesPhi();
			break;
		}
	} while (nAdjust != 0 && m_lpGLCanvas->m_bProcessing);
	m_lpGLCanvas->m_lpMeshCenter->m_fnClearBarrier();
	--i;
	LDLTS1.clear();
	LDLTS0.clear();
	m_lpGLCanvas->m_lpMeshCenter->m_fnRepresentElectricField();
	if (m_lpGLCanvas->m_bMovieMode)
	{
	    m_lpGLCanvas->m_wxstrMovieDir.clear();
		m_lpGLCanvas->m_bMovieMode = false;
	}
	if (m_lpGLCanvas->m_bRecordMode)
    {
        ofs.close();
        m_lpGLCanvas->m_wxstrRecordFile.clear();
        m_lpGLCanvas->m_bRecordMode = false;
    }
}


void CAlgorithmExecutor::m_fnGlobalAdjust()
{
	CDirectionAdjustor DA(m_lpGLCanvas->m_lpMeshCenter->m_mshSurface);
	CMeshBase::CPolyhedron::Vertex_handle vhMaximum, vhMinimum;
	CMeshBase::CLDLTSv LDLTS0, LDLTS1;
	CMeshBase::FT ftEnergy, ftMinEnergy;
	double dblBottom, dblTop;
	int iTestTimes, nAdjust, iFile, iFolder, iGlobalIndex, iLocalIndex;
	wxString wxstrFileName, wxstrFolderName;
	std::ofstream ofs, ofs_r;
	std::vector<CMeshBase::FT>::iterator iterEnergy;
    CLICPainter LICPainter(m_lpGLCanvas);
	iTestTimes = 0;
	iFile = 0;
	iFolder = 0;
	iGlobalIndex = 0;
	iLocalIndex = 0;
	if (m_lpGLCanvas->m_bMovieMode)
	{
        LICPainter.m_fnVideoTestSize(2.0, 640, 720);
        wxstrFolderName = wxString::Format("\\%d%d", iFolder / 10, iFolder % 10);
        wxFileName::Mkdir(m_lpGLCanvas->m_wxstrMovieDir + wxstrFolderName);
        ++iFolder;
	}
	if (m_lpGLCanvas->m_bRecordMode)
    {
        ofs_r.open(m_lpGLCanvas->m_wxstrRecordFile.ToAscii());
    }
	DA.m_fnInitDirectionMatrix(LDLTS0);
	DA.m_fnInitDensityMatrix(LDLTS1);
	LDLTS0.decompose();
	LDLTS1.decompose();
	do
	{
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
		DA.m_fnGenerateForces();
		if (m_lpGLCanvas->m_bMovieMode)
		{
            LICPainter.m_fnInitBuffer();
            LICPainter.m_fnGenerateViewCoordinates();
            m_lpGLCanvas->m_lpMeshCenter->m_fnFillVertexCrossField();
            LICPainter.m_fnSetLocations(true);
            LICPainter.m_fnDetectSeam();
			LICPainter.m_fnAssignOriginalColor(m_lpGLCanvas->m_bmpOriginal);
            LICPainter.m_fnCrossTrace(200);
            LICPainter.m_fnEnhanceBitmap(false);
			LICPainter.m_fnAppendSingularities(2.0);
			LICPainter.m_fnWriteBitmap(m_lpGLCanvas->m_bmpLIC, m_lpGLCanvas->m_iBmpX, m_lpGLCanvas->m_iBmpY);
			wxstrFileName = wxString::Format("\\DirectionField%d%d%d.bmp", (iFile / 100), ((iFile / 10) % 10), (iFile % 10));
			m_lpGLCanvas->m_bmpLIC.SaveFile(m_lpGLCanvas->m_wxstrMovieDir + wxstrFolderName + wxstrFileName, wxBITMAP_TYPE_BMP);
			++iFile;
		}
		if (m_lpGLCanvas->m_bRecordMode)
        {
            ofs_r << iGlobalIndex << '\t' << iLocalIndex << '\t' << DA.m_fnNumberOfSingularities() << '\t' << DA.m_fnEnergy(m_lpGLCanvas->m_lpMeshCenter->m_ftChargeRadius) << "\t\n";
            ++iLocalIndex;
        }
		switch (m_iAdjustType)
		{
		case 0:
			nAdjust = DA.m_fnAdjustSingularitiesTheta();
			break;
		case 1:
			nAdjust = DA.m_fnAdjustSingularitiesPhi();
			break;
		}
	} while (nAdjust != 0 && m_lpGLCanvas->m_bProcessing);
	m_lpGLCanvas->m_lpMeshCenter->m_fnClearBarrier();
	m_lpGLCanvas->m_lpMeshCenter->m_fnMarkSingularity();
	DA.m_fnFillPoisson(LDLTS1);
	LDLTS1.solve();
	DA.m_fnAssignDensity(LDLTS1);
	ftMinEnergy = DA.m_fnEnergy(m_lpGLCanvas->m_lpMeshCenter->m_ftChargeRadius);
	m_lpGLCanvas->m_lpMeshCenter->m_fnBackupField();
	do
	{
	    if (m_lpGLCanvas->m_bMovieMode)
        {
            wxstrFolderName = wxString::Format("\\%d%d", iFolder / 10, iFolder % 10);
            wxFileName::Mkdir(m_lpGLCanvas->m_wxstrMovieDir + wxstrFolderName);
            ++iFolder;
            DA.m_fnModifyDensity(m_lpGLCanvas->m_lpMeshCenter->m_ftChargeRadius, dblBottom, dblTop);
            wxstrFileName = _("\\conf_coef.obj");
            ofs.open((m_lpGLCanvas->m_wxstrMovieDir + wxstrFolderName + wxstrFileName).ToAscii());
            m_lpGLCanvas->m_lpMeshCenter->m_fnExportMountain(ofs, dblBottom, dblTop, 0.0);
            ofs.close();
            wxstrFileName = _("\\previous_field.txt");
            ofs.open((m_lpGLCanvas->m_wxstrMovieDir + wxstrFolderName + wxstrFileName).ToAscii());
            m_lpGLCanvas->m_lpMeshCenter->m_fnExportCrossField(ofs);
            ofs.close();
            wxstrFileName = _("\\previous_singularities.txt");
            ofs.open((m_lpGLCanvas->m_wxstrMovieDir + wxstrFolderName + wxstrFileName).ToAscii());
            m_lpGLCanvas->m_lpMeshCenter->m_fnExportSingularities(ofs);
            ofs.close();
        }
		DA.m_fnFindPoles(m_lpGLCanvas->m_lpMeshCenter->m_ftChargeRadius, vhMinimum, vhMaximum);
		if (vhMaximum != NULL && vhMinimum != NULL)
		{
			DA.m_fnMoveSingularity(vhMinimum, vhMaximum);
			if (m_lpGLCanvas->m_bMovieMode)
            {
                wxstrFileName = _("\\crevasse.obj");
                ofs.open((m_lpGLCanvas->m_wxstrMovieDir + wxstrFolderName + wxstrFileName).ToAscii());
                m_lpGLCanvas->m_lpMeshCenter->m_fnExportCrevasse(ofs);
                ofs.close();
                wxstrFolderName = wxString::Format("\\%d%d", iFolder / 10, iFolder % 10);
                wxFileName::Mkdir(m_lpGLCanvas->m_wxstrMovieDir + wxstrFolderName);
                ++iFolder;
                iFile = 0;
            }
            if (m_lpGLCanvas->m_bRecordMode)
            {
                ++iGlobalIndex;
                iLocalIndex = 0;
            }
			m_lpGLCanvas->m_lpMeshCenter->m_fnClearFreeEdge();
		}
		do
		{
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
			DA.m_fnGenerateForces();
			if (m_lpGLCanvas->m_bMovieMode)
			{
                LICPainter.m_fnInitBuffer();
                LICPainter.m_fnGenerateViewCoordinates();
                m_lpGLCanvas->m_lpMeshCenter->m_fnFillVertexCrossField();
                LICPainter.m_fnSetLocations(true);
                LICPainter.m_fnDetectSeam();
                LICPainter.m_fnAssignOriginalColor(m_lpGLCanvas->m_bmpOriginal);
                LICPainter.m_fnCrossTrace(200);
                LICPainter.m_fnEnhanceBitmap(false);
                LICPainter.m_fnAppendSingularities(2.0);
                LICPainter.m_fnWriteBitmap(m_lpGLCanvas->m_bmpLIC, m_lpGLCanvas->m_iBmpX, m_lpGLCanvas->m_iBmpY);
                wxstrFileName = wxString::Format("\\DirectionField%d%d%d.bmp", (iFile / 100), ((iFile / 10) % 10), (iFile % 10));
                m_lpGLCanvas->m_bmpLIC.SaveFile(m_lpGLCanvas->m_wxstrMovieDir + wxstrFolderName + wxstrFileName, wxBITMAP_TYPE_BMP);
				ofs.close();
				++iFile;
			}
            if (m_lpGLCanvas->m_bRecordMode)
            {
                ofs_r << iGlobalIndex << '\t' << iLocalIndex << '\t' << DA.m_fnNumberOfSingularities() << '\t' << DA.m_fnEnergy(m_lpGLCanvas->m_lpMeshCenter->m_ftChargeRadius) << "\t\n";
                ++iLocalIndex;
            }
			switch (m_iAdjustType)
			{
			case 0:
				nAdjust = DA.m_fnAdjustSingularitiesTheta();
				break;
			case 1:
				nAdjust = DA.m_fnAdjustSingularitiesPhi();
				break;
			}
		} while (nAdjust != 0 && m_lpGLCanvas->m_bProcessing);
		m_lpGLCanvas->m_lpMeshCenter->m_fnClearBarrier();
		ftEnergy = DA.m_fnEnergy(m_lpGLCanvas->m_lpMeshCenter->m_ftChargeRadius);
		if (ftEnergy < ftMinEnergy)
		{
			iTestTimes = 0;
			ftMinEnergy = ftEnergy;
			m_lpGLCanvas->m_lpMeshCenter->m_fnBackupField();
		}
		else
		{
			++iTestTimes;
			if (iTestTimes > m_nIterateTimes)
			{
				m_lpGLCanvas->m_bProcessing = false;
			}
		}
	} while (m_lpGLCanvas->m_bProcessing);
	m_lpGLCanvas->m_lpMeshCenter->m_fnRestoreField();
	m_lpGLCanvas->m_lpMeshCenter->m_fnGenerateMoment();
	m_lpGLCanvas->m_lpMeshCenter->m_fnMarkSingularity();
	DA.m_fnFillPoisson(LDLTS1);
	LDLTS1.solve();
	DA.m_fnAssignDensity(LDLTS1);
	DA.m_fnGenerateNablaField();
	m_lpGLCanvas->m_lpMeshCenter->m_fnRepresentElectricField();
	if (m_lpGLCanvas->m_bMovieMode)
	{
        m_lpGLCanvas->m_wxstrMovieDir.clear();
		m_lpGLCanvas->m_bMovieMode = false;
	}
	if (m_lpGLCanvas->m_bRecordMode)
    {
        ofs_r.close();
        m_lpGLCanvas->m_wxstrRecordFile.clear();
        m_lpGLCanvas->m_bRecordMode = false;
    }
}

void CAlgorithmExecutor::m_fnCountSingularities()
{
	int i;
	std::ofstream ofs;
	std::ifstream ifs;
	ofs.open(m_wxstrOutputData.ToAscii());
	CDirectionGenerator DG(m_lpGLCanvas->m_lpMeshCenter->m_mshSurface);
	CDirectionAdjustor DA(m_lpGLCanvas->m_lpMeshCenter->m_mshSurface);
	CMeshBase::CPolyhedron::Vertex_handle vhMaximum, vhMinimum;
	CMeshBase::CLDLTSv LDLTS0, LDLTS1;
	CMeshBase::FT ftEnergy, ftMinEnergy;
	bool bStepContinue;
	int iTestTimes, nAdjust;
	iTestTimes = 0;
	DA.m_fnInitDirectionMatrix(LDLTS0);
	DA.m_fnInitDensityMatrix(LDLTS1);
	LDLTS0.decompose();
	LDLTS1.decompose();
	for (i = 0; i <= m_nTestTimes && m_lpGLCanvas->m_bProcessing; ++i)
	{
		ifs.open(m_wxstrInputField.ToAscii());
		m_lpGLCanvas->m_lpMeshCenter->m_fnImportDirectionField(ifs, false, false);
		ifs.close();
		m_nIterateTimes = 3;
		m_lpGLCanvas->m_lpMeshCenter->m_ftChargeRadius = m_lpGLCanvas->m_lpMeshCenter->m_ftGeometricDimension * exp(CMeshBase::FT(m_dblStartLevel - m_dblDecrement * i) * log(10.0));
		m_lpGLCanvas->m_lpMeshCenter->m_fnClearFreeEdge();
		m_lpGLCanvas->m_lpMeshCenter->m_fnDetectSeamType();
		m_lpGLCanvas->m_lpMeshCenter->m_fnMarkSingularity();
		do
		{
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
			DA.m_fnGenerateForces();

			switch (m_iAdjustType)
			{
			case 0:
				nAdjust = DA.m_fnAdjustSingularitiesTheta();
				break;
			case 1:
				nAdjust = DA.m_fnAdjustSingularitiesPhi();
				break;
			}
		} while (nAdjust != 0 && m_lpGLCanvas->m_bProcessing);
		m_lpGLCanvas->m_lpMeshCenter->m_fnClearBarrier();
		m_lpGLCanvas->m_lpMeshCenter->m_fnMarkSingularity();
		DA.m_fnFillPoisson(LDLTS1);
		LDLTS1.solve();
		DA.m_fnAssignDensity(LDLTS1);
		ftMinEnergy = DA.m_fnEnergy(m_lpGLCanvas->m_lpMeshCenter->m_ftChargeRadius);
		m_lpGLCanvas->m_lpMeshCenter->m_fnBackupField();
		bStepContinue = true;
		do
		{
			DA.m_fnFindPoles(m_lpGLCanvas->m_lpMeshCenter->m_ftChargeRadius, vhMinimum, vhMaximum);
			if (vhMaximum != NULL && vhMinimum != NULL)
			{
				DA.m_fnMoveSingularity(vhMinimum, vhMaximum);
			}
			m_lpGLCanvas->m_lpMeshCenter->m_fnClearFreeEdge();
			do
			{
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
				DA.m_fnGenerateForces();

				switch (m_iAdjustType)
				{
				case 0:
					nAdjust = DA.m_fnAdjustSingularitiesTheta();
					break;
				case 1:
					nAdjust = DA.m_fnAdjustSingularitiesPhi();
					break;
				}
			} while (nAdjust != 0 && m_lpGLCanvas->m_bProcessing);
			m_lpGLCanvas->m_lpMeshCenter->m_fnClearBarrier();
			ftEnergy = DA.m_fnEnergy(m_lpGLCanvas->m_lpMeshCenter->m_ftChargeRadius);
			if (ftEnergy < ftMinEnergy)
			{
				iTestTimes = 0;
				ftMinEnergy = ftEnergy;
				m_lpGLCanvas->m_lpMeshCenter->m_fnBackupField();
			}
			else
			{
				++iTestTimes;
				if (iTestTimes > m_nIterateTimes)
				{
					bStepContinue = false;
				}
			}
		} while (bStepContinue && m_lpGLCanvas->m_bProcessing);
		m_lpGLCanvas->m_lpMeshCenter->m_fnRestoreField();
		m_lpGLCanvas->m_lpMeshCenter->m_fnGenerateMoment();
		m_lpGLCanvas->m_lpMeshCenter->m_fnMarkSingularity();
		DA.m_fnFillPoisson(LDLTS1);
		LDLTS1.solve();
		DA.m_fnAssignDensity(LDLTS1);
		DA.m_fnGenerateNablaField();
		m_lpGLCanvas->m_lpMeshCenter->m_fnRepresentElectricField();
		ofs << m_dblStartLevel - m_dblDecrement * i << '\t' << DA.m_fnNumberOfSingularities() << '\t' << DA.m_fnEnergy(m_lpGLCanvas->m_lpMeshCenter->m_ftChargeRadius) << '\n';
	}
	ofs.close();
	m_lpGLCanvas->m_bProcessing = false;
}

void CAlgorithmExecutor::m_fnCoarseDecurl()
{
	CDirectionAdjustor DA(m_lpGLCanvas->m_lpMeshCenter->m_mshSurface);
	CMeshBase::CLDLTSv LDLTS;
	CMeshBase::FT ftEnergy, ftMinEnergy, ftMaxElectricStrength, ftMaxCurlStrength;
	int iTestTimes, iFile, iIndex;
	bool bMaxCurlInitialized;
	bMaxCurlInitialized = false;
	std::ofstream ofs, ofs_r;
	wxString wxstrFileName;
    CLICPainter LICPainter(m_lpGLCanvas);
	iTestTimes = 0;
	DA.m_fnInitDirectionMatrix(LDLTS);
	LDLTS.decompose();
	DA.m_fnGenerateMoment();
	DA.m_fnFillMoment(LDLTS);
	LDLTS.solve();
	DA.m_fnAdjustDirection(LDLTS);
	DA.m_fnGenerateMoment();
	DA.m_fnGenerateNablaField();
	ftMinEnergy = ftEnergy = DA.m_fnCurlEnergy();
	m_lpGLCanvas->m_lpMeshCenter->m_fnBackupField();
	iFile = 0;
	iIndex = 0;
	if (m_lpGLCanvas->m_bMovieMode)
    {
        LICPainter.m_fnVideoTestSize(2.0, 320, 720);
    }
    if (m_lpGLCanvas->m_bRecordMode)
    {
        ofs_r.open(m_lpGLCanvas->m_wxstrRecordFile.ToAscii());
    }
	do
	{
		if (m_lpGLCanvas->m_bMovieMode)
        {
            LICPainter.m_fnInitBuffer();
            LICPainter.m_fnGenerateViewCoordinates();
            m_lpGLCanvas->m_lpMeshCenter->m_fnFillVertexCrossField();
            LICPainter.m_fnSetLocations(true);
            LICPainter.m_fnDetectSeam();
			LICPainter.m_fnAssignOriginalColor(m_lpGLCanvas->m_bmpOriginal);
            LICPainter.m_fnCrossTrace(200);
            LICPainter.m_fnEnhanceBitmap(false);
			LICPainter.m_fnAppendSingularities(2.0);
			LICPainter.m_fnWriteBitmap(m_lpGLCanvas->m_bmpLIC, m_lpGLCanvas->m_iBmpX, m_lpGLCanvas->m_iBmpY);
			wxstrFileName = wxString::Format("\\DirectionField%d%d%d.bmp", (iFile / 100), ((iFile / 10) % 10), (iFile % 10));
			m_lpGLCanvas->m_bmpLIC.SaveFile(m_lpGLCanvas->m_wxstrMovieDir + wxstrFileName, wxBITMAP_TYPE_BMP);


			m_lpGLCanvas->m_lpMeshCenter->m_fnFillVertexNablaTheta();
			ftMaxElectricStrength = m_lpGLCanvas->m_lpMeshCenter->m_ftMaxElectricField;
            LICPainter.m_fnInitBuffer();
            LICPainter.m_fnGenerateViewCoordinates();
            m_lpGLCanvas->m_lpMeshCenter->m_fnFillVertexElectricField();
            if (m_lpGLCanvas->m_lpMeshCenter->m_ftMaxElectricField > ftMaxElectricStrength)
            {
                ftMaxElectricStrength = m_lpGLCanvas->m_lpMeshCenter->m_ftMaxElectricField;
            }
            LICPainter.m_fnSetLocations(false);
            LICPainter.m_fnDetectSeam();
			LICPainter.m_fnAssignOriginalColor(m_lpGLCanvas->m_bmpOriginal);
            LICPainter.m_fnLinearTrace(200);
            LICPainter.m_fnEnhanceBitmap(true, ftMaxElectricStrength);
			LICPainter.m_fnAppendSingularities(2.0);
			LICPainter.m_fnWriteBitmap(m_lpGLCanvas->m_bmpLIC, m_lpGLCanvas->m_iBmpX, m_lpGLCanvas->m_iBmpY);
			wxstrFileName = wxString::Format("\\ElectricField%d%d%d.bmp", (iFile / 100), ((iFile / 10) % 10), (iFile % 10));
			m_lpGLCanvas->m_bmpLIC.SaveFile(m_lpGLCanvas->m_wxstrMovieDir + wxstrFileName, wxBITMAP_TYPE_BMP);

            LICPainter.m_fnInitBuffer();
            LICPainter.m_fnGenerateViewCoordinates();
            m_lpGLCanvas->m_lpMeshCenter->m_fnFillVertexNablaTheta();
            LICPainter.m_fnSetLocations(false);
            LICPainter.m_fnDetectSeam();
			LICPainter.m_fnAssignOriginalColor(m_lpGLCanvas->m_bmpOriginal);
            LICPainter.m_fnLinearTrace(200);
            LICPainter.m_fnEnhanceBitmap(true, ftMaxElectricStrength);
			LICPainter.m_fnAppendSingularities(2.0);
			LICPainter.m_fnWriteBitmap(m_lpGLCanvas->m_bmpLIC, m_lpGLCanvas->m_iBmpX, m_lpGLCanvas->m_iBmpY);
			wxstrFileName = wxString::Format("\\NablaTheta%d%d%d.bmp", (iFile / 100), ((iFile / 10) % 10), (iFile % 10));
			m_lpGLCanvas->m_bmpLIC.SaveFile(m_lpGLCanvas->m_wxstrMovieDir + wxstrFileName, wxBITMAP_TYPE_BMP);

            LICPainter.m_fnInitBuffer();
            LICPainter.m_fnGenerateViewCoordinates();
            m_lpGLCanvas->m_lpMeshCenter->m_fnFillVertexCurlField();
            if (!bMaxCurlInitialized)
            {
                ftMaxCurlStrength = m_lpGLCanvas->m_lpMeshCenter->m_ftMaxCurlField;
                bMaxCurlInitialized = true;
            }
            LICPainter.m_fnSetLocations(false);
            LICPainter.m_fnDetectSeam();
			LICPainter.m_fnAssignOriginalColor(m_lpGLCanvas->m_bmpOriginal);
            LICPainter.m_fnLinearTrace(200);
            LICPainter.m_fnEnhanceBitmap(true, ftMaxCurlStrength);
			LICPainter.m_fnAppendSingularities(2.0);
			LICPainter.m_fnWriteBitmap(m_lpGLCanvas->m_bmpLIC, m_lpGLCanvas->m_iBmpX, m_lpGLCanvas->m_iBmpY);
			wxstrFileName = wxString::Format("\\CurlField%d%d%d.bmp", (iFile / 100), ((iFile / 10) % 10), (iFile % 10));
			m_lpGLCanvas->m_bmpLIC.SaveFile(m_lpGLCanvas->m_wxstrMovieDir + wxstrFileName, wxBITMAP_TYPE_BMP);

			wxstrFileName = wxString::Format("\\CrossField%d%d%d.txt", (iFile / 100), ((iFile / 10) % 10), (iFile % 10));
            ofs.open((m_lpGLCanvas->m_wxstrMovieDir + wxstrFileName).ToAscii());
            m_lpGLCanvas->m_lpMeshCenter->m_fnExportDirectionField(ofs, false);
            ofs.close();
			wxstrFileName = wxString::Format("\\Singularities%d%d%d.obj", (iFile / 100), ((iFile / 10) % 10), (iFile % 10));
            ofs.open((m_lpGLCanvas->m_wxstrMovieDir + wxstrFileName).ToAscii());
            m_lpGLCanvas->m_lpMeshCenter->m_fnExportSingularities(ofs);
            ofs.close();
        }
		DA.m_fnMarkCurlPath();
		if (m_lpGLCanvas->m_bMovieMode)
        {
			wxstrFileName = wxString::Format("\\Crevasse%d%d%d.obj", (iFile / 100), ((iFile / 10) % 10), (iFile % 10));
            ofs.open((m_lpGLCanvas->m_wxstrMovieDir + wxstrFileName).ToAscii());
            m_lpGLCanvas->m_lpMeshCenter->m_fnExportCrevasse(ofs);
            ofs.close();
            ++iFile;
        }
		if (m_lpGLCanvas->m_bRecordMode)
        {
            ofs_r << iIndex << '\t' << ftEnergy << "\t\n";
            ++iIndex;
        }
		DA.m_fnReleaseCurl();
		DA.m_fnGenerateMoment();
		DA.m_fnFillMoment(LDLTS);
		LDLTS.solve();
		DA.m_fnAdjustDirection(LDLTS);
		DA.m_fnGenerateMoment();
		DA.m_fnGenerateNablaField();
		ftEnergy = DA.m_fnCurlEnergy();
		if (ftEnergy < ftMinEnergy)
		{
			iTestTimes = 0;
			ftMinEnergy = ftEnergy;
			m_lpGLCanvas->m_lpMeshCenter->m_fnBackupField();
		}
		else
		{
			++iTestTimes;
			if (iTestTimes > m_nIterateTimes)
			{
				m_lpGLCanvas->m_bProcessing = false;
			}
		}
	} while (m_lpGLCanvas->m_bProcessing && m_lpGLCanvas->m_bProcessing);
	LDLTS.clear();
	m_lpGLCanvas->m_lpMeshCenter->m_fnRestoreField();
	m_lpGLCanvas->m_lpMeshCenter->m_fnGenerateMoment();
	DA.m_fnGenerateNablaField();
	m_lpGLCanvas->m_lpMeshCenter->m_fnRepresentElectricField();
	if (m_lpGLCanvas->m_bMovieMode)
    {
        m_lpGLCanvas->m_bMovieMode = false;
        m_lpGLCanvas->m_wxstrMovieDir.clear();
    }
    if (m_lpGLCanvas->m_bRecordMode)
    {
        ofs_r.clear();
        m_lpGLCanvas->m_bRecordMode = false;
        m_lpGLCanvas->m_wxstrRecordFile.clear();
    }
}


void CAlgorithmExecutor::m_fnFineDecurl()
{
	CDirectionAdjustor DA(m_lpGLCanvas->m_lpMeshCenter->m_mshSurface);
	CMeshBase::CLDLTSv LDLTS0, LDLTS1;
	wxString wxstrFileName;
	std::ofstream ofs;
	CMeshBase::FT ftEnergy, ftMinEnergy, ftMaxElectricStrength, ftMaxCurlStrength;
	int iTestTimes, iFile;
	ftMaxCurlStrength = 0.5165;
    CLICPainter LICPainter(m_lpGLCanvas);
	iTestTimes = 0;
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
	ftMinEnergy = DA.m_fnCurlEnergy();
	m_lpGLCanvas->m_lpMeshCenter->m_fnBackupField();
	iFile = 0;
	if (m_lpGLCanvas->m_bMovieMode)
    {
        LICPainter.m_fnVideoTestSize(2.0, 600, 600);
    }
    if (m_lpGLCanvas->m_bRecordMode)
    {
        ofs.open(m_lpGLCanvas->m_wxstrRecordFile.ToAscii());
    }
	do
	{
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
		ftEnergy = DA.m_fnCurlEnergy();
		if (ftEnergy < ftMinEnergy)
		{
			iTestTimes = 0;
			ftMinEnergy = ftEnergy;
			m_lpGLCanvas->m_lpMeshCenter->m_fnBackupField();
		}
		else
		{
			++iTestTimes;
			if (iTestTimes > m_nIterateTimes)
			{
				m_lpGLCanvas->m_bProcessing = false;
			}
		}
		if (m_lpGLCanvas->m_bMovieMode)
		{
            LICPainter.m_fnInitBuffer();
            LICPainter.m_fnGenerateViewCoordinates();
            m_lpGLCanvas->m_lpMeshCenter->m_fnFillVertexCrossField();
            LICPainter.m_fnSetLocations(true);
            LICPainter.m_fnDetectSeam();
			LICPainter.m_fnAssignOriginalColor(m_lpGLCanvas->m_bmpOriginal);
            LICPainter.m_fnCrossTrace(200);
            LICPainter.m_fnEnhanceBitmap(false);
			LICPainter.m_fnAppendSingularities(2.0);
			LICPainter.m_fnWriteBitmap(m_lpGLCanvas->m_bmpLIC, m_lpGLCanvas->m_iBmpX, m_lpGLCanvas->m_iBmpY);
			wxstrFileName = wxString::Format("\\DirectionField%d%d%d.bmp", (iFile / 100), ((iFile / 10) % 10), (iFile % 10));
			m_lpGLCanvas->m_bmpLIC.SaveFile(m_lpGLCanvas->m_wxstrMovieDir + wxstrFileName, wxBITMAP_TYPE_BMP);


			m_lpGLCanvas->m_lpMeshCenter->m_fnFillVertexNablaTheta();
			ftMaxElectricStrength = m_lpGLCanvas->m_lpMeshCenter->m_ftMaxElectricField;
            LICPainter.m_fnInitBuffer();
            LICPainter.m_fnGenerateViewCoordinates();
            m_lpGLCanvas->m_lpMeshCenter->m_fnFillVertexElectricField();
            if (m_lpGLCanvas->m_lpMeshCenter->m_ftMaxElectricField > ftMaxElectricStrength)
            {
                ftMaxElectricStrength = m_lpGLCanvas->m_lpMeshCenter->m_ftMaxElectricField;
            }
            LICPainter.m_fnSetLocations(false);
            LICPainter.m_fnDetectSeam();
			LICPainter.m_fnAssignOriginalColor(m_lpGLCanvas->m_bmpOriginal);
            LICPainter.m_fnLinearTrace(200);
            LICPainter.m_fnEnhanceBitmap(true, ftMaxElectricStrength);
			LICPainter.m_fnAppendSingularities(2.0);
			LICPainter.m_fnWriteBitmap(m_lpGLCanvas->m_bmpLIC, m_lpGLCanvas->m_iBmpX, m_lpGLCanvas->m_iBmpY);
			wxstrFileName = wxString::Format("\\ElectricField%d%d%d.bmp", (iFile / 100), ((iFile / 10) % 10), (iFile % 10));
			m_lpGLCanvas->m_bmpLIC.SaveFile(m_lpGLCanvas->m_wxstrMovieDir + wxstrFileName, wxBITMAP_TYPE_BMP);

            LICPainter.m_fnInitBuffer();
            LICPainter.m_fnGenerateViewCoordinates();
            m_lpGLCanvas->m_lpMeshCenter->m_fnFillVertexNablaTheta();
            LICPainter.m_fnSetLocations(false);
            LICPainter.m_fnDetectSeam();
			LICPainter.m_fnAssignOriginalColor(m_lpGLCanvas->m_bmpOriginal);
            LICPainter.m_fnLinearTrace(200);
            LICPainter.m_fnEnhanceBitmap(true, ftMaxElectricStrength);
			LICPainter.m_fnAppendSingularities(2.0);
			LICPainter.m_fnWriteBitmap(m_lpGLCanvas->m_bmpLIC, m_lpGLCanvas->m_iBmpX, m_lpGLCanvas->m_iBmpY);
			wxstrFileName = wxString::Format("\\NablaTheta%d%d%d.bmp", (iFile / 100), ((iFile / 10) % 10), (iFile % 10));
			m_lpGLCanvas->m_bmpLIC.SaveFile(m_lpGLCanvas->m_wxstrMovieDir + wxstrFileName, wxBITMAP_TYPE_BMP);

            LICPainter.m_fnInitBuffer();
            LICPainter.m_fnGenerateViewCoordinates();
            m_lpGLCanvas->m_lpMeshCenter->m_fnFillVertexCurlField();
            LICPainter.m_fnSetLocations(false);
            LICPainter.m_fnDetectSeam();
			LICPainter.m_fnAssignOriginalColor(m_lpGLCanvas->m_bmpOriginal);
            LICPainter.m_fnLinearTrace(200);
            LICPainter.m_fnEnhanceBitmap(true, ftMaxCurlStrength);
			LICPainter.m_fnAppendSingularities(2.0);
			LICPainter.m_fnWriteBitmap(m_lpGLCanvas->m_bmpLIC, m_lpGLCanvas->m_iBmpX, m_lpGLCanvas->m_iBmpY);
			wxstrFileName = wxString::Format("\\CurlField%d%d%d.bmp", (iFile / 100), ((iFile / 10) % 10), (iFile % 10));
			m_lpGLCanvas->m_bmpLIC.SaveFile(m_lpGLCanvas->m_wxstrMovieDir + wxstrFileName, wxBITMAP_TYPE_BMP);

			LICPainter.m_fnInitBuffer();
            LICPainter.m_fnGenerateViewCoordinates();
            m_lpGLCanvas->m_lpMeshCenter->m_fnFillVertexCurlField();
            LICPainter.m_fnSetLocations(false);
            LICPainter.m_fnDetectSeam();
			LICPainter.m_fnAssignOriginalColor(m_lpGLCanvas->m_bmpOriginal);
            LICPainter.m_fnNoFieldRender(ftMaxCurlStrength);
			LICPainter.m_fnAppendSingularities(2.0);
			LICPainter.m_fnWriteBitmap(m_lpGLCanvas->m_bmpLIC, m_lpGLCanvas->m_iBmpX, m_lpGLCanvas->m_iBmpY);
			wxstrFileName = wxString::Format("\\PureCurlField%d%d%d.bmp", (iFile / 100), ((iFile / 10) % 10), (iFile % 10));
			m_lpGLCanvas->m_bmpLIC.SaveFile(m_lpGLCanvas->m_wxstrMovieDir + wxstrFileName, wxBITMAP_TYPE_BMP);

			wxstrFileName = wxString::Format("\\CrossField%d%d%d.txt", (iFile / 100), ((iFile / 10) % 10), (iFile % 10));
            ofs.open((m_lpGLCanvas->m_wxstrMovieDir + wxstrFileName).ToAscii());
            m_lpGLCanvas->m_lpMeshCenter->m_fnExportDirectionField(ofs, false);
            ofs.close();
			wxstrFileName = wxString::Format("\\Singularities%d%d%d.txt", (iFile / 100), ((iFile / 10) % 10), (iFile % 10));
            ofs.open((m_lpGLCanvas->m_wxstrMovieDir + wxstrFileName).ToAscii());
            m_lpGLCanvas->m_lpMeshCenter->m_fnExportSingularities(ofs);
            ofs.close();
		}
		if (m_lpGLCanvas->m_bRecordMode)
        {
            ofs << iFile << '\t' << ftEnergy << "\t\n";
        }
		++iFile;
	} while (m_lpGLCanvas->m_bProcessing);
	m_lpGLCanvas->m_lpMeshCenter->m_fnRestoreField();
	m_lpGLCanvas->m_lpMeshCenter->m_fnGenerateMoment();
	m_lpGLCanvas->m_lpMeshCenter->m_fnMarkSingularity();
	DA.m_fnFillPoisson(LDLTS1);
	LDLTS1.solve();
	DA.m_fnAssignDensity(LDLTS1);
	LDLTS1.clear();
	LDLTS0.clear();
	DA.m_fnGenerateNablaField();
	m_lpGLCanvas->m_lpMeshCenter->m_fnRepresentElectricField();
	if (m_lpGLCanvas->m_bMovieMode)
	{
		m_lpGLCanvas->m_bMovieMode = false;
		m_lpGLCanvas->m_wxstrMovieDir.clear();
	}
	if (m_lpGLCanvas->m_bRecordMode)
    {
        m_lpGLCanvas->m_bRecordMode = false;
        m_lpGLCanvas->m_wxstrRecordFile.clear();
        ofs.close();
    }
}


void CAlgorithmExecutor::m_fnConstrainedDecurl()
{
	CDirectionAdjustor DA(m_lpGLCanvas->m_lpMeshCenter->m_mshSurface);
	CMeshBase::CLDLTSv LDLTS0, LDLTS1;
	wxString wxstrFileName;
	std::ofstream ofs;
	CMeshBase::FT ftEnergy, ftMinEnergy, ftMaxElectricStrength, ftMaxCurlStrength;
	int iTestTimes, iFile;
	bool bMaxCurlInitialized;
	bMaxCurlInitialized = false;
    CLICPainter LICPainter(m_lpGLCanvas);
	iTestTimes = 0;
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
	ftMinEnergy = ftEnergy = DA.m_fnCurlEnergy();
	m_lpGLCanvas->m_lpMeshCenter->m_fnBackupField();
	iFile = 0;
	if (m_lpGLCanvas->m_bMovieMode)
    {
        LICPainter.m_fnVideoTestSize(2.0, 320, 720);
    }
    if (m_lpGLCanvas->m_bRecordMode)
    {
        ofs.open(m_lpGLCanvas->m_wxstrRecordFile.ToAscii());
    }
	do
	{
		DA.m_fnGenerateForces();
		DA.m_fnConstrainedDecurl();//DA.m_fnCurlAdjust();
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
		ftEnergy = DA.m_fnCurlEnergy();
		if (ftEnergy < ftMinEnergy)
		{
			iTestTimes = 0;
			ftMinEnergy = ftEnergy;
			m_lpGLCanvas->m_lpMeshCenter->m_fnBackupField();
		}
		else
		{
			++iTestTimes;
			if (iTestTimes > m_nIterateTimes)
			{
				m_lpGLCanvas->m_bProcessing = false;
			}
		}
		if (m_lpGLCanvas->m_bMovieMode)
		{
            LICPainter.m_fnInitBuffer();
            LICPainter.m_fnGenerateViewCoordinates();
            m_lpGLCanvas->m_lpMeshCenter->m_fnFillVertexCrossField();
            LICPainter.m_fnSetLocations(true);
            LICPainter.m_fnDetectSeam();
			LICPainter.m_fnAssignOriginalColor(m_lpGLCanvas->m_bmpOriginal);
            LICPainter.m_fnCrossTrace(200);
            LICPainter.m_fnEnhanceBitmap(false);
			LICPainter.m_fnAppendSingularities(2.0);
			LICPainter.m_fnWriteBitmap(m_lpGLCanvas->m_bmpLIC, m_lpGLCanvas->m_iBmpX, m_lpGLCanvas->m_iBmpY);
			wxstrFileName = wxString::Format("\\DirectionField%d%d%d.bmp", (iFile / 100), ((iFile / 10) % 10), (iFile % 10));
			m_lpGLCanvas->m_bmpLIC.SaveFile(m_lpGLCanvas->m_wxstrMovieDir + wxstrFileName, wxBITMAP_TYPE_BMP);


			m_lpGLCanvas->m_lpMeshCenter->m_fnFillVertexNablaTheta();
			ftMaxElectricStrength = m_lpGLCanvas->m_lpMeshCenter->m_ftMaxElectricField;
            LICPainter.m_fnInitBuffer();
            LICPainter.m_fnGenerateViewCoordinates();
            m_lpGLCanvas->m_lpMeshCenter->m_fnFillVertexElectricField();
            if (m_lpGLCanvas->m_lpMeshCenter->m_ftMaxElectricField > ftMaxElectricStrength)
            {
                ftMaxElectricStrength = m_lpGLCanvas->m_lpMeshCenter->m_ftMaxElectricField;
            }
            LICPainter.m_fnSetLocations(false);
            LICPainter.m_fnDetectSeam();
			LICPainter.m_fnAssignOriginalColor(m_lpGLCanvas->m_bmpOriginal);
            LICPainter.m_fnLinearTrace(200);
            LICPainter.m_fnEnhanceBitmap(true, ftMaxElectricStrength);
			LICPainter.m_fnAppendSingularities(2.0);
			LICPainter.m_fnWriteBitmap(m_lpGLCanvas->m_bmpLIC, m_lpGLCanvas->m_iBmpX, m_lpGLCanvas->m_iBmpY);
			wxstrFileName = wxString::Format("\\ElectricField%d%d%d.bmp", (iFile / 100), ((iFile / 10) % 10), (iFile % 10));
			m_lpGLCanvas->m_bmpLIC.SaveFile(m_lpGLCanvas->m_wxstrMovieDir + wxstrFileName, wxBITMAP_TYPE_BMP);

            LICPainter.m_fnInitBuffer();
            LICPainter.m_fnGenerateViewCoordinates();
            m_lpGLCanvas->m_lpMeshCenter->m_fnFillVertexNablaTheta();
            LICPainter.m_fnSetLocations(false);
            LICPainter.m_fnDetectSeam();
			LICPainter.m_fnAssignOriginalColor(m_lpGLCanvas->m_bmpOriginal);
            LICPainter.m_fnLinearTrace(200);
            LICPainter.m_fnEnhanceBitmap(true, ftMaxElectricStrength);
			LICPainter.m_fnAppendSingularities(2.0);
			LICPainter.m_fnWriteBitmap(m_lpGLCanvas->m_bmpLIC, m_lpGLCanvas->m_iBmpX, m_lpGLCanvas->m_iBmpY);
			wxstrFileName = wxString::Format("\\NablaTheta%d%d%d.bmp", (iFile / 100), ((iFile / 10) % 10), (iFile % 10));
			m_lpGLCanvas->m_bmpLIC.SaveFile(m_lpGLCanvas->m_wxstrMovieDir + wxstrFileName, wxBITMAP_TYPE_BMP);

            LICPainter.m_fnInitBuffer();
            LICPainter.m_fnGenerateViewCoordinates();
            m_lpGLCanvas->m_lpMeshCenter->m_fnFillVertexCurlField();
            if (!bMaxCurlInitialized)
            {
                ftMaxCurlStrength = m_lpGLCanvas->m_lpMeshCenter->m_ftMaxCurlField;
                bMaxCurlInitialized = true;
            }
            LICPainter.m_fnSetLocations(false);
            LICPainter.m_fnDetectSeam();
			LICPainter.m_fnAssignOriginalColor(m_lpGLCanvas->m_bmpOriginal);
            LICPainter.m_fnLinearTrace(200);
            LICPainter.m_fnEnhanceBitmap(true, ftMaxCurlStrength);
			LICPainter.m_fnAppendSingularities(2.0);
			LICPainter.m_fnWriteBitmap(m_lpGLCanvas->m_bmpLIC, m_lpGLCanvas->m_iBmpX, m_lpGLCanvas->m_iBmpY);
			wxstrFileName = wxString::Format("\\CurlField%d%d%d.bmp", (iFile / 100), ((iFile / 10) % 10), (iFile % 10));
			m_lpGLCanvas->m_bmpLIC.SaveFile(m_lpGLCanvas->m_wxstrMovieDir + wxstrFileName, wxBITMAP_TYPE_BMP);

			wxstrFileName = wxString::Format("\\CrossField%d%d%d.txt", (iFile / 100), ((iFile / 10) % 10), (iFile % 10));
            ofs.open((m_lpGLCanvas->m_wxstrMovieDir + wxstrFileName).ToAscii());
            m_lpGLCanvas->m_lpMeshCenter->m_fnExportDirectionField(ofs, false);
            ofs.close();
			wxstrFileName = wxString::Format("\\Singularities%d%d%d.txt", (iFile / 100), ((iFile / 10) % 10), (iFile % 10));
            ofs.open((m_lpGLCanvas->m_wxstrMovieDir + wxstrFileName).ToAscii());
            m_lpGLCanvas->m_lpMeshCenter->m_fnExportSingularities(ofs);
            ofs.close();
		}
		if (m_lpGLCanvas->m_bRecordMode)
        {
            ofs << iFile << '\t' << ftEnergy << "\t\n";
        }
		++iFile;
	} while (m_lpGLCanvas->m_bProcessing);
	m_lpGLCanvas->m_lpMeshCenter->m_fnRestoreField();
	m_lpGLCanvas->m_lpMeshCenter->m_fnGenerateMoment();
	m_lpGLCanvas->m_lpMeshCenter->m_fnMarkSingularity();
	DA.m_fnFillPoisson(LDLTS1);
	LDLTS1.solve();
	DA.m_fnAssignDensity(LDLTS1);
	LDLTS1.clear();
	LDLTS0.clear();
	DA.m_fnGenerateNablaField();
	m_lpGLCanvas->m_lpMeshCenter->m_fnRepresentElectricField();
	if (m_lpGLCanvas->m_bMovieMode)
	{
		m_lpGLCanvas->m_bMovieMode = false;
		m_lpGLCanvas->m_wxstrMovieDir.clear();
	}
	if (m_lpGLCanvas->m_bRecordMode)
    {
        m_lpGLCanvas->m_bRecordMode = false;
        m_lpGLCanvas->m_wxstrRecordFile.clear();
        ofs.close();
    }
}

void CAlgorithmExecutor::m_fnOptimizeMesh()
{
	int i;
    CMeshOptimizer Optimizer(m_lpGLCanvas->m_lpMeshCenter->m_mshSurface, m_lpGLCanvas->m_lpMeshCenter->m_mshRemeshedSurface);
	for (i = 0; i < m_nIterateTimes && m_lpGLCanvas->m_bProcessing; ++i)
    {
        Optimizer.m_fnFillInfo();
        Optimizer.m_fnGenerateAdjustment();
        Optimizer.m_fnAdjustVertices();
        Optimizer.m_fnProjectToOriginal();
    }
    ++(m_lpGLCanvas->m_lpMeshCenter->m_psProcessStatus);
}

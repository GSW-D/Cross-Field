#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include "CProcessStatusManager.h"
//included by "CDialogs.h"
#ifndef CPROCESSSTATUSRESERVER_H
#define CPROCESSSTATUSRESERVER_H


class CMeshInfoReserver: public CProcessStatusManager
{
public:
    struct CDataConvertor
    {
        union
        {
            struct
            {
                int iLoWord, iHiWord;
            };
            struct
            {
                double dblData;
            };
        };
    };
	void m_fnReadMesh(std::ifstream &ifs);
	void m_fnWriteMesh(std::ofstream &ofs);
	int m_fnImportParaObj(std::ifstream &ifs);
	void m_fnWriteNewMesh(std::ofstream &ofs);
	void m_fnReadNewMesh(std::ifstream &ifs);
	void m_fnExportFixedDirection(std::ofstream &ofs);
	void m_fnExportHardFeature(std::ofstream &ofs);
	void m_fnExportFixedMark(std::ofstream &ofs);
	void m_fnImportFixedDirection(std::ifstream &ifs);
	void m_fnImportDirectionField(std::ifstream &ifs, bool bWithIndex, bool bFrameField);
	void m_fnExportDirectionField(std::ofstream &ofs, bool bWithIndex);
	void m_fnExportParaDif(std::ofstream& ofs);
	void m_fnImportIntCones(std::ifstream& ifs);
	void m_fnExportMountain(std::ofstream &ofs, double dblMinClr, double dblMaxClr, double dblHeightScl = 1.0);
	void m_fnExportQuadSingularities(std::ofstream &ofs);
	void m_fnExportElectricityField(std::ofstream &ofs);
	void m_fnExportCurlField(std::ofstream &ofs);
	void m_fnExportCrossField(std::ofstream &ofs);
	void m_fnExportSingularities(std::ofstream &ofs);
	void m_fnExportCrevasse(std::ofstream &ofs);
	void m_fnExportTexturePara(std::ofstream &ofs, bool bNormalizeTexture);
	void m_fnExportSelected(std::ofstream &ofs, bool bNormalizeTexture);
	void m_fnExportPartPara(std::ofstream &ofs, bool bNormalizeTexture);
	void m_fnExportAngleError(std::ofstream &ofs, int iPara, double dblMin, double dblMax);
	void m_fnExportConfCoefColor(std::ofstream &ofs);
	void m_fnExportTriangularPara(std::ofstream &ofs);
    void m_fnReadProcessStatus(std::ifstream &ifs);
    void m_fnWriteProcessStatus(std::ofstream &ofs);
	void m_fnExportObjFile(std::ofstream &ofs);
	void m_fnExportQuadObjFile(std::ofstream &ofs);
	void m_fnExportQuadCylinder(std::ofstream &ofs);
protected:
    void m_fnWriteOriginalCoordinates(std::ofstream &ofsWriter);
    void m_fnReadOriginalCoordinates(std::ifstream &ifsReader, std::vector<double> &vecCoordinates);
    void m_fnWriteOriginalFacetLinkage(std::ofstream &ofsWriter);
    void m_fnReadOriginalFacetLinkage(std::ifstream &ifsReader, std::vector<int> &vecIndex);
    void m_fnWriteChartDir(std::ofstream &ofsWriter);
    void m_fnReadChartDir(std::ifstream &ifsReader);
    void m_fnWriteDensity(std::ofstream &ofsWriter);
	void m_fnReadDensity(std::ifstream &ifsReader);
	void m_fnWriteFreeEdges(std::ofstream &ofsWriter);
	void m_fnReadFreeEdges(std::ifstream &ifsReader);
	void m_fnWriteGlobalParameter(std::ofstream &ofsWriter);
	void m_fnReadGlobalParameter(std::ifstream &ifsWriter);
};
#endif //CPROCESSSTATUSRESERVER_H

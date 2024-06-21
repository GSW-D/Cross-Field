#include "CMeshInfoInitializer.h"
//included by "CMeshCenter.h"
#ifndef CMESHINFOREPRESENTOR_H
#define CMESHINFOREPRESENTOR_H
class CMeshInfoRepresentor: public CMeshInfoInitializer
{
public:
	void m_fnRepresentConformalCoefficient();
	void m_fnRepresentGaussCurvature();
	void m_fnRepresentAngleError(int iDir, CMeshBase::FT*lpftErr);
	void m_fnRepresentParaCurl(int iDir, CMeshBase::FT* lpftCurl);
	void m_fnRepresentLocalPara();
	void m_fnRepresentGlobalPara();
	void m_fnRepresentElectricField();
	void m_fnReportSingularities(int &iPositive, int &iNegative);
	void m_fnColorMap(double dblValue, double dblMaxValue, double dblMinValue, double *lpdblColor);
protected:
	double m_dblMeanLocalParaScl;
};
#endif
#include <float.h>
#include "MeshIO.h"
#include "CConstrainedLeastSquareSolver.hpp"
//included by "CParaConvertor.h"
#define ftHalfPi  1.570796326794896619231321691639751442e+00
#define ftTwoPi 6.283185307179586476925286766559005768e+00
#define ftThirdPi 1.047197551196597746154214461093167628e+00
#define ftPi 3.141592653589793238462643383279502884e+00
#define ftSixthPi 5.235987755982988730771072305465838140e-01
#ifndef CMESHBASE_H
#define CMESHBASH_H
namespace CMeshBase
{
    typedef CGAL::Cartesian<double> Kernel;
    typedef Kernel::Vector_3 Vector_3;
    typedef Kernel::Point_3 Point_3;
    typedef Kernel::FT FT;
    struct CMeshItems
    {
		struct Facet
		{
			union
			{
				struct
				{
					FT ftChartDir;
					FT ftFrameDir;
				};
				struct
				{
					FT lpftChebshevDir[2];
				};
			};
			FT ftChartDirCopy;
			FT ftCurvatureChart;
			FT ftEccentricity;
			FT ftTwArea;
			Point_3 ptIncent;
			FT ftRadius;
			Vector_3 vtNorm, vtPrincipalAxis, vtAuxiliaryAxis;
			Vector_3 vtNablaTheta, vtNablaPhi, vtDifNabla;
			FT lpftElectricFieldSum[2], lpftElectricFieldDif[2];
			FT lpftLocalParaDir[2], lpftLocalParaScl[2];
			FT ftLocalParaScl;
			FT lpftAngleErr[2], lpftCurl[2];
			int lpiTraceInd[4];
			Vector_3 lpvtTrace[4];
			signed long iDistOrd;
			unsigned bDirFixed;
			unsigned bUnif;
			unsigned bSelected;
			double lpdblColor[3];
			int iIndex;
			int iTempIndex;
			int iEquationIndex;
		};
		struct Halfedge
		{
			unsigned iLocalInd;
			FT ftDot, ftAngle, ftSqLen, ftLen, ftLaplaceCoef, ftPolarAxisDif, ftDirectionalCoef;
			Vector_3 vtVector;
			unsigned bFreeEdge;
			signed iSeamType;
			signed iSeamTypeCopy;
			unsigned bInStack;
			unsigned bSelected;
			unsigned bSharp;
			FT ftMom;
			FT ftCurl;
			FT ftForce, ftForceBarrier;
			FT *lpftGlobalPara;
			int *lpiGlobalIndex;
			CGAL::Polyhedron<Kernel, CMeshItems>::Halfedge_handle hhCrevassePrev, hhCrevasseNext;
			Vector_3 lpvtViewField[2];
			double lpdblColor[9];
			int iCrevasseId;
			signed iAccSeamType;
			int iTempIndex;
		};
		struct Vertex
		{
			signed iSourceId;
			signed iBordId;
			signed nNumPorts;
			signed iDegreeBias;
			unsigned nCrevasseSize;
			unsigned bDistConf;
			unsigned iSelected;
			CGAL::Polyhedron<Kernel, CMeshItems>::Halfedge_handle hhCome;
			FT ftAngleDefect;
			FT ftDist;
			FT ftAngleBias;
			FT ftDensity, ftViewDensity;
			Point_3 ptMountainPoint;
			FT ftVoronoiArea;
			FT ftGaussCurvature;
			Vector_3 vtCurlForce, vtThetaForce, vtPhiForce, vtViewForce;
			FT ftSqCurlForceLen;
			signed int iHeapPos;
			FT ftDensityAdjustment;
			double lpdblColor[3];
			double lpdblViewCoord[3];
			int iTempIndex;
			int iEquationIndex;
			int iViewX, iViewY;
			FT ftDepth, ftFieldStrength;
		};
    };
    typedef CGAL::Polyhedron<Kernel, CMeshItems> CPolyhedron;
	struct CSamplerItems
	{
		struct Facet
		{
			unsigned nDepth;
			FT ftArea, ftLogArea;
		};
		struct Halfedge
		{
			Vector_3 vtVector;
			FT ftSqLen, ftLen, ftAngle, ftLogLenRatio;
			bool bCrevasse;
			unsigned nDepth;
		};
		struct Vertex
		{
			CPolyhedron::Halfedge_handle hhBase;
			Vector_3 vtAdjustment;
			int iFreedom;
			int iBorderId;
			bool bVisible;
			Point_3 ptViewCoord;
			unsigned nDepth;
			int iViewX, iViewY;
			FT ftDepth;
			int iMark;
		};
	};
    typedef CGAL::Polyhedron<Kernel, CSamplerItems> CSamplerPolyhedron;
	typedef CLeastSquareSolver<FT> CLstSqrSv;
	typedef CConstrainedLeastSquareSolver<FT> CConsLstSqrSv;
	typedef CMinimalNormSolver<FT> CMinNormSv;
	typedef CLDLTSolver<FT> CLDLTSv;
	typedef CLDLTSolver<std::complex<CMeshBase::FT> > CCmplxSv;
};

#endif // CMESHBASE_H

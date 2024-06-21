#include "CMeshBase.h"
//included by "CTracer.h"
#ifndef CPARACONVERTOR_H
#define CPARACONVERTOR_H
class CParaConvertor
{
    public:
        CParaConvertor();
        ~CParaConvertor();
		static void m_fnParameterConvert(int iSeamType, CMeshBase::FT *lpftInput, CMeshBase::FT *lpftOutput, CMeshBase::FT *lpftOffset = NULL);
        static void m_fnParameterInverse(int iSeamType, CMeshBase::FT *lpftInput, CMeshBase::FT *lpftOutput, CMeshBase::FT *lpftOffset = NULL);
        static void m_fnAdjustmentConvert(int iSeamType, CMeshBase::FT *lpftInput, CMeshBase::FT *lpftOutput);
        static void m_fnExtractOffset(int iSeamType, CMeshBase::FT *lpftOriginal, CMeshBase::FT *lpftImage, CMeshBase::FT *lpftOffset);
        static void m_fnSolveFixedPoint(int iSeamType, CMeshBase::FT *lpftOffset, CMeshBase::FT *lpftSolution);
    protected:
    private:
};

#endif // CPARACONVERTOR_H

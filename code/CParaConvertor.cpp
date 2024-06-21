#include "CParaConvertor.h"

CParaConvertor::CParaConvertor()
{
    //ctor
}

CParaConvertor::~CParaConvertor()
{
    //dtor
}


void CParaConvertor::m_fnParameterConvert(int iSeamType, CMeshBase::FT *lpftInput, CMeshBase::FT *lpftOutput, CMeshBase::FT *lpftOffset)//, CMeshBase::FT ftCoef)
{
    int lpiSign[2], iRotate, i;
	CMeshBase::FT lpftBuffer[2];
	if (lpftInput != NULL)
	{
		for(i = 0; i < 2; ++i)
		{
			lpiSign[i] = (((iSeamType + 2 + i) % 4) / 2) * 2 - 1;
		}
		iRotate = (iSeamType + 2) % 2;
		for(i = 0; i < 2; ++i)
		{
			lpftBuffer[i] = lpftInput[(iRotate + i) % 2] * lpiSign[i];
		}
	}
	else
	{
		memset(lpftBuffer, 0, 2 * sizeof(CMeshBase::FT));
	}
    if(lpftOffset != NULL)
    {
        for(i = 0; i < 2; ++i)
        {
            lpftBuffer[i] += lpftOffset[i];
        }
    }
	memcpy(lpftOutput, lpftBuffer, 2 * sizeof(CMeshBase::FT));
}

void CParaConvertor::m_fnParameterInverse(int iSeamType, CMeshBase::FT *lpftInput, CMeshBase::FT *lpftOutput, CMeshBase::FT *lpftOffset)
{
    int lpiSign[2], iRotate, i;
	CMeshBase::FT lpftBuffer[2];
	if (lpftInput != NULL)
	{
		for(i = 0; i < 2; ++i)
		{
			lpiSign[i] = (((iSeamType + 2 + i) % 4) / 2) * 2 - 1;
		}
		iRotate = (iSeamType + 2) % 2;
		for(i = 0; i < 2; ++i)
		{
			lpftBuffer[(iRotate + i) % 2] = lpftInput[i] * lpiSign[i];
		}
	}
	else
	{
		memset(lpftBuffer, 0, 2 * sizeof(CMeshBase::FT));
	}
    if(lpftOffset != NULL)
    {
        for(i = 0; i < 2; ++i)
        {
            lpftBuffer[(iRotate + i) % 2] -= lpftOffset[i] * lpiSign[i];
        }
    }
	memcpy(lpftOutput, lpftBuffer, 2 * sizeof(CMeshBase::FT));
}

void CParaConvertor::m_fnAdjustmentConvert(int iSeamType, CMeshBase::FT *lpftInput, CMeshBase::FT *lpftOutput)
{
    int lpiSign[2], iRotate, i;
	CMeshBase::FT lpftBuffer[2];
	if (lpftInput != NULL)
	{
		for(i = 0; i < 2; ++i)
		{
			lpiSign[i] = (((iSeamType + 4 + i) % 4) / 2) * 2 - 1;
		}
		iRotate = (iSeamType + 2) % 2;
		for(i = 0; i < 2; ++i)
		{
			lpftBuffer[i] = lpftInput[(iRotate + i) % 2] * lpiSign[i];
		}
	}
	else
	{
		memset(lpftBuffer, 0, 2 * sizeof(CMeshBase::FT));
	}
	memcpy(lpftOutput, lpftBuffer, 2 * sizeof(CMeshBase::FT));
}




void CParaConvertor::m_fnExtractOffset(int iSeamType, CMeshBase::FT *lpftOriginal, CMeshBase::FT *lpftImage, CMeshBase::FT *lpftOffset)
{
    int lpiSign[2], iRotate, i;
	if (lpftOriginal != NULL)
	{
		for(i = 0; i < 2; ++i)
		{
			lpiSign[i] = (((iSeamType + 2 + i) % 4) / 2) * 2 - 1;
		}
		iRotate = (iSeamType + 2) % 2;
		for(i = 0; i < 2; ++i)
		{
			lpftOffset[i] = lpftImage[i] - lpftOriginal[(iRotate + i) % 2] * lpiSign[i];
		}
	}
	else
	{
		memcpy(lpftOffset, lpftImage, 2 * sizeof(CMeshBase::FT));
	}
}

void CParaConvertor::m_fnSolveFixedPoint(int iSeamType, CMeshBase::FT *lpftOffset, CMeshBase::FT *lpftSolution)
{
    int i;
    switch(iSeamType)
    {
    case -2:
        for(i = 0; i < 2; ++i)
        {
            lpftSolution[i] = lpftOffset[i] / 2;
        }
        break;
    case -1:
        for(i = 0; i < 2; ++i)
        {
            lpftSolution[i] = (lpftOffset[0] + lpftOffset[1] * (i * 2 - 1)) / 2;
        }
        break;
    case 0:
        for(i = 0; i < 2; ++i)
        {
            lpftSolution[i] = 0;
        }
        break;
    case 1:
        for(i = 0; i < 2; ++i)
        {
            lpftSolution[i] = (lpftOffset[1] + lpftOffset[0] * (1 - i * 2)) / 2;
        }
        break;
    }
}

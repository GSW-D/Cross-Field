#include "CDistanceHeap.h"
//included by "CMeshInfoInitializer.h"
#ifndef CCREVASSEPROCESSOR_H
#define CCREVASSEPROCESSOR_H
struct CCrevasseProcessor
{
	struct CCrevasse: public CParaConvertor
	{
		int iCountStatus, lpiEqInd[2], lpiIntHead[2], lpiIntGroupInd[2], lpiIntGroupSize[2], lpiIntMark[2];
		CMeshBase::FT lpftParaRange[2], lpftGroupRange[2];
		CMeshBase::FT lpftIntRange[2];
		CMeshBase::FT lpftOffset[2];
		std::vector<CMeshBase::CPolyhedron::Halfedge_handle> m_vechhElements;
		CCrevasse* m_lpPrev;
		CCrevasse* m_lpNext;
		CCrevasse* m_lpOpposite;
		CCrevasse* prev(){ return m_lpPrev; }
		CCrevasse* next(){ return m_lpNext; }
		CCrevasse* opposite(){ return m_lpOpposite; }
		CMeshBase::CPolyhedron::Halfedge_handle first_element(){ return *(m_vechhElements.begin()); }
		CMeshBase::CPolyhedron::Halfedge_handle last_element(){ return *(m_vechhElements.rbegin()); }
		bool is_border(){ return m_vechhElements[0]->is_border(); }
		int seamtype(){ return m_vechhElements[0]->iSeamType; }
		void parameter_convert(CMeshBase::FT *lpftInput, CMeshBase::FT *lpftOutput)
		{
			m_fnParameterConvert(seamtype(), lpftInput, lpftOutput, lpftOffset);
		}
		void parameter_inverse(CMeshBase::FT *lpftInput, CMeshBase::FT *lpftOutput)
		{
			m_fnParameterInverse(seamtype(), lpftInput, lpftOutput, lpftOffset);
		}
		int combine_inner_maps_next(CMeshBase::FT *lpftInput, CMeshBase::FT *lpftOutput, int iAccSeamType = 0)
		{
			CCrevasse *lpCrevasseRotate;
			memcpy(lpftOutput, lpftInput, 2 * sizeof(CMeshBase::FT));
			if (is_border())
			{
				lpCrevasseRotate = this;
				do
				{
					iAccSeamType = (iAccSeamType + lpCrevasseRotate->seamtype() + 6) % 4 - 2;
					lpCrevasseRotate->parameter_convert(lpftOutput, lpftOutput);
					lpCrevasseRotate = lpCrevasseRotate->opposite()->prev();
				} while (lpCrevasseRotate != this);
			}
			return iAccSeamType;
		}
		int combine_inner_maps_prev(CMeshBase::FT *lpftInput, CMeshBase::FT *lpftOutput, int iAccSeamType = 0)
		{
			CCrevasse *lpCrevasseRotate;
			memcpy(lpftOutput, lpftInput, 2 * sizeof(CMeshBase::FT));
			if (is_border())
			{
				lpCrevasseRotate = this;
				do
				{
					iAccSeamType = (iAccSeamType + lpCrevasseRotate->seamtype() + 6) % 4 - 2;
					lpCrevasseRotate->parameter_convert(lpftOutput, lpftOutput);
					lpCrevasseRotate = lpCrevasseRotate->opposite()->next();
				} while (lpCrevasseRotate != this);
			}
			return iAccSeamType;
		}
	};
};
#endif
#include "CComplexOptimizer.h"
CComplexOptimizer::CComplexOptimizer(CMeshBase::CPolyhedron &mshSurface) : m_mshSurface(mshSurface)
{
	;
}

void CComplexOptimizer::m_fnGenerateFreeMoment()
{
	CMeshBase::CPolyhedron::Edge_iterator ei;
	CMeshBase::CPolyhedron::Facet_iterator fi;
	CMeshBase::FT ftAngleDif;
	for (ei = m_mshSurface.edges_begin(); ei != m_mshSurface.edges_end(); ++ei)
	{
		if (ei->is_border_edge())
		{
			ei->ftMom = 0;
			ei->opposite()->ftMom = 0;
		}
		else
		{
			ftAngleDif = ei->opposite()->facet()->ftChartDir - ei->facet()->ftChartDir + ei->ftPolarAxisDif;
			ei->ftMom = ftAngleDif;
			ei->ftMom = ei->ftMom - CMeshBase::FT(int(ei->ftMom / ftTwoPi + 15.5) - 15) * ftTwoPi;
			ei->opposite()->ftMom = -ei->ftMom;
		}
	}
}

void CComplexOptimizer::m_fnFillComplexMatrix(CMeshBase::CCmplxSv &CMS)
{
	CMeshBase::CPolyhedron::Facet_iterator fi;
	std::queue<CMeshBase::CPolyhedron::Facet_handle> qufh;
	CMeshBase::CPolyhedron::Facet_handle fh0, fh;
	int iIndex;
	MatElem<std::complex <CMeshBase::FT> > meDiag, meTri;
	CMeshBase::CPolyhedron::Halfedge_handle hh0, hh;
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
	{
		fi->iEquationIndex = -1;
	}
	iIndex = 0;
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
	{
		if (fi->iEquationIndex == -1)
		{
			fi->iEquationIndex = iIndex;
			qufh.push(&(*fi));
			while (!qufh.empty())
			{
				fh0 = qufh.front();
				qufh.pop();
				hh0 = fh0->halfedge();
				hh = hh0;
				do
				{
					if (!hh->opposite()->is_border() && hh->ftLaplaceCoef <= 0.0)
					{
						fh = hh->opposite()->facet();
						if (fh->iEquationIndex == -1)
						{
							fh->ftChartDir -= hh->ftMom;
							hh->opposite()->next()->ftMom += hh->ftMom;
							hh->opposite()->next()->opposite()->ftMom = -hh->opposite()->next()->ftMom;
							hh->opposite()->prev()->ftMom += hh->ftMom;
							hh->opposite()->prev()->opposite()->ftMom = -hh->opposite()->prev()->ftMom;
							hh->ftMom = 0;
							hh->opposite()->ftMom = 0;
							fh->iEquationIndex = iIndex;
							qufh.push(fh);
						}
					}
					hh = hh->next();
				} while (hh != hh0);
			}
			++iIndex;
		}
	}
	CMS.init(iIndex);
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
	{
		if (fi->iEquationIndex != -1)
		{
			meDiag.row = fi->iEquationIndex;
			meDiag.col = fi->iEquationIndex;
			meDiag.data = 0;
			hh0 = fi->halfedge();
			hh = hh0;
			do
			{
				if (!hh->opposite()->is_border() && hh->ftLaplaceCoef > 0.0)
				{
					meDiag.data += 1.0 / hh->ftLaplaceCoef;
					if (hh->opposite()->facet()->iEquationIndex != -1)
					{
						meTri.row = fi->iEquationIndex;
						meTri.col = hh->opposite()->facet()->iEquationIndex;
						meTri.data = -std::complex<CMeshBase::FT>(cos(hh->ftMom * 4), sin(hh->ftMom * 4)) / (hh->ftLaplaceCoef * sqrt(fi->ftTwArea * hh->opposite()->facet()->ftTwArea));
						CMS.add_data(meTri.row, meTri.col, meTri.data);
					}
				}
				hh = hh->next();
			} while (hh != hh0);
			meDiag.data /= fi->ftTwArea;
			CMS.add_data(meDiag.row, meDiag.col, meDiag.data);
		}
	}
}

void CComplexOptimizer::m_fnInitVector(CMeshBase::CCmplxSv &CMS)
{
	int i;
	std::complex<CMeshBase::FT> cmpOne, *lpcmpCur;
	cmpOne = std::complex<CMeshBase::FT>(1.0);
	lpcmpCur = CMS.m_b;
	for (i = 0; i < CMS.m_dim; ++i)
	{
		*lpcmpCur = cmpOne;
		++lpcmpCur;
	}
}

void CComplexOptimizer::m_fnNormalizeVector(CMeshBase::CCmplxSv &CMS)
{
	CMeshBase::FT ftSumNorm;
	std::complex<CMeshBase::FT> *lpcmpX, *lpcmpB;
	int i;
	ftSumNorm = 0;
	lpcmpX = CMS.m_x;
	lpcmpB = CMS.m_b;
	for (i = 0; i < CMS.m_dim; ++i)
	{
		ftSumNorm += lpcmpX->real() * lpcmpX->real() + lpcmpX->imag() * lpcmpX->imag();
		++lpcmpX;
	}
	ftSumNorm = sqrt(ftSumNorm / CMS.m_dim);
	lpcmpB = CMS.m_b;
	lpcmpX = CMS.m_x;
	for (i = 0; i < CMS.m_dim; ++i)
	{
		*lpcmpB = *lpcmpX / ftSumNorm;
		++lpcmpB;
		++lpcmpX;
	}
}

void CComplexOptimizer::m_fnFillDirection(CMeshBase::CCmplxSv &CMS)
{
	CMeshBase::CPolyhedron::Facet_iterator fi;
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
	{
		if (fi->iEquationIndex != -1)
		{
			fi->ftChartDir += atan2(CMS.m_x[fi->iEquationIndex].imag(), CMS.m_x[fi->iEquationIndex].real()) / 4;
			fi->ftChartDir = fi->ftChartDir - CMeshBase::FT(int(fi->ftChartDir / ftTwoPi + 15.5) - 15) * ftTwoPi;
		}
	}
}
#include <string.h>
#include <stdlib.h>
#include <fstream>
#include "CPolyhedron.h"
//included by "CMeshBase.h"
#ifndef MESHIO_H
#define MESHIO_H
struct StringSeparator
{
	char lpstrData[256], *lpcLocalStart, *lpcLocalEnd;
	StringSeparator(){ lpstrData[0] = '\0'; lpcLocalEnd = lpstrData; }
	StringSeparator& operator>>(unsigned int &iInteger)
	{
		int nLen;
		char lpstrBuffer[256];
		lpcLocalStart = lpcLocalEnd;
		while ((*lpcLocalStart < '0' || *lpcLocalStart > '9') && (*lpcLocalStart != '\0'))
		{
			++lpcLocalStart;
		}
		lpcLocalEnd = lpcLocalStart;
		while ((*lpcLocalEnd >= '0' && *lpcLocalEnd <= '9') && (*lpcLocalEnd != '\0'))
		{
			++lpcLocalEnd;
		}
		nLen = int(lpcLocalEnd - lpcLocalStart);
		memcpy(lpstrBuffer, lpcLocalStart, nLen * sizeof(char));
		lpstrBuffer[nLen] = '\0';
		iInteger = atoi(lpstrBuffer);
		return *this;
	}
	StringSeparator& operator>>(double &dblReal)
	{
		int nLen;
		char lpstrBuffer[256];
		lpcLocalStart = lpcLocalEnd;
		while ((*lpcLocalStart < '0' || *lpcLocalStart > '9') && (*lpcLocalStart != '.' && *lpcLocalStart != '-' && *lpcLocalStart != 'e' && *lpcLocalStart != 'E' && *lpcLocalStart != '\0'))
		{
			++lpcLocalStart;
		}
		lpcLocalEnd = lpcLocalStart;
		while ((*lpcLocalEnd >= '0' && *lpcLocalEnd <= '9' || *lpcLocalEnd == '.' || *lpcLocalEnd == '-' || *lpcLocalEnd == 'e' || *lpcLocalEnd == 'E') && (*lpcLocalEnd != '\0'))
		{
			++lpcLocalEnd;
		}
		nLen = int(lpcLocalEnd - lpcLocalStart);
		memcpy(lpstrBuffer, lpcLocalStart, nLen * sizeof(char));
		lpstrBuffer[nLen] = '\0';
		dblReal = atof(lpstrBuffer);
		return *this;
	}
	void clear(){ lpstrData[0] = '\0'; lpcLocalEnd = lpstrData; }
};

template <typename PolyhedronTraits, typename PolyhedronItems>
static std::istream& operator>>(std::istream& istr, CGAL::Polyhedron<PolyhedronTraits, PolyhedronItems> &surface)
{
	StringSeparator strsep;
	typename PolyhedronTraits::FT ftX, ftY, ftZ;
	unsigned int nOldNumV, nNumV, nNumE, nNumF, iLineAcc, iReadStatus, iVertexInd, nFacetSize, iFacetVertex;
	unsigned int *lpiIndex;
	typename CGAL::Polyhedron<PolyhedronTraits, PolyhedronItems>::Vertex_iterator vi;
//	CGAL::Polyhedron<PolyhedronTraits, PolyhedronItems>::Halfedge_handle hh0, hh;
	bool bLegal;
	istr.getline(strsep.lpstrData, 255);
	if (strcmp(strsep.lpstrData, "OFF") == 0)
	{
		nOldNumV = surface.m_lstVertices.size();
		strsep.clear();
		iReadStatus = 0;
		while (!istr.eof())
		{
			istr.getline(strsep.lpstrData, 255);
			if (strlen(strsep.lpstrData) > 0)
			{
				if (strsep.lpstrData[0] != '#')
				{
					switch (iReadStatus)
					{
					case 0:
					{
						strsep >> nNumV >> nNumF >> nNumE;
						strsep.clear();
						iLineAcc = 0;
						++iReadStatus;
					}
						break;
					case 1:
					{
						strsep >> ftX >> ftY >> ftZ;
						strsep.clear();
						surface.append_vertex(ftX, ftY, ftZ);
						++iLineAcc;
						if (iLineAcc == nNumV)
						{
							surface.begin_facets();
							iLineAcc = 0;
							++iReadStatus;
						}
					}
						break;
					case 2:
					{
						strsep >> nFacetSize;
						if (nFacetSize > 0)
						{
							lpiIndex = new unsigned int[nFacetSize];
							bLegal = true;
							for (iFacetVertex = 0; iFacetVertex < nFacetSize; ++iFacetVertex)
							{
								strsep >> iVertexInd;
								if (iVertexInd >= 0 && iVertexInd + nOldNumV < surface.m_lstVertices.size())
								{
									lpiIndex[iFacetVertex] = iVertexInd + nOldNumV;
								}
								else
								{
									bLegal = false;
									break;
								}
							}
							if (bLegal)
							{
								surface.append_facet(lpiIndex, nFacetSize);
							}
							delete[]lpiIndex;
							strsep.clear();
							++iLineAcc;
							if (iLineAcc == nNumF)
							{
								iLineAcc = 0;
								++iReadStatus;
								surface.end_facets();
							}
						}
					}
						break;
					}
				}
			}
		}
	}
	surface.refresh_degree();
	return istr;
}

template <typename PolyhedronTraits, typename PolyhedronItems>
static std::ostream& operator<<(std::ostream& ostr, CGAL::Polyhedron<PolyhedronTraits, PolyhedronItems> &surface)
{
	typename CGAL::Polyhedron<PolyhedronTraits, PolyhedronItems>::Vertex_iterator vi;
	typename CGAL::Polyhedron<PolyhedronTraits, PolyhedronItems>::Facet_iterator fi;
	typename CGAL::Polyhedron<PolyhedronTraits, PolyhedronItems>::Halfedge_handle hh0, hh;
	unsigned int iIndex;
	iIndex = 0;
	for (vi = surface.m_lstVertices.begin(); vi != surface.m_lstVertices.end(); ++vi)
	{
		vi->iIndex = iIndex;
		++iIndex;
	}
	ostr << "OFF\n";
	ostr << surface.m_lstVertices.size() << ' ' << surface.m_lstFacets.size() << ' ' << surface.m_lstHalfedges.size() / 2 << " \n";
	for (vi = surface.m_lstVertices.begin(); vi != surface.m_lstVertices.end(); ++vi)
	{
		ostr << vi->point();
	}
	for (fi = surface.m_lstFacets.begin(); fi != surface.m_lstFacets.end(); ++fi)
	{
		ostr << fi->m_nSize << ' ';
		hh0 = fi->m_lpHE;
		hh = hh0;
		do
		{
			ostr << hh->m_lpVertex->iIndex << ' ';
			hh = hh->m_lpNext;
		} while (hh != hh0);
		ostr << '\n';
	}
	return ostr;
}
#endif

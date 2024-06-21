#include <list>
#include <stdarg.h>
#include "CGAL.h"
namespace CGAL
{
	struct Polyhedron_items
	{
		struct Facet{};
		struct Halfedge{};
		struct Vertex{};
	};
	template <typename PolyhedronTraits, typename PolyhedronItems = Polyhedron_items>
	struct Polyhedron
	{
		struct HalfedgeDS_facet_base;
		struct HalfedgeDS_halfedge_base;
		struct HalfedgeDS_vertex_base;
		typedef HalfedgeDS_vertex_base *Vertex_handle;
		typedef HalfedgeDS_halfedge_base *Halfedge_handle;
		typedef HalfedgeDS_facet_base *Facet_handle;
		struct HalfedgeDS_facet_base : PolyhedronItems::Facet
		{
			Halfedge_handle m_lpHE;
			unsigned int m_nSize;
			HalfedgeDS_facet_base(){ m_lpHE = NULL; m_nSize = 0; }
			Halfedge_handle halfedge(){ return m_lpHE; }
			unsigned int size() { return m_nSize; }
		};
		struct HalfedgeDS_halfedge_base : PolyhedronItems::Halfedge
		{
			Halfedge_handle m_lpPrev, m_lpNext, m_lpOpposite;
			Vertex_handle m_lpVertex;
			Facet_handle m_lpFacet;
			HalfedgeDS_halfedge_base(){ m_lpPrev = m_lpNext = m_lpOpposite = NULL; m_lpVertex = NULL; m_lpFacet = NULL; }
			Halfedge_handle next(){ return m_lpNext; }
			Halfedge_handle prev(){ return m_lpPrev; }
			Halfedge_handle opposite(){ return m_lpOpposite; }
			Facet_handle facet(){ return m_lpFacet; }
			Vertex_handle vertex(){ return m_lpVertex; }
			bool is_border(){ return m_lpFacet == NULL; }
			bool is_border_edge(){ return (is_border() || opposite()->is_border()); }
		};
		struct HalfedgeDS_vertex_base : PolyhedronItems::Vertex
		{
			Halfedge_handle m_lpHE;
			typename PolyhedronTraits::Point_3 ptPosition;
			unsigned int m_nDegree;
			unsigned int iIndex;
			HalfedgeDS_vertex_base(){ m_lpHE = NULL; m_nDegree = 0; }
			HalfedgeDS_vertex_base(typename PolyhedronTraits::FT _ftX, typename PolyhedronTraits::FT _ftY, typename PolyhedronTraits::FT _ftZ) : ptPosition(_ftX, _ftY, _ftZ){ m_lpHE = NULL; m_nDegree = 0; }
			Halfedge_handle halfedge(){ return m_lpHE; }
			typename PolyhedronTraits::Point_3& point(){ return ptPosition; };
			unsigned int degree(){ return m_nDegree; }
		};
		typedef typename PolyhedronTraits::FT FT;
		typedef typename PolyhedronTraits::Vector_3 Vector_3;
		typedef typename PolyhedronTraits::Point_3 Point_3;
		struct Vertex_iterator : std::list<HalfedgeDS_vertex_base>::iterator
		{
			Vertex_iterator(typename std::list<HalfedgeDS_vertex_base>::iterator vi) : std::list<HalfedgeDS_vertex_base>::iterator(vi)
			{
			}
			Vertex_iterator(){}
			Vertex_handle operator()()
			{
				return &(*this);
			}
		};
		struct Halfedge_iterator : std::list<HalfedgeDS_halfedge_base>::iterator
		{
			Halfedge_iterator(typename std::list<HalfedgeDS_halfedge_base>::iterator hi) : std::list<HalfedgeDS_halfedge_base>::iterator(hi)
			{
			}
			Halfedge_iterator(){}
			Halfedge_handle operator()()
			{
				return &(*this);
			}
		};
		struct Facet_iterator : std::list<HalfedgeDS_facet_base>::iterator
		{
			Facet_iterator(typename std::list<HalfedgeDS_facet_base>::iterator fi) : std::list<HalfedgeDS_facet_base>::iterator(fi)
			{
			}
			Facet_iterator(){}
			Facet_handle operator()()
			{
				return &(this);
			}
		};
		std::list<HalfedgeDS_vertex_base> m_lstVertices;
		std::list<HalfedgeDS_halfedge_base> m_lstHalfedges;
		std::list<HalfedgeDS_facet_base> m_lstFacets;
		HalfedgeDS_vertex_base **m_lplpVertexOrder;
		struct Edge_iterator : public Halfedge_iterator
		{
			Edge_iterator(){}
			Edge_iterator(Halfedge_iterator hi) : Halfedge_iterator(hi){}
			Edge_iterator& operator++()
			{
				Halfedge_iterator::operator++();
				Halfedge_iterator::operator++();
				return *this;
			}
//			Halfedge_handle operator()()
//			{
//				return &(*Halfedge_iterator);
//			}
		};
		Polyhedron(){ m_lplpVertexOrder = NULL; }
		Vertex_iterator vertices_begin(){ return Vertex_iterator(m_lstVertices.begin()); }
		Vertex_iterator vertices_end(){ return Vertex_iterator(m_lstVertices.end()); }
		Halfedge_iterator halfedges_begin(){ return Halfedge_iterator(m_lstHalfedges.begin()); }
		Halfedge_iterator halfedges_end(){ return Halfedge_iterator(m_lstHalfedges.end()); }
		Facet_iterator facets_begin(){ return Facet_iterator(m_lstFacets.begin()); }
		Facet_iterator facets_end(){ return Facet_iterator(m_lstFacets.end()); }
		Edge_iterator edges_begin(){ return Edge_iterator(halfedges_begin()); }
		Edge_iterator edges_end(){ return Edge_iterator(halfedges_end()); }
		unsigned int size_of_vertices(){ return int(m_lstVertices.size()); }
		unsigned int size_of_facets() { return int(m_lstFacets.size()); };
		unsigned int size_of_edges(){ return int(m_lstHalfedges.size()) / 2; }
		unsigned int size_of_halfedges(){ return int(m_lstHalfedges.size()); }
		void clear(){ m_lstVertices.clear(); m_lstFacets.clear(); m_lstHalfedges.clear(); };
		void begin_facets()
		{
			unsigned int i, n;
			Vertex_iterator vi;
			if (m_lplpVertexOrder != NULL)
			{
				delete[]m_lplpVertexOrder;
			}
			n = int(m_lstVertices.size());
			if (n > 0)
			{
				m_lplpVertexOrder = new HalfedgeDS_vertex_base*[m_lstVertices.size()];
				i = 0;
				for (vi = m_lstVertices.begin(); vi != m_lstVertices.end(); ++vi)
				{
					m_lplpVertexOrder[i] = &(*vi);
					++i;
				}
			}
		}
		void end_facets(bool bErase = true)
		{
			delete[]m_lplpVertexOrder;
			m_lplpVertexOrder = NULL;
			Vertex_iterator vi;
			if (bErase)
			{
				for (vi = m_lstVertices.begin(); vi != m_lstVertices.end(); ++vi)
				{
					if (vi->m_lpHE == NULL)
					{
						vi = m_lstVertices.erase(vi);
						--vi;
					}
				}
			}
		}
		void append_vertex(FT ftX, FT ftY, FT ftZ);
		void append_facet(HalfedgeDS_vertex_base **lplpFacetVertices, unsigned int nFacetSize);
		void append_facet(unsigned int *lpiIndex, unsigned int nFacetSize);
		void append_facet(unsigned int nFacetSize, ...);
		void refresh_degree();
		Halfedge_handle flip_edge(Halfedge_handle hh);
		Halfedge_handle split_edge(Halfedge_handle hh);
		Halfedge_handle split_facet(Halfedge_handle hh0, Halfedge_handle hh1);
		Halfedge_handle split_loop(Halfedge_handle hh0, Halfedge_handle hh1, Halfedge_handle hh2);
	};
	template <typename Kernel, typename PolyherdonItems>
	void Polyhedron<Kernel, PolyherdonItems>::append_vertex(typename Kernel::FT ftX, typename Kernel::FT ftY, typename Kernel::FT ftZ)
	{
		m_lstVertices.push_back(HalfedgeDS_vertex_base(ftX, ftY, ftZ));
	}

	template <typename Kernel, typename PolyherdonItems>
	void Polyhedron<Kernel, PolyherdonItems>::append_facet(HalfedgeDS_vertex_base **lplpFacetVertices, unsigned int nFacetSize)
	{
		unsigned int i, j;
		HalfedgeDS_vertex_base *lpvtTo, *lpvtFrom;
		HalfedgeDS_halfedge_base newHalfedge;
		HalfedgeDS_facet_base newFacet, *lpFacet;
		Halfedge_iterator hi0, hi1;
		Halfedge_handle hh0, hh;
		Facet_iterator fi;
		HalfedgeDS_halfedge_base **lplpFacetHalfedges, **lplpBorderHalfedges;
		bool bLegal;
		if (nFacetSize > 0)
		{
			bLegal = true;
			for (i = 0; i < nFacetSize; ++i)
			{
				j = (i + 1) % nFacetSize;
				lpvtTo = lplpFacetVertices[j];
				lpvtFrom = lplpFacetVertices[i];
				hh = NULL;
				if (lpvtTo->m_lpHE != NULL)
				{
					bLegal = false;
					hh0 = lpvtTo->m_lpHE;
					hh = hh0;
					//Test whether there are border edges to append facet
					do
					{
						if (hh->m_lpFacet == NULL)
						{
							bLegal = true;
						}
						hh = hh->m_lpNext->m_lpOpposite;
					} while (hh != hh0);
					//Test whether the edge linking the 2 points have conflict facet
					do
					{
						if (hh->m_lpOpposite->m_lpVertex == lpvtFrom)
						{
							if (hh->m_lpFacet != NULL)
							{
								bLegal = false;
							}
							break;
						}
						hh = hh->m_lpNext->m_lpOpposite;
					} while (hh != hh0);
				}
				if (!bLegal)
				{
					break;
				}
			}
		}
		else
		{
			bLegal = false;
		}
		if (bLegal)
		{
			lplpFacetHalfedges = new HalfedgeDS_halfedge_base*[nFacetSize];
			lplpBorderHalfedges = new HalfedgeDS_halfedge_base*[nFacetSize];
			newFacet.m_nSize = nFacetSize;
			m_lstFacets.push_back(newFacet);
			fi = m_lstFacets.end();
			--fi;
			lpFacet = &(*fi);
			newHalfedge.m_lpFacet = NULL;
			newHalfedge.m_lpNext = NULL;
			newHalfedge.m_lpPrev = NULL;
			for (i = 0; i < nFacetSize; ++i)
			{
				j = (i + 1) % nFacetSize;
				lpvtTo = lplpFacetVertices[j];
				lpvtFrom = lplpFacetVertices[i];
				hh = NULL;
				//Store the borders to append facet
				if (lpvtTo->m_lpHE != NULL)
				{
					hh0 = lpvtTo->m_lpHE;
					hh = hh0;
					//Find an edge linking the 2 points
					do
					{
						if (hh->m_lpOpposite->m_lpVertex == lpvtFrom)
						{
							break;
						}
						hh = hh->m_lpNext->m_lpOpposite;
					} while (hh != hh0);
					if (hh->m_lpOpposite->m_lpVertex != lpvtFrom)//If there is no edge linking the 2 points, a border edge is required to find
					{
						hh = hh0;
						while (hh->m_lpFacet != NULL)
						{
							hh = hh->m_lpNext->m_lpOpposite;
						}
					}
					lplpBorderHalfedges[i] = hh;
				}
				else
				{
					lplpBorderHalfedges[i] = NULL;
				}
				//Store halfedges around the new facet
				if (lpvtTo->m_lpHE != NULL && hh->m_lpOpposite->m_lpVertex == lpvtFrom)//Store an existing edge linking the 2 points
				{
					lplpFacetHalfedges[i] = hh;
				}
				else//Create a pair of halfedges linking the 2 points.
				{
					newHalfedge.m_lpVertex = lpvtTo;
					m_lstHalfedges.push_back(newHalfedge);
					newHalfedge.m_lpVertex = lpvtFrom;
					m_lstHalfedges.push_back(newHalfedge);
					hi1 = m_lstHalfedges.end();
					--hi1;
					hi0 = hi1;
					--hi0;
					hi0->m_lpOpposite = &(*hi1);
					hi1->m_lpOpposite = &(*hi0);
					if (lpvtTo->m_lpHE == NULL)
					{
						lpvtTo->m_lpHE = &(*hi0);
					}
					lplpFacetHalfedges[i] = &(*hi0);
				}
				lplpFacetHalfedges[i]->m_lpFacet = lpFacet;
			}
			for (i = 0; i < nFacetSize; ++i)
			{
				j = (i + 1) % nFacetSize;
				//9 cases to construct relations between new halfedges
				if (lplpBorderHalfedges[i] == NULL)
				{
					if (lplpBorderHalfedges[j] == NULL)
					{
						lplpFacetHalfedges[i]->m_lpOpposite->m_lpPrev = lplpFacetHalfedges[j]->m_lpOpposite;
						lplpFacetHalfedges[j]->m_lpOpposite->m_lpNext = lplpFacetHalfedges[i]->m_lpOpposite;
						lplpFacetHalfedges[i]->m_lpNext = lplpFacetHalfedges[j];
						lplpFacetHalfedges[j]->m_lpPrev = lplpFacetHalfedges[i];
					}
					else if (lplpBorderHalfedges[j] == lplpFacetHalfedges[j])
					{
						lplpFacetHalfedges[i]->m_lpOpposite->m_lpPrev = lplpFacetHalfedges[j]->m_lpPrev;
						lplpFacetHalfedges[j]->m_lpPrev->m_lpNext = lplpFacetHalfedges[i]->m_lpOpposite;
						lplpFacetHalfedges[i]->m_lpNext = lplpFacetHalfedges[j];
						lplpFacetHalfedges[j]->m_lpPrev = lplpFacetHalfedges[i];
					}
					else
					{
						lplpFacetHalfedges[i]->m_lpOpposite->m_lpPrev = lplpFacetHalfedges[j]->m_lpOpposite;
						lplpFacetHalfedges[j]->m_lpOpposite->m_lpNext = lplpFacetHalfedges[i]->m_lpOpposite;
						lplpFacetHalfedges[i]->m_lpNext = lplpFacetHalfedges[j];
						lplpFacetHalfedges[j]->m_lpPrev = lplpFacetHalfedges[i];
					}
				}
				else if (lplpBorderHalfedges[i] == lplpFacetHalfedges[i])
				{
					if (lplpBorderHalfedges[j] == NULL)
					{
						lplpFacetHalfedges[j]->m_lpOpposite->m_lpNext = lplpFacetHalfedges[i]->m_lpNext;
						lplpFacetHalfedges[i]->m_lpNext->m_lpPrev = lplpFacetHalfedges[j]->m_lpOpposite;
						lplpFacetHalfedges[j]->m_lpPrev = lplpFacetHalfedges[i];
						lplpFacetHalfedges[i]->m_lpNext = lplpFacetHalfedges[j];
					}
					else if (lplpBorderHalfedges[j] == lplpFacetHalfedges[j])
					{
						if (lplpFacetHalfedges[i]->m_lpNext != lplpFacetHalfedges[j])
						{
							hh0 = lplpFacetHalfedges[i]->opposite()->m_lpPrev;
							while (hh0->m_lpFacet != NULL)
							{
								hh0 = hh0->opposite()->prev();
							}
							hh = hh0->m_lpNext;
							hh0->m_lpNext = lplpFacetHalfedges[i]->m_lpNext;
							lplpFacetHalfedges[i]->m_lpNext->m_lpPrev = hh0;
							hh->m_lpPrev = lplpFacetHalfedges[j]->m_lpPrev;
							lplpFacetHalfedges[j]->m_lpPrev->m_lpNext = hh;
							lplpFacetHalfedges[j]->m_lpPrev = lplpFacetHalfedges[i];
							lplpFacetHalfedges[i]->m_lpNext = lplpFacetHalfedges[j];
						}
					}
					else
					{
						lplpFacetHalfedges[j]->m_lpOpposite->m_lpNext = lplpFacetHalfedges[i]->m_lpNext;
						lplpFacetHalfedges[i]->m_lpNext->m_lpPrev = lplpFacetHalfedges[j]->m_lpOpposite;
						lplpFacetHalfedges[j]->m_lpPrev = lplpFacetHalfedges[i];
						lplpFacetHalfedges[i]->m_lpNext = lplpFacetHalfedges[j];
					}
				}
				else
				{
					if (lplpBorderHalfedges[j] == NULL)
					{
						lplpFacetHalfedges[j]->m_lpOpposite->m_lpNext = lplpBorderHalfedges[i]->m_lpNext;
						lplpBorderHalfedges[i]->m_lpNext->m_lpPrev = lplpFacetHalfedges[j]->m_lpOpposite;
						lplpFacetHalfedges[i]->m_lpOpposite->m_lpPrev = lplpBorderHalfedges[i];
						lplpBorderHalfedges[i]->m_lpNext = lplpFacetHalfedges[i]->m_lpOpposite;
						lplpFacetHalfedges[i]->m_lpNext = lplpFacetHalfedges[j];
						lplpFacetHalfedges[j]->m_lpPrev = lplpFacetHalfedges[i];
					}
					else if (lplpBorderHalfedges[j] == lplpFacetHalfedges[j])
					{
						lplpFacetHalfedges[j]->m_lpPrev->m_lpNext = lplpFacetHalfedges[i]->m_lpOpposite;
						lplpFacetHalfedges[i]->m_lpOpposite->m_lpPrev = lplpFacetHalfedges[j]->m_lpPrev;
						lplpFacetHalfedges[i]->m_lpNext = lplpFacetHalfedges[j];
						lplpFacetHalfedges[j]->m_lpPrev = lplpFacetHalfedges[i];
					}
					else
					{
						lplpFacetHalfedges[j]->m_lpOpposite->m_lpNext = lplpBorderHalfedges[i]->m_lpNext;
						lplpBorderHalfedges[i]->m_lpNext->m_lpPrev = lplpFacetHalfedges[j]->m_lpOpposite;
						lplpFacetHalfedges[i]->m_lpOpposite->m_lpPrev = lplpBorderHalfedges[i];
						lplpBorderHalfedges[i]->m_lpNext = lplpFacetHalfedges[i]->m_lpOpposite;
						lplpFacetHalfedges[i]->m_lpNext = lplpFacetHalfedges[j];
						lplpFacetHalfedges[j]->m_lpPrev = lplpFacetHalfedges[i];
					}
				}
			}
			fi->m_lpHE = lplpFacetHalfedges[nFacetSize - 1];
			while (fi->m_lpHE->m_lpVertex != lplpFacetVertices[0])
			{
				fi->m_lpHE = fi->m_lpHE->m_lpNext;
			}
			delete[]lplpBorderHalfedges;
			delete[]lplpFacetHalfedges;
		}
	}

	template <typename Kernel, typename PolyherdonItems>
	void Polyhedron<Kernel, PolyherdonItems>::append_facet(unsigned int *lpiIndex, unsigned int nFacetSize)
	{
		unsigned int i;
		HalfedgeDS_vertex_base **lplpFacetVertices;
		if (m_lplpVertexOrder != NULL && nFacetSize > 0)
		{
			lplpFacetVertices = new HalfedgeDS_vertex_base*[nFacetSize];
			for (i = 0; i < nFacetSize; ++i)
			{
				lplpFacetVertices[i] = m_lplpVertexOrder[lpiIndex[i]];
			}
			append_facet(lplpFacetVertices, nFacetSize);
			delete[]lplpFacetVertices;
		}
	}

	template <typename Kernel, typename PolyherdonItems>
	void Polyhedron<Kernel, PolyherdonItems>::append_facet(unsigned int nFacetSize, ...)
	{
		unsigned int i;
		va_list arg_ptr;
		HalfedgeDS_vertex_base **lplpFacetVertices;
		if (m_lplpVertexOrder != NULL && nFacetSize > 0)
		{
			lplpFacetVertices = new HalfedgeDS_vertex_base*[nFacetSize];
			va_start(arg_ptr, nFacetSize);
			for (i = 0; i < nFacetSize; ++i)
			{
				lplpFacetVertices[i] = m_lplpVertexOrder[va_arg(arg_ptr, unsigned int)];
			}
			va_end(arg_ptr);
			append_facet(lplpFacetVertices, nFacetSize);
			delete[]lplpFacetVertices;
		}
	}

	template <typename Kernel, typename PolyhedronItems>
	void Polyhedron<Kernel, PolyhedronItems>::refresh_degree()
	{
		Vertex_iterator vi;
		Halfedge_handle hh0, hh;
		if (m_lstVertices.size() > 0)
		{
			for (vi = m_lstVertices.begin(); vi != m_lstVertices.end(); ++vi)
			{
				vi->m_nDegree = 0;
				hh0 = vi->m_lpHE;
				if (hh0 != NULL)
				{
					hh = hh0;
					do
					{
						++vi->m_nDegree;
						hh = hh->m_lpNext->m_lpOpposite;
					} while (hh != hh0);
				}
				else
				{
					vi = m_lstVertices.erase(vi);
					--vi;
				}
			}
		}
	}


	template <typename Kernel, typename PolyhedronItems>
	typename Polyhedron<Kernel, PolyhedronItems>::Halfedge_handle Polyhedron<Kernel, PolyhedronItems>::flip_edge(Halfedge_handle hh)
	{
		Halfedge_handle hh0, hh1, hh2, hh3, hh4, hh5;
		hh4 = NULL;
		if (!hh->is_border_edge() && hh->prev() == hh->next()->next() && hh->opposite()->prev() == hh->opposite()->next()->next())
		{
			hh0 = hh->next();
			hh1 = hh->prev();
			hh2 = hh->opposite()->next();
			hh3 = hh->opposite()->prev();
			hh4 = hh;
			hh5 = hh->opposite();
			hh0->m_lpPrev = hh3;
			hh0->m_lpNext = hh5;
			hh1->m_lpPrev = hh4;
			hh1->m_lpNext = hh2;
			hh2->m_lpPrev = hh1;
			hh2->m_lpNext = hh4;
			hh3->m_lpPrev = hh5;
			hh3->m_lpNext = hh0;
			hh4->m_lpPrev = hh2;
			hh4->m_lpNext = hh1;
			hh5->m_lpPrev = hh0;
			hh5->m_lpNext = hh3;
			hh4->m_lpFacet->m_lpHE = hh4;
			hh5->facet()->m_lpHE = hh5;
			hh4->vertex()->m_lpHE = hh3;
			hh5->vertex()->m_lpHE = hh1;
			hh4->m_lpVertex = hh0->vertex();
			hh5->m_lpVertex = hh2->vertex();
			hh0->m_lpFacet = hh5->m_lpFacet;
			hh2->m_lpFacet = hh4->m_lpFacet;
			++(hh0->m_lpVertex->m_nDegree);
			--(hh1->m_lpVertex->m_nDegree);
			++(hh2->m_lpVertex->m_nDegree);
			--(hh3->m_lpVertex->m_nDegree);
		}
		return hh4;
	}

	template <typename Kernel, typename PolyhedronItems>
	typename Polyhedron<Kernel, PolyhedronItems>::Halfedge_handle Polyhedron<Kernel, PolyhedronItems>::split_edge(Halfedge_handle hh)
	{
		HalfedgeDS_vertex_base newVertex;
		HalfedgeDS_halfedge_base newHalfedge;
		Vertex_handle vh;
		Halfedge_handle hh0, hh1;
		Halfedge_iterator hi;
		newVertex.m_lpHE = NULL;
		newVertex.ptPosition = typename Kernel::Point_3((hh->m_lpVertex->ptPosition.x() + hh->m_lpPrev->m_lpVertex->ptPosition.x()) / 2, (hh->m_lpVertex->ptPosition.y() + hh->m_lpPrev->m_lpVertex->ptPosition.y()) / 2, (hh->m_lpVertex->ptPosition.z() + hh->m_lpPrev->m_lpVertex->ptPosition.z()) / 2);
		m_lstVertices.push_back(newVertex);
		vh = &(*(m_lstVertices.rbegin()));
		newHalfedge.m_lpOpposite = NULL;
		newHalfedge.m_lpNext = NULL;
		newHalfedge.m_lpPrev = NULL;
		newHalfedge.m_lpVertex = NULL;
		newHalfedge.m_lpFacet = NULL;
		m_lstHalfedges.push_back(newHalfedge);
		hh0 = &(*(m_lstHalfedges.rbegin()));
		m_lstHalfedges.push_back(newHalfedge);
		hh1 = &(*(m_lstHalfedges.rbegin()));
		hh0->m_lpOpposite = hh1;
		hh1->m_lpOpposite = hh0;
		hh->m_lpOpposite->m_lpNext->m_lpPrev = hh1;
		hh1->m_lpNext = hh->m_lpOpposite->m_lpNext;
		hh1->m_lpPrev = hh->m_lpOpposite;
		hh->m_lpOpposite->m_lpNext = hh1;
		hh0->m_lpPrev = hh->m_lpPrev;
		hh->m_lpPrev->m_lpNext = hh0;
		hh->m_lpPrev = hh0;
		hh0->m_lpNext = hh;
		hh0->m_lpVertex = vh;
		hh->m_lpOpposite->m_lpVertex = vh;
		hh1->m_lpVertex = hh0->m_lpPrev->m_lpVertex;
		hh1->m_lpVertex->m_lpHE = hh1;
		hh0->m_lpVertex->m_lpHE = hh0;
		vh->m_lpHE = hh0;
		hh1->m_lpFacet = hh->opposite()->m_lpFacet;
		hh0->m_lpFacet = hh->m_lpFacet;
		if (hh0->m_lpFacet != NULL)
		{
			++(hh0->m_lpFacet->m_nSize);
		}
		if (hh1->m_lpFacet != NULL)
		{
			++(hh1->m_lpFacet->m_nSize);
		}
		vh->m_nDegree = 2;
		return hh0;
	}

	template <typename Kernel, typename PolyhedronItems>
	typename Polyhedron<Kernel, PolyhedronItems>::Halfedge_handle Polyhedron<Kernel, PolyhedronItems>::split_facet(Halfedge_handle hh0, Halfedge_handle hh1)
	{
		HalfedgeDS_facet_base newFacet;
		HalfedgeDS_halfedge_base newHalfedge;
		Facet_handle fh0, fh1;
		Halfedge_handle hh2, hh3, hh;
		if (hh0->m_lpFacet == hh1->m_lpFacet && hh0->m_lpFacet != NULL && hh0->m_lpNext != hh1 && hh1->m_lpNext != hh0)
		{
			newFacet.m_lpHE = NULL;
			m_lstFacets.push_back(newFacet);
			fh0 = hh0->m_lpFacet;
			fh1 = &(*(m_lstFacets.rbegin()));
			newHalfedge.m_lpOpposite = NULL;
			newHalfedge.m_lpPrev = NULL;
			newHalfedge.m_lpNext = NULL;
			newHalfedge.m_lpFacet = NULL;
			newHalfedge.m_lpVertex = NULL;
			m_lstHalfedges.push_back(newHalfedge);
			hh2 = &(*(m_lstHalfedges.rbegin()));
			m_lstHalfedges.push_back(newHalfedge);
			hh3 = &(*(m_lstHalfedges.rbegin()));
			hh2->m_lpOpposite = hh3;
			hh3->m_lpOpposite = hh2;
			hh0->m_lpNext->m_lpPrev = hh3;
			hh3->m_lpNext = hh0->m_lpNext;
			hh0->m_lpNext = hh2;
			hh2->m_lpPrev = hh0;
			hh1->m_lpNext->m_lpPrev = hh2;
			hh2->m_lpNext = hh1->m_lpNext;
			hh1->m_lpNext = hh3;
			hh3->m_lpPrev = hh1;
			hh2->m_lpVertex = hh1->m_lpVertex;
			++(hh1->m_lpVertex->m_nDegree);
			hh3->m_lpVertex = hh0->m_lpVertex;
			++(hh0->m_lpVertex->m_nDegree);
			hh2->m_lpFacet = fh0;
			fh0->m_lpHE = hh2;
			fh1->m_nSize = 0;
			hh = hh1;
			do
			{
				++(fh1->m_nSize);
				hh->m_lpFacet = fh1;
				hh = hh->m_lpNext;
			} while (hh != hh1);
			fh1->m_lpHE = hh3;
			fh0->m_nSize = fh0->m_nSize + 2 - fh1->m_nSize;
		}
		return hh2;
	}

	template <typename Kernel, typename PolyhedronItems>
	typename Polyhedron<Kernel, PolyhedronItems>::Halfedge_handle Polyhedron<Kernel, PolyhedronItems>::split_loop(Halfedge_handle hh0, Halfedge_handle hh1, Halfedge_handle hh2)
	{
		HalfedgeDS_facet_base newFacet;
		HalfedgeDS_halfedge_base newHalfedge;
		HalfedgeDS_vertex_base newVertex;
		Facet_handle lpfh[2];
		Halfedge_handle lphh[12];
		Vertex_handle lpvh[6];
		int i;
		if (hh0->m_lpVertex == hh1->m_lpPrev->m_lpVertex && hh1->m_lpVertex == hh2->m_lpPrev->m_lpVertex && hh2->m_lpVertex == hh0->m_lpPrev->m_lpVertex)
		{
			newFacet.m_lpHE = NULL;
			for (i = 0; i < 2; ++i)
			{
				m_lstFacets.push_back(newFacet);
				lpfh[i] = &(*(m_lstFacets.rbegin()));
			}
			lphh[0] = hh0;
			lphh[2] = hh1;
			lphh[4] = hh2;
			for (i = 0; i < 3; ++i)
			{
				lphh[i * 2 + 1] = lphh[i * 2]->m_lpOpposite;
			}
			newHalfedge.m_lpOpposite = NULL;
			newHalfedge.m_lpPrev = NULL;
			newHalfedge.m_lpNext = NULL;
			newHalfedge.m_lpFacet = NULL;
			newHalfedge.m_lpVertex = NULL;
			for (i = 6; i < 12; ++i)
			{
				m_lstHalfedges.push_back(newHalfedge);
				lphh[i] = &(*(m_lstHalfedges.rbegin()));
			}
			for (i = 0; i < 3; ++i)
			{
				lpvh[i] = lphh[i * 2]->m_lpVertex;
			}
			newVertex.m_lpHE = NULL;
			for (i = 3; i < 6; ++i)
			{
				newVertex.ptPosition = lphh[(i - 3) * 2 + 1]->m_lpVertex->ptPosition;
				m_lstVertices.push_back(newVertex);
				lpvh[i] = &(*(m_lstVertices).rbegin());
			}
			for (i = 0; i < 6; i += 2)
			{
				lphh[i]->m_lpPrev->m_lpNext = lphh[i + 6];
				lphh[i]->m_lpNext->m_lpPrev = lphh[i + 6];
				lphh[i]->m_lpFacet->m_lpHE = lphh[i + 6];
				lphh[i]->m_lpVertex->m_lpHE = lphh[i + 6];
				lphh[i + 6]->m_lpPrev = lphh[i]->m_lpPrev;
				lphh[i + 6]->m_lpNext = lphh[i]->m_lpNext;
				lphh[i + 6]->m_lpOpposite = lphh[i + 7];
				lphh[i + 6]->m_lpFacet = lphh[i]->m_lpFacet;
				lphh[i + 6]->m_lpVertex = lpvh[i];
			}
			for (i = 1; i < 6; i += 2)
			{
				lphh[i + 6]->m_lpPrev = lphh[(i + 2) % 6 + 6];
				lphh[i + 6]->m_lpNext = lphh[(i + 4) % 6 + 6];
				lphh[i + 6]->m_lpOpposite = lphh[i + 5];
				lphh[i + 6]->m_lpFacet = lpfh[1];
				lphh[i + 6]->m_lpVertex = lpvh[(i / 2 + 2) % 3];
			}
			for (i = 1; i < 6; i += 2)
			{
				lphh[i]->m_lpVertex = lpvh[i / 2 + 3];
				lpvh[i + 3]->m_lpHE = lphh[i * 2 + 1];
			}
			for (i = 0; i < 6; i += 2)
			{
				lphh[i]->m_lpPrev = lphh[(i + 4) % 6];
				lphh[i]->m_lpNext = lphh[(i + 2) % 6];
				lphh[i]->m_lpFacet = lpfh[0];
				lphh[i]->m_lpVertex = lpvh[(i / 2 + 1) % 3 + 3];
			}
			lpfh[0]->m_lpHE = lphh[0];
			lpfh[1]->m_lpHE = lphh[7];
		}
		return lphh[6];
	}
}

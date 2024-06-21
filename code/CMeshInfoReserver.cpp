#include "CMeshInfoReserver.h"

void CMeshInfoReserver::m_fnReadMesh(std::ifstream &ifs)
{
	CMeshBase::CPolyhedron::Halfedge_iterator hi;
	ifs >> m_mshSurface;
}

void CMeshInfoReserver::m_fnWriteMesh(std::ofstream &ofs)
{
	ofs.precision(16);
	ofs << m_mshSurface;
}

int CMeshInfoReserver::m_fnImportParaObj(std::ifstream &ifs)
{
	std::string line, comp_index, key_word;
	std::istringstream line_stream;
	double component;
	int i, ind, nVertices, nFacets, nParas, iContainsPara;
	std::vector<double> coordinate, parameter;
	std::vector<int> facet_vertex, facet_para;
	size_t slash;
	CMeshBase::CPolyhedron::Facet_iterator fi;
	CMeshBase::CPolyhedron::Halfedge_iterator hi;
	CMeshBase::CPolyhedron::Halfedge_handle hh0, hh;
	nVertices = 0;
	nFacets = 0;
	nParas = 0;
	while (!ifs.eof())
	{
		std::getline(ifs, line);
		line_stream.str(line);
		line_stream >> key_word;
		if (strcmp(key_word.c_str(), "v") == 0)
		{
			for (i = 0; i < 3; ++i)
			{
				line_stream >> component;
				coordinate.push_back(component);
			}
			++nVertices;
		}
		else if (strcmp(key_word.c_str(), "vt") == 0)
		{
			for (i = 0; i < 2; ++i)
			{
				line_stream >> component;
				parameter.push_back(component);
			}
			++nParas;
		}
		else if (strcmp(key_word.c_str(), "f") == 0)
		{
			for (i = 0; i < 3; ++i)
			{
				line_stream >> comp_index;
				slash = comp_index.find('/');
				ind = atoi(comp_index.substr(0, slash).c_str());
				facet_vertex.push_back(ind - 1);
				if (slash != std::string::npos)
				{
					comp_index = comp_index.substr(slash + 1);
					slash = comp_index.find('/');
					ind = atoi(comp_index.substr(0, slash).c_str());
					facet_para.push_back(ind - 1);
				}
			}
			++nFacets;
		}
		line_stream.clear();
		key_word.clear();
	}
	for (i = 0; i < nVertices; ++i)
	{
		m_mshSurface.append_vertex(coordinate[i * 3], coordinate[i * 3 + 1], coordinate[i * 3 + 2]);
	}
	m_mshSurface.begin_facets();
	for (i = 0; i < nFacets; ++i)
	{
		m_mshSurface.append_facet(3, facet_vertex[i * 3], facet_vertex[i * 3 + 1], facet_vertex[i * 3 + 2]);
	}
	m_mshSurface.end_facets();
	m_mshSurface.refresh_degree();
	if (parameter.size() > 0 && facet_para.size() == facet_vertex.size())
	{
		if (m_lpftGlobalPara != NULL)
		{
			delete[]m_lpftGlobalPara;
		}
		m_lpftGlobalPara = new CMeshBase::FT[parameter.size()];
		for (i = 0; i < parameter.size(); ++i)
		{
			m_lpftGlobalPara[i] = parameter[i];
		}
		i = 0;
		for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
		{
			hh0 = fi->halfedge();
			hh = hh0;
			do
			{
				hh->lpftGlobalPara = m_lpftGlobalPara + facet_para[i] * 2;
				++i;
				hh = hh->next();
			} while (hh != hh0);
		}
		for (hi = m_mshSurface.halfedges_begin(); hi != m_mshSurface.halfedges_end(); ++hi)
		{
			if (hi->is_border())
			{
				hi->lpftGlobalPara = hi->opposite()->prev()->lpftGlobalPara;
			}
		}
		iContainsPara = 1;
	}
	else
	{
		iContainsPara = 0;
	}
	return iContainsPara;
}

void CMeshInfoReserver::m_fnWriteNewMesh(std::ofstream &ofs)
{
	ofs.precision(16);
	ofs << m_mshRemeshedSurface;
}


void CMeshInfoReserver::m_fnReadNewMesh(std::ifstream &ifs)
{
	m_mshRemeshedSurface.clear();
	ifs >> m_mshRemeshedSurface;
}

void CMeshInfoReserver::m_fnExportFixedDirection(std::ofstream &ofs)
{
	CMeshBase::CPolyhedron::Facet_iterator fi;
	int nConstraints, iIndex;
	CMeshBase::Vector_3 vtDir;
	nConstraints = 0;
	iIndex = 0;
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
	{
		fi->iIndex = iIndex;
		if (fi->bDirFixed)
		{
			++nConstraints;
		}
		++iIndex;
	}
	ofs << nConstraints << '\n';
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
	{
		if (fi->bDirFixed)
		{
			vtDir = fi->vtPrincipalAxis * cos(fi->ftChartDir) + fi->vtAuxiliaryAxis * sin(fi->ftChartDir);
			ofs << fi->iIndex << '\t' << vtDir.x() << '\t' << vtDir.y() << '\t' << vtDir.z() << "\t\n";
		}
	}
}

void CMeshInfoReserver::m_fnExportHardFeature(std::ofstream &ofs)
{
	CMeshBase::CPolyhedron::Facet_iterator fi;
	CMeshBase::CPolyhedron::Halfedge_handle hh0, hh;
	int iNumFeatures, iFacetInd;
	iNumFeatures = 0;
	iFacetInd = 0;
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
	{
		hh0 = fi->halfedge();
		hh = hh0;
		do
		{
			if (hh->opposite()->is_border())
			{
				++iNumFeatures;
			}
			hh = hh->next();
		} while (hh != hh0);
		fi->iIndex = iFacetInd;
		++iFacetInd;
	}
	ofs << iNumFeatures << '\n';
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
	{
		hh0 = fi->halfedge();
		hh = hh0;
		do
		{
			if (hh->opposite()->is_border())
			{
				ofs << fi->iIndex << '\t' << hh->prev()->iLocalInd << "\t\n";
			}
			hh = hh->next();
		} while (hh != hh0);
		fi->iIndex = iFacetInd;
	}
}

void CMeshInfoReserver::m_fnExportFixedMark(std::ofstream &ofs)
{
	CMeshBase::CPolyhedron::Facet_iterator fi;
	CMeshBase::Vector_3 vtDir;
	CMeshBase::CPolyhedron::Halfedge_handle hh0, hh;
	bool bInnerFacet;
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
	{
		if (fi->bDirFixed)
		{
			bInnerFacet = true;
			hh0 = fi->halfedge();
			hh = hh0;
			do
			{
				if (hh->opposite()->is_border())
				{
					bInnerFacet = false;
				}
				hh = hh->next();
			} while (hh != hh0);
			if (bInnerFacet)
			{
				ofs << fi->ptIncent.x() << '\t' << fi->ptIncent.y() << '\t' << fi->ptIncent.z() << '\t';
				vtDir = fi->vtPrincipalAxis * cos(fi->ftChartDir) + fi->vtAuxiliaryAxis * sin(fi->ftChartDir);
				ofs << vtDir.x() << '\t' << vtDir.y() << '\t' << vtDir.z() << "\t";
				vtDir = -fi->vtPrincipalAxis * sin(fi->ftChartDir) + fi->vtAuxiliaryAxis * cos(fi->ftChartDir);
				ofs << vtDir.x() << '\t' << vtDir.y() << '\t' << vtDir.z() << "\t\n";
			}
		}
	}
}
//{
//	CMeshBase::CPolyhedron::Vertex_iterator vi;
//	CMeshBase::CPolyhedron::Halfedge_handle hh0, hh;
//	CMeshBase::Vector_3 lpvtDir[2];
//	int iNumSum;
//	for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
//	{
//		iNumSum = 0;
//		lpvtDir[0] = CMeshBase::Vector_3(0, 0, 0);
//		lpvtDir[1] = CMeshBase::Vector_3(0, 0, 0);
//		hh0 = vi->halfedge();
//		hh = hh0;
//		do
//		{
//			if (hh->facet()->bDirFixed)
//			{
//				lpvtDir[0] = lpvtDir[0] + (hh->facet()->vtPrincipalAxis) * cos(hh->facet()->ftChartDir) + (hh->facet()->vtAuxiliaryAxis) * sin(hh->facet()->ftChartDir);
//				lpvtDir[1] = lpvtDir[1] - (hh->facet()->vtPrincipalAxis) * sin(hh->facet()->ftChartDir) + (hh->facet()->vtAuxiliaryAxis) * cos(hh->facet()->ftChartDir);
//				++iNumSum;
//			}
//			hh = hh->next()->opposite();
//		} while (hh != hh0);
//		if (iNumSum == 4)
//		{
//			ofs << vi->point().x() << '\t' << vi->point().y() << '\t' << vi->point().z() << '\t';
//			ofs << lpvtDir[0].x() << '\t' << lpvtDir[0].y() << '\t' << lpvtDir[0].z() << '\t';
//			ofs << lpvtDir[1].x() << '\t' << lpvtDir[1].y() << '\t' << lpvtDir[1].z() << "\t\n";
//		}
//	}
//}

void CMeshInfoReserver::m_fnImportFixedDirection(std::ifstream &ifs)
{
	int iIndex, i, n;
	CMeshBase::CPolyhedron::Facet_iterator fi;
	CMeshBase::CPolyhedron::Facet_handle *lpfhBuffer;
	CMeshBase::FT ftDir, ftX, ftY, ftZ;
	CMeshBase::Vector_3 vtDir;
	lpfhBuffer = new CMeshBase::CPolyhedron::Facet_handle[m_mshSurface.size_of_facets()];
	i = 0;
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
	{
		lpfhBuffer[i] = &(*fi);
		++i;
	}
	//m_lstFixedFacets.clear();
	//m_lstFixedDirections.clear();
	ifs >> n;
	for (i = 0; i < n; ++i)
	{
		ifs >> iIndex >> ftX >> ftY >> ftZ;
		vtDir = CMeshBase::Vector_3(ftX, ftY, ftZ);
		ftDir = atan2(lpfhBuffer[iIndex]->vtAuxiliaryAxis * vtDir, lpfhBuffer[iIndex]->vtPrincipalAxis * vtDir);
		lpfhBuffer[iIndex]->ftChartDir = ftDir;
		lpfhBuffer[iIndex]->bDirFixed = true;
		//m_lstFixedFacets.push_back(lpfhBuffer[iIndex]);
		//m_lstFixedDirections.push_back(ftDir);
	}
	delete[]lpfhBuffer;
}

void CMeshInfoReserver::m_fnImportDirectionField(std::ifstream &ifs, bool bWithIndex, bool bFrameField)
{
	CMeshBase::CPolyhedron::Facet_iterator fi;
	CMeshBase::FT ftX, ftY, ftZ;
	CMeshBase::Vector_3 lpvtDir[2];
	int iIndex;
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
	{
		if (bWithIndex)
		{
			ifs >> iIndex;
		}
		ifs >> ftX >> ftY >> ftZ;
		lpvtDir[0] = CMeshBase::Vector_3(ftX, ftY, ftZ);
		fi->ftChartDir = atan2(fi->vtAuxiliaryAxis * lpvtDir[0], fi->vtPrincipalAxis * lpvtDir[0]);
		if (bFrameField)
		{
			ifs >> ftX >> ftY >> ftZ;
			lpvtDir[1] = CMeshBase::Vector_3(ftX, ftY, ftZ);
			if (CGAL::cross_product(lpvtDir[0], lpvtDir[1]) * fi->vtNorm < 0)
			{
				fi->ftFrameDir = atan2(-fi->vtAuxiliaryAxis * lpvtDir[1], -fi->vtPrincipalAxis * lpvtDir[1]);
			}
			else
			{
				fi->ftFrameDir = atan2(fi->vtAuxiliaryAxis * lpvtDir[1], fi->vtPrincipalAxis * lpvtDir[1]);
			}
		}
		else
		{
			fi->ftFrameDir = atan2(fi->vtPrincipalAxis * lpvtDir[0], -fi->vtAuxiliaryAxis * lpvtDir[0]);
		}
		//ifs >> ftX >> ftY >> ftZ;
	}
}

void CMeshInfoReserver::m_fnExportDirectionField(std::ofstream &ofs, bool bWithIndex)
{
	CMeshBase::CPolyhedron::Facet_iterator fi;
	CMeshBase::Vector_3 vtDir;
	int i;
	i = 0;
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
	{
		if (bWithIndex)
		{
			ofs << i << '\t';
		}
		vtDir = fi->vtPrincipalAxis * cos(fi->ftChartDir) + fi->vtAuxiliaryAxis * sin(fi->ftChartDir);
		ofs << vtDir.x() << '\t' << vtDir.y() << '\t' << vtDir.z() << "\t\n";
		++i;
	}
}

void CMeshInfoReserver::m_fnExportParaDif(std::ofstream& ofs)
{
	CMeshBase::CPolyhedron::Facet_iterator fi;
	CMeshBase::CPolyhedron::Edge_iterator ei;
	CMeshBase::Vector_3 lpvtDirection[4];
	CMeshBase::FT lpftLocalPara[4], lpftErr[2];
	CMeshBase::CPolyhedron::Facet_handle fh;
	int iSeamType, iInd;
	iInd = 0;
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
	{
		fi->iIndex = iInd;
		++iInd;
	}
	for (ei = m_mshSurface.edges_begin(); ei != m_mshSurface.edges_end(); ++ei)
	{
		if (!ei->is_border_edge())
		{
			fh = ei->facet();
			lpvtDirection[0] = (cos(fh->ftChartDir) * fh->vtPrincipalAxis + sin(fh->ftChartDir) * fh->vtAuxiliaryAxis) / fh->ftRadius;
			lpvtDirection[1] = ((-sin(fh->ftChartDir)) * fh->vtPrincipalAxis + cos(fh->ftChartDir) * fh->vtAuxiliaryAxis) / fh->ftRadius;
			lpftLocalPara[0] = fh->ftLocalParaScl * (lpvtDirection[0] * ei->vtVector);
			lpftLocalPara[1] = fh->ftLocalParaScl * (lpvtDirection[1] * ei->vtVector);
			fh = ei->opposite()->facet();
			iSeamType = ei->iSeamType;
			lpvtDirection[2] = (cos(fh->ftChartDir - iSeamType * ftHalfPi) * fh->vtPrincipalAxis + sin(fh->ftChartDir - iSeamType * ftHalfPi) * fh->vtAuxiliaryAxis) / fh->ftRadius;
			lpvtDirection[3] = ((-sin(fh->ftChartDir - iSeamType * ftHalfPi)) * fh->vtPrincipalAxis + cos(fh->ftChartDir - iSeamType * ftHalfPi) * fh->vtAuxiliaryAxis) / fh->ftRadius;
			lpftLocalPara[2] = fh->ftLocalParaScl * (lpvtDirection[2] * ei->vtVector);
			lpftLocalPara[3] = fh->ftLocalParaScl * (lpvtDirection[3] * ei->vtVector);
			lpftErr[0] = fabs((lpftLocalPara[2] - lpftLocalPara[0]) / (lpftLocalPara[2] + lpftLocalPara[0]));
			lpftErr[1] = fabs((lpftLocalPara[3] - lpftLocalPara[1]) / (lpftLocalPara[3] + lpftLocalPara[1]));
			ofs << ei->facet()->iIndex << '\t' << ei->opposite()->facet()->iIndex << '\t' << lpftErr[0] << '\t' << lpftErr[1] << "\t\n";
		}
	}
}

void CMeshInfoReserver::m_fnImportIntCones(std::ifstream& ifs)
{
	int iIndex;
	CMeshBase::CPolyhedron::Vertex_iterator vi;
	CMeshBase::CPolyhedron::Vertex_handle* lpvhBuffer;
	CMeshBase::FT ftConeDefect;
	CMeshBase::CPolyhedron::Vertex_handle vh;
	lpvhBuffer = new CMeshBase::CPolyhedron::Vertex_handle[m_mshSurface.size_of_vertices()];
	iIndex = 0;
	for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
	{
		vi->nNumPorts = 0;
		lpvhBuffer[iIndex] = &(*vi);
		++iIndex;
	}
	while (!ifs.eof())
	{
		iIndex = -1;
		ifs >> iIndex;
		if (iIndex != -1)
		{
			ifs >> ftConeDefect;
			vh = lpvhBuffer[iIndex - 1];
			if (vh->iBordId == -1)
			{
				vh->nNumPorts = 4 - (int(ftConeDefect / ftHalfPi + 15.5) - 15);
			}
			else
			{
				vh->nNumPorts = 2 - (int(ftConeDefect / ftHalfPi + 15.5) - 15);
			}
		}
	}
	delete[]lpvhBuffer;
}

void CMeshInfoReserver::m_fnExportMountain(std::ofstream &ofs, double dblMinClr, double dblMaxClr, double dblHeightScl)
{
	CMeshBase::CPolyhedron::Vertex_iterator vi;
	CMeshBase::CPolyhedron::Halfedge_handle hh0, hh;
	CMeshBase::CPolyhedron::Facet_iterator fi;
	CMeshBase::Vector_3 vtNormal, vtRotGrad;
	int iInd, i;
	iInd = 0;
	for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
	{
		vtNormal = CMeshBase::Vector_3(0, 0, 0);
		hh0 = vi->halfedge();
		hh = hh0;
		do
		{
			if (!hh->is_border())
			{
				vtNormal = vtNormal + hh->facet()->vtNorm;
			}
			hh = hh->next()->opposite();
		} while (hh != hh0);
		vtNormal = vtNormal / sqrt(CGAL::squared_length(vtNormal));
		vi->ptMountainPoint = vi->point() + vtNormal * vi->ftViewDensity * dblHeightScl;
		ofs << "v " << vi->ptMountainPoint.x() << ' ' << vi->ptMountainPoint.y() << ' ' << vi->ptMountainPoint.z() << " \n";
		vi->iTempIndex = iInd;
		++iInd;
	}
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
	{
		vtRotGrad = CMeshBase::Vector_3(0, 0, 0);
		hh0 = fi->halfedge();
		hh = hh0;
		do
		{
			vtRotGrad = vtRotGrad + hh->prev()->vtVector * ((hh->vertex()->ftViewDensity - dblMinClr) / (dblMaxClr - dblMinClr));
			hh = hh->next();
		} while (hh != hh0);
		vtRotGrad = vtRotGrad / fi->ftTwArea;
		hh = hh0;
		do
		{
			ofs << "vt " << vtRotGrad * (hh->vertex()->point() - fi->ptIncent) * 0.5 + 0.5 << ' ' << (hh->vertex()->ftViewDensity - dblMinClr) / (dblMaxClr - dblMinClr) << " \n";
			hh = hh->next();
		} while (hh != hh0);
	}
	iInd = 0;
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
	{
		i = 0;
		ofs << "f ";
		hh0 = fi->halfedge();
		hh = hh0;
		do
		{
			ofs << hh->vertex()->iTempIndex + 1 << '/' << iInd * 3 + i + 1 << ' ';
			++i;
			hh = hh->next();
		} while (hh != hh0);
		ofs << '\n';
		++iInd;
	}
}

void CMeshInfoReserver::m_fnExportQuadSingularities(std::ofstream &ofs)
{
	CMeshBase::CSamplerPolyhedron::Vertex_iterator vi;
	ofs.precision(16);
	for (vi = m_mshRemeshedSurface.vertices_begin(); vi != m_mshRemeshedSurface.vertices_end(); ++vi)
	{
		if (vi->degree() != 4 && vi->iBorderId == -1)
		{
			ofs << vi->point().x() << '\t' << vi->point().y() << '\t' << vi->point().z() << '\t' << vi->degree() <<"\t\n";
		}
	}
}

void CMeshInfoReserver::m_fnExportElectricityField(std::ofstream &ofs)
{
	CMeshBase::CPolyhedron::Facet_iterator fi;
	CMeshBase::Vector_3 vtElectricityField;
	int iIndex;
	iIndex = 0;
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
	{
		vtElectricityField = fi->vtPrincipalAxis * fi->lpftElectricFieldSum[0] + fi->vtAuxiliaryAxis * sin(fi->lpftElectricFieldSum[1]);
		ofs << (iIndex + 1) << ' ' << vtElectricityField.x() << ' ' << vtElectricityField.y() << ' ' << vtElectricityField.z() << '\n';
		++iIndex;
	}
}


void CMeshInfoReserver::m_fnExportCurlField(std::ofstream &ofs)
{
	CMeshBase::CPolyhedron::Facet_iterator fi;
	CMeshBase::Vector_3 vtCurlField;
	int iIndex;
	iIndex = 0;
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
	{
		vtCurlField = fi->vtPrincipalAxis * fi->lpftElectricFieldDif[0] + fi->vtAuxiliaryAxis * sin(fi->lpftElectricFieldDif[1]);
		ofs << (iIndex + 1) << ' ' << vtCurlField.x() << ' ' << vtCurlField.y() << ' ' << vtCurlField.z() << '\n';
		++iIndex;
	}
}

void CMeshInfoReserver::m_fnExportCrossField(std::ofstream &ofs)
{
	CMeshBase::CPolyhedron::Facet_iterator fi;
	CMeshBase::Vector_3 vtCrossField;
	int iIndex;
	iIndex = 0;
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
	{
		vtCrossField = fi->vtPrincipalAxis * cos(fi->ftChartDir) + fi->vtAuxiliaryAxis * sin(fi->ftChartDir);
		ofs << vtCrossField.x() << ' ' << vtCrossField.y() << ' ' << vtCrossField.z() << '\n';
		++iIndex;
	}
}

void CMeshInfoReserver::m_fnExportSingularities(std::ofstream &ofs)
{
	CMeshBase::CPolyhedron::Vertex_iterator vi;
//	CMeshBase::CPolyhedron::Halfedge_handle hh0, hh;
//	int iSumInd;
	for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
	{
		if (vi->nNumPorts != 0)
		{
			//if (vi->iBordId == -1)
			//{
			//	iSumInd = 4 - vi->nNumPorts;
			//}
			//else
			//{
			//	iSumInd = 2 - vi->nNumPorts;
			//}
			//hh0 = vi->halfedge();
			//hh = hh0;
			//do
			//{
			//	if (hh->prev()->vertex()->nNumPorts != 0)
			//	{
			//		if (hh->prev()->vertex()->iBordId == -1)
			//		{
			//			iSumInd += (4 - hh->prev()->vertex()->nNumPorts);
			//		}
			//		else
			//		{
			//			iSumInd += (2 - hh->prev()->vertex()->nNumPorts);
			//		}
			//	}
			//	hh = hh->next()->opposite();
			//} while (hh != hh0);
			//if (iSumInd != 0)
			//{
				ofs << vi->point().x() << '\t' << vi->point().y() << '\t' << vi->point().z() << '\t' << vi->nNumPorts << "\t\n";
			//}
			//ofs << (vi->iIndex + 1) << ' ' << (4 - int(vi->nNumPorts)) << '\n';
		}
	}
}

void CMeshInfoReserver::m_fnExportCrevasse(std::ofstream &ofs)
{
    CMeshBase::CPolyhedron::Vertex_iterator vi;
    CMeshBase::CPolyhedron::Edge_iterator ei;
    int iIndex;
    for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
    {
        vi->iTempIndex = -1;
    }
    iIndex = 0;
    for (ei = m_mshSurface.edges_begin(); ei != m_mshSurface.edges_end(); ++ei)
    {
        if (ei->bFreeEdge && !ei->is_border_edge())
        {
            if (ei->prev()->vertex()->iTempIndex == -1)
            {
                ofs << "v " << ei->prev()->vertex()->point().x() << ' ' << ei->prev()->vertex()->point().y() << ' ' << ei->prev()->vertex()->point().z() << " \n";
                ei->prev()->vertex()->iTempIndex = iIndex;
                ++iIndex;
            }
            if (ei->vertex()->iTempIndex == -1)
            {
                ofs << "v " << ei->vertex()->point().x() << ' ' << ei->vertex()->point().y() << ' ' << ei->vertex()->point().z() << " \n";
                ei->vertex()->iTempIndex = iIndex;
                ++iIndex;
            }
        }
    }
    for (ei = m_mshSurface.edges_begin(); ei != m_mshSurface.edges_end(); ++ei)
    {
        if (ei->bFreeEdge && !ei->is_border_edge())
        {
            ofs << "l " << ei->prev()->vertex()->iTempIndex + 1 << ' ' << ei->vertex()->iTempIndex + 1 << '\n';
        }
    }
}

void CMeshInfoReserver::m_fnWriteProcessStatus(std::ofstream &ofs)
{
    ofs << int(m_psProcessStatus)<<'\n';
	ofs << "4\n";
    if(m_psProcessStatus > PS_EMPTY)
    {
        m_fnWriteOriginalCoordinates(ofs);
        m_fnWriteOriginalFacetLinkage(ofs);
        ofs << "NumberOfHalfEdges: " << (m_mshSurface.size_of_halfedges()) << '\n';
    }
    if(m_psProcessStatus > PS_GEOMETRY_FILLED)
    {
        m_fnWriteChartDir(ofs);
    }
    if(m_psProcessStatus > PS_DIRECTION_SCATTER)
    {
        m_fnWriteDensity(ofs);
    }
	if (m_psProcessStatus > PS_TRACK_PATH)
	{
		m_fnWriteFreeEdges(ofs);
	}
	if (m_psProcessStatus > PS_CREVASSE_GLOBAL)
	{
		m_fnWriteGlobalParameter(ofs);
	}
}

void CMeshInfoReserver::m_fnReadProcessStatus(std::ifstream &ifs)
{
    int iProcessStatus;
    std::vector<double> vecCoordinates;
    std::vector<int> vecIndex;
    char lpstrBuffer[64];
    int nNumHalfEdges, nFieldDegree;
	CDirectionAdjustor DA(m_mshSurface);
	CParameterGenerator PG(m_mshSurface, m_vecCrevasses);
	CParameterAdjustor PA(m_mshSurface, m_vecCrevasses);
	CSampler Sampler(m_mshSurface, m_vecCrevasses, m_mshRemeshedSurface, m_vecSamplePoints, m_vecSampleBlocks, m_vecFacetLinkages);
	CMeshOptimizer Optimizer(m_mshSurface, m_mshRemeshedSurface);
    ifs >> iProcessStatus;
	ifs >> nFieldDegree;
    m_psProcessStatus = PROCESS_STATUS(iProcessStatus);
    if(m_psProcessStatus > PS_EMPTY)
    {
        m_fnReadOriginalCoordinates(ifs, vecCoordinates);
        m_fnReadOriginalFacetLinkage(ifs, vecIndex);
        ifs >> lpstrBuffer;
        if(strcmp(lpstrBuffer, "NumberOfHalfEdges:") == 0)
        {
            ifs >> nNumHalfEdges;
            CTriPolyhedronModifier PolyhedronModifier(vecCoordinates, vecIndex, m_mshSurface);
			PolyhedronModifier.GenerateMesh();
            vecCoordinates.clear();
            vecIndex.clear();
			m_mshSurface.refresh_degree();
        }
    }
    if(m_psProcessStatus > PS_GEOMETRY_BLANK)
    {
        m_fnFillGeometricInfo();
    }
    if(m_psProcessStatus > PS_GEOMETRY_FILLED)
    {
        m_fnReadChartDir(ifs);
        m_fnDetectSeamType();
        m_fnMarkSingularity();
    }
    if(m_psProcessStatus > PS_DIRECTION_SCATTER)
    {
        m_fnReadDensity(ifs);
    }
	if (m_psProcessStatus > PS_TRACK_PATH)
	{
		m_fnReadFreeEdges(ifs);
	}
	if (m_psProcessStatus > PS_DIRECTION_DECOMPOSE)
	{
		PG.m_fnMarkCrevasse();
		PG.m_fnExtractCrevasses();
		m_fnIsAligned();
		PG.m_fnUniformizeChart();
	}
	if (m_psProcessStatus > PS_DIRECTION_UNIFORM)
	{
		DA.m_fnGenerateMoment();
	}
	if (m_psProcessStatus > PS_CREVASSE_GLOBAL)
	{
		m_fnReadGlobalParameter(ifs);
		m_fnRepresentGlobalPara();
	}
	if (m_psProcessStatus > PS_PARAMETER_INT)
	{
		Sampler.m_fnInitSampleBlocks();
	}
	if (m_psProcessStatus > PS_NEW_POINTS)
	{
		Sampler.m_fnGenerateLinkage();
	}
	if (m_psProcessStatus > PS_NEW_LINKAGE)
	{
		Sampler.m_fnGenerateFacets();
	}
	if (m_psProcessStatus > PS_NEW_SCATTER)
	{
		Sampler.m_fnCorrectFacetDirection();
	}
	if (m_psProcessStatus > PS_NEW_UNIFORM)
	{
		Sampler.m_fnExtractNewMesh();
		m_mshRemeshedSurface.refresh_degree();
		m_fnCountNewMeshSize();
		m_fnMarkNewMeshBorder();
		m_fnMarkNewMeshCrevasse();
		Optimizer.m_fnMarkFreedom();
	}
    ifs.close();
}

void CMeshInfoReserver::m_fnWriteOriginalCoordinates(std::ofstream &ofsWriter)
{
    CMeshBase::CPolyhedron::Vertex_iterator vi;
    CDataConvertor DataConvertor;
    ofsWriter << "Coordinates:\n";
    ofsWriter << (m_mshSurface.size_of_vertices()) << '\n';
    for(vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
    {
        DataConvertor.dblData = vi->point().x();
        ofsWriter << DataConvertor.iLoWord << '\t' << DataConvertor.iHiWord << '\t';
        DataConvertor.dblData = vi->point().y();
        ofsWriter << DataConvertor.iLoWord << '\t' << DataConvertor.iHiWord << '\t';
        DataConvertor.dblData = vi->point().z();
        ofsWriter << DataConvertor.iLoWord << '\t' << DataConvertor.iHiWord << '\t';
        ofsWriter << '\n';
    }
}

void CMeshInfoReserver::m_fnReadOriginalCoordinates(std::ifstream &ifsReader, std::vector<double> &vecCoordinates)
{
    int nNumCoordinates, i;
    char lpstrBuffer[64];
    CDataConvertor DataConvertor;
    ifsReader >> lpstrBuffer;
    if(strcmp(lpstrBuffer, "Coordinates:") == 0)
    {
        ifsReader >> nNumCoordinates;
        nNumCoordinates *= 3;
        for(i = 0; i < nNumCoordinates; ++i)
        {
            ifsReader >> DataConvertor.iLoWord >> DataConvertor.iHiWord;
            vecCoordinates.push_back(DataConvertor.dblData);
        }
    }
}

void CMeshInfoReserver::m_fnWriteOriginalFacetLinkage(std::ofstream &ofsWriter)
{
    CMeshBase::CPolyhedron::Facet_iterator fi;
    CMeshBase::CPolyhedron::Halfedge_handle hh0, hh;
    ofsWriter << "FacetLinkages:\n";
    ofsWriter << (m_mshSurface.size_of_facets()) << '\n';
    for(fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
    {
        hh0 = fi->halfedge();
        hh = hh0;
        do
        {
            ofsWriter << hh->vertex()->iIndex << '\t';
            hh = hh->next();
        }while(hh != hh0);
        ofsWriter << '\n';
    }
}

void CMeshInfoReserver::m_fnReadOriginalFacetLinkage(std::ifstream &ifsReader, std::vector<int> &vecIndex)
{
    int nNumIndex, i, iIndex;
    char lpstrBuffer[64];
    ifsReader >> lpstrBuffer;
    if(strcmp(lpstrBuffer, "FacetLinkages:") == 0)
    {
        ifsReader >> nNumIndex;
        nNumIndex *= 3;
        for(i = 0; i < nNumIndex; ++i)
        {
            ifsReader >> iIndex;
            vecIndex.push_back(iIndex);
        }
    }
}


void CMeshInfoReserver::m_fnWriteChartDir(std::ofstream &ofsWriter)
{
    CMeshBase::CPolyhedron::Facet_iterator fi;
    CDataConvertor DataConvertor;
    ofsWriter << "ChartDirection:\n";
    for(fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
    {
        DataConvertor.dblData = fi->ftChartDir;
        ofsWriter << DataConvertor.iLoWord << '\t' << DataConvertor.iHiWord <<'\n';
    }
}

void CMeshInfoReserver::m_fnReadChartDir(std::ifstream &ifsReader)
{
    char lpstrBuffer[64];
    CMeshBase::CPolyhedron::Facet_iterator fi;
    CDataConvertor DataConvertor;
    ifsReader >> lpstrBuffer;
    if(strcmp(lpstrBuffer, "ChartDirection:") == 0)
    {
        for(fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
        {
            ifsReader >> DataConvertor.iLoWord >> DataConvertor.iHiWord;
            fi->ftChartDir = DataConvertor.dblData;
        }
    }
}

void CMeshInfoReserver::m_fnWriteDensity(std::ofstream &ofsWriter)
{
    CMeshBase::CPolyhedron::Vertex_iterator vi;
    CDataConvertor DataConvertor;
    ofsWriter << "Density:\n";
    for(vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
    {
        DataConvertor.dblData = vi->ftDensity;
        ofsWriter << DataConvertor.iLoWord << '\t' << DataConvertor.iHiWord<< '\n';
    }
}

void CMeshInfoReserver::m_fnReadDensity(std::ifstream &ifsReader)
{
    char lpstrBuffer[64];
    CMeshBase::CPolyhedron::Vertex_iterator vi;
    CDataConvertor DataConvertor;
    ifsReader >> lpstrBuffer;
    if(strcmp(lpstrBuffer, "Density:") == 0)
    {
        for(vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
        {
            ifsReader >> DataConvertor.iLoWord >> DataConvertor.iHiWord;
            vi->ftDensity = DataConvertor.dblData;
        }
    }
}

void CMeshInfoReserver::m_fnWriteFreeEdges(std::ofstream &ofsWriter)
{
	int iTempIndex;
	CMeshBase::CPolyhedron::Edge_iterator ei;
	ofsWriter << "FreeEdges:\n";
	iTempIndex = 0;
	for (ei = m_mshSurface.edges_begin(); ei != m_mshSurface.edges_end(); ++ei)
	{
		if (ei->bFreeEdge)
		{
			ofsWriter << iTempIndex << '\n';
		}
		++iTempIndex;
	}
	ofsWriter << -1 << '\n';
}

void CMeshInfoReserver::m_fnReadFreeEdges(std::ifstream &ifsReader)
{
	CMeshBase::CPolyhedron::Halfedge_handle *lphhHalfedgeBuffer, *lphhHalfedgeFiller;
	CMeshBase::CPolyhedron::Edge_iterator ei;
	int iEdgeInd;
	char lpstrBuffer[64];
	ifsReader >> lpstrBuffer;
	if (strcmp(lpstrBuffer, "FreeEdges:") == 0)
	{
		lphhHalfedgeBuffer = new CMeshBase::CPolyhedron::Halfedge_handle[m_mshSurface.size_of_edges()];
		lphhHalfedgeFiller = lphhHalfedgeBuffer;
		for (ei = m_mshSurface.edges_begin(); ei != m_mshSurface.edges_end(); ++ei)
		{
			*lphhHalfedgeFiller = &(*ei);
			++lphhHalfedgeFiller;
		}
		ifsReader >> iEdgeInd;
		while (iEdgeInd != -1)
		{
			lphhHalfedgeBuffer[iEdgeInd]->bFreeEdge = true;
			lphhHalfedgeBuffer[iEdgeInd]->opposite()->bFreeEdge = true;
			ifsReader >> iEdgeInd;
		}
		delete[]lphhHalfedgeBuffer;
	}
}

void CMeshInfoReserver::m_fnWriteGlobalParameter(std::ofstream &ofsWriter)
{
	CDataConvertor DataConvertor;
	int i;
	ofsWriter << "GlobalParameters:\n";
	ofsWriter << m_nGlobalParameterSize << '\n';
	for (i = 0; i < m_nGlobalParameterSize; ++i)
	{
		DataConvertor.dblData = m_lpftGlobalPara[i * 2];
		ofsWriter << DataConvertor.iLoWord << '\t' << DataConvertor.iHiWord << '\t';
		DataConvertor.dblData = m_lpftGlobalPara[i * 2 + 1];
		ofsWriter << DataConvertor.iLoWord << '\t' << DataConvertor.iHiWord << '\n';
	}
}

void CMeshInfoReserver::m_fnReadGlobalParameter(std::ifstream &ifsReader)
{
	char lpstrBuffer[64];
	int i;
	CDataConvertor DataConvertor;
	ifsReader >> lpstrBuffer;
	if (strcmp(lpstrBuffer, "GlobalParameters:") == 0)
	{
		ifsReader >> m_nGlobalParameterSize;
		for (i = 0; i < m_nGlobalParameterSize; ++i)
		{
			ifsReader >> DataConvertor.iLoWord >> DataConvertor.iHiWord;
			m_lpftGlobalPara[i * 2] = DataConvertor.dblData;
			ifsReader >> DataConvertor.iLoWord >> DataConvertor.iHiWord;
			m_lpftGlobalPara[i * 2 + 1] = DataConvertor.dblData;
		}
	}
}

void CMeshInfoReserver::m_fnExportTexturePara(std::ofstream &ofs, bool bNormalizeTexture)
{
	CMeshBase::CPolyhedron::Vertex_iterator vi;
	CMeshBase::CPolyhedron::Halfedge_handle hh0, hh, hhStart, hhEnd;
	CMeshBase::CPolyhedron::Facet_iterator fi;
	int iIndex;
	CMeshBase::FT lpftOffset[2], ftCoef;
	if (bNormalizeTexture)
	{
		lpftOffset[0] = DBL_MAX;
		lpftOffset[1] = DBL_MAX;
		ftCoef = 0;
		for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
		{
			hh0 = fi->halfedge();
			hh = hh0;
			do
			{
				if (hh->lpftGlobalPara[0] < lpftOffset[0])
				{
					lpftOffset[0] = hh->lpftGlobalPara[0];
				}
				if (hh->lpftGlobalPara[1] < lpftOffset[1])
				{
					lpftOffset[1] = hh->lpftGlobalPara[1];
				}
				if (hh->lpftGlobalPara[0] - lpftOffset[0] > ftCoef)
				{
					ftCoef = hh->lpftGlobalPara[0] - lpftOffset[0];
				}
				if (hh->lpftGlobalPara[1] - lpftOffset[1] > ftCoef)
				{
					ftCoef = hh->lpftGlobalPara[1] - lpftOffset[1];
				}
				hh = hh->next();
			} while (hh != hh0);
		}
	}
	iIndex = 0;
	for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
	{
		vi->iTempIndex = iIndex;
		ofs << "v " << vi->point().x() << ' ' << vi->point().y() << ' ' << vi->point().z() << " \n";
		++iIndex;
	}
	iIndex = 0;
	for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
	{
		if (vi->iBordId == -1)
		{
			if (vi->nCrevasseSize < 2)
			{
				hh0 = vi->halfedge();
				hh = hh0;
				do
				{
					hh->iTempIndex = iIndex;
					hh = hh->next()->opposite();
				} while (hh != hh0);
				if (bNormalizeTexture)
				{
					ofs << "vt " << (hh0->lpftGlobalPara[0] - lpftOffset[0]) / ftCoef << ' ' << (hh0->lpftGlobalPara[1] - lpftOffset[1]) / ftCoef << " \n";
				}
				else
				{
					ofs << "vt " << hh0->lpftGlobalPara[0] << ' ' << hh0->lpftGlobalPara[1] << " \n";
				}
				++iIndex;
			}
			else
			{
				hh0 = vi->halfedge()->hhCrevassePrev;
				hhStart = hh0;
				do
				{
					hhEnd = hhStart->hhCrevasseNext->opposite();
					hh = hhStart;
					while (hh != hhEnd)
					{
						hh->iTempIndex = iIndex;
						hh = hh->next()->opposite();
					}
					if (bNormalizeTexture)
					{
						ofs << "vt " << (hhStart->lpftGlobalPara[0] - lpftOffset[0]) / ftCoef << ' ' << (hhStart->lpftGlobalPara[1] - lpftOffset[1]) / ftCoef << " \n";
					}
					else
					{
						ofs << "vt " << hhStart->lpftGlobalPara[0] << ' ' << hhStart->lpftGlobalPara[1] << " \n";
					}
					hhStart = hhEnd;
					++iIndex;
				} while (hhStart != hh0);
			}
		}
		else
		{
			hh0 = vi->halfedge();
			while (!hh0->is_border())
			{
				hh0 = hh0->next()->opposite();
			}
			hhStart = hh0->next()->opposite();
			do
			{
				hhEnd = hhStart->hhCrevasseNext->opposite();
				hh = hhStart;
				while (hh != hhEnd)
				{
					hh->iTempIndex = iIndex;
					hh = hh->next()->opposite();
				}
				if (bNormalizeTexture)
				{
					ofs << "vt " << (hhStart->lpftGlobalPara[0] - lpftOffset[0]) / ftCoef << ' ' << (hhStart->lpftGlobalPara[1] - lpftOffset[1]) / ftCoef << " \n";
				}
				else
				{
					ofs << "vt " << hhStart->lpftGlobalPara[0] << ' ' << hhStart->lpftGlobalPara[1] << " \n";
				}
				hhStart = hhEnd;
				++iIndex;
			} while (hhStart != hh0);
		}
	}
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
	{
		ofs << "f ";
		hh0 = fi->halfedge();
		hh = hh0;
		do
		{
			ofs << hh->vertex()->iTempIndex + 1 << '/' << hh->iTempIndex + 1 << ' ';
			hh = hh->next();
		} while (hh != hh0);
		ofs << '\n';
	}
}

void CMeshInfoReserver::m_fnExportSelected(std::ofstream &ofs, bool bNormalizeTexture)
{
	CMeshBase::CPolyhedron::Vertex_iterator vi;
	CMeshBase::CPolyhedron::Halfedge_handle hh0, hh, hhStart, hhEnd;
	CMeshBase::CPolyhedron::Facet_iterator fi;
	bool bSelected;
	int iIndex;
	CMeshBase::FT lpftOffset[2], ftCoef;
	if (bNormalizeTexture)
	{
		lpftOffset[0] = DBL_MAX;
		lpftOffset[1] = DBL_MAX;
		ftCoef = 0;
		for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
		{
			if (fi->bSelected)
			{
				hh0 = fi->halfedge();
				hh = hh0;
				do
				{
					if (hh->lpftGlobalPara[0] < lpftOffset[0])
					{
						lpftOffset[0] = hh->lpftGlobalPara[0];
					}
					if (hh->lpftGlobalPara[1] < lpftOffset[1])
					{
						lpftOffset[1] = hh->lpftGlobalPara[1];
					}
					if (hh->lpftGlobalPara[0] - lpftOffset[0] > ftCoef)
					{
						ftCoef = hh->lpftGlobalPara[0] - lpftOffset[0];
					}
					if (hh->lpftGlobalPara[1] - lpftOffset[1] > ftCoef)
					{
						ftCoef = hh->lpftGlobalPara[1] - lpftOffset[1];
					}
					hh = hh->next();
				} while (hh != hh0);
			}
		}
	}
	for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
	{
		vi->iTempIndex = -1;
		hh0 = vi->halfedge();
		hh = hh0;
		do
		{
			hh->iTempIndex = -1;
			hh = hh->next()->opposite();
		} while (hh != hh0);
	}
	iIndex = 0;
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
	{
		if (fi->bSelected)
		{
			hh0 = fi->halfedge();
			hh = hh0;
			do
			{
				if (hh->vertex()->iTempIndex == -1)
				{
					hh->vertex()->iTempIndex = iIndex;
					ofs << "v " << hh->vertex()->point().x() << ' ' << hh->vertex()->point().y() << ' ' << hh->vertex()->point().z() << " \n";
					++iIndex;
				}
				hh = hh->next();
			} while (hh != hh0);
		}
	}
	iIndex = 0;
	for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
	{
		if (vi->iTempIndex != -1)
		{
			if (vi->nCrevasseSize < 2)
			{
				hh0 = vi->halfedge();
				hh = hh0;
				do
				{
					hh->iTempIndex = iIndex;
					hh = hh->next()->opposite();
				} while (hh != hh0);
				if (bNormalizeTexture)
				{
					ofs << "vt " << (hh0->lpftGlobalPara[0] - lpftOffset[0]) / ftCoef << ' ' << (hh0->lpftGlobalPara[1] - lpftOffset[1]) / ftCoef << " \n";
				}
				else
				{
					ofs << "vt " << hh0->lpftGlobalPara[0] << ' ' << hh0->lpftGlobalPara[1] << " \n";
				}
				++iIndex;
			}
			else
			{
				hh0 = vi->halfedge()->hhCrevassePrev;
				hhStart = hh0;
				do
				{
					hhEnd = hhStart->hhCrevasseNext->opposite();
					bSelected = false;
					hh = hhStart;
					while (hh != hhEnd)
					{
						if (!hh->is_border() && hh->facet()->bSelected)
						{
							bSelected = true;
						}
						hh = hh->next()->opposite();
					}
					if (bSelected)
					{
						hh = hhStart;
						while (hh != hhEnd)
						{
							hh->iTempIndex = iIndex;
							hh = hh->next()->opposite();
						}
						if (bNormalizeTexture)
						{
							ofs << "vt " << (hhStart->lpftGlobalPara[0] - lpftOffset[0]) / ftCoef << ' ' << (hhStart->lpftGlobalPara[1] - lpftOffset[1]) / ftCoef << " \n";
						}
						else
						{
							ofs << "vt " << hhStart->lpftGlobalPara[0] << ' ' << hhStart->lpftGlobalPara[1] << " \n";
						}
						++iIndex;
					}
					hhStart = hhEnd;
				} while (hhStart != hh0);
			}
		}
	}
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
	{
		if (fi->bSelected)
		{
			ofs << "f ";
			hh0 = fi->halfedge();
			hh = hh0;
			do
			{
				ofs << hh->vertex()->iTempIndex + 1 << '/' << hh->iTempIndex + 1 << ' ';
				hh = hh->next();
			} while (hh != hh0);
			ofs << '\n';
		}
	}
}

void CMeshInfoReserver::m_fnExportPartPara(std::ofstream &ofs, bool bNormalizeTexture)
{
	CMeshBase::CPolyhedron::Vertex_iterator vi;
	CMeshBase::CPolyhedron::Halfedge_handle hh0, hh, hhStart, hhEnd;
	CMeshBase::CPolyhedron::Facet_iterator fi;
	bool bSelected;
	int iIndex;
	CMeshBase::FT lpftOffset[2], ftCoef;
	if (bNormalizeTexture)
	{
		lpftOffset[0] = DBL_MAX;
		lpftOffset[1] = DBL_MAX;
		ftCoef = 0;
		for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
		{
			if (fi->bSelected)
			{
				hh0 = fi->halfedge();
				hh = hh0;
				do
				{
					if (hh->lpftGlobalPara[0] < lpftOffset[0])
					{
						lpftOffset[0] = hh->lpftGlobalPara[0];
					}
					if (hh->lpftGlobalPara[1] < lpftOffset[1])
					{
						lpftOffset[1] = hh->lpftGlobalPara[1];
					}
					if (hh->lpftGlobalPara[0] - lpftOffset[0] > ftCoef)
					{
						ftCoef = hh->lpftGlobalPara[0] - lpftOffset[0];
					}
					if (hh->lpftGlobalPara[1] - lpftOffset[1] > ftCoef)
					{
						ftCoef = hh->lpftGlobalPara[1] - lpftOffset[1];
					}
					hh = hh->next();
				} while (hh != hh0);
			}
		}
	}
	iIndex = 0;
	for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
	{
		vi->iTempIndex = iIndex;
		ofs << "v " << vi->point().x() << ' ' << vi->point().y() << ' ' << vi->point().z() << " \n";
		++iIndex;
	}
	iIndex = 0;
	for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
	{
		if (vi->nCrevasseSize < 2)
		{
			hh0 = vi->halfedge();
			hh = hh0;
			bSelected = false;
			do
			{
				if (!hh->is_border() && hh->facet()->bSelected)
				{
					bSelected = true;
				}
				hh = hh->next()->opposite();
			} while (hh != hh0);
			if (bSelected)
			{
				do
				{
					hh->iTempIndex = iIndex;
					hh = hh->next()->opposite();
				} while (hh != hh0);
				if (bNormalizeTexture)
				{
					ofs << "vt " << (hh0->lpftGlobalPara[0] - lpftOffset[0]) / ftCoef << ' ' << (hh0->lpftGlobalPara[1] - lpftOffset[1]) / ftCoef << " \n";
				}
				else
				{
					ofs << "vt " << hh0->lpftGlobalPara[0] << ' ' << hh0->lpftGlobalPara[1] << " \n";
				}
				++iIndex;
			}
		}
		else
		{
			hh0 = vi->halfedge()->hhCrevassePrev;
			hhStart = hh0;
			do
			{
				hhEnd = hhStart->hhCrevasseNext->opposite();
				bSelected = false;
				hh = hhStart;
				while (hh != hhEnd)
				{
					if (!hh->is_border() && hh->facet()->bSelected)
					{
						bSelected = true;
					}
					hh = hh->next()->opposite();
				}
				if (bSelected)
				{
					hh = hhStart;
					while (hh != hhEnd)
					{
						hh->iTempIndex = iIndex;
						hh = hh->next()->opposite();
					}
					if (bNormalizeTexture)
					{
						ofs << "vt " << (hhStart->lpftGlobalPara[0] - lpftOffset[0]) / ftCoef << ' ' << (hhStart->lpftGlobalPara[1] - lpftOffset[1]) / ftCoef << " \n";
					}
					else
					{
						ofs << "vt " << hhStart->lpftGlobalPara[0] << ' ' << hhStart->lpftGlobalPara[1] << " \n";
					}
					++iIndex;
				}
				hhStart = hhEnd;
			} while (hhStart != hh0);
		}
	}
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
	{
		if (fi->bSelected)
		{
			ofs << "f ";
			hh0 = fi->halfedge();
			hh = hh0;
			do
			{
				ofs << hh->vertex()->iTempIndex + 1 << '/' << hh->iTempIndex + 1 << ' ';
				hh = hh->next();
			} while (hh != hh0);
			ofs << '\n';
		}
	}
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
	{
		if (!fi->bSelected)
		{
			ofs << "f ";
			hh0 = fi->halfedge();
			hh = hh0;
			do
			{
				ofs << hh->vertex()->iTempIndex + 1 << ' ';
				hh = hh->next();
			} while (hh != hh0);
			ofs << '\n';
		}
	}
}

void CMeshInfoReserver::m_fnExportAngleError(std::ofstream &ofs, int iPara, double dblMin, double dblMax)
{
	CMeshBase::CPolyhedron::Facet_iterator fi;
	CMeshBase::CPolyhedron::Vertex_iterator vi;
	CMeshBase::CPolyhedron::Halfedge_handle hh0, hh;
	double dblValue;
	int iIndex;
	iIndex = 0;
	for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
	{
		ofs << "v " << vi->point().x() << ' ' << vi->point().y() << ' ' << vi->point().z() << " \n";
		vi->iTempIndex = iIndex;
		++iIndex;
	}
	iIndex = 0;
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
	{
		if (iPara == 0 || iPara == 1)
		{
			dblValue = (fi->lpftAngleErr[iPara] - dblMin) / (dblMax - dblMin);
		}
		else
		{
			dblValue = (fi->lpftAngleErr[0] + fi->lpftAngleErr[1] - dblMin * 2) / ((dblMax - dblMin) * 2);
		}
		if (dblValue > 1.0)
		{
			dblValue = 1.0;
		}
		if (dblValue < 0.0)
		{
			dblValue = 0.0;
		}
		ofs << "vt " << 0.5 << " " << dblValue << '\n';
		fi->iTempIndex = iIndex;
		++iIndex;
	}
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
	{
		ofs << "f ";
		hh0 = fi->halfedge();
		hh = hh0;
		do
		{
			ofs << (hh->vertex()->iTempIndex + 1) << '/' << fi->iTempIndex + 1 << " ";
			hh = hh->next();
		} while (hh != hh0);
		ofs << '\n';
	}
}

void CMeshInfoReserver::m_fnExportConfCoefColor(std::ofstream &ofs)
{
	CMeshBase::CPolyhedron::Vertex_iterator vi;
	CMeshBase::CPolyhedron::Facet_iterator fi;
	CMeshBase::CPolyhedron::Halfedge_handle hh0, hh;
	int iInd;
	CMeshBase::FT ftMax, ftMin, ftAve;
	ftMax = -DBL_MAX;
	ftMin = DBL_MAX;
	ftAve = 0;
	for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
	{
		if (vi->ftDensity > ftMax)
		{
			ftMax = vi->ftDensity;
		}
		if (vi->ftDensity < ftMin)
		{
			ftMin = vi->ftDensity;
		}
		ftAve += vi->ftDensity;
	}
	ftAve /= m_mshSurface.size_of_vertices();
	for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
	{
		ofs << "v " << vi->point().x() << ' ' << vi->point().y() << ' ' << vi->point().z() << " \n";
	}
	iInd = 0;
	for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
	{
		ofs << "vt 0.5 " << (vi->ftDensity - ftMin) / (ftMax - ftMin) << " \n";
		vi->iIndex = iInd;
		++iInd;
	}
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
	{
		ofs << "f ";
		hh0 = fi->halfedge();
		hh = hh0;
		do
		{
			ofs << hh->vertex()->iIndex + 1 << '/' << hh->vertex()->iIndex + 1 << ' ';
			hh = hh->next();
		} while (hh != hh0);
		ofs << '\n';
	}
	ofs << "# " << ftMin << ' ' << ftAve << ' ' << ftMax << " \n";
}

void CMeshInfoReserver::m_fnExportTriangularPara(std::ofstream &ofs)
{
	CMeshBase::CPolyhedron::Vertex_iterator vi;
	CMeshBase::CPolyhedron::Facet_iterator fi;
	CMeshBase::CPolyhedron::Halfedge_handle hh0, hh;
	int iIndex, i;
	iIndex = 0;
	for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
	{
		ofs << "v " << vi->point().x() << ' ' << vi->point().y() << ' ' << vi->point().z() << '\n';
		vi->iIndex = iIndex;
		++iIndex;
	}
	ofs << "vt " << 0 << ' ' << 0 << '\n';
	ofs << "vt " << sqrt(0.375) + sqrt(0.125) << ' ' << sqrt(0.375) - sqrt(0.125) << '\n';
	ofs << "vt " << sqrt(0.375) - sqrt(0.125) << ' ' << sqrt(0.375) + sqrt(0.125) << '\n';
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
	{
		ofs << "f ";
		hh0 = fi->halfedge();
		hh = hh0;
		for (i = 0; i < 3; ++i)
		{
			ofs << hh->vertex()->iIndex + 1 << '/' << (i + 1);
			if (i == 2)
			{
				ofs << '\n';
			}
			else
			{
				ofs << ' ';
			}
			hh = hh->next();
		}
	}
}

void CMeshInfoReserver::m_fnExportObjFile(std::ofstream &ofs)
{
	CMeshBase::CPolyhedron::Vertex_iterator vi;
	CMeshBase::CPolyhedron::Facet_iterator fi;
	CMeshBase::CPolyhedron::Halfedge_handle hh0, hh;
	int iIndex;
	iIndex = 0;
	for (vi = m_mshSurface.vertices_begin(); vi != m_mshSurface.vertices_end(); ++vi)
	{
		vi->iIndex = iIndex;
		ofs << "v " << vi->point().x() << ' ' << vi->point().y() << ' ' << vi->point().z() << '\n';
		++iIndex;
	}
	for (fi = m_mshSurface.facets_begin(); fi != m_mshSurface.facets_end(); ++fi)
	{
		ofs << "f ";
		hh0 = fi->halfedge();
		hh = hh0;
		do
		{
			ofs << hh->vertex()->iIndex + 1;
			if (hh->next() == hh0)
			{
				ofs << '\n';
			}
			else
			{
				ofs << ' ';
			}
			hh = hh->next();
		} while (hh != hh0);
	}
}

void CMeshInfoReserver::m_fnExportQuadObjFile(std::ofstream &ofs)
{
	CMeshBase::CSamplerPolyhedron::Vertex_iterator vi;
	CMeshBase::CSamplerPolyhedron::Facet_iterator fi;
	CMeshBase::CSamplerPolyhedron::Halfedge_handle hh0, hh, hhStart;
	int iIndex, i;
	iIndex = 0;
	for (vi = m_mshRemeshedSurface.vertices_begin(); vi != m_mshRemeshedSurface.vertices_end(); ++vi)
	{
		vi->iIndex = iIndex;
		ofs << "v " << vi->point().x() << ' ' << vi->point().y() << ' ' << vi->point().z() << '\n';
		++iIndex;
	}
	ofs << "vt " << 0 << ' ' << 0 << '\n';
	ofs << "vt " << 1.0 << ' ' << 0 << '\n';
	ofs << "vt " << 1.0 << ' ' << 1.0 << '\n';
	ofs << "vt " << 0 << ' ' << 1.0 << '\n';
	for (fi = m_mshRemeshedSurface.facets_begin(); fi != m_mshRemeshedSurface.facets_end(); ++fi)
	{
		ofs << "f ";
		hhStart = NULL;
		hh0 = fi->halfedge();
		hh = hh0;
		do
		{
			if (hh->vertex()->iMark == 1)
			{
				hhStart = hh;
			}
			hh = hh->next();
		} while (hh != hh0);
		if (hhStart != NULL)
		{
			hh = hhStart;
		}
		for (i = 0; i < 4; ++i)
		{
			ofs << hh->vertex()->iIndex + 1 << '/' << (i + 1);
			if (i == 3)
			{
				ofs << '\n';
			}
			else
			{
				ofs << ' ';
			}
			hh = hh->next();
		}
	}
}


void CMeshInfoReserver::m_fnExportQuadCylinder(std::ofstream &ofs)
{
	CMeshBase::CSamplerPolyhedron::Vertex_iterator vi;
	CMeshBase::CSamplerPolyhedron::Edge_iterator ei;
	int iIndex;
	iIndex = 0;
	for (vi = m_mshRemeshedSurface.vertices_begin(); vi != m_mshRemeshedSurface.vertices_end(); ++vi)
	{
		ofs << "v " << vi->point().x() << ' ' << vi->point().y() << ' ' << vi->point().z() << " \n";
		vi->iIndex = iIndex;
		++iIndex;
	}
	ofs << "g border_edges\n";
	for (ei = m_mshRemeshedSurface.edges_begin(); ei != m_mshRemeshedSurface.edges_end(); ++ei)
	{
		if (ei->is_border_edge())
		{
			ofs << "l " << ei->prev()->vertex()->iIndex + 1 << ' ' << ei->vertex()->iIndex + 1 << " \n";
		}
	}
	ofs << "g inner_edges\n";
	for (ei = m_mshRemeshedSurface.edges_begin(); ei != m_mshRemeshedSurface.edges_end(); ++ei)
	{
		if (!ei->is_border_edge())
		{
			ofs << "l " << ei->prev()->vertex()->iIndex + 1 << ' ' << ei->vertex()->iIndex + 1 << " \n";
		}
	}
}
/*
    
Copyright (c) 2024 Centre for Water Systems,
                   University of Exeter

This file is part of cafloodpro.

cafloodpro is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>
*/



//! \file GDAL_Functions.cpp
//!  Contains the structure and function for licensing options
//! \author Mike Gibson, University of Exeter, 
//! contact: m.j.gibson [at] exeter.ac.uk
//! \date 2019-12

#ifdef CAAPI_GDAL_OPTION


#include "GDAL_Functions.hpp"



int saveRasterAs(CA::AsciiGrid<CA::Real> grid, const std::string driverName, std::string geoCS, std::string fileName)
{


	GDALAllRegister();

	GDALDataset* dataSet;
	char **papszOptions = NULL;

	
	//GDALDriverManager *driverMng = new GDALDriverManager();
	GDALDriverManager* driverMng = GetGDALDriverManager();

	GDALDriver* driver;
	driver = driverMng->GetDriverByName(driverName.c_str());

	dataSet = driver->Create(fileName.c_str(), grid.ncols, grid.nrows, 1, GDT_CFloat64, papszOptions);


	double adfGeoTransform[6] = { grid.xllcorner, grid.cellsize,0,grid.yllcorner + (grid.nrows * grid.cellsize), 0, -grid.cellsize };
	OGRSpatialReference oSRS;
	char* pszSRS_WKT = NULL;
	GDALRasterBand* poBand;
	dataSet->SetGeoTransform(adfGeoTransform);

	//oSRS.SetUTM(11, TRUE);
	oSRS.SetWellKnownGeogCS(geoCS.c_str());
	//oSRS.importFromEPSG
	oSRS.exportToWkt(&pszSRS_WKT);

	//double gridData[grid.ncols * grid.nrows];
	//for(int r = 0; r < grid.nrows; )



	dataSet->SetProjection(pszSRS_WKT);
	CPLFree(pszSRS_WKT);
	poBand = dataSet->GetRasterBand(1);
	//poBand->RasterIO(GF_Write, 0, 0, grid.ncols, grid.nrows,(void*)grid.data.[0], grid.ncols, grid.nrows, GDT_Float64, 0, 0);
	/* Once we're done, close properly the dataset */
	GDALClose((GDALDatasetH)dataSet);




	return 0;
}









#endif
//##############################################################################



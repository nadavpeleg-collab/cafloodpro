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

//! \file SpatialTemporalDense.cpp
//! contact: m.j.gibson [at] exeter.ac.uk
//! \date 2024-05


#include"ca2D.hpp"
#include"ArgsData.hpp"
#include"SpatialTemporalDense.hpp"
#include"Utilities.hpp"
#include<iostream>
#include<fstream>

#ifdef CAAPI_NETCDF_OPTION
#include <netcdf.h>

//#include "GDAL_Functions.hpp"
#include "gdal_priv.h"

#endif


// -------------------------//
// Include the CA 2D functions //
// -------------------------//
//#include CA_2D_INCLUDE(writeSTIndices)
#include CA_2D_INCLUDE(addRainSTD)
#include CA_2D_INCLUDE(removeNegatives)


int initSpatialTemporalDenseFromCSV(const std::string& dir,
    const std::vector<std::string>& filenames, SpatialTemporalDense& std, bool userNetCDF)
{

  // Default value.
  //std.l   = 0;
  std.end = 0;

  

  // If the vector is empty we don't do anything.
  if(filenames.size()==0)
    return 0;

#ifndef CAAPI_NETCDF_OPTION

  // However if the vector is not empty we nee to check that there is
  // the righ amount of elements.
  if(filenames.size()==1)
  {
    std::cerr<<"Spatial Temporal variable setup must contains at least two files"<<std::endl;
    return 1;
  }

  // Load the CSV file with the indices.
  std::string csvfile (dir+filenames[0]);
  std::ifstream ifile(csvfile.c_str());
  
  if(!ifile)
  {
    std::cerr<<"Error opening CSV file: "<<csvfile<<std::endl;
    return 1;
  }

  // Parse the file line by line until the end of file 
  // and retrieve the tokens of each line.
  while(!ifile.eof())
  {
    // If true the token was identified;
    bool found_tok = false;

    std::vector<std::string> tokens( CA::getLineTokens(ifile, ',') );
    
    // If the tokens vector is empty we reached the eof or an
    // empty line... continue.
    if(tokens.empty())
      continue;       

    if(CA::compareCaseInsensitive("Name",tokens[0],true))
    {
      std::string str;
      READ_TOKEN(found_tok,str,tokens[1],tokens[0]);
      
      std.name = CA::trimToken(str," \t\r");
    }



    if(CA::compareCaseInsensitive("Times",tokens[0],true))
    {
      // Get the first token string without the "Time" substring.
      std::string str = tokens[0].substr(5, tokens[0].size() - 5);



      // Read the sequence of values.
      for (size_t i=1; i<tokens.size(); ++i)
      {
	    CA::Real value;
	    READ_TOKEN(found_tok,value,tokens[i],tokens[0]);

        std.times.push_back(value);

	    // Check if it is the last end
	    std.end = std::max(value, std.end);
      }//end of looping times values list
    }// end of TIMES token 
    
    // If the token was not identified stop!
    if(!found_tok)
    {
      std::cerr<<"Element '"<<CA::trimToken(tokens[0])<<"' not identified"<<std::endl; \
      return 1;
    }
  }

  // Read the rest of the remaining file names
  for (size_t i = 1; i < filenames.size(); ++i)
  {


      try
      {
          // ATTENTION. Load only the header here .. not the actual data
          //CA::readAsciiGridHeader(st.map,dir+filenames[0]);
          CA::AsciiGrid<CA::Real> inputGrid;
          inputGrid.readAsciiGridHeader(dir + filenames[i]);


          std.maps.push_back(inputGrid);

          std.map_filenames.push_back(dir + filenames[i]);


      }
      catch (std::runtime_error & e)
      {
          std::cerr << e.what() << std::endl;
          return 1;
      }
  }// end of looping remaining file names

    // check sizes match here....

#else

/* Handle errors by printing an error message and exiting with a
 * non-zero status. */
#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}

    /* This will be the netCDF ID for the file and data variable. */
    int ncid, varid;

    /* Loop indexes, and error handling. */
    int retval;

    std::string netCDFfile(dir + filenames[0]);

    /* Open the file. NC_NOWRITE tells netCDF we want read-only access
     * to the file.*/
    if ((retval = nc_open(netCDFfile.c_str(), NC_NOWRITE, &ncid)))
        ERR(retval);

    std::string timeName(filenames[2]);

    /* Get the varid of the data variable, based on its name. */
    if ((retval = nc_inq_varid(ncid, timeName.c_str(), &varid)))
        ERR(retval);

    nc_type rh_type;
    int rh_ndims;
    int  rh_dimids[NC_MAX_VAR_DIMS];
    int rh_natts;
    size_t latlength;
    int latid;

    /* Check details of the time variable. */
    if (retval = nc_inq_var(ncid, varid, 0, &rh_type, &rh_ndims, rh_dimids,&rh_natts))
        ERR(retval);


    if (retval = nc_inq_dimid(ncid, timeName.c_str(), &latid))
        ERR(retval);
    //if (status != NC_NOERR) handle_error(status);
    if (retval = nc_inq_dimlen(ncid, latid, &latlength))
        ERR(retval);
    //if (status != NC_NOERR) handle_error(status);

    //std::cout << "dim len: " << latlength << std::endl;


    if (rh_type == 5) {

        std::vector<float> data_in;

        data_in.resize(latlength);

        if ((retval = nc_get_var_float(ncid, varid, data_in.data())))
            ERR(retval);

        for (int x = 0; x < latlength; x++) {
            //std::cout << "data[" << x << "]: " << data_in[x] << std::endl;
            std.times.push_back(data_in[x]);
            // Check if it is the last end
            std.end = std::max((CA::Real)data_in[x], std.end);
        }
    }


    if (rh_type == 6) {

        std::vector<double> data_in;

        data_in.resize(latlength);

        if ((retval = nc_get_var_double(ncid, varid, data_in.data())))
            ERR(retval);

        for (int x = 0; x < latlength; x++) {
            //std::cout << "data[" << x << "]: " << data_in[x] << std::endl;
            std.times.push_back(data_in[x]);
            // Check if it is the last end
            std.end = std::max((CA::Real)data_in[x], std.end);
        }

    }

    /* Close the file, freeing all resources. */
    if ((retval = nc_close(ncid)))
        ERR(retval);




    //----------------------------------------------------------


    GDALAllRegister();

    GDALDataset* dataSet;

    std::ostringstream netCDFfileGDALOpenoss;
    netCDFfileGDALOpenoss << "NETCDF:" << netCDFfile << ":" << filenames[1];
    std::string netCDFfileGDALOpen;



    dataSet = (GDALDataset*)GDALOpen(netCDFfile.c_str(), GA_ReadOnly);

    int rasterCount = dataSet->GetRasterCount();

    // ulx, xres, xskew, uly, yskew, yres 
    double adfGeoTransform[6];
    dataSet->GetGeoTransform(adfGeoTransform);

    for (int currentRaster = 1; currentRaster < rasterCount + 1; ++currentRaster) {

        CA::AsciiGrid<CA::Real> inputGrid;

        inputGrid.ncols = dataSet->GetRasterXSize();
        inputGrid.nrows = dataSet->GetRasterYSize();
        inputGrid.xllcorner = adfGeoTransform[0];
        // check for that negative?
        inputGrid.yllcorner = adfGeoTransform[3] + (dataSet->GetRasterYSize() * adfGeoTransform[5]);

        // check for negative and that they are both the same value?
        inputGrid.cellsize = adfGeoTransform[1];

        GDALRasterBand* poBand;

        // might be zero in c++?
        poBand = dataSet->GetRasterBand(currentRaster);

        inputGrid.nodata = poBand->GetNoDataValue();

        // Allocate data and set the default value to nodata.
        inputGrid.data.clear();
        inputGrid.data.resize(inputGrid.ncols * inputGrid.nrows, inputGrid.nodata);


        int   nXSize = poBand->GetXSize();
        int   nYSize = poBand->GetYSize();


        poBand->RasterIO(GF_Read, 0, 0, nXSize, nYSize,
            static_cast<void*>(inputGrid.data.data()), nXSize, nYSize, GDT_Float64,
            0, 0);


        std.maps.push_back(inputGrid);

        std::ostringstream oss;
        oss << netCDFfile << currentRaster;
        //std::cout << oss.str();

        std.map_filenames.push_back(oss.str());

    }



    //delete(poBand);
    delete(dataSet);

    // check sizes match here....

#endif

  

  return 0;
}


int loadSpatialTemporalDenseData(CA::Grid& GRID, SpatialTemporalDense& std, CA::Box domain, bool loadGrids)
{

    int maxCount = std.map_filenames.size();
    if(!loadGrids)
        maxCount = std.maps.size();

    // Read the rest of the remaining file names
    for (size_t i = 0; i < maxCount; ++i)
    {


        try
        {
            // ATTENTION. Load only the actual data
            //CA::readAsciiGrid(st.map,dir+filenames[0]);
            
            
            if (loadGrids) {
                std.maps[i].readAsciiGrid(std.map_filenames[i]);
            }

            CA::CellBuffReal TMP(GRID);

            CA::BoxList  fulldomain;
            CA::Box      fullbox = GRID.box();
            fulldomain.add(fullbox);

            TMP.fill(fulldomain, 0.0);
            
            // Insert the DEM data into the elevation cell buffer. The first
            // parameter is the region of the buffer to write the data which is
            // the real domain, i.e. the original non extended domain.
            TMP.insertData(domain, &std.maps[i].data[0], std.maps[i].ncols, std.maps[i].nrows);


            CA::Execute::function(fulldomain, removeNegatives, GRID, TMP, std.maps[i].nodata);


            CA::BoxList  realdomain;
            realdomain.add(domain);

            CA::Real rainMax = 0.0;

            TMP.sequentialOp(realdomain, rainMax, CA::Seq::Max);

            std.rain_maxes.push_back(rainMax);

            //----

            CA::Real rainTotal = 0.0;

            TMP.sequentialOp(realdomain, rainTotal, CA::Seq::Add);

            std.rain_totals.push_back(rainTotal);



        }
        catch (std::runtime_error & e)
        {
            std::cerr << e.what() << std::endl;
            return 1;
        }


    }// end of looping remaining file names


    return 0;

}// end of loadSpatialTemporalDenseData







STRainDense::STRainDense(CA::Grid&  GRID, const SpatialTemporalDense& std):
  _grid(GRID),
  _std(std),
  _data(),
    RAIN_STD(GRID)
{
  //_data.values.resize(_st.seqs.size());
  //_data.steps.resize(_st.seqs.size());
    _data.index = 0;
    //_data.RAIN_STD(GRID);

      // Create the full (extended) computational domain of CA grid. 
    CA::BoxList  fulldomain;
    CA::Box      fullbox = _grid.box();
    fulldomain.add(fullbox);
    // Se the default value of the elevation to be nodata.
    RAIN_STD.fill(fulldomain, 0.0);

    _data.volume = 0.0;

    //_data.volumes.resize(std.rain_totals.size());

    _data.firstUploadDone = false;


}

STRainDense::~STRainDense()
{
}


void STRainDense::addDomain(CA::BoxList& compdomain)
{


	// NOTE: this needs improving as, currently adding the entire domain for every spatial area.

    if(_std.map_filenames.size()>0)
    {
        // Ok easy here ... we add all the full box domain to be sure.
        compdomain.add(_grid.box()); //???????
    }
  
}


void STRainDense::analyseArea(CA::CellBuffReal& TMP, CA::CellBuffState& MASK, CA::BoxList&  domain)
{

}


CA::Real STRainDense::potentialVA(CA::Real t, CA::Real period_time_dt)
{
    CA::Real potential_va = 0.0;

    // Use the maximum amount of rain that is going to fall in the next
    // period to calculate the possible water depth.  Use the critical
    // velocity with this water depth to compute the potential velocity.


    // Get the index.
    size_t index = _data.index;

    // If the index is larger than the available rain/time, do
    // nothing.
    if (index >= _std.times.size())
        return potential_va;

    // Check if the simulation time now is equal or higher than the
    // time when this rain intensity ends. If it is the case,
    // increase the index to the next rain intensity.
    if (t >=_std.times[_data.index])
        index++;

    // If the index is larger than the available rain/time, do
    // nothing.
    if (index >= _std.times.size())
        return potential_va;
    
    // Get the value of the rain (transformed in metres from mm) for next period.
    CA::Real rain = (_std.rain_maxes[index] * 0.001) * (period_time_dt / 3600.0);

    // Compute the potential velocity.
    potential_va = std::max(potential_va, std::sqrt(static_cast<CA::Real>(9.81)* (rain)));
    

    return potential_va;
}





CA::Real STRainDense::volume()
{
  CA::Real rain_volume = 0.0;

  rain_volume = _data.volume;
  
  return rain_volume;
}





CA::Real STRainDense::endTime()
{
  
  // return the last time in the time series, if there are any
    CA::Real t_end = 0.0;

    // assuming the times are in order, return the last time
    if(_std.times.size()>0)
        t_end = _std.times[_std.times.size() - 1];

    return t_end;
}


void STRainDense::prepare(CA::Real t, CA::Real period_time_dt, CA::Real next_dt)
{

    // Get the index.
    size_t index = _data.index;

    // If the index is larger than the available rain/time, do
    // nothing.
    if (index >= _std.times.size())
        return;

    _data.volume = 0.0;

    // Check if the simulation time now is equal or higher than the
    // time when this rain intensity ends. If it is the case,
    // increase the index to the next rain intensity.
    if (t >= _std.times[_data.index])
        index++;

    // If the index is larger than the available rain/time, do
    // nothing.
    if (index >= _std.times.size())
        return;

    // Get the rain (transformed in metres from mm) for the next period.
    CA::Real period_rain = (_std.rain_totals[index] * 0.001) * (period_time_dt / 3600.0);

    // Add the volume. (total rain depths, multiple by a cell area)
    _data.volume = period_rain * _grid.area();

    if (_data.firstUploadDone == false && _data.index == 0) {


        CA::BoxList  fulldomain;
        CA::Box      fullbox = _grid.box();
        fulldomain.add(fullbox);

        RAIN_STD.fill(fulldomain, 0.0);

        CA::Box      realbox(_grid.box().x() + 1, _grid.box().y() + 1, _grid.box().w() - 2, _grid.box().h() - 2);
        CA::BoxList  realdomain;
        realdomain.add(realbox);


        // Insert the DEM data into the elevation cell buffer. The first
        // parameter is the region of the buffer to write the data which is
        // the real domain, i.e. the original non extended domain.
        RAIN_STD.insertData(realbox, &_std.maps[index].data[0], _std.maps[index].ncols, _std.maps[index].nrows);

        _data.currentNoData = _std.maps[index].nodata;


        //std::ofstream outfile;
        //outfile.open("temp.asc");
        //RAIN_STD.dump(outfile);
        //outfile.close();

        _data.firstUploadDone = true;

    }

    if (index != _data.index) {
        CA::BoxList  fulldomain;
        CA::Box      fullbox = _grid.box();
        fulldomain.add(fullbox);

        RAIN_STD.fill(fulldomain, 0.0);

        CA::Box      realbox(_grid.box().x() + 1, _grid.box().y() + 1, _grid.box().w() - 2, _grid.box().h() - 2);
        CA::BoxList  realdomain;
        realdomain.add(realbox);


        // Insert the DEM data into the elevation cell buffer. The first
        // parameter is the region of the buffer to write the data which is
        // the real domain, i.e. the original non extended domain.
        RAIN_STD.insertData(realbox, &_std.maps[index].data[0], _std.maps[index].ncols, _std.maps[index].nrows);

        _data.currentNoData = _std.maps[index].nodata;

        _data.firstUploadDone = true;

    }



    // Update index.
    _data.index = index;


}// end of prepare


void STRainDense::add(CA::CellBuffReal& WD, CA::CellBuffState& MASK, CA::Real t, CA::Real period_time_dt)
{


    // Check if the simulation time now is equal or higher than the
    // time when this rain intensity ends. If it is the case,
    // increase the index to the next rain intensity.
    if (t > _std.end)
        return;

    CA::BoxList  fulldomain;
    CA::Box      fullbox = _grid.box();
    fulldomain.add(fullbox);

    CA::Execute::function(fulldomain, addRainSTD, _grid, WD, MASK, RAIN_STD, period_time_dt, _data.currentNoData);



}// end of add function
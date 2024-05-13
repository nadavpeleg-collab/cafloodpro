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

//! \file preProc.cpp
//! Perform the pre processing of the data for a CA 2D model.
//! contact: a.s.chen [at] exeter.ac.uk
//! \date 2024-05


#include"ca2D.hpp"
#include"ArgsData.hpp"
#include"Setup.hpp"
#include"SpatialTemporal.hpp"
#include"SpatialTemporalDense.hpp"





// -------------------------//
// Include the CA 2D functions //
// -------------------------//
#include CA_2D_INCLUDE(geoRatio)
#include CA_2D_INCLUDE(addHeight)


//! Perform the pre processing of the data for a CA 2D model. It
//! mainly save the elevation file.
//! \attention The GRID size add an extra set of cells in each directions.
int preProc(const ArgsData& ad, const Setup& setup, const std::string& ele_file,
	const std::vector<SpatialTemporal>& sts,
    SpatialTemporalDense& _std)
{
  
  if(setup.output_computation)
  {
    std::cout<<"Pre-processing : "<<setup.sim_name<< std::endl;
    std::cout<<"------------------------------------------" << std::endl; 
  }

  // ----  Timer ----
  
  // Get starting time.
  CA::Clock total_timer;

  // Check if pre-process data is already there
  if( CA::Grid::exist(ad.data_dir,setup.preproc_name+"_Grid","0") &&
      CA::CellBuffReal::existData(ad.data_dir,setup.preproc_name+"_ELV","0") )
  {

    if(setup.output_console)
      std::cout<<"Pre-proc data already exist"<< std::endl;    
  }
  else
  {

    // -- ELEVATION FILE(s) --
    
	  CA::AsciiGrid<CA::Real> eg;

	  if (setup.loadTerrain)
	  {
		  //AsciiGrid<CA::Real> eg;
		  eg.readAsciiGrid(ele_file);
		  //eg.copyData(setup.eg);
	  }
	  else
	  {
		  //eg = setup.eg;
		  eg.copyData(setup.eg);
	  }

	  
      //eg.writeAsciiGrid(ad.output_dir + ad.sdir + "elevationTest.tif");
	  


    
    // -- EXTRA CELLS --
    
    // Given the way the boundary cells are computed (see later), we
    // need to add an extra set of cells in each direction.
    CA::AsciiGrid<CA::Real> grid;
    grid.ncols     = eg.ncols+2;
    grid.nrows     = eg.nrows+2;
    grid.xllcorner = eg.xllcorner-eg.cellsize; // ATTENTION HERE!
    grid.yllcorner = eg.yllcorner-eg.cellsize; // ATTENTION HERE!
    grid.cellsize  = eg.cellsize;
    grid.nodata    = eg.nodata;
    
    // ----  CA GRID ----
    
    // Create the square regular CA grid of the extended size of the DEM.
    // The internal implementation could be different than a square regular grid.
    CA::Grid  GRID(grid.ncols, grid.nrows, grid.cellsize, grid.xllcorner, grid.yllcorner, ad.args.active());
    
    // Set if to print debug information on CA function.
    GRID.setCAPrint(false);
    
    // Set the place where the direct I/O data is going to be saved.
    GRID.setDataDir(ad.data_dir);
    
    // Save the data of the Grid
    if( !GRID.save(setup.preproc_name+"_Grid","0") )
    {
      std::cerr<<"Error while saving the GRID information"<<std::endl;
      return 1;
    }
    
    if(setup.output_console)
      std::cout<<"Saved Grid information"<< std::endl;
    
  
    // Create the full (extended) computational domain of CA grid. 
    CA::BoxList  fulldomain;
    CA::Box      fullbox = GRID.box();
    fulldomain.add(fullbox);
    
    // Create a borders object that contains all the borders and
    // corners of the grid.
    CA::Borders borders;
    
    borders.addSegment(CA::Top);
    borders.addSegment(CA::Bottom);
    borders.addSegment(CA::Right);
    borders.addSegment(CA::Left);
    
    borders.addCorner(CA::TopLeft);
    borders.addCorner(CA::TopRight);
    borders.addCorner(CA::BottomLeft);
    borders.addCorner(CA::BottomRight);
    
    // Create the real computational domain of CA grid, i.e. the
    // original DEM size not the extended one.
    CA::BoxList  realdomain;
    CA::Box      realbox(GRID.box().x()+1,GRID.box().y()+1,GRID.box().w()-2,GRID.box().h()-2);
    realdomain.add(realbox);
    
    // -- INITIALISE ELEVATION ---
    
    // Create the elevation cell buffer.
    // It contains a "real" value in each cell of the grid.
    CA::CellBuffReal  ELV(GRID);
    
	//std::cout << "Grid created" << std::endl;

    // Set the border of the elevation buffer to be no data. 
    ELV.bordersValue(borders,grid.nodata);

	//std::cout << "Borders done" << std::endl;
    
    // Se the default value of the elevation to be nodata.
    ELV.fill(fulldomain, grid.nodata);

	//exit(1);

	//std::cout << "Fill done" << std::endl;
    
    // Insert the DEM data into the elevation cell buffer. The first
    // parameter is the region of the buffer to write the data which is
    // the real domain, i.e. the original non extended domain.
    ELV.insertData(realbox, &eg.data[0], eg.ncols, eg.nrows);

	
	//ELV.syncSubGrids(ad.output_dir);

//#ifdef CA2D_MPI
//	std::stringstream subGridIdStream;
//	subGridIdStream << GRID.caGrid().rank;
//	std::string outFileName = ad.output_dir + "temp" + subGridIdStream.str() + ".txt";
//	std::cout << "Writing: " + outFileName;
//	std::ofstream outFile(outFileName);
//	ELV.dump(outFile);
//#endif // CA2D_MPI
    

    if(setup.output_console)
      std::cout<<"Saved Elevation data"<< std::endl;

    //-------------- spatio temporal Dense -------------

    if (_std.map_filenames.size() > 1 || setup.loadSpatialRainDense == false) {

#ifndef CAAPI_GDAL_OPTION
        if (loadSpatialTemporalDenseData(GRID, _std, realbox, setup.loadSpatialRainDense) != 0)
#else
        if (loadSpatialTemporalDenseData(GRID, _std, realbox, !setup.rainNetCDF) != 0)
#endif
        {
            std::cerr << "Error while Loading Spatiotemporal Dense data" << std::endl;
            return 1;
        }


        if (setup.output_console)
            std::cout << "Loaded Spatiotemporal Dense data" << std::endl;
    }



    //-------------- spatio temporal spare -------------

    // The Spatial/Temporal variables vector must have at least three
    // elements. If there are more than three elements. The extra one
    // are ignored.
    if(sts.size()<3)
    {
      std::cerr<<"Error while processing spatial temporal data"<<std::endl;
      return 1;    
    }
    
    // Save the spatial temporal roughness.
    if(saveSpatialTemporalMap(GRID,sts[0],setup.preproc_name+"_ROUGH","0", setup.loadRoughnesses)!=0)
    {
      std::cerr<<"Error while saving the Roughness data"<<std::endl;
      return 1;
    }

    // Save the spatial temporal infiltration.
    if(saveSpatialTemporalMap(GRID,sts[2],setup.preproc_name+"_INF","0", setup.loadSpatialInfiltrations)!=0)
    {
      std::cerr<<"Error while saving the Infiltration data"<<std::endl;
      return 1;
    }
    if(setup.output_console)
      std::cout<<"Saved Spatial Temporal data"<< std::endl;

    /*
    if (setup.output_console)
        std::cout << "loading Spatial Temporal Dense data" << std::endl;


    if (setup.output_console)
        std::cout << "loaded Spatial Temporal Dense data" << std::endl;
        */
	//sts[0].map

	if (setup.addHeights) {


		// load the map and data
		//sts[3].map.readAsciiGrid(ad.working_dir + ad.sdir + setup.heighAdjustmentMap_files[1]);


		//------------------------------------------

		//std::cout << "Making height Adjustments...["<< setup.addedHeight.size() <<"]" << std::endl;

		CA::TableReal ADDED_HEIGHTS(GRID, setup.addedHeight.size());

		ADDED_HEIGHTS.update(0, setup.addedHeight.size(), &setup.addedHeight[0], setup.addedHeight.size());

		//ADDED_HEIGHTS.dump(std::cout);

		CA::CellBuffState BUILDING_MASK(GRID);
		//BUILDING_MASK.insertData(realbox, &sts[0].map.data[0], sts[0].map.ncols, sts[0].map.nrows);
		BUILDING_MASK.insertData(realbox, &(setup.heightAdjustmentMap->data[0]), setup.heightAdjustmentMap->ncols, setup.heightAdjustmentMap->nrows);
		CA::Execute::function(fulldomain, addHeight, GRID, ELV, BUILDING_MASK, ADDED_HEIGHTS);
		//////////////////////////////////////////////////////////

		//ELV.dump(std::cout);

		//delete now unneeded map
		delete(setup.heightAdjustmentMap);
	}


	// Save the data of the Elevation
	if (!ELV.saveData(setup.preproc_name + "_ELV", "0"))
	{
		std::cerr << "Error while saving the Elevation data" << std::endl;
		return 1;
	}

    // Compute the Geographic ratio if the coordinate system is a
    // Geographic one, i.e. lat/long
    if(setup.latlong)
    {
      if(setup.output_console)
	std::cout<<"Prepare Geographic Coordinate System"<< std::endl;

      // Create the Geographic ratio buffer.
      CA::CellBuffReal  RATIO(GRID);

      // The earth radio.
      CA::Real R = 6371;
      // Compute the dy in meters using Equirectangular
      // approximation. This should be constant at all latitude
      // longitude.
      CA::Real dy = std::abs( GRID.length()*(PI/180.0)*R*1000 );
      // Compute the the lenght of a side in radian.
      CA::Real l  = GRID.length()*(PI/180.0);
      
      // Compute the geographic ratio.
      // ATTENTION This work only on a square grid.
      CA::Execute::function(fulldomain, geoRatio, GRID, RATIO, dy, R, l);

      // Save the data of the Elevation
      if( !RATIO.saveData(setup.preproc_name+"_RATIO","0") )
      {
	std::cerr<<"Error while saving the Geographic Coordinate System data"<<std::endl;
	return 1;
      }
            
      if(setup.output_console)
	std::cout<<"Saved Geographic Coordinate System"<< std::endl;
    }
  }





















  // ---- TIME OUTPUT ----

  if(setup.output_computation)
  {
    std::cout<<"-----------------" << std::endl; 
    std::cout<<"Total run time taken (s) = " << total_timer.millisecond()/1000.0 << std::endl;
    std::cout<<"-----------------" << std::endl; 
  }

  
  return 0;
}
//---------------------------------------------------------------------------------------------------------------------------



//---------------------------------------------------------------------------------------------------------------------------
int Test_CADDIES()
{

	//    CA::Grid  GRID(grid.ncols, grid.nrows, grid.cellsize, grid.xllcorner, grid.yllcorner, ad.args.active());
	std::cout << "Making a GRID!!!!" << std::endl;
	CA::Grid GRID(100, 100, 1, 31900, 17700);
	std::cout << "Made a GRID!!!!" << std::endl;


	std::cout << "Making a Cellbuff!!!!" << std::endl;
	CA::CellBuffReal  ELV(GRID);
	std::cout << "Made a Cellbuff!!!!" << std::endl;






	return 0;
}// end of Test_CADDIES
//---------------------------------------------------------------------------------------------------------------------------

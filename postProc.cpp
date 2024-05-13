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

//! \file postProc.cpp
//! Perform the post processing of the data for a CA 2D model.
//! contact: a.s.chen [at] exeter.ac.uk
//! \date 2024-05


#include"ca2D.hpp"
#include"Masks.hpp"
#include"ArgsData.hpp"
#include"Setup.hpp"
#include"TimePlot.hpp"
#include"RasterGrid.hpp"


// -------------------------//
// Include the CA 2D functions //
// -------------------------//
#include CA_2D_INCLUDE(zeroedWD)
#include CA_2D_INCLUDE(zeroedVA)
#include CA_2D_INCLUDE(flux2velocity)


//! Inline function that given the actual time, the time of the
//! nearest action and the period of the possible new action, find the
//! possible NEW simulation time of possible new nearest action.
//! \attention An action is: the simulation finished, output to
//! console, raster, timeplot, or rain event change, etc..
inline CA::Real nextTimeNearestAction(CA::Real t, CA::Real t_nearest, CA::Real period)
{
  return std::min(t_nearest, t + period - std::fmod(t,period));
}


// -------------------------//
// Include the CA 2D functions //
// -------------------------//
#include CA_2D_INCLUDE(makeWL)

//#include CA_2D_INCLUDE(checkFloodNoDataNeighbourhood)



//! Perform the post processing of the data for a CA 2D model. 
int postProc(const ArgsData& ad, const Setup& setup, CA::AsciiGrid<CA::Real>& eg,
	     const std::vector<TimePlot>& tps, std::vector<RasterGrid>& rgs,
		 CA::CellBuffReal* timeGrid)
{
  

	//std::cout << eg.getString() << std::endl;

  if(setup.output_computation)
  {
    std::cout<<"Post-processing : "<<setup.sim_name<< std::endl;
    std::cout<<"------------------------------------------" << std::endl; 
  }

  // ----  Timer ----
  
  // Get starting time.
  CA::Clock total_timer;



  




  // ----  CA GRID ----
  
  // Load the CA Grid from the DataDir. 
  // ATTENTION this should have an extra set of cells in each
  // direction.  The internal implementation could be different than a
  // square regular grid.
  CA::Grid  GRID(ad.data_dir,setup.preproc_name+"_Grid","0", ad.args.active());

  if(setup.output_console)
    std::cout<<"Loaded Grid data"<< std::endl;

  // Set if to print debug information on CA function.
  GRID.setCAPrint(false);

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

  

  


  // CREATE ASCII GRID TEMPORARY BUFFERS
  CA::AsciiGrid<CA::Real>& agtmp1 = eg;
  CA::AsciiGrid<CA::Real>  agtmp2 = eg;

  // -- INITIALISE ELEVATION ---

  // Create the elevation cell buffer.
  // It contains a "real" value in each cell of the grid.
  CA::CellBuffReal  ELV(GRID);

  // Create the water depth buffer.
  // It contains a "real" value in each cell of the grid.
  CA::CellBuffReal  WD(GRID);

  // Set the border of the elevation buffer to be no data. 
  ELV.bordersValue(borders,agtmp1.nodata);

  // Se the default value of the elevation to be nodata.
  ELV.fill(fulldomain, agtmp1.nodata);

  // Load the data not from the DEM file but from the pre-processed
  // file.
  if(!ELV.loadData(setup.preproc_name+"_ELV","0") )
  {
    std::cerr<<"Error while loading the Elevation pre-processed file"<<std::endl;
    return 1;
  }

  if(setup.output_console)
    std::cout<<"Loaded Elevation data"<< std::endl;

  // ----  CELL BUFFERS ----
    
  // Create two temporary Cell buffer
  CA::CellBuffReal TMP1(GRID);
  CA::CellBuffReal TMP2(GRID);

  //CA::EdgeBuffReal TMP3(GRID);
  CA::EdgeBuffReal* TMP3_EDGE1;
  CA::EdgeBuffReal* TMP3_EDGE2;

  // Create the MASK cell buffer. The mask is usefull to check wich
  // cell has data and nodata and which cell has neighbourhood with
  // data.
  CA::CellBuffState MASK(GRID);  
  
  // ----  SCALAR VALUES ----

  CA::Unsigned iter      = 0;		        // The actual iteration number.
  CA::Real     t         = setup.time_start;    // The actual time in seconds
  CA::Real     t_nearest = setup.time_end;      // The next time step
  CA::Real     nodata    = agtmp1.nodata; // capture terrain no data value, to display actual no data cells

  // set output no_data value to zero
  agtmp1.nodata = 0;

  // -- CREATE FULL MASK ---
  
  CA::createCellMask(fulldomain,GRID,ELV,MASK,nodata);  

  // ----  INIT RASTER GRID ----

  // List of raster grid data
  std::vector<RGData> rgdatas(rgs.size());

  // Peak buffers.
  RGPeak rgpeak;

  // Allocate temporary Ascii file data to ouptut raster.
  if(rgs.size()>0)
  {
    agtmp1.data.resize(agtmp1.ncols * agtmp1.nrows, nodata);
    agtmp2.data.resize(agtmp2.ncols * agtmp2.nrows, nodata);
  }

  for(size_t i = 0; i<rgs.size(); ++i)
  {
    std::string filename = ad.output_dir+ad.sdir+setup.short_name+"_"+setup.rastergrid_files[i]; 
    initRGData(filename, GRID, nodata, rgs[i], rgdatas[i], rgpeak);
  }

  // ---- CONTAINER DATA TO REMOVE ---

  // Create the container with the data to remove.
  typedef std::vector< std::pair<std::string, std::string> > RemoveID;
  RemoveID removeIDsCB;
  RemoveID removeIDsEB;


  // Indicate that the water depth buffer was loaded.
  bool WDloaded  = false;

  // ------------------------- MAIN LOOP -------------------------------
  while(t<=setup.time_end)
  {
    // Need to realod WD for this time step.
    WDloaded  = false;

    // -------  OUTPUTS --------    

    // Retrieve the string of the time.
    std::string strtime;
    CA::toString(strtime,t);
    
    // ----  RASTER GRID ----
    for(size_t i = 0; i<rgdatas.size(); ++i)
    {
      // Check if it is time to process raster grid!
      // Or the final extend need to be processed
      if(t >= rgdatas[i].time_next || (rgs[i].final && t==setup.time_end) )
      {
	if(!WDloaded)
	{
	  // The water depth buffer is practically always needed (for WD/WL and VEL).
	  // Reset the buffer
	  WD.fill(fulldomain,nodata);
	  
	  // Load the water depth data.
	  if(!WD.loadData(setup.short_name+"_WD",strtime) )
	  {
	    std::cerr<<"Missing water depth data: "<<strtime<<std::endl;
	    continue;
	  }
	  
	  // Set the water depth to zero if it less than tollerance.
	  CA::State boundary = (setup.rast_boundary) ? 1 : 0;
	  CA::Execute::function(fulldomain, zeroedWD, GRID, WD, MASK, setup.rast_wd_tol, boundary);
	
	  // Add the ID to remove.
	  removeIDsCB.push_back(std::make_pair(setup.short_name+"_WD",strtime));

	  WDloaded=true;
	}
	// Perform the action depending on the type of variable.
	switch(rgs[i].pv)
	{
	case PV::VEL:
	  {	  
	    std::string filenameV;
	    std::string filenameA;

		std::string filenameEV;
	    
	    // Create the name of the files if it is going to output a
	    // velocity field or simply rasters.
	    if(setup.rast_vel_as_vect == true)
	    {
	      //filenameV = removeExtension(rgdatas[i].filename) + "_" + strtime + ".csv";
	      //filenameA = removeExtension(rgdatas[i].filename) + "_" + strtime + ".csvt";
			// note the file name contains no extension in the Dll framework
			filenameV = rgdatas[i].filename + "_" + strtime + ".csv";
			filenameA = rgdatas[i].filename + "_" + strtime + ".csvt";
	    }
	    else
	    {
	      //filenameV = removeExtension(rgdatas[i].filename) + "_V_" + strtime + ".asc";
	      //filenameA = removeExtension(rgdatas[i].filename) + "_A_" + strtime + ".asc";
			//filenameV = rgdatas[i].filename + "_V_" + strtime + ".asc";
			//filenameA = rgdatas[i].filename + "_A_" + strtime + ".asc";

			filenameV = removeExtension(rgdatas[i].filename) + "_V_" + strtime;
			filenameA = removeExtension(rgdatas[i].filename) + "_A_" + strtime;
			filenameEV = removeExtension(rgdatas[i].filename) + "_EV_" + strtime;
	    }
	    
	    // Reset the buffer
	    TMP1.fill(fulldomain,nodata);
	    TMP2.fill(fulldomain,nodata);

		// declare the necassary memory
		TMP3_EDGE1 = new CA::EdgeBuffReal(GRID);
		TMP3_EDGE2 = new CA::EdgeBuffReal(GRID);


		TMP3_EDGE1->clear(nodata);
		TMP3_EDGE2->clear(nodata);

		// problems with fill command in hex version
		//TMP3_EDGE1->fill(realdomain, nodata);
		//TMP3_EDGE1->fill(fulldomain, nodata);
		//TMP3_EDGE2->fill(fulldomain, nodata);
	    
	    // Load the velocity data on TMP1.
	    if(! TMP1.loadData(setup.short_name+"_V",strtime) )
	    {
	      std::cerr<<"Missing velocity data to create file: "<<filenameV<<std::endl;
	      continue;
	    }
	    
	    // Load the angle data on TMP2.
	    if(! TMP2.loadData(setup.short_name+"_A",strtime) )
	    {
	      std::cerr<<"Missing angle data to create file: "<<filenameA<<std::endl;
	      continue;
	    }


		


		// Load the angle data on TMP3.
		//if (!TMP3_EDGE1->loadData(setup.short_name + "_EV", strtime))
		//{
		//	std::cerr << "Missing edge velocity data to create file: " << filenameEV << std::endl;
		//	continue;
		//}

		// calculate the velocities from 
		CA::Execute::function(fulldomain, zeroedVA, GRID, TMP1, TMP2, WD, MASK, setup.rast_wd_tol);
	    
	    
	    // Retrieve the velocity and angle data
	    TMP1.retrieveData(realbox, &agtmp1.data[0], agtmp1.ncols, agtmp1.nrows);	  
	    TMP2.retrieveData(realbox, &agtmp2.data[0], agtmp2.ncols, agtmp2.nrows);


		//std::ofstream outputStream;
		//outputStream.open(filenameEV);
		//if (setup.output_console)
		//	std::cout << "Write Edge fluxes Grid: " << filenameEV << std::endl;
		//
		//TMP3_EDGE1->dump(outputStream, ",");
		//outputStream.close();


		CA::Real dt = rgs[i]._timeStepList[strtime];

		
		// Set the V and A to zero if water depth is less than tollerance.
		CA::Execute::function(fulldomain, flux2velocity, GRID,
			WD, ELV, 
			(*TMP3_EDGE1), (*TMP3_EDGE2), 
			MASK,
			dt
			//setup.rast_wd_tol
			);

		//std::ofstream outputStream2;
		//outputStream2.open(filenameEV);
		//if (setup.output_console)
		//	std::cout << "Write Edge velocity Grid: " << filenameEV << std::endl;

		//TMP3_EDGE2->dump(outputStream2, ",");
		//outputStream2.close();


		//removeIDsEB.push_back(std::make_pair(setup.short_name + "_EV", strtime));

		delete(TMP3_EDGE1);
		delete(TMP3_EDGE2);

	    
	    // check if we need to output a velocity field
	    if(setup.rast_vel_as_vect)
	    {
	      if(setup.output_console)
		std::cout<<"Write Raster Grid: "<<filenameV<<std::endl;
	      
	      // Create the CSV file
	      //std::ofstream file( filenameV.c_str() ); 
	      FILE* fout = fopen(filenameV.c_str(), "w");	      

	      // Set  manipulators 
	      //file.setf(std::ios::fixed, std::ios::floatfield);
	      //file.precision(6); 
	      
	      // Write the header
	      //file<<"X, Y, Speed, Angle_RAD, Angle_DEG, Angle_QGIS"<<std::endl;
	      fprintf(fout,"X, Y, Speed, Angle_RAD, Angle_DEG, Angle_QGIS\n");
	     
	      // Loop thourhg the grid points.
	      for(CA::Unsigned j_reg=realbox.y(), j_mem=0; j_reg<realbox.h()+realbox.y(); ++j_reg, ++j_mem)
	      {
		for(CA::Unsigned i_reg=realbox.x(), i_mem=0; i_reg<realbox.w()+realbox.x(); ++i_reg, ++i_mem)
		{
		  // Create the point and find the coordinates.
		  CA::Point p(i_reg,j_reg);
		  p.setCoo(GRID);
		  
		  // Retrieve the speed and angle values (in radians),
		  // the compute the angle value in degrees.
		  CA::Real V  = agtmp1.data[j_mem * agtmp1.ncols + i_mem];
		  CA::Real AR = agtmp2.data[j_mem * agtmp2.ncols + i_mem];
		  CA::Real AD = AR*180/PI;
		  
		  // Compute the angle needed by QGis.
		  CA::Real AQ = -AR*180/PI + 90; 		
		  
		  // Write the results if the veclocity is more than zero.
		  if(V>0)
		  {
		    //file<<p.coo().x()<<","<<p.coo().y()<<","<<V<<","<<AR<<","<<AD<<","<<AQ<<","<<std::endl;
		    fprintf(fout,"%.12f,%.12f,%.6f,%.6f,%.6f,%.6f,\n", p.coo().x(),p.coo().y(),V,AR,AD,AQ);
		  }
		}
	      }        
	      
	      // Close the file.
	      //file.close();	      
	      fclose(fout);
	    }
	    // NOPE, simply output the rasters.
	    else
	    {
	      if(setup.output_console)
		std::cout<<"Write Raster Grid: "<<filenameV<<" "<<filenameA<<std::endl;
	      
	      // Write the  data.
	      //CA::writeAsciiGrid(agtmp1,filenameV,setup.rast_places);	      
	      //CA::writeAsciiGrid(agtmp2,filenameA,setup.rast_places);	      	  
		  agtmp1.writeAsciiGrid(filenameV, setup.rast_places);
		  agtmp2.writeAsciiGrid(filenameA, setup.rast_places);
	    }
	    
	    // Add the ID to remove.
	    removeIDsCB.push_back(std::make_pair(setup.short_name+"_V",strtime));
	    removeIDsCB.push_back(std::make_pair(setup.short_name+"_A",strtime));
	  }
	  break;
	  
	case PV::WL:
	  {
	    // Create the name of the file.
	    //std::string filename = removeExtension(rgdatas[i].filename) + "_" + strtime + ".asc";	
		  //std::string filename = rgdatas[i].filename + "_" + strtime + ".asc";
		  std::string filename = removeExtension(rgdatas[i].filename) + "_" + strtime;

	    // Reset the buffer
	    TMP1.fill(fulldomain,nodata);

	    // Make the water depth data and elevation into the water level.
	    CA::Execute::function(fulldomain, makeWL, GRID, TMP1, WD, ELV, MASK);
	    
	    // Retrieve the data
	    TMP1.retrieveData(realbox, &agtmp1.data[0], agtmp1.ncols, agtmp1.nrows);	  
	    	    
	    if(setup.output_console)
	      std::cout<<"Write Raster Grid: "<<filename<<std::endl;
	    
	    // Write the data.
	    //CA::writeAsciiGrid(agtmp1,filename,setup.rast_places);	      
		agtmp1.writeAsciiGrid(filename, setup.rast_places);
	  }
	  break;
	case PV::WD:
	  {
	    // Create the name of the file.
	    //std::string filename = removeExtension(rgdatas[i].filename) + "_" + strtime + ".asc";	
		  //std::string filename = rgdatas[i].filename + "_" + strtime + ".asc";
		  std::string filename = removeExtension(rgdatas[i].filename) + "_" + strtime;

	    // Water depth already loaded.
	    
	    // Retrieve the data
	    WD.retrieveData(realbox, &agtmp1.data[0], agtmp1.ncols, agtmp1.nrows);	  
	    	    
	    if(setup.output_console)
	      std::cout<<"Write Raster Grid: "<<filename<<std::endl;
	    
	    // Write the data.
	    //CA::writeAsciiGrid(agtmp1,filename,setup.rast_places);    
		agtmp1.writeAsciiGrid(filename, setup.rast_places);
	  }
	  break;
	default:
	  break;
	}
	// Update the next time to process a raster grid.
	rgdatas[i].time_next += rgs[i].period;
      }

      // Its time to find the possible next time of interest.
      t_nearest = nextTimeNearestAction(t,t_nearest, rgs[i].period);

    }

    //Finish after t reach end
    if(t == setup.time_end)
      break;

    // Set the nearest important time as the enxt time and set the
    // possible nearest important time as the end of the simulation.
    t         = t_nearest;
    t_nearest = setup.time_end;   
  }

  // ----  PEAK RASTER GRID ----

  // Need to realod WD only once for PEAK.
  WDloaded  = false;
  
  for(size_t i = 0; i<rgdatas.size(); ++i)
  {
    // Check if the peak values need to be saved.
    if(rgs[i].peak == true)
    {

      if(!WDloaded)
      {
	// The water depth buffer is practically always needed (for WD/WL and VEL).
	// Reset the buffer
	WD.fill(fulldomain,nodata);
	
	// Load the water depth data.
	if(! WD.loadData(setup.short_name+"_WD","PEAK") )
	{
	  std::cerr<<"Missing water depth data: "<<"PEAK"<<std::endl;
	  continue;
	}

	// Set the water depth to zero if it less than tollerance.
	CA::State boundary = (setup.rast_boundary) ? 1 : 0;
	CA::Execute::function(fulldomain, zeroedWD, GRID, WD, MASK, setup.rast_wd_tol, boundary);

	// Add the ID to remove.
	removeIDsCB.push_back(std::make_pair(setup.short_name+"_WD","PEAK"));

	// Only once.
	WDloaded = true;
      }
	
      switch(rgs[i].pv)
      {
      case PV::VEL:
	{
	  // Create the name of the file.
	  //std::string filenameV = removeExtension(rgdatas[i].filename) + "_V_PEAK.asc";	
		//std::string filenameV = rgdatas[i].filename + "_V_PEAK.asc";
		std::string filenameV = removeExtension(rgdatas[i].filename) + "_V_PEAK";
	  
	  // Reset the buffer
	  TMP1.fill(fulldomain,nodata);
	  
	  // Load the data on TMP1.
	  if(! TMP1.loadData(setup.short_name+"_V","PEAK") )
	  {
	    std::cerr<<"Missing velocity data to create file: "<<filenameV<<std::endl;
	    continue;
	  }

	  // Set the V and A to zero if water depth is less than tollerance.
	  CA::Execute::function(fulldomain, zeroedVA, GRID, TMP1, TMP2, WD, MASK, setup.rast_wd_tol);
	  	  
	  // Retrieve the data
	  TMP1.retrieveData(realbox, &agtmp1.data[0], agtmp1.ncols, agtmp1.nrows);	  
	  
	  if(setup.output_console)
	    std::cout<<"Write Raster Grid: "<<filenameV<<std::endl;
	  
	  // Write the data.
	  //CA::writeAsciiGrid(agtmp1, filenameV, setup.rast_places);
	  agtmp1.writeAsciiGrid(filenameV, setup.rast_places);

	  // Add the ID to remove.
	  removeIDsCB.push_back(std::make_pair(setup.short_name+"_V","PEAK"));
	}
	break;

      case PV::WL:
	{
	  // Create the name of the file.
	  //std::string filename = removeExtension(rgdatas[i].filename) + "_PEAK.asc";	
		//std::string filename = rgdatas[i].filename + "_PEAK.asc";
		std::string filename = removeExtension(rgdatas[i].filename) + "_WLPEAK";
	  
	  // Reset the buffer
	  TMP1.fill(fulldomain,nodata);
	  	  

	  // Make the water depth data and elevation into the water level.
	  CA::Execute::function(fulldomain, makeWL, GRID, TMP1, WD, ELV, MASK);
	  
	  // Retrieve the data
	  TMP1.retrieveData(realbox, &agtmp1.data[0], agtmp1.ncols, agtmp1.nrows);	  
	  
	  if(setup.output_console)
	    std::cout<<"Write Raster Grid: "<<filename<<std::endl;
	  
	  // Write the data.
	  //CA::writeAsciiGrid(agtmp1,filename,setup.rast_places);	 
	  agtmp1.writeAsciiGrid(filename, setup.rast_places);
	}
	break;
	
      case PV::WD:
	{
	  // Create the name of the file.
	  //std::string filename = removeExtension(rgdatas[i].filename) + "_PEAK.asc";	
		//std::string filename = rgdatas[i].filename + "_PEAK.asc";
		std::string filename = removeExtension(rgdatas[i].filename) + "_PEAK";

	  // Water depth already loaded.
	  
	  // Retrieve the data
	  WD.retrieveData(realbox, &agtmp1.data[0], agtmp1.ncols, agtmp1.nrows);	  
	  
	  if(setup.output_console)
	    std::cout<<"Write Raster Grid: "<<filename<<std::endl;
	  
	  // Write the data.
	  //CA::writeAsciiGrid(agtmp1,filename,setup.rast_places);      
	  agtmp1.writeAsciiGrid(filename, setup.rast_places);
	}
	break;
      default:
	break;      
      }      
    }

  }

  

 



  //-------------------------------------------------

  if(setup.remove_data)
  {
    for(size_t i=0; i<removeIDsCB.size(); i++)
    {
      std::pair <std::string,std::string>& ID = removeIDsCB[i];; 
      CA::CellBuffReal::removeData(ad.data_dir,ID.first,ID.second);
    }    
    for(size_t i=0; i<removeIDsEB.size(); i++)
    {
      std::pair <std::string,std::string>& ID = removeIDsEB[i];; 
      CA::EdgeBuffReal::removeData(ad.data_dir,ID.first,ID.second);
    }    
  }

  if(setup.remove_prec_data)
  {
    // Remove Elevation data.
    CA::CellBuffReal::removeData(ad.data_dir,setup.preproc_name+"_ELV","0");
    
    // Remove Grid data.
    CA::Grid::remove(ad.data_dir,setup.preproc_name+"_Grid","0");

    if(CA::CellBuffReal::existData(ad.data_dir,setup.preproc_name+"_RATIO","0"))
      CA::CellBuffReal::removeData(ad.data_dir,setup.preproc_name+"_RATIO","0");

    if(CA::CellBuffState::existData(ad.data_dir,setup.preproc_name+"_ROUGH","0"))
      CA::CellBuffState::removeData(ad.data_dir,setup.preproc_name+"_ROUGH","0");

    if(CA::CellBuffState::existData(ad.data_dir,setup.preproc_name+"_INF","0"))
      CA::CellBuffState::removeData(ad.data_dir,setup.preproc_name+"_INF","0");
  }

  if(setup.output_console)
    std::cout<<"Cleaned data"<< std::endl;


  if(setup.output_computation)
  {
    std::cout<<"-----------------" << std::endl; 
    std::cout<<"Total run time taken (s) = " << total_timer.millisecond()/1000.0 << std::endl;
    std::cout<<"-----------------" << std::endl; 
  }

  
  return 0;
}

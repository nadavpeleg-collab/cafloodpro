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

//! \file terrainInfo.cpp
//! Return the terrain info of the CA 2D model to simulate.
//! contact: a.s.chen [at] exeter.ac.uk
//! \date 2024-05

#include"ca2D.hpp"
#include"Masks.hpp"
#include"ArgsData.hpp"
#include"Setup.hpp"

#include CA_2D_INCLUDE(setBoundaryEle)
#include CA_2D_INCLUDE(computeArea)
#include CA_2D_INCLUDE(computeAreaGEO)
#include CA_2D_INCLUDE(computeDataCells)
#include CA_2D_INCLUDE(computeSlope)
#include CA_2D_INCLUDE(computeSlopeGEO)

// Return the terrrain info of the CA2D model to simulate. 
// \attention The preProc function should be called before this one.
// \warning If "Remove Pre-proc data is true, this function remove them."
int terrainInfo(const ArgsData& ad,Setup& setup, const CA::AsciiGrid<CA::Real>& eg)
{
  
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

  // Print the Grid Information.
  if(setup.output_console)
  {
    std::cout<<"-----------------" << std::endl; 
    GRID.printInfo(std::cout);
    std::cout<<"-----------------" << std::endl; 
  }

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

  
  // -- INITIALISE ELEVATION ---
  
  // Create the elevation cell buffer.
  // It contains a "real" value in each cell of the grid.
  CA::CellBuffReal  ELV(GRID);

  // Set the border of the elevation buffer to be no data. 
  ELV.bordersValue(borders,eg.nodata);

  // Se the default value of the elevation to be nodata.
  ELV.fill(fulldomain, eg.nodata);

  // Load the data not from the DEM file but from the pre-processed
  // file.
  if(!ELV.loadData(setup.preproc_name+"_ELV","0") )
  {
    std::cerr<<"Error while loading the Elevation pre-processed file"<<std::endl;
    return 1;
  }

  // Highest elevation
  CA::Real     high_elv = 90000;

  // Find highest elevation
  ELV.sequentialOp(fulldomain, high_elv, CA::Seq::Max);

  if(setup.output_console)
  {
    std::cout<<"Loaded Elevation data"<< std::endl;
    std::cout<<"Highest elevation = "<< high_elv<<std::endl;
  }


  // ----  CELL BUFFERS ----
    
  // Create temporary buffers
  CA::CellBuffReal TMP1(GRID);

  // Create the MASK cell buffer. The mask is usefull to check wich
  // cell has data and nodata and which cell has neighbourhood with
  // data.
  CA::CellBuffState MASK(GRID);
   
  // If the we are using the Geographic coordinate system create a GEO
  // RATIO buffer and load it from the pre/processed data. 
  cpp11::shared_ptr<CA::CellBuffReal> PRATIO; 


  if(setup.latlong)
  {
    PRATIO.reset( new CA::CellBuffReal(GRID) );
    
    // Load the data not from the DEM file but from the pre-processed
    // file.
    if(!(*PRATIO).loadData(setup.preproc_name+"_RATIO","0") )
    {
      std::cerr<<"Error while loading the Geographic Coordinate System pre-processed data"<<std::endl;
      return 1;
    }    
    if(setup.output_console)
    {
      std::cout<<"Loaded Geographic Coordinate System data"<< std::endl;
    }
  }
  

  // ----  SCALAR VALUES ----

  CA::Real     nodata = eg.nodata;

  // The total number of cells and the total area of the domain.
  CA::Real    total_cells = 0.0;
  CA::Real    total_area  = 0.0;

  // The average slope in percent.
  CA::Real    avg_slope   = 0.0;

  // The erth radious
  CA::Real R         = 6371; // Km
  
  // The average dx and dy of the grid.  
  CA::Real dy        = GRID.length();
  CA::Real dx        = GRID.length();

  // These two values need to be calculated using the Equirectangular
  // approximation if we are using Geographic coordinate system
  if(setup.latlong)
  {   
    dy = std::abs( GRID.length()*(PI/180.0)*R*1000 );

    // Find the middle point of the grid.
    CA::Real midcoo =  GRID.yCoo()+GRID.yNum()/2*GRID.length();

    dx = std::abs( (GRID.length()*(PI/180.0))* std::cos(midcoo*(PI/180.0)) * R * 1000 );
    
    if(setup.output_console)
    {
      std::cout<<"--------------------------------------------------------" << std::endl; 
      std::cout<<"Using Geographic Coordinate System"<< std::endl;
      std::cout<<"-----------------" << std::endl; 
      std::cout<<"DX : "<<dx<<std::endl;
      std::cout<<"DY : "<<dy<<std::endl;
      std::cout<<"--------------------------------------------------------" << std::endl; 
    }    
  }

  // -- CREATE FULL MASK ---
  
  CA::createCellMask(fulldomain,GRID,ELV,MASK,nodata);
  
  //***********************//
  //* BOUNDARY DISCUSSION *//
  //***********************//

  // The MASK is used to check if a cell with nodata has a neighbour
  // with data (Bit 31 set to true in the MASK). This kind of cells
  // are called boundary cells and they are ignored by the computation
  // (not visited), i.e. the rain is not added into them, the outflow
  // are not computed. However, these boundary cells are used as
  // neighbour cell thus the bouandary cell can have inflow (outflow
  // from the data cell). Thus the elevation value in these bounday
  // cell is changed with the Boundary Elevation value. If this value
  // is very high the global baundary is a CLOSED one, if it is VERY
  // negative the global boundary is an OPEN one.

  // ATTENTION The water that finish in the boundary cell (open
  // boundary case) is not removed (it stay in the WD buffer). 

  // Set the boundary cell elevation to the given boundary_elv value.
  CA::Execute::function(fulldomain, setBoundaryEle, GRID, ELV, MASK, setup.boundary_elv);

  //CA_DUMP_BUFF(ELV,0);  
 
  // -- SLOPE, TOTAL AREA and TOTAL CELLS  ---

  // Compute the slope.
  TMP1.clear();
  if(setup.latlong)
    CA::Execute::function(fullbox, computeSlopeGEO, GRID, TMP1, ELV, MASK, (*PRATIO),dy);    
  else
    CA::Execute::function(fullbox, computeSlope, GRID, TMP1, ELV, MASK);    
  TMP1.sequentialOp(fullbox, avg_slope, CA::Seq::Add);
    
  // Compute the total area,
  TMP1.clear();
  if(setup.latlong)
    CA::Execute::function(fullbox, computeAreaGEO, GRID, TMP1, MASK,(*PRATIO),dy);    
  else
    CA::Execute::function(fullbox, computeArea, GRID, TMP1, MASK);    
  TMP1.sequentialOp(fullbox, total_area, CA::Seq::Add);
  
  // Compute the total number of data cells.
  TMP1.clear();
  CA::Execute::function(fullbox, computeDataCells, GRID, TMP1, MASK);    
  TMP1.sequentialOp(fullbox, total_cells, CA::Seq::Add);
    
  avg_slope/=total_cells;
  
  if(setup.output_console)
  {
    std::cout<<"--------------------------------------------------------" << std::endl;   
    std::cout<<"Total number of data cells : "<<total_cells<<std::endl;
    std::cout<<"Total area of data cells   : "<<total_area<<std::endl;
    std::cout<<"Average area of data cell  : "<<total_area/total_cells<<std::endl;
    std::cout<<"Average % slope            : "<<avg_slope*100<<std::endl;
    std::cout<<"--------------------------------------------------------" << std::endl;   
  }

  // Remove pre-proc data if requested

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

  

  // ---- TIME OUTPUT ----

  if(setup.output_console && setup.output_computation)
  {
    std::cout<<"-----------------" << std::endl; 
    std::cout<<"Total run time taken (s) = " << total_timer.millisecond()/1000.0 << std::endl;
    std::cout<<"-----------------" << std::endl; 
  }
  
  return 0;
}





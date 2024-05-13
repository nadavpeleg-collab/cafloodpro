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

//! \file CADDIES2D_2.cpp
//! Perform the CADDIES2D flood modelling
//! contact: a.s.chen [at] exeter.ac.uk
//! \date 2024-05


#include"ca2D.hpp"
#include"Masks.hpp"
#include"ArgsData.hpp"
#include"Setup.hpp"
#include"Rain.hpp"
#include"Inflow.hpp"
#include"WaterLevel.hpp"
#include"TimePlot.hpp"
#include"RasterGrid.hpp"
#include"TSPlot.hpp"
#include"SpatialTemporal.hpp"
#include<ctime>
#include"SpatialTemporalDense.hpp"



//-----------------------------//
// Include the CA 2D functions //
//-----------------------------//
#include CA_2D_INCLUDE(geoRatio)
#include CA_2D_INCLUDE(addHeight)

#include CA_2D_INCLUDE(outflowWCA2Dv1)
#include CA_2D_INCLUDE(waterdepthWCA2Dv1)
#include CA_2D_INCLUDE(velocityWCA2Dv1)
#include CA_2D_INCLUDE(infiltration)

#include CA_2D_INCLUDE(setBoundaryEle)
#include CA_2D_INCLUDE(removeUpstr)
#include CA_2D_INCLUDE(updatePEAKC)
#include CA_2D_INCLUDE(updatePEAKE)

#include CA_2D_INCLUDE(outflowWCA2Dv1GEO)
#include CA_2D_INCLUDE(waterdepthWCA2Dv1GEO)
#include CA_2D_INCLUDE(velocityWCA2Dv1GEO)
#include CA_2D_INCLUDE(infiltrationGEO)

#include CA_2D_INCLUDE(outflowWCA2Dv2)

#include CA_2D_INCLUDE(computeArea)
#include CA_2D_INCLUDE(computeDataCells)
#include CA_2D_INCLUDE(computeSlope)

// LISFLOOD like ca
#include CA_2D_INCLUDE(outflowDiffusive)
#include CA_2D_INCLUDE(velocityDiffusive)
#include CA_2D_INCLUDE(velocityDiffusivev2)

#include CA_2D_INCLUDE(outflowInertial)
#include CA_2D_INCLUDE(velocityInertial)

// Generic waterdepth 
#include CA_2D_INCLUDE(waterdepth)

// time grid outouts
#include CA_2D_INCLUDE(timeToThreshold)

// banding system
#include CA_2D_INCLUDE(banding)

// check for neighbouring no data on flood
#include CA_2D_INCLUDE(checkFloodNoDataNeighbourhood)

#include CA_2D_INCLUDE(zeroedWD)
#include CA_2D_INCLUDE(zeroedVA)
#include CA_2D_INCLUDE(makeWL)


// WCA2Dv2 model base on:
// CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
//                                                                                            
// Paper name : A WEIGHTED CELLULAR AUTOMATA 2D 
//              INUNDATION MODEL FOR RAPID FLOOD ANALYSIS 
// Authors    : M. Guidolin, A. S. Chen, B. Ghimire, E. Keedwell, 
//              S. Djordjević, D.A. Savić 
//
// CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


/* Round a number to the given digits.*/
inline CA::Real fround(CA::Real n, unsigned d)
{
  return std::floor(n * std::pow(static_cast<CA::Real>(10.0), static_cast<CA::Real>(d)) 
		    + static_cast<CA::Real>(.5)) 
    / std::pow(static_cast<CA::Real>(10.0), static_cast<CA::Real>(d) );
}


/* Extend a domain */
inline CA::Box extendBox(CA::Box extent, CA::Box& fullbox, CA::Unsigned lines, const char* borders)
{
  // Check that the new box is not too big by using full.
  // Check which border need to be extend by using the borders array.
  // !!!!! ATTENTION!!!! THIS WORK ONLY FOR a four neighbour square-grid

  // Width and heigh line multiplier which depend on many directions need to be extend.
  CA::Unsigned wmul = 0;
  CA::Unsigned hmul = 0;

  
  //std::cout<<"BORDERS  = ";
  //for(int i=0; i<=caNeighbours; i++)
  //  std::cout<<(int)borders[i];
  //std::cout<<std::endl;
  

  // ATTENTION There is a problem here when the lines is more than 1
  // and the X and Y coordinate reach the fullbox. The W and H
  // continue to grow also if they shouldn't.
  
  // Extend reft
  if(borders[3]>0)
  {
    if(extent.x()>fullbox.x()+lines) 
      extent.setX(extent.x()-lines);
    else
      extent.setX(fullbox.x());
    
    wmul++;
  }

  // Extend up
  if(borders[2]>0)
  {
    if(extent.y()>fullbox.y()+lines) 
      extent.setY(extent.y()-lines);	
    else
      extent.setY(fullbox.y());
    
    hmul++;
  }

  // Extend right
  if(borders[1]>0)
    wmul++;

  // Extend down
  if(borders[4]>0)
    hmul++;

  // Extend width and heigh depending on the active directions 
  extent.setW(std::min((extent.w()+lines*wmul),fullbox.w())); 
  extent.setH(std::min((extent.h()+lines*hmul),fullbox.h()));   
  return extent;
}

//! Output to console the information about the simulation.
void outputConsole(CA::Unsigned iter, CA::Unsigned oiter, CA::Real t, CA::Real dt, 
		   CA::Real avgodt, CA::Real minodt, CA::Real maxodt, 
		   CA::Real vamax, CA::Real upstr_elv, const CA::BoxList& domain,
		   const Setup& setup)
{
  std::cout<<"-----"<<std::endl;
  std::cout<<"Total iterations = "<<iter<<" Simulation time (MIN) = "<<t/60
	   <<" Last DT = "<<dt<<std::endl;
  std::cout<<"Last iterations  = "<<oiter <<" Average DT ="<<avgodt/oiter
	   <<" Min DT = "<<minodt<<" Max DT = "<<maxodt<<std::endl;	
  std::cout<<"UPSTRELV = "<<upstr_elv<<std::endl;
  std::cout<<"VAMAX    = "<<vamax<<std::endl;    
  if(setup.expand_domain)
  {
    CA::Box B(domain.extent());
    std::cout<<"DOMAIN   = ("<<B.x()<<","<<B.y()<<"):("<<B.w()<<","<<B.h()<<")"<<std::endl;
  }
  std::cout<<"-----"<<std::endl;

}


//! Compute the next time step as a fraction of the period time step.
//! Check that dt is between  min max.
void computeDT(CA::Real& dt, CA::Unsigned& dtfrac, CA::Real dtn1, const Setup& setup)
{
  CA::Unsigned dtfracmax = static_cast<CA::Unsigned>(setup.time_maxdt / setup.time_mindt);
  
  // Find the fraction of time_maxdt second that is just less than dtn1.
  // If dt is smaller than dtn1 we need to decrease the
  // fraction and stop as soon as the result is higher than dtn1.
  // If dt is bigger than dtn1 we need to increase the
  // fraction and stop as soon as the restuts is lower tha dtn1.
  if(dt<=dtn1)
  {
    for(; dtfrac>=1;dtfrac--)
    {
      CA::Real tmpdt = setup.time_maxdt/static_cast<CA::Real>(dtfrac);
      if(tmpdt>=dtn1)
	break;
    }
    if(dtfrac==1)
      dt = setup.time_maxdt;
    else
    {
      dtfrac+=1;
      dt = setup.time_maxdt/static_cast<CA::Real>(dtfrac);
    }
  }
  else
  {
    for(; dtfrac<=dtfracmax;dtfrac++)
    {
      CA::Real tmpdt = setup.time_maxdt/static_cast<CA::Real>(dtfrac);
      if(tmpdt<=dtn1)
	break;
    }
    dt = setup.time_maxdt/static_cast<CA::Real>(dtfrac);
  }

  // Check that dt is between  min max.
  dt = std::min(std::max(dt,setup.time_mindt), setup.time_maxdt);  
}

// Calculates log2 of number.  
// Needed for Visual Studio
double Log2( double n )  
{  
    // log(n)/log(2) is log2.  
    return log( n ) / log( 2.0f );  
}

// Find the number of Bits needed to store value
CA::Unsigned getNumberBits(CA::Unsigned v)
{
  return std::ceil( Log2( static_cast<double>(v) ) );
}

// find the wetted area of grid for initial processing domain
CA::Box findWettedArea(CA::AsciiGrid<CA::Real> wetGrid, CA::Grid& GRID)
{

	//CA::Box      realbox(GRID.box().x() + 1, GRID.box().y() + 1, GRID.box().w() - 2, GRID.box().h() - 2);

	int firstX = 0;
	int firstY = 0;

	int lastX = 0;
	int lastY = 0;

	bool foundFirstWetCell = false;

	for (int col = 0; col < wetGrid.ncols; ++col) {

		bool wetCellInRow = false;

		for (int row = 0; row < wetGrid.nrows; ++row) {

			// check if cell is wet
			if (wetGrid.data[col + (row*wetGrid.ncols)] > 0.0) {


				if (!foundFirstWetCell) {
					foundFirstWetCell = true;
					firstX = col;
					firstY = row;
				}
				else {
					if (col < firstX) {
						firstX = col;
					}
				}

				if (col > lastX) {
					lastX = col;
				}

				if (row > lastY) {
					lastY = row;
				}


			}// end of wet cell

		}// end of looping rows
	}// end of looping columns


	 //if (!foundFirstWetCell) {
	 //return NULL;
	 //}
	 //else {


	 //std::cout << "First X,Y: " << firstX << "," << firstY << std::endl;
	 //std::cout << "Last  X,Y: " << lastX  << "," << lastY  << std::endl;


	return CA::Box(GRID.box().x() + 1 + firstX, GRID.box().y() + 1 + firstY, lastX - firstX + 1, lastY - firstY + 1);

	//}

}// end of findWettedArea




//---------------------------------------------------------------------------------------------------------------------------









int CADDIES2D(const ArgsData& ad,Setup& setup, CA::AsciiGrid<CA::Real>& eg,
	      const std::vector<SpatialTemporal>& sts, 
	      const std::vector<RainEvent>& res, const std::vector<WLEvent>& wles, 
	      const std::vector<IEvent>& ies, 
	      const std::vector<TimePlot>& tps, std::vector<RasterGrid>& rgs,
		  CA::CellBuffReal**  timeGrid, double& thresholdDepth,
		  double& currentTime, CA::State& resultcode,
		  bool& abortIt,
          SpatialTemporalDense& _std,
			bool preProcessing
		  )
{
	if(preProcessing){
		// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// start - pre processing
		// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

		if (setup.output_computation)
		{
			std::cout << "Pre-processing : " << setup.sim_name << std::endl;
			std::cout << "------------------------------------------" << std::endl;
		}

		if (setup.output_console)
		{
			time_t t = time(0);   // get time now
			struct tm * now = localtime(&t);
			std::cout << "Simulation : " << setup.sim_name << std::endl;
			std::cout << "Model      : " << setup.model_type << std::endl;
			std::cout << "Date Start : " << (now->tm_year + 1900) << "-" << (now->tm_mon + 1) << '-' << now->tm_mday
				<< " " << now->tm_hour << ":" << now->tm_min << ":" << now->tm_sec << std::endl;
			std::cout << "------------------------------------------" << std::endl;
		}

		// -- ELEVATION FILE(s) --

		//CA::AsciiGrid<CA::Real> eg;

		if (setup.loadTerrain)
		{
			//AsciiGrid<CA::Real> eg;
            eg.readAsciiGrid(setup.elevation_ASCII);
			//eg.copyData(setup.eg);
		}
		else
		{
			//eg = setup.eg;
            eg.copyData(setup.eg);
		}
	}
  CA::Real outputNoData = -1;

  // Check the model.
  switch(setup.model_type)
  {
  case MODEL::WCA2Dv1 :
  case MODEL::WCA2Dv2 :
    // OK!
    break;
  case MODEL::Diffusive:
  case MODEL::Inertial :
    std::cout<<"-------------------------------------------------------------"<<std::endl;
    std::cout<<"ATTENTION: The model "<<setup.model_type<< " is EXPERIMENTAL!"<<std::endl;
    std::cout<<"-------------------------------------------------------------"<<std::endl;
    break;
  default:
    std::cerr<<"Error the simulation does not support the model: "<<setup.model_type<<std::endl;
    return 1;
  }

  if(setup.output_console)
  {
    time_t t = time(0);   // get time now
    struct tm * now = localtime( & t );
    std::cout<<"Simulation : "<<setup.sim_name<< std::endl;
    std::cout<<"Model      : "<<setup.model_type<< std::endl;
    std::cout<<"Date Start : "<<(now->tm_year + 1900)<<"-"<<(now->tm_mon + 1)<<'-'<<now->tm_mday
	     <<" "<<now->tm_hour<<":"<<now->tm_min<<":"<<now->tm_sec << std::endl;
    std::cout<<"------------------------------------------" << std::endl; 
  }
  
  // ----  Timer ----
  
  // Get starting time.
  CA::Clock total_timer;

  // ----  CA GRID ----

CA::Grid GRID;

if(preProcessing)
{
// -- EXTRA CELLS --

	// Given the way the boundary cells are computed (see later), we
	// need to add an extra set of cells in each direction.
	CA::AsciiGrid<CA::Real> grid;
	grid.ncols = eg.ncols + 2;
	grid.nrows = eg.nrows + 2;
	grid.xllcorner = eg.xllcorner - eg.cellsize; // ATTENTION HERE!
	grid.yllcorner = eg.yllcorner - eg.cellsize; // ATTENTION HERE!
	grid.cellsize = eg.cellsize;
	grid.nodata = eg.nodata;

	// ----  CA GRID ----

	// Create the square regular CA grid of the extended size of the DEM.
	// The internal implementation could be different than a square regular grid.
	//CA::Grid  GRID(grid.ncols, grid.nrows, grid.cellsize, grid.xllcorner, grid.yllcorner, ad.args.active());
	GRID = CA::Grid(grid.ncols, grid.nrows, grid.cellsize, grid.xllcorner, grid.yllcorner, ad.args.active());

	// Create two temporary Cell buffer
	//CA::CellBuffReal TMP1(GRID);
	//CA::CellBuffReal TMP2(GRID);

	// Set the place where the direct I/O data is going to be saved.
	GRID.setDataDir(ad.data_dir);


}
else{
  
  // Load the CA Grid from the DataDir. 
  // ATTENTION this should have an extra set of cells in each
  // direction.  The internal implementation could be different than a
  // square regular grid.
  //CA::Grid  GRID(ad.data_dir,setup.preproc_name+"_Grid","0", ad.args.active());
    GRID = CA::Grid(ad.data_dir, setup.preproc_name + "_Grid", "0", ad.args.active());

  if(setup.output_console)
    std::cout<<"Loaded Grid data"<< std::endl;
}

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
  // Create the computational domain, i.e the area where the water is
  // present and moving.
  CA::BoxList  compdomain;
  //CA::Box newDomainBox = CA::Box(compdomain.extent());
  //std::cout << "rank[" << GRID.caGrid().rank << "]newDomainBox.y():" << newDomainBox.y() << std::endl;
  //std::cout << "rank[" << GRID.caGrid().rank << "]newDomainBox.x():" << newDomainBox.x() << std::endl;
  //std::cout << "rank[" << GRID.caGrid().rank << "]newDomainBox.w():" << newDomainBox.w() << std::endl;
  //std::cout << "rank[" << GRID.caGrid().rank << "]newDomainBox.h():" << newDomainBox.h() << std::endl;
  // Create the computational domain used for sequential
  // compuattaion. SequentialOp does not work with multiple boxes.
  CA::BoxList  seqdomain(fullbox);
  
  // -- INITIALISE ELEVATION ---
  
  // Create the elevation cell buffer.
  // It contains a "real" value in each cell of the grid.
  CA::CellBuffReal  ELV(GRID);
  std::cout << "-----------------!" << std::endl;
  // create a pointer to a cell buffer for time output grid
  // in order to minimise memory usuage, using a pointer, which is only declared if appropriated option is selected
  //CA::CellBuffReal* timeGrid;
  //CA::CellBuffReal**  timeGrid;

  if (setup.timeGridOut){
	  
	  (*timeGrid) = new CA::CellBuffReal(GRID);

	  //(*timeGrid)->fill(fullbox, eg.nodata);
	  (*timeGrid)->fill(fullbox, outputNoData);


	  //std::cout << "blahhhhhhhhhhhhhhhhhhhh" << std::endl;
  }

  std::cout << "-----------------!!" << std::endl;

  // Set the border of the elevation buffer to be no data. 
  ELV.bordersValue(borders,eg.nodata);

  // Se the default value of the elevation to be nodata.
  ELV.fill(fulldomain, eg.nodata);

  std::cout << "-----------------!!!" << std::endl;

if(preProcessing)
{
	// Insert the DEM data into the elevation cell buffer. The first
	// parameter is the region of the buffer to write the data which is
	// the real domain, i.e. the original non extended domain.
	ELV.insertData(realbox, &eg.data[0], eg.ncols, eg.nrows);

	// The Spatial/Temporal variables vector must have at least three
	// elements. If there are more than three elements. The extra one
	// are ignored.
	if (sts.size()<3)
	{
		std::cerr << "Error while processing spatial temporal data" << std::endl;
		return 1;
	}

	if (setup.addHeights) {

		

		//std::cout << "Making height Adjustments...["<< setup.addedHeight.size() <<"]" << std::endl;

		CA::TableReal ADDED_HEIGHTS(GRID,setup.addedHeight.size());

		ADDED_HEIGHTS.update(0, setup.addedHeight.size(), &setup.addedHeight[0], setup.addedHeight.size());

		//ADDED_HEIGHTS.dump(std::cout);

		CA::CellBuffState BUILDING_MASK(GRID);
		//BUILDING_MASK.insertData(realbox, &sts[0].map.data[0], sts[0].map.ncols, sts[0].map.nrows);
		BUILDING_MASK.insertData(realbox, &(setup.heightAdjustmentMap->data[0]), setup.heightAdjustmentMap->ncols, setup.heightAdjustmentMap->nrows);
		CA::Execute::function(fulldomain, addHeight, GRID, ELV, BUILDING_MASK, ADDED_HEIGHTS);//////////////////////////////////////////////////////////

		//ELV.dump(std::cout);

		//delete now unneeded map
		delete(setup.heightAdjustmentMap);
	}

	// Compute the Geographic ratio if the coordinate system is a
	// Geographic one, i.e. lat/long
	if (setup.latlong)
	{
		if (setup.output_console)
			std::cout << "Prepare Geographic Coordinate System" << std::endl;

		// Create the Geographic ratio buffer.
		CA::CellBuffReal  RATIO(GRID);

		// The earth radio.
		CA::Real R = 6371;
		// Compute the dy in meters using Equirectangular
		// approximation. This should be constant at all latitude
		// longitude.
		CA::Real dy = std::abs(GRID.length()*(PI / 180.0)*R * 1000);
		// Compute the the lenght of a side in radian.
		CA::Real l = GRID.length()*(PI / 180.0);

		// Compute the geographic ratio.
		// ATTENTION This work only on a square grid.
		CA::Execute::function(fulldomain, geoRatio, GRID, RATIO, dy, R, l);

		
	}// end of if setup latlong


}
else
{
  // Load the data not from the DEM file but from the pre-processed
  // file.
  if(!ELV.loadData(setup.preproc_name+"_ELV","0") )
  {
    std::cerr<<"Error while loading the Elevation pre-processed file"<<std::endl;
    return 1;
  }
}
	

	

  std::cout << "-----------------!!!!" << std::endl;

  // Highest elevation
  CA::Real     high_elv = 90000;

  // Find highest elevation
  ELV.sequentialOp(fulldomain, high_elv, CA::Seq::Max);

  if(setup.output_console)
  {
    std::cout<<"Loaded Elevation data"<< std::endl;
    std::cout<<"Highest elevation = "<< high_elv<<std::endl;
  }


  // Decomposed the domain if requested.

  // ATTENTION Upstream reduction and Expand Domain are disable for
  // safety reasons.
  CA::BoxList dd;
  if(setup.decompose_domain)
  {
    CA::domainDecomposition(GRID.box(),ELV,eg.nodata,setup.decompose_eff,200,dd);
    setup.ignore_upstream = false;
    setup.expand_domain   = false;
  }

  // ----  CELL BUFFERS ----
    
  // Create  the water depth cell buffer.
  CA::CellBuffReal WD(GRID);
   
  // Create the MASK cell buffer. The mask is usefull to check wich
  // cell has data and nodata and which cell has neighbourhood with
  // data.

  // Attention The MASK is allocated in a pointer and linked into
  // another pointer due to the need of spatial temporal
  // variables. See the discossion later.
  cpp11::shared_ptr<CA::CellBuffState> PMASK1 (new CA::CellBuffState(GRID));
  cpp11::shared_ptr<CA::CellBuffState> PMASK2 = PMASK1;
  
  // Create the velocity cell buffer with the speed of the velocity
  CA::CellBuffReal V(GRID);

  // Create the velocity cell buffer with the angle of the velocity
  CA::CellBuffReal A(GRID);

  // If the we are using the Geographic coordinate system create a GEO
  // RATIO buffer and load it from the pre/processed data. 
  cpp11::shared_ptr<CA::CellBuffReal> PRATIO; 

  // If the we are using the WCA2Dv2 model create the buffer with the
  // dt value.
  cpp11::shared_ptr<CA::CellBuffReal> PDT; 
  

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

  // Allocate the buffer for the time step
  // Only WCA2Dv1 does not need this.
  if(setup.model_type != MODEL::WCA2Dv1)
  {
    PDT.reset( new CA::CellBuffReal(GRID) );
    (*PDT).fill(fulldomain, setup.time_updatedt);
  }

  // ----  EDGES BUFFERS ----

  // Create the outflow edge buffer(s). Many models uses a double buffer
  // technique to save time.
  CA::EdgeBuffReal OUTF1(GRID);
  CA::EdgeBuffReal OUTF2(GRID);

  CA::EdgeBuffReal *POUTF1 = &OUTF1;
  CA::EdgeBuffReal *POUTF2 = &OUTF2;

  // Allocate the buffer for total flow for  an edge for an update step.
  // Only WCA2Dv1 does need this.
  cpp11::shared_ptr<CA::EdgeBuffReal> PTOT; 
  if(setup.model_type == MODEL::WCA2Dv1)
  {
    PTOT.reset( new CA::EdgeBuffReal(GRID) );
    (*PTOT).clear();
  }

  // ---- ALARMS ----

  // Create the alarm(s) checked during the outflux computation.
  // Alarm 0: indicates when there is going to be an outflux outside of the computational domain of main cell.
  // Alarm 1: indicates that the outflux was on the edges toward neighbour cell 1.
  // Alarm 2: indicates that the outflux was on the edges toward neighbour cell 2.
  //...
  CA::Alarms  OUTFALARMS(GRID,caNeighbours+1);

  // Create the alarm(s) checked during the velocity computation
  // Alarm 0: indicates when there is still water movement over the elevation threshould.
  CA::Alarms  VELALARMS(GRID,1);


  // ----  SCALAR VALUES ----

  CA::Unsigned iter   = 0;		   // The actual iteration number.
  CA::Real     t      = setup.time_start;  // The actual time in seconds
  CA::Real     dt     = setup.time_maxdt;  // Starting delta time.
  CA::Real     dtn1   = 0.0;
  CA::Real     nodata = eg.nodata;


  currentTime = t;

  // The simulation time when the events will not produce any further water.
  CA::Real     t_end_events = setup.time_start; 

  // The level of water that can be ignored.
  CA::Real     ignore_wd  = setup.ignore_wd;

  // The water difference between cell that can be ignored.
  CA::Real     tol_delwl = setup.tolerance;

  // The minimum level of water that is used for computing the
  // velocity. 
  CA::Real     tol_va    = std::min(setup.ignore_wd,static_cast<CA::Real>(0.001));

  // The tollerance of the slope used to compute the dt. The slope
  // between two points must be higher than this value.
  // The tollerance is passed in percentile.
  CA::Real     tol_slope = setup.tol_slope/100;

  // This is the period of when the velocity is computed.
  CA::Real     period_time_dt   = setup.time_updatedt;         // Alias updatedt
  CA::Real     time_dt          = t + period_time_dt;          // The next simulation time when to update dt.
  CA::Unsigned iter_dt          = 0;			       // The number of iteration before next update dt.
  
  // This is the time from the start of an update dt.
  CA::Real     start_updatedt   = 0.0;
  
  // This is the previous update step dt.
  CA::Real     previous_dt      = dt;

  // The fraction of time_maxdt second used to compute the next the dt
  CA::Unsigned dtfrac    = 1;		  
  
  // The parameter of the time step used by various model.
  // The WCA2Dv2 use it only for the potential va and it is fixed.  
  CA::Real     alpha  = setup.time_alpha;

  // The default roughness
  CA::Real     rough = setup.roughness_global;

  CA::Unsigned oiter  = 0;                 // Output step num of iteration.
  CA::Real     minodt = setup.time_maxdt;  // Output step minimum dt.
  CA::Real     maxodt = 0.0;               // Output step maximum dt.
  CA::Real     avgodt = 0.0;	           // Average dt;
  CA::Real     time_output  = t + setup.output_period;	   // The time of the next output.

  CA::Real     rain_volume   = 0.0;
  CA::Real     inflow_volume = 0.0;
  CA::Real     wd_volume     = 0.0;
  CA::Real     inf_volume    = 0.0;

  // Maximum velocity.
  CA::Real     vamax=0.0;

  // Possible dt for WCA2Dv2 model. This is used to compute the
  // possible dt create by an event and created by the outflow.
  CA::Real     possible_dt;

  // Upstream elevation value. This value is the elevation where the
  // water "probably cannot reach anymore".
  CA::Real     upstr_elv = high_elv;

  // The potential velocity that an event has created. This is used by
  // the WCA2Dv1 model togher with alpha to compute a potential time
  // step created by an event.
  CA::Real     potential_va;

  // If true the outflow computation will check the box cell for a
  // possible expanding domain.
  bool         checkbox = setup.expand_domain;         

  // Variable which indicates if the PEAK values need to be updated.
  // The PEAK values can be updated on every update step or on every
  // steps.
  bool UpdatePEAK    = false;

  // Variable which indicates if the raster habe been saved in the lat
  // iteration.
  bool RGwritten     = false;

  // If true perform the infiltration step.
  bool useInfiltration  = false;
  
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
  
  CA::createCellMask(fulldomain,GRID,ELV,(*PMASK1),nodata);
  
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
  CA::Execute::function(fulldomain, setBoundaryEle, GRID, ELV, (*PMASK1), setup.boundary_elv);

  //CA_DUMP_BUFF(ELV,0);  

  //**********************************//
  //* SPATIAL TEMPORAL VAR DICUSSION *//
  //**********************************//

  // The spatial/temporal variables, like roughness, rain and
  // infiltration are usually composed of two 'items': 1) a
  // grid map which contains various indices, and 2) a table where
  // each index is assigned a value/time sequence. These two items can
  // be build by the CAAPI tool spatialIndexingAG.

  // Given that the number of different value/time sequences are known
  // the idea is to marge the three different maps in a single buffer
  // in order to save memory. Furthermore, to save even further
  // memory, these maps are included in the MASK buffer. This BUFFER
  // have various bites that are not used. So the idea is to use these
  // free bites.

  // ATTENTION, if the number of value/time sequences in the various
  // maps exceed the number of space available in the free bites of the
  // MASK. Then the mask get duplicated. The first MASK contains the
  // value/time sequences for the spatial roughness and it is used
  // into the CA function that computes the fluxes, the second MASK
  // contains the value/time sequence for the spatial rain and
  // infiltration and it is used into the CA function that update the
  // water depth. 

  // WARNING For the grid map, the nodata value should be considere as
  // index zero. Thus there cannot be a zero index into the CVS.  The
  // index zero for the roughness and infiltration is allocated to the
  // default values given into the Setup file.

  // This variable indicate the number of bits that are free (not
  // used) in a state variable of a cell of the mask. The State type
  // should be 32 bit (at least). The mask uses the first caNeighbours
  // + 1 bits to indicate the status of the neighbours cell and the
  // main cell. Furthermore, it uses the bit 30 and 31. So we need to
  // subtract these used bits from the number of bits available.
  int  mask_free_bits = 32 - (caNeighbours+1) - 2;

  // The Spatial/Temporal variables vector must have at least three
  // elements. If there are more than three elements. The extra one
  // are ignored.
  if(sts.size()<3)
  {
    std::cerr<<"Error while processing spatial temporal data"<<std::endl;
    return 1;    
  }

  // This array will contain in positions:
  // 0 and 1: The starting and Ending position of roughness bits in the MASK.
  // 2 and 3: The starting and Ending position of rain bits in the MASK.
  // 4 and 5: The starting and Ending position of infiltration bits in the MASK.
  std::vector<CA::State> bits_pos(6);

  // OK let start by processing the number of bits needed to store the
  // number of indices (value/time sequences) of each spatial/temporal
  // variable.
  int rough_bits = getNumberBits(sts[0].seqs.size());  
  int rain_bits  = getNumberBits(sts[1].seqs.size());  
  int inf_bits   = getNumberBits(sts[2].seqs.size());  
  
  // First we check if rain_bits and inf_bits can be stored inside a
  // single MASK by cheking the free bits. If this is not the case
  // (dooub it) Then there are too many elements in rain or infiltration.
  if((rain_bits+inf_bits)>mask_free_bits)
  {
    std::cerr<<"Error while processing spatial temporal data"<<std::endl;
    std::cerr<<"Too many values in the Spatial/Temporal rain or infiltration"<<std::endl;
    return 1;    
  }

  // Now we check if the rain_bits/inf_bits and rough bits can be
  // stored in the free_bits of the mask.
  if((rough_bits+rain_bits+inf_bits)<=mask_free_bits)
  {
    // YAY! We can use only one MASK buffer.
    // Set up the starting and stopping bits for each spatial/temporal
    // variable considering only one MASK.
    bits_pos[0] = (caNeighbours+1); 
    bits_pos[1] = bits_pos[0]+rough_bits;
    bits_pos[2] = bits_pos[1];	                  // The starting position of rain is the end of rough.
    bits_pos[3] = bits_pos[2]+rain_bits;
    bits_pos[4] = bits_pos[3];	                  // The starting position of inf is the end of rain.
    bits_pos[5] = bits_pos[4]+inf_bits;
  }
  else
  {
    // NAY! We need to spend more memory!!    

    // Allocate a new buffer into the second pointer of the mask and
    // copy the data from the original mask.
    PMASK2.reset(new CA::CellBuffState(GRID));
    (*PMASK2).copy((*PMASK1));  

    // Set up the starting and stopping bits for each spatial/temporal
    // variable considering rough for the first mask and rain/inf for
    // the second mask.
    bits_pos[0] = (caNeighbours+1); 
    bits_pos[1] = bits_pos[0]+rough_bits;
    bits_pos[2] = (caNeighbours+1); 
    bits_pos[3] = bits_pos[2]+rain_bits;
    bits_pos[4] = bits_pos[3];	                  // The starting position of inf is the end of rain.
    bits_pos[5] = bits_pos[4]+inf_bits;
  }
  
  // ----  TABLES ----

  //! CA function table with the spatial value for the roughness.
  //! \attention This should contains the inverse of the roughness, i.e. 1/n
  //! \warning For the Inertial model this contain the n^2 value.
  CA::TableReal ROUGH(GRID,sts[0].seqs.size());

  //! CA function table with the spatial value for the rain.
  CA::TableReal RAIN(GRID,sts[1].seqs.size());  

  //! CA function table with the spatial value for the infiltration.
  //! This contain an infiltration value for each update step, no an
  //! infiltration rate.
  CA::TableReal INF(GRID,sts[2].seqs.size());  
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// spatio temporal load
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  // ----  INIT SPATIAL TEMPORAL  ----

	if(preProcessing){
	    // Initialise the roughness indices into the first MASK buffer pointer.
	    loadAndPopulateSpatialTemporalMap(GRID, sts[0], (*PMASK1), bits_pos[0], bits_pos[1], setup.loadRoughnesses);
	
	
	    // Initialise the rain indices into the second MASK buffer pointer.
	    loadAndPopulateSpatialTemporalMap(GRID, sts[1], (*PMASK2), bits_pos[2], bits_pos[3],setup.loadSpatialRain);

	    // Initialise the infiltration indices into the third MASK buffer pointer.
	    loadAndPopulateSpatialTemporalMap(GRID, sts[2], (*PMASK2), bits_pos[4], bits_pos[5], setup.loadSpatialInfiltrations);

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


    }
    else{
      // Initialise the roughness indices into the first MASK buffer pointer.
      populateSpatialTemporalIndices(GRID,sts[0],(*PMASK1),bits_pos[0], bits_pos[1], setup.preproc_name+"_ROUGH","0");
  
      // need to check if file load is necessary here for spatial rain......

      // Initialise the rain indices into the second MASK buffer pointer.
      populateSpatialTemporalIndices(GRID,sts[1],(*PMASK2),bits_pos[2], bits_pos[3]);

      // Initialise the infiltration indices into the second MASK buffer pointer.
      populateSpatialTemporalIndices(GRID,sts[2],(*PMASK2),bits_pos[4], bits_pos[5],setup.preproc_name+"_INF","0");
    }
  
  // Initialise the object that manage the sequences of the
  // spatial/temporal roughness.
  STRough strough(GRID,sts[0]);

  // Initialise the object that manage the sequences of the
  // spatial/temporal rain.
  STRain  strain(GRID,sts[1]);

  // Add the area with spatial/temporal rainf in the computational domain.
  strain.addDomain(compdomain);

  // Analyse the area where a water level event will happen. WD
  // is used as temporary buffer.
  if(setup.check_vols)
    strain.analyseArea(WD,(*PMASK1),fulldomain);

  // Initialise the object that manage the sequences of the
  // spatial/temporal infiltration.
  STInf  stinf(GRID,sts[2]);

  // Check if the infiltration computation is needed.
  useInfiltration = stinf.needed();

  if(setup.output_console)
  {
    std::cout<<"--------------------------------------------------------" << std::endl; 
    std::cout<<"Loaded Spatial Temporal data"<< std::endl;
    std::cout<<"-----------------" << std::endl; 
    std::cout<<"Number bits for roughness    : "<<rough_bits<<std::endl;
    std::cout<<"Number bits for rain         : "<<rain_bits<<std::endl;
    std::cout<<"Number bits for infiltration : "<<inf_bits<<std::endl;
    if(PMASK1.get()!=PMASK2.get())
      std::cout<<"New buffer allocated         : yes"<<std::endl;
    else
      std::cout<<"New buffer allocated         : no"<<std::endl;
    if(useInfiltration)
      std::cout<<"Infiltration computation     : yes"<< std::endl;
    else
      std::cout<<"Infiltration computation     : no"<< std::endl;
    std::cout<<"--------------------------------------------------------" << std::endl; 
  }
  
  // ----  INIT WATER LEVEL EVENT  ----

  // Initialise the object that manage the water level events.
  WaterLevelManager wl_manager(GRID,wles);
  
  // Add the area with water level in the computational domain.
  wl_manager.addDomain(compdomain);

  // Get the elevation information.
  wl_manager.getElevation(ELV);

  // Analyse the area where a water level event will happen. WD
  // is used as temporary buffer.
  if(setup.check_vols)
    wl_manager.analyseArea(WD,(*PMASK1),fulldomain, PRATIO.get(),dy);

  // ----  INIT RAIN EVENT  ----

  // Initialise the object that manage the rain.
  RainManager rain_manager(GRID,res);

  // Add the area with rain in the computational domain.
  rain_manager.addDomain(compdomain);

  // Analyse the area where it will rain to use for volume cheking. WD
  // is used as temporary buffer.
  if(setup.check_vols)
    rain_manager.analyseArea(WD,(*PMASK1),fulldomain, PRATIO.get(),dy);

  // ---- INIT SPATIAL TEMPORAL RAIN DENSE EVENT ----
  
  // Create the spatial Temporal Dense rain veriable
  STRainDense STD_manager(GRID, _std);

  STD_manager.addDomain(compdomain);

  // Analyse the area where it will rain to use for volume cheking. WD
  // is used as temporary buffer.
  if (setup.check_vols)
      STD_manager.analyseArea(WD, (*PMASK1), fulldomain);


  // ----  INIT INFLOW EVENT  ----

  // Initialise the object that manage the inflow.
  InflowManager inflow_manager(GRID,ies);

  // Add the area with the inflow in the computational domain.
  inflow_manager.addDomain(compdomain);

  //CA::Box newDomainBox = CA::Box(compdomain.extent());
  //std::cout << "rank[" << GRID.caGrid().rank << "]newDomainBox.y():" << newDomainBox.y() << std::endl;
  //std::cout << "rank[" << GRID.caGrid().rank << "]newDomainBox.x():" << newDomainBox.x() << std::endl;
  //std::cout << "rank[" << GRID.caGrid().rank << "]newDomainBox.w():" << newDomainBox.w() << std::endl;
  //std::cout << "rank[" << GRID.caGrid().rank << "]newDomainBox.h():" << newDomainBox.h() << std::endl;

  // Analyse the area where it will inflow to use for volume
  // cheking. WD is used as temporary buffer.
  if(setup.check_vols)
    inflow_manager.analyseArea(WD,(*PMASK1),fulldomain, PRATIO.get(),dy);
  
  // ----  INIT TIME PLOTS AND RASTER GRID ----

  std::string basefilename = ad.output_dir+ad.sdir+setup.short_name;

  // Initialise the object that manage the time plots.
  TPManager tp_manager(GRID,ELV,tps,basefilename,setup.timeplot_files);
  
  // Initialise the object that manage the time plots.
  RGManager rg_manager(GRID,rgs,basefilename,setup.rastergrid_files);    

  // List of raster grid data
  std::vector<RGData> rgdatas(rgs.size());
 
  // Peak buffers
  RGPeak rgpeak;

  for(size_t i = 0; i<rgs.size(); ++i)
  {
    std::string filename = ad.output_dir+ad.sdir+setup.short_name+"_"+setup.rastergrid_files[i]; 
    initRGData(filename, GRID, nodata, rgs[i], rgdatas[i], rgpeak);
  }

  // Time Step plot manager
  TSPlot tsplot(basefilename+"_ts.csv",setup.ts_plot);
 
  // -- INITIALISE  ---

  // Clear the buffer(s) to zero (borders included).
  OUTF1.clear();
  OUTF2.clear();
  A.clear();
  V.clear();
  WD.clear();


  // if hot start add the initial water depths
  if (setup.hotStart) {

	  // assuming the water detphs grid is the same size and position as the elevation

	  // Se the default value of the elevation to be nodata.
	  //ELV.fill(fulldomain, grid.nodata);

	  // Insert the DEM data into the elevation cell buffer. The first
	  // parameter is the region of the buffer to write the data which is
	  // the real domain, i.e. the original non extended domain.
	  WD.insertData(realbox, &(setup.startingDepths->data[0]), setup.startingDepths->ncols, setup.startingDepths->nrows);


	  // add the domain of wetted area to the compuational domain
	  CA::Box initialWettedBox = findWettedArea((*setup.startingDepths), GRID);
	  //if (initialWettedBox != NULL) {
	  compdomain.add(initialWettedBox);
	  //}


      delete(setup.startingDepths);


  }// end of if hotstart


  if (setup.expand_domain && setup.output_computation)
  {
	  CA::Box B(compdomain.extent());
	  std::cout << "DOMAIN   = (" << B.x() << "," << B.y() << "):(" << B.w() << "," << B.h() << ")" << std::endl;
  }

  // If there is not request to expand domain. Set the computtational
  // and extend domain to full domain.
  if(!setup.expand_domain)
  {
    compdomain.clear();
    compdomain.add(fullbox);
  }
  else
  {
    // Set the sequential domain to the compdomain.
    seqdomain = compdomain;
  }

  //newDomainBox = CA::Box(compdomain.extent());
  //std::cout << "rank[" << GRID.caGrid().rank << "]newDomainBox.y():" << newDomainBox.y() << std::endl;
  //std::cout << "rank[" << GRID.caGrid().rank << "]newDomainBox.x():" << newDomainBox.x() << std::endl;
  //std::cout << "rank[" << GRID.caGrid().rank << "]newDomainBox.w():" << newDomainBox.w() << std::endl;
  //std::cout << "rank[" << GRID.caGrid().rank << "]newDomainBox.h():" << newDomainBox.h() << std::endl;

  if(setup.decompose_domain)
  {
    compdomain = dd;
    if(setup.output_console)
    {
      // Print the decomposed domain
      std::cout<<"------------------------------------------" << std::endl; 
      std::cout<<"Domain decomposed:"<<std::endl;
      for(CA::BoxList::ConstIter b = compdomain.begin(); b != compdomain.end(); b++)
      {
	if((*b).e()!=0.0)
	{
	  std::cout<<"("<<(*b).x()<<","<<(*b).y()<<"):("<<(*b).w()<<","<<(*b).h()<<"):("<<(*b).e()<<")"
		   <<std::endl;
	}
      }
      std::cout<<"------------------------------------------" << std::endl; 
    }
  }

	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// core
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  // -- CALCULATE POSSIBLE INITIAL DT ---

  if (setup.hotStart) {

	  // --- UPDATE VA  ---

	  switch (setup.model_type)
	  {
	  case MODEL::WCA2Dv1:

		  // Clear the Velocity and angle.
		  V.clear();
		  A.clear();

		  // Compute the velocity using the averege outflux for the WCA2Dv1 model.
		  // Attention the tollerance is different here. 
		  // Furthermore, this version uses the spacial/temporal roughness
		  // ROUGH where the indices are saved in special bits in the MASK
		  // which position is indicate by STPOS table.
		  // Check if there is water movement over the upstream elevation threshould.
		  if (setup.latlong)
		  {
			  CA::Execute::function(compdomain, velocityWCA2Dv1GEO, GRID, V, A, WD, ELV, (*PTOT), (*PMASK1),
				  VELALARMS,
				  ROUGH, bits_pos[0], bits_pos[1],
				  tol_va, period_time_dt, upstr_elv,
				  (*PRATIO), dy);
		  }
		  else
		  {
			  CA::Execute::function(compdomain, velocityWCA2Dv1, GRID, V, A, WD, ELV, (*PTOT), (*PMASK1), VELALARMS,
				  ROUGH, bits_pos[0], bits_pos[1],
				  tol_va, period_time_dt, upstr_elv);
		  }

		  // CLear the total outflux. 
		  (*PTOT).clear();

		  break;

	  case MODEL::WCA2Dv2:
		  // new velocity system

		  //std::cout << "blah" << std::endl;

		  // Clear the Velocity and angle.
		  V.clear();
		  A.clear();

		  if (setup.latlong)
		  {
			  std::cout << "GEO not implemented" << std::endl;
		  }
		  else
		  {
			  // Compute the velocity using the last outflux (OUTF2)
			  // Compute dt using Hunter formual
			  CA::Execute::function(compdomain, velocityDiffusivev2, GRID, V, A, (*PDT),
			  //CA::Execute::function(compdomain, velocityDiffusive, GRID, V, A, (*PDT),
				  WD, ELV, (*POUTF2), (*PMASK1), VELALARMS,
				  ROUGH, bits_pos[0], bits_pos[1],
				  tol_va, tol_slope, dt, upstr_elv);
		  }


		  break;
	  case MODEL::Diffusive:

		  // Clear the Velocity and angle.
		  V.clear();
		  A.clear();

		  if (setup.latlong)
		  {
			  std::cout << "GEO not implemented" << std::endl;
		  }
		  else
		  {
			  // Compute the velocity using the last outflux (OUTF2)
			  // Compute dt using Hunter formual
			  CA::Execute::function(compdomain, velocityDiffusive, GRID, V, A, (*PDT),
				  WD, ELV, (*POUTF2), (*PMASK1), VELALARMS,
				  ROUGH, bits_pos[0], bits_pos[1],
				  tol_va, tol_slope, dt, upstr_elv);
		  }
		  break;

	  case MODEL::Inertial:
		  if (setup.latlong)
		  {
			  std::cout << "GEO not implemented" << std::endl;
		  }
		  else
		  {
			  // Compute the velocity using the last outflux (OUTF2)
			  // Compute dt using celery
			  CA::Execute::function(compdomain, velocityInertial, GRID, V, A, (*PDT),
				  WD, ELV, (*POUTF2), (*PMASK1), VELALARMS,
				  tol_va, alpha, dt, upstr_elv);
		  }
		  break;
	  }


	  // Retrieve the maximum velocity 
	  V.sequentialOp(seqdomain, vamax, CA::Seq::MaxAbs);



  }//end of if hotstart

  // Find the possible velocity caused by the events.  
  potential_va = 0.0;
  potential_va = std::max( potential_va, rain_manager  .potentialVA(t,period_time_dt) );
  potential_va = std::max( potential_va, inflow_manager.potentialVA(t,period_time_dt) );
  //int _rank;
  //MPI_Comm_rank(MPI_COMM_WORLD, &_rank);
  //std::cout << "Rank[" << _rank << "]: potential_va: " << potential_va << std::endl;
  potential_va = std::max( potential_va, wl_manager    .potentialVA(t,period_time_dt) );
  potential_va = std::max( potential_va, strain        .potentialVA(t,period_time_dt) );
  potential_va = std::max( potential_va, STD_manager   .potentialVA(t, period_time_dt));

  //if(ad.info)
  //  std::cout << "potential_va: " << potential_va << std::endl;
  
  
  switch(setup.model_type)
  {
  case MODEL::WCA2Dv1:
  case MODEL::WCA2Dv2:
  case MODEL::Diffusive:    
  case MODEL::Inertial:    
    // Compute the possible next time step using the critical velocity equations.
    // Use the smaller between dx/dy
    // I Don't like using alpha. But at the moment this is the
    // simplest way to find the potential inpact of events.
    dtn1 = std::min(setup.time_maxdt,alpha*std::min(dx,dy)/potential_va);    
    break;
  }    
  //int _rank;
  //MPI_Comm_rank(MPI_COMM_WORLD, &_rank);
  //std::cout << "Rank[" << _rank << "]: dtn1: " << dtn1 << std::endl;

  if (setup.hotStart) {

	  // Find the maximum velocity
	  CA::Real grid_max_va = vamax;


	  switch (setup.model_type)
	  {
	  case MODEL::WCA2Dv1:
		  //case MODEL::WCA2Dv2:
		  // Compute the possible next dt from the grid velocity and from
		  // the potential velocity fron an event.
		  // Use the minimum between dx and dy.
		  //dtn1 = setup.time_maxdt;
		  //dtn1 = std::min(dtn1, alpha*std::min(dx, dy) / potential_va);
		  dtn1 = std::min(dtn1, alpha*std::min(dx, dy) / grid_max_va);
		  break;

	  case MODEL::WCA2Dv2:

		  dtn1 = setup.time_maxdt;
		  // Retrieve the possible dt using the WCA2Dv2 diffusive formula.
		  // This is very similar to the LISFLOOD-FP diffusive formula.
		  (*PDT).sequentialOp(seqdomain, possible_dt, CA::Seq::Min);

		  // I Don't like using alpha. But at the moment this is the
		  // simplest way to find the potential inpact of events.
		  dtn1 = std::min(dtn1, alpha*std::min(dx, dy) / potential_va);

		  // Furthermore we are using alpha to keep a minimum time step
		  // depending on the velocity like in WCA2Dv1
		  dtn1 = std::min(dtn1, alpha*std::min(dx, dy) / grid_max_va);

		  //std::cerr<<"@@ PDT = "<<possible_dt <<" NDT = "<<dtn1<<" ST = "<<t/60<<" @@"<<std::endl;

		  // Use the possible dt and the dtn1 to find the time step.
		  dtn1 = std::min(dtn1, possible_dt);

		  // Reset the PDT.
		  //A.copy((*PDT)); 
		  (*PDT).fill(fulldomain, setup.time_updatedt);
		  break;

	  case MODEL::Diffusive:
	  case MODEL::Inertial:
		  // Retrieve the possible dt computed.
		  // Alpha is already computed in possible dt.
		  dtn1 = setup.time_maxdt;
		  (*PDT).sequentialOp(seqdomain, possible_dt, CA::Seq::Min);

		  /// Compute potential inpact of events.
		  dtn1 = std::min(dtn1, alpha*std::min(dx, dy) / potential_va);

		  //std::cerr<<"@@ PDT = "<<possible_dt <<" NDT = "<<dtn1<<" ST = "<<t/60<<" @@"<<std::endl;

		  // Use the possible dt and the dtn1 to find the time step.
		  dtn1 = std::min(dtn1, possible_dt);

		  break;
	  }
  }// end of if hotstart
  
  // Compute the next time step as a fraction fo the period time step.
  // Check that dt is between  min max.
  computeDT(dt,dtfrac,dtn1,setup);

  // Update the number of iteration before the next update dt.
  iter_dt = floor(setup.time_updatedt/dt + 0.5);

  // -- PREPARE THE EVENTS MANAGERS WITH THE NEW DT ---

  // Get the amount of events that would happen in each area for each dt
  // for the next period.
  rain_manager  .prepare(t,period_time_dt,dt);
  inflow_manager.prepare(t,period_time_dt,dt);
  wl_manager    .prepare(t,period_time_dt,dt);
  STD_manager   .prepare(t,period_time_dt,dt);

  // -- UPDATE  SPATIAL/TEMPORAL VARIABLE---
  
  // Update spatial roughness
  strough.update(ROUGH,t,period_time_dt,dt, rough, setup.model_type);

  // Update spatial rain
  strain.update(RAIN,t,period_time_dt,dt);

  // Update spatial infiltration
  stinf.update(INF,t,period_time_dt,dt);

  // -- PREMATURE END --

  // get the time when the events will not add any further water.
  t_end_events = std::max(t_end_events, rain_manager.endTime());
  t_end_events = std::max(t_end_events, inflow_manager.endTime());
  t_end_events = std::max(t_end_events, wl_manager.endTime());
  t_end_events = std::max(t_end_events, strain.endTime());
  t_end_events = std::max(t_end_events, STD_manager.endTime());

  if(setup.output_console)
  {
    std::cout<<"The events will end at "<<t_end_events<<" (s) simulation time"<< std::endl;
    std::cout<<"------------------------------------------" << std::endl; 
  }

	CA::State boundary = (setup.rast_boundary) ? 1 : 0;
  // ------------------------- MAIN LOOP -------------------------------

  if(setup.output_console && setup.output_computation)
  {
    std::cout<<"-----------------" << std::endl; 
    std::cout<<"Initialisation time taken (s) = " << total_timer.millisecond()/1000.0 << std::endl;
    std::cout<<"-----------------" << std::endl;     
  }
  if(setup.output_console)
  {
    std::cout<<"Start main loop" << std::endl;
    std::cout<<"-----------------" << std::endl; 
  }

  /*
if(preProcessing)
{
	// CREATE ASCII GRID TEMPORARY BUFFERS
	CA::AsciiGrid<CA::Real>& agtmp1 = setup.eg;
	CA::AsciiGrid<CA::Real>  agtmp2 = setup.eg;

	// Create the MASK cell buffer. The mask is usefull to check wich
	// cell has data and nodata and which cell has neighbourhood with
	// data.
	CA::CellBuffState MASK(GRID);
	CA::createCellMask(fulldomain, GRID, ELV, MASK, nodata);
}*/
    
  while(iter<setup.time_maxiters && t<setup.time_end)
  {

#ifdef CA2D_MPI

	  
	  WD.syncSubGrids();






#endif // CA2D_MPI
		
    // Set this to false. This will be set to the righ value during an
    // update step or before the update itself.
    UpdatePEAK = false;

    // The raster have not been written this itertaion.
    RGwritten     = false;

    // If there is tehe request to expand the domain.
    // Deactivate Box alarm(s) and set them.
    if(setup.expand_domain)
    {
      OUTFALARMS.deactivateAll();
      OUTFALARMS.set();
    }

    // --- CONSOLE OUTPUT ---

    // Check if it is time to output to console.
    if(setup.output_console && t>= time_output)
    {
      outputConsole(iter,oiter,t,dt,avgodt,minodt,maxodt,vamax,upstr_elv,compdomain,setup);

      oiter  = 0;
      avgodt = 0.0;
      minodt = setup.time_maxdt;  // Output step minimum dt.
      maxodt = 0.0;               // Output step maximum dt.
      
      // Compute the next output time.
      time_output += setup.output_period; 

      if(setup.check_vols == true)
      {
	// Compute the total volume of water that is in the water
	// depth (included the boundary cell).
	WD.sequentialOp(fulldomain, wd_volume, CA::Seq::Add);
	wd_volume *= dx*dy; // NEED TO BE IMPROVED!

	std::cout<<"Volume check:"<<std::endl;
	std::cout<<"RAIN = "<<rain_volume<<" INFLOW = "<<inflow_volume<<" INFILT = "<<-inf_volume
		 <<" WD = "<<wd_volume<<std::endl;	
	std::cout<<"-----------------" << std::endl; 
      }

      if(setup.output_console && setup.output_computation)
      {
	std::cout<<"Partial run time taken (s) = " << total_timer.millisecond()/1000.0 << std::endl;
	std::cout<<"-----------------" << std::endl; 
      }
    }

    // --- SIMULATION TIME ---

    // Set the new time step.
    t+=dt;

    // Round the time step to be of 0.01 second precision and then
    // check if it is a multiple of 60 (with a 0.01 precision). If it
    // is the case we need to reset the time of the simulation. If we
    // don't do this the floating point error will creep into the time
    // of simulation.
    CA::Real tround = fround(t, 2);
    if(std::fmod(tround,period_time_dt) < static_cast<CA::Real>(0.01))
    {
      t = tround;
    }

    // Compute output step information.
    avgodt+=dt;
    
    if(dt>maxodt) maxodt=dt;
    if(dt<minodt) minodt=dt;    

    // --- COMPUTE OUTFLUX ---

    switch(setup.model_type)
    {
    case MODEL::WCA2Dv1:
      // Compute outflow using WCA2Dv1.
      // This version check if there is an outflow in the border of the box.
      // Uses the firt pointer to mask. The one that should contains the
      // spatial roughness indices.
      // Furthermore, this version uses the spacial/temporal roughness
      // ROUGH where the indices are saved in special bits in the MASK
      // which position is indicate by STPOS table.
      if(setup.latlong)
      {
	CA::Execute::function(compdomain, outflowWCA2Dv1GEO, GRID, (*POUTF1), (*POUTF2),
			      ELV, WD, (*PMASK1), OUTFALARMS, 
			      ROUGH, bits_pos[0], bits_pos[1], ignore_wd, tol_delwl, dt,
			      (*PRATIO), dy);
      }
      else
      {
	CA::Execute::function(compdomain, outflowWCA2Dv1, GRID, (*POUTF1), (*POUTF2),
			      ELV, WD, (*PMASK1), OUTFALARMS, 
			      ROUGH, bits_pos[0], bits_pos[1], ignore_wd, tol_delwl,dt);
      }
      break;
    case MODEL::WCA2Dv2:
      // Compute outflow using WCA2Dv2.
      if(setup.latlong)
      {
	std::cout<<"GEO not implemented"<<std::endl;	
      }
      else
      { 

	// This save a division operation for each cell.
	CA::Real ratio_dt = dt/previous_dt; 
	CA::State expand = setup.expand_domain;

	CA::Execute::function(compdomain, outflowWCA2Dv2, GRID, (*POUTF1), (*POUTF2),
			      ELV, WD, (*PMASK1), OUTFALARMS, 
			      ROUGH, bits_pos[0], bits_pos[1], ignore_wd, tol_delwl,dt,ratio_dt,expand);       


      }     
      break;
    case MODEL::Diffusive:
      if(setup.latlong)
      {
	std::cout<<"GEO not implemented"<<std::endl;	
      }
      else
      { 
	// Compute outflow using  diffusive formula
	CA::Execute::function(compdomain, outflowDiffusive, GRID, (*POUTF1),
			      ELV, WD, (*PMASK1), OUTFALARMS, 
			      ROUGH, bits_pos[0], bits_pos[1], ignore_wd, tol_delwl,dt);
      }     
      break;
    case MODEL::Inertial:
      if(setup.latlong)
      {
	std::cout<<"GEO not implemented"<<std::endl;	
      }
      else
      {	
	// Compute the outflow in OUTF1 using the inertial formula. OUTF2 contains the previous step
	// ATTENTIO! No border alarm!
	CA::Execute::function(compdomain, outflowInertial, GRID, (*POUTF1), (*POUTF2),
			      ELV, WD, (*PMASK1), OUTFALARMS,
			      ROUGH, bits_pos[0], bits_pos[1], ignore_wd, dt, previous_dt);
      }
      break;

    }    

    // If there is a request to expand the domain.
    // Get the alarms states.
    if(setup.expand_domain)
    {
      OUTFALARMS.get();
    
      // Check if the box alarm is active, that mean there is some
      // outflux on the border of the computational domain.
      if(OUTFALARMS.isActivated(0))
      {
	// !!!!! ATTENTION!!!! THIS has been tested only on a four neighbour square-grid.

	// Set the computational domain to be the extend version and
	// create the new extended one.
	CA::Box extent(compdomain.extent());
	compdomain.clear();
	compdomain.add(extendBox(extent, fullbox,1 , OUTFALARMS.isActivatedArray()));
	// Set the sequential domain to the compdomain.
	seqdomain = compdomain;
      }  
    }
    
    // --- UPDATE WL AND WD  ---
    

    // Update the water depth with the outflux and store the total
    // amount of outflux
    switch(setup.model_type)
    {
    case MODEL::WCA2Dv1:
      if(setup.latlong)
      {
	CA::Execute::function(compdomain, waterdepthWCA2Dv1GEO, GRID, WD, (*POUTF1), (*PTOT), (*PMASK2), 
			      RAIN, bits_pos[2], bits_pos[3], dt, period_time_dt,
			      (*PRATIO), dy);
      }
      else
      {
	CA::Execute::function(compdomain, waterdepthWCA2Dv1, GRID, WD, (*POUTF1), (*PTOT), (*PMASK2), 
			      RAIN, bits_pos[2], bits_pos[3], dt, period_time_dt);      
      }

      // Swap the double buffer
      // Now POUTF1 is zeroed while POUTF2 contains the previous flux.
      std::swap(POUTF1,POUTF2);
      
      break;

    case MODEL::WCA2Dv2:
    case MODEL::Diffusive:
    case MODEL::Inertial:
      if(setup.latlong)
      {
	std::cout<<"GEO not implemented"<<std::endl;	
      }
      else
      {
	// Generic water depth, use OUTF1, erase OUTF2.
	CA::Execute::function(compdomain, waterdepth, GRID, WD, (*POUTF1),  (*POUTF2), (*PMASK2), 
			      RAIN, bits_pos[2], bits_pos[3], dt, period_time_dt);      
      }

      // Swap the double buffer
      // Now POUTF1 is zeroed while POUTF2 contains the previous flux.
      std::swap(POUTF1,POUTF2);
      break;
    }

	//--------------------------------------------------------------------------------

	// time grid outputs

	if (setup.timeGridOut){

		// launch the kernel that calculates the time outputs
		CA::Execute::function(compdomain, timeToThreshold, GRID,
			ELV, WD, 
			(*PMASK2), OUTFALARMS,
			**timeGrid,
			ignore_wd, setup.thresholdDepth,
			t,
			//eg.nodata
			outputNoData
		);


	}





	//--------------------------------------------------------------------------------
    // 
    start_updatedt+=dt;
    
    // --- EXTRA LATERAL EVENT(s) ---

    // Add the eventual rain events.
    rain_manager.add(WD,(*PMASK1),t,dt);

    // Add the eventual inflow events.
    inflow_manager.add(WD,(*PMASK1),t,dt,PRATIO.get(),dy);

    // Add the eventual water level events.
    wl_manager.add(WD,ELV,(*PMASK1),t,dt);

    // Add the eventual Spatial temporal Dense Rain events
    STD_manager.add(WD, (*PMASK1),t, dt);

    // --- COMPUTE NEXT DT, I.E. PERIOD STEP ---

    // Update previous dt
    previous_dt = dt;
    
    // Check if the dt need to be re-computed.
    if(t>=time_dt || --iter_dt==0)
    {   
      // Reset the start of updatedt.
      start_updatedt = 0.0;

      /*
if(preProcessing)
{
			// Compute the Inflitration if needed. It uses the CellBuffer A (with
			// velocity angle) as temporary buffer where to store the volume
			// of water removed.
			if (useInfiltration)
			{
				// Clear the angle.
				A.clear();

				if (setup.latlong)
					CA::Execute::function(fulldomain, infiltrationGEO, GRID, WD, (*PMASK2), A,
						INF, bits_pos[4], bits_pos[5],
						(*PRATIO), dy);
				else
					CA::Execute::function(fulldomain, infiltration, GRID, WD, (*PMASK2), A,
						INF, bits_pos[4], bits_pos[5]);


				//std::cout << "Infiltration rates:" << std::endl;
				//INF.dump(std::cout);


				if (setup.check_vols)
				{
					// Retrieve the volume removed through infiltration
					CA::Real vol = 0;
					A.sequentialOp(seqdomain, vol, CA::Seq::Add);
					inf_volume += vol;
				}
			}
}// end of preprocssing, outputs; i.e. no temporay files*/
      // Deactivate the alarms checked during the velocity
      // calculation if necessary.
      if(setup.ignore_upstream)
      {
	VELALARMS.deactivateAll();
	VELALARMS.set();
      }

      // Lets make sure there are not any rounding errors.
      t = time_dt;

      // Plot the time steps if requested.
      tsplot.output(t,dt);

      // The peak value need to be updated
      UpdatePEAK = true;

      // Update the total volume from the events for the last period.
      rain_volume   += rain_manager.volume();
      inflow_volume += inflow_manager.volume();
      // NO water level at the moment.
      rain_volume   += strain.volume();

      rain_volume += STD_manager.volume();

      // --- UPDATE VA  ---

      switch(setup.model_type)
      {
      case MODEL::WCA2Dv1:

	// Clear the Velocity and angle.
	V.clear();
	A.clear();

	// Compute the velocity using the averege outflux for the WCA2Dv1 model.
	// Attention the tollerance is different here. 
	// Furthermore, this version uses the spacial/temporal roughness
	// ROUGH where the indices are saved in special bits in the MASK
	// which position is indicate by STPOS table.
	// Check if there is water movement over the upstream elevation threshould.
	if(setup.latlong)
	{
	  CA::Execute::function(compdomain, velocityWCA2Dv1GEO, GRID, V, A, WD, ELV, (*PTOT), (*PMASK1), 
				VELALARMS,
				ROUGH, bits_pos[0], bits_pos[1],
				tol_va, period_time_dt,upstr_elv,
				(*PRATIO), dy);
	}
	else
	{
	  CA::Execute::function(compdomain, velocityWCA2Dv1, GRID, V, A, WD, ELV, (*PTOT), (*PMASK1), VELALARMS,
				ROUGH, bits_pos[0], bits_pos[1],
				tol_va, period_time_dt,upstr_elv);
	}

	// CLear the total outflux. 
	(*PTOT).clear();

	break;

      case MODEL::WCA2Dv2:
	// new velocity system

	//std::cout << "blah" << std::endl;

	// Clear the Velocity and angle.
	V.clear();
	A.clear();

	if (setup.latlong)
	{
		std::cout << "GEO not implemented" << std::endl;
	}
	else
	{
		// Compute the velocity using the last outflux (OUTF2)
		// Compute dt using Hunter formual
		//CA::Execute::function(compdomain, velocityDiffusivev2, GRID, V, A, (*PDT),
		CA::Execute::function(compdomain, velocityDiffusive, GRID, V, A, (*PDT),
			WD, ELV, (*POUTF2), (*PMASK1), VELALARMS,
			ROUGH, bits_pos[0], bits_pos[1],
			tol_va, tol_slope, dt, upstr_elv);
	}


	break;
      case MODEL::Diffusive:

	// Clear the Velocity and angle.
	V.clear();
	A.clear();

	if(setup.latlong)
	{
	  std::cout<<"GEO not implemented"<<std::endl;	
	}	
	else
	{
	  // Compute the velocity using the last outflux (OUTF2)
	  // Compute dt using Hunter formual
 	  CA::Execute::function(compdomain, velocityDiffusive, GRID, V, A, (*PDT), 
				WD, ELV, (*POUTF2), (*PMASK1), VELALARMS,
				ROUGH, bits_pos[0], bits_pos[1],
				tol_va, tol_slope,dt,upstr_elv);
	}
	break;

      case MODEL::Inertial:
	if(setup.latlong)
	{
	  std::cout<<"GEO not implemented"<<std::endl;	
	}	
	else
	{
	  // Compute the velocity using the last outflux (OUTF2)
	  // Compute dt using celery
	  CA::Execute::function(compdomain, velocityInertial, GRID, V, A, (*PDT), 
				WD, ELV, (*POUTF2), (*PMASK1), VELALARMS,
				tol_va, alpha, dt,upstr_elv);
	}
	break;
      }
	                
      // Retrieve the maximum velocity 
      V.sequentialOp(seqdomain, vamax,CA::Seq::MaxAbs);                

      // Find the maximum velocity
      //CA::Real grid_max_va = vamax*caNeighbours;
	  CA::Real grid_max_va = vamax;

	  //std::cout << "vamax       : " << vamax << std::endl;

      // Find the possible velocity caused by the events.
      potential_va = 0.0;
      potential_va = std::max(potential_va, rain_manager  .potentialVA(t,period_time_dt) );
      potential_va = std::max(potential_va, inflow_manager.potentialVA(t,period_time_dt) );
      potential_va = std::max(potential_va, wl_manager    .potentialVA(t,period_time_dt) );
      potential_va = std::max(potential_va, strain        .potentialVA(t,period_time_dt) );
      potential_va = std::max(potential_va, STD_manager.potentialVA(t, period_time_dt));

      //if (ad.info)
      //    std::cout << "potential_va: " << potential_va << std::endl;

      switch(setup.model_type)
      {
      case MODEL::WCA2Dv1:
	//case MODEL::WCA2Dv2:
	// Compute the possible next dt from the grid velocity and from
	// the potential velocity fron an event.
	// Use the minimum between dx and dy.
	dtn1 = setup.time_maxdt;     
	dtn1 = std::min(dtn1,alpha*std::min(dx,dy)/potential_va);
	dtn1 = std::min(dtn1,alpha*std::min(dx,dy)/grid_max_va);
	break;

      case MODEL::WCA2Dv2:

	dtn1 = setup.time_maxdt;     
	// Retrieve the possible dt using the WCA2Dv2 diffusive formula.
	// This is very similar to the LISFLOOD-FP diffusive formula.
	(*PDT).sequentialOp(seqdomain, possible_dt,CA::Seq::Min);                
	
	// I Don't like using alpha. But at the moment this is the
	// simplest way to find the potential inpact of events.
	dtn1 = std::min(dtn1,alpha*std::min(dx,dy)/potential_va);

	// Furthermore we are using alpha to keep a minimum time step
	// depending on the velocity like in WCA2Dv1
	dtn1 = std::min(dtn1,alpha*std::min(dx,dy)/grid_max_va);

	//std::cerr<<"@@ PDT = "<<possible_dt <<" NDT = "<<dtn1<<" ST = "<<t/60<<" @@"<<std::endl;

	// Use the possible dt and the dtn1 to find the time step.

	//std::cout << "dtn1       : " << dtn1 << std::endl;
	//std::cout << "possible_dt: " << possible_dt << std::endl;


	dtn1 = std::min(dtn1,possible_dt);		

	// Reset the PDT.
	//A.copy((*PDT)); 
	(*PDT).fill(fulldomain, setup.time_updatedt);          
	break;

      case MODEL::Diffusive:
      case MODEL::Inertial:
	// Retrieve the possible dt computed.
	// Alpha is already computed in possible dt.
	dtn1 = setup.time_maxdt;     
	(*PDT).sequentialOp(seqdomain, possible_dt,CA::Seq::Min);                

	/// Compute potential inpact of events.
	dtn1 = std::min(dtn1,alpha*std::min(dx,dy)/potential_va);

	//std::cerr<<"@@ PDT = "<<possible_dt <<" NDT = "<<dtn1<<" ST = "<<t/60<<" @@"<<std::endl;

	// Use the possible dt and the dtn1 to find the time step.
	dtn1 = std::min(dtn1,possible_dt);		

	break;
      }    
                  	     
      // Compute the next time step as a fraction fo the period time step.
      // Check that dt is between  min max.
      computeDT(dt,dtfrac,dtn1,setup);
      
      // Update the number of iteration before the next update dt.
      iter_dt = floor(setup.time_updatedt/dt + 0.5);
      
      // When the dt need to be recomputed.
      time_dt += period_time_dt;
    
      // Prepare the Events manager for the next update .
      rain_manager  .prepare(t,period_time_dt,dt);
      inflow_manager.prepare(t,period_time_dt,dt);
      wl_manager    .prepare(t,period_time_dt,dt);
      STD_manager   .prepare(t,period_time_dt,dt);

      // Update spatial roughness
      strough.update(ROUGH,t,period_time_dt,dt,rough,setup.model_type);

      // Update spatial rain
      strain.update(RAIN,t,period_time_dt,dt);

      // Update spatial infiltration
      stinf.update(INF,t,period_time_dt,dt);

      if(setup.ignore_upstream)
      {
	// Check the ALARMS computed during the velocity step.
	VELALARMS.get();
	
	// If the alarms one is not active, then there was not any cell
	// which had the water level over the upstream elevation
	// threshold and some flux at the same time. So we can remove
	// the cell from the computation and then lower the upstream
	// threshold.
	// ATTENTION This action is performed only when all the events
	// finished to add water to the domain.
	if(!VELALARMS.isActivated(0) && t > t_end_events)
	{
	  CA::Execute::function(fulldomain, removeUpstr, GRID, (*PMASK1), ELV, upstr_elv);
	  // Do it gain if PMASK2 differ from PMASK1
	  if(PMASK1.get()!=PMASK2.get())	   
	    CA::Execute::function(fulldomain, removeUpstr, GRID, (*PMASK2), ELV, upstr_elv);
	  upstr_elv -= setup.upstream_reduction;
	}
      }



	  // Compute the Inflitration if needed. It uses the CellBuffer A (with
	  // velocity angle) as temporary buffer where to store the volume
	  // of water removed.
	  if (useInfiltration)
	  {
		  // Clear the angle.
		  A.clear();

		  if (setup.latlong)
			  CA::Execute::function(fulldomain, infiltrationGEO, GRID, WD, (*PMASK2), A,
				  INF, bits_pos[4], bits_pos[5],
				  (*PRATIO), dy);
		  else
			  CA::Execute::function(fulldomain, infiltration, GRID, WD, (*PMASK2), A,
				  INF, bits_pos[4], bits_pos[5]);

		  if (setup.check_vols)
		  {
			  // Retrieve the volume removed through infiltration
			  CA::Real vol = 0;
			  A.sequentialOp(seqdomain, vol, CA::Seq::Add);
			  inf_volume += vol;
		  }
	  }


    } // COMPUTE NEXT DT.

    // -------  OUTPUTS --------
         
    // Output time plots.
    tp_manager.output(t, iter, WD, V, setup.output_console);

    // Check if we need to update the peak at every time step.
    if(setup.update_peak_dt)
      UpdatePEAK = true;
    
    // Update the peak
    if(UpdatePEAK)
      rg_manager.updatePeak(compdomain,WD,V,(*PMASK1));


    if (!preProcessing) {
        // Output raster grid. Keep track if the raster have been written.
        RGwritten = rg_manager.output(true, t, previous_dt,
            WD, V, A, (*POUTF2),
            setup.short_name, setup.output_console,
            (iter >= setup.time_maxiters - 1 || t >= setup.time_end));

        // Output raster grid. Keep track if the raster have been written.
        //RGwritten = rg_manager.output(false, t, WD, V, A, setup.short_name, setup.output_console,(iter >= setup.time_maxiters - 1 || t >= setup.time_end));

        //rg_manager.writeOutput(true, t, WD, V, A, setup.output_console, (iter >= setup.time_maxiters - 1 || t >= setup.time_end),
        //	GRID, setup, fulldomain, realbox, eg, (*PMASK1), nodata, ELV, ad
        //	);


    }
//if(preProcessing)
    else
{

bool VAsaved = false;
bool VAPEAKsaved = false;
bool WDsaved = false;
bool WDPEAKsaved = false;


// Create two temporary Cell buffer
    CA::CellBuffReal TMP1(GRID);
    CA::CellBuffReal TMP2(GRID);

    // CREATE ASCII GRID TEMPORARY BUFFERS
    CA::AsciiGrid<CA::Real>& agtmp1 = setup.eg;
    CA::AsciiGrid<CA::Real>  agtmp2 = setup.eg;


for (size_t i = 0; i<rgs.size(); ++i)
{
// Check if it is time to plot or the output is forced!
if (t >= rgdatas[i].time_next )
{
	if (setup.output_console)
	{
		std::cout << "Write Raster Grid (MIN " << t / 60 << "): ";
	}

	// Retrieve the string of the time.
	std::string strtime;
	CA::toString(strtime, std::floor(t + 0.5));

	std::string filename;

	

	// Save the buffer using direct I/O where the main ID is the
	// buffer name and the subID is the timestep.
	switch (rgs[i].pv)
	{
	case PV::VEL:
		if (!VAsaved)
		{


			std::string filenameV;
			std::string filenameA;

			// Create the name of the files if it is going to output a
			// velocity field or simply rasters.
			if (setup.rast_vel_as_vect == true)
			{
				//filenameV = removeExtension(rgdatas[i].filename) + "_" + strtime + ".csv";
				//filenameA = removeExtension(rgdatas[i].filename) + "_" + strtime + ".csvt";
				// note the file name contains no extension in the Dll framework
				filenameV = rgdatas[i].filename + ".csv";
				filenameA = rgdatas[i].filename + ".csvt";
			}
			else
			{
				//filenameV = removeExtension(rgdatas[i].filename) + "_V_" + strtime + ".asc";
				//filenameA = removeExtension(rgdatas[i].filename) + "_A_" + strtime + ".asc";
				//filenameV = rgdatas[i].filename + "_V_" + strtime + ".asc";
				//filenameA = rgdatas[i].filename + "_A_" + strtime + ".asc";

				filenameV = removeExtension(rgdatas[i].filename) + "_V_" + strtime;
				filenameA = removeExtension(rgdatas[i].filename) + "_A_" + strtime;
			}

			TMP1.copy(V);
			TMP2.copy(A);


			// Set the V and A to zero if water depth is less than tollerance.
			CA::Execute::function(fulldomain, zeroedVA, GRID, TMP1, TMP2, WD, (*PMASK1), setup.rast_wd_tol);
			//CA::Execute::function(fulldomain, zeroedVA, GRID, V, A, WD, MASK, setup.rast_wd_tol);

			// Retrieve the velocity and angle data
			TMP1.retrieveData(realbox, &agtmp1.data[0], agtmp1.ncols, agtmp1.nrows);
			TMP2.retrieveData(realbox, &agtmp2.data[0], agtmp2.ncols, agtmp2.nrows);

			// check if we need to output a velocity field
			if (setup.rast_vel_as_vect)
			{
				if (setup.output_console)
					std::cout << "Write Raster Grid: " << filenameV << std::endl;

				// Create the CSV file
				//std::ofstream file( filenameV.c_str() ); 
				FILE* fout = fopen(filenameV.c_str(), "w");

				// Set  manipulators 
				//file.setf(std::ios::fixed, std::ios::floatfield);
				//file.precision(6); 

				// Write the header
				//file<<"X, Y, Speed, Angle_RAD, Angle_DEG, Angle_QGIS"<<std::endl;
				fprintf(fout, "X, Y, Speed, Angle_RAD, Angle_DEG, Angle_QGIS\n");

				// Loop thourhg the grid points.
				for (CA::Unsigned j_reg = realbox.y(), j_mem = 0; j_reg<realbox.h() + realbox.y(); ++j_reg, ++j_mem)
				{
					for (CA::Unsigned i_reg = realbox.x(), i_mem = 0; i_reg<realbox.w() + realbox.x(); ++i_reg, ++i_mem)
					{
						// Create the point and find the coordinates.
						CA::Point p(i_reg, j_reg);
						p.setCoo(GRID);

						// Retrieve the speed and angle values (in radians),
						// the compute the angle value in degrees.
						CA::Real V = agtmp1.data[j_mem * agtmp1.ncols + i_mem];
						CA::Real AR = agtmp2.data[j_mem * agtmp2.ncols + i_mem];
						CA::Real AD = AR * 180 / PI;

						// Compute the angle needed by QGis.
						CA::Real AQ = -AR * 180 / PI + 90;

						// Write the results if the veclocity is more than zero.
						if (V>0)
						{
							//file<<p.coo().x()<<","<<p.coo().y()<<","<<V<<","<<AR<<","<<AD<<","<<AQ<<","<<std::endl;
							fprintf(fout, "%.12f,%.12f,%.6f,%.6f,%.6f,%.6f,\n", p.coo().x(), p.coo().y(), V, AR, AD, AQ);
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
				if (setup.output_console)
					std::cout << "Write Raster Grid: " << filenameV << " " << filenameA << std::endl;

				// Write the  data.
				//CA::writeAsciiGrid(agtmp1,filenameV,setup.rast_places);	      
				//CA::writeAsciiGrid(agtmp2,filenameA,setup.rast_places);	      	  
				agtmp1.writeAsciiGrid(filenameV, setup.rast_places);
				agtmp2.writeAsciiGrid(filenameA, setup.rast_places);
			}





		}
		// ATTENTION! The break is removed since in order to
		// post-process VA we need WD. 
		break;
	case PV::WL:

		///TMP1.copy(WD);

		// Create the name of the file.
		//std::string filename = removeExtension(rgdatas[i].filename) + "_" + strtime + ".asc";	
		//std::string filename = rgdatas[i].filename + "_" + strtime + ".asc";
		filename = removeExtension(rgdatas[i].filename) + "_" + strtime;

		// Reset the buffer
		TMP1.fill(fulldomain, nodata);

		TMP2.copy(WD);

		// Set the water depth to zero if it less than tollerance.

		CA::Execute::function(fulldomain, zeroedWD, GRID, TMP2, (*PMASK1), setup.rast_wd_tol, boundary);

		// Make the water depth data and elevation into the water level.
		CA::Execute::function(fulldomain, makeWL, GRID, TMP1, WD, ELV, (*PMASK1));
		//CA::Execute::function(fulldomain, makeWL, GRID, TMP1, TMP2, ELV, MASK);//??

		// Retrieve the data
		TMP1.retrieveData(realbox, &agtmp1.data[0], agtmp1.ncols, agtmp1.nrows);

		if (setup.output_console)
			std::cout << "Write Raster Grid: " << filename << std::endl;

		// Write the data.
		//CA::writeAsciiGrid(agtmp1,filename,setup.rast_places);	      
		agtmp1.writeAsciiGrid(filename, setup.rast_places);

		break;
	case PV::WD:


		// Create the name of the file.
		//std::string filename = removeExtension(rgdatas[i].filename) + "_" + strtime + ".asc";	
		//std::string filename = rgdatas[i].filename + "_" + strtime + ".asc";
		filename = removeExtension(rgdatas[i].filename) + "_" + strtime;

		// Water depth already loaded.
		TMP1.copy(WD);

		// Set the water depth to zero if it less than tollerance.
		
		CA::Execute::function(fulldomain, zeroedWD, GRID, TMP1, (*PMASK1), setup.rast_wd_tol, boundary);


		// Retrieve the data
		//WD.retrieveData(realbox, &agtmp1.data[0], agtmp1.ncols, agtmp1.nrows);
		TMP1.retrieveData(realbox, &agtmp1.data[0], agtmp1.ncols, agtmp1.nrows);

		if (setup.output_console)
			std::cout << "Write Raster Grid: " << filename << std::endl;

		// Write the data.
		//CA::writeAsciiGrid(agtmp1,filename,setup.rast_places);    
		agtmp1.writeAsciiGrid(filename, setup.rast_places);

		break;
	default:
		break;
	}
	

	// Update the next time to save a raster grid.
	rgdatas[i].time_next += rgs[i].period;
}

			
		}
		
}// end of preprocessing// end of processing, outputing files, instead of temporarly files
    // ---- END OF ITERATION ----

    // Increase time step (and output time step)
    iter++;    
    oiter++;	      
	currentTime = t;


	if (abortIt){
		return 0;
	}


  }  // end of main loop here..?!

  // --- END OF MAIN LOOP ---###########################################################################################################################################

#ifdef CA2D_MPI
  //std::stringstream subGridIdStream;
  //subGridIdStream << GRID.caGrid().rank;
  //std::string outFileName = ad.output_dir + "tempWD" + subGridIdStream.str() + ".txt";
  //std::cout << "Writing: " + outFileName;
  //std::ofstream outFile(outFileName);


  //WD.clear();
  //CA::Execute::function(fulldomain, banding, GRID,WD, (*PMASK1));

  //CA::Box newDomainBox = CA::Box(compdomain.extent());
  //std::cout << "rank[" << GRID.caGrid().rank << "]newDomainBox.y():" << newDomainBox.y() << std::endl;
  //std::cout << "rank[" << GRID.caGrid().rank << "]newDomainBox.x():" << newDomainBox.x() << std::endl;
  //std::cout << "rank[" << GRID.caGrid().rank << "]newDomainBox.w():" << newDomainBox.w() << std::endl;
  //std::cout << "rank[" << GRID.caGrid().rank << "]newDomainBox.h():" << newDomainBox.h() << std::endl;


  //CA::Execute::function(compdomain, computeArea, GRID, WD, (*PMASK2));

  //WD.dump(outFile);
#endif // CA2D_MPI



  // --- OUTPUT PEAK & FINAL ---  
  
	   // Check if the raster have not been written in the last
	   // iteration. If not, we need to make sure that we save the PEAK and
	   // the FINAL extend (if requested).
	//if (!RGwritten)
	//{
		// Make sure to output the last peack value.  
		//rg_manager.updatePeak(compdomain, WD, V, (*PMASK1));
		//rg_manager.output(false, t, WD, V, A, setup.short_name, setup.output_console, true);
		//rg_manager.outputPeak(t, WD, V, setup.short_name, setup.output_console);
	//}

if(!preProcessing){
  // Check if the raster have not been written in the last
  // iteration. If not, we need to make sure that we save the PEAK and
  // the FINAL extend (if requested).
  if(!RGwritten)
  {
    // Make sure to output the last peack value.  
    rg_manager.updatePeak(compdomain,WD,V,(*PMASK1));
    rg_manager.output(true, t, previous_dt,
		WD, V, A, (*POUTF2),
		setup.short_name, setup.output_console, true);
    rg_manager.outputPeak(t,WD, V, setup.short_name, setup.output_console);
  }
}
  // --- CONSOLE OUTPUT ---
  
  // Check if it is time to output to console.
  if(setup.output_console && t>= time_output)
  {    
    outputConsole(iter,oiter,t,dt,avgodt,minodt,maxodt,vamax,upstr_elv,compdomain,setup);

    if(setup.check_vols == true)
    {
      // Compute the total volume of water that is in the water
      // depth (included the boundary cell).
      WD.sequentialOp(fulldomain, wd_volume, CA::Seq::Add);
      wd_volume *= dx*dy;

      std::cout<<"Volume check:"<<std::endl;
      std::cout<<"RAIN = "<<rain_volume<<" INFLOW = "<<inflow_volume<<" INFILT = "<<-inf_volume
	       <<" WD = "<<wd_volume<<std::endl;	      
      std::cout<<"-----------------" << std::endl;       
    }
  }

  // ---- TIME OUTPUT ----

  if(setup.output_console && setup.output_computation)
  {
    std::cout<<"-----------------" << std::endl; 
    std::cout<<"Total run time taken (s) = " << total_timer.millisecond()/1000.0 << std::endl;
    std::cout<<"-----------------" << std::endl; 
  }

  if(setup.output_console)
  {
    time_t t = time(0);   // get time now
    struct tm * now = localtime( & t );
    std::cout<<"Simulation : "<<setup.sim_name<< std::endl;
    std::cout<<"Model      : "<<setup.model_type<< std::endl;
    std::cout<<"Date End   : "<<(now->tm_year + 1900)<<"-"<<(now->tm_mon + 1)<<'-'<<now->tm_mday
	     <<" "<<now->tm_hour<<":"<<now->tm_min<<":"<<now->tm_sec << std::endl;
    std::cout<<"------------------------------------------" << std::endl; 
  }


  // ------------- time output grid ----------------

  // save the time output grids if called for
  if ((*timeGrid) != NULL){
  //if (timeGrid != NULL){


	  //std::cout << "blah!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
	  //std::cout << eg.getString() << std::endl;

	  //timeGrid->saveData()

	  // Create the name of the file.
	  //std::string filename = removeExtension(rgdatas[i].filename) + "_PEAK.asc";	
	  //std::string filename = rgdatas[0].filename + "_TIME.asc"; // assumption!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  //std::string filename = ad.output_dir + ad.sdir + setup.short_name + "_" + setup.timeGridName + "_TIME.asc";
	  std::string filename = ad.output_dir + ad.sdir + setup.short_name + "_" + setup.timeGridName + "_TIME";

	  // Water depth already loaded.
	  //CA::AsciiGrid<CA::Real> timeGridOutput = eg;
	  CA::AsciiGrid<CA::Real> timeGridOutput;


	  

	  timeGridOutput.ncols = eg.ncols;
	  timeGridOutput.nrows = eg.nrows;
	  timeGridOutput.xllcorner = eg.xllcorner;
	  timeGridOutput.yllcorner = eg.yllcorner;
	  //timeGridOutput.nodata = eg.nodata;
	  timeGridOutput.nodata = outputNoData;


	  timeGridOutput.cellsize = eg.cellsize;
	  //timeGridOutput.data.resize(eg.data.size(), eg.nodata);
	  //timeGridOutput.data.resize(timeGridOutput.ncols*timeGridOutput.nrows, eg.nodata);
	  timeGridOutput.data.resize(timeGridOutput.ncols*timeGridOutput.nrows, outputNoData);




	  //std::cout << eg.getString() << std::endl;

	  // Retrieve the data
	  (*timeGrid)->retrieveData(realbox, &timeGridOutput.data[0], timeGridOutput.ncols, timeGridOutput.nrows);

	  //std::cout << "blah!!!!!!!!!!!!!!!!!!"<< std::endl;

	  if (setup.output_console)
		  std::cout << "Write Time Raster Grid: " << filename << std::endl;

	  // Write the data.
	  //CA::writeAsciiGrid(agtmp1,filename,setup.rast_places);      
	  //agtmp1.writeAsciiGrid(filename, setup.rast_places);
	  timeGridOutput.writeAsciiGrid(filename, setup.rast_places);


  }// end of if time grid outputs

// checking for flood neighbouring no data cells ###################################################################################

  // Create  the water depth cell buffer.
  CA::CellBuffState FLAGS(GRID);

  FLAGS.fill(fullbox, 0);

  // set the 29th bit in borders, so as single kernel can indentify the difference between near no data inside the grid and at the edge
  //ELV.bordersValue(borders, agtmp1.nodata);

  CA::State borderMask = 0;
  borderMask = caWriteBitsState(1, borderMask, 29, 30);
  //(*PMASK1).bordersValue(borders, borderMask, CA::Bdr::Operator::Add);
  (*PMASK1).bordersValue(borders, borderMask);

  // full domain include borders
  CA::Execute::function(fulldomain, checkFloodNoDataNeighbourhood, GRID, 
  //CA::Execute::function(realdomain, checkFloodNoDataNeighbourhood, GRID,
	  ELV,
	  WD
	  //*(rgpeak.WD)
	  , (*PMASK1), 
	  FLAGS,
	  eg.nodata
	  );

  // Find largest code (0 -default, 1 water on a no_Data cell, 2 water neighbouring no_data)
  CA::State code = -1;

  // Find highest elevation
  FLAGS.sequentialOp(fulldomain, code, CA::Seq::Max);


  std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
  std::cout << "CODE: " << code << std::endl;
  std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;

  resultcode = code;


  //##################################################################################################################################

  // ------------- Badnded ascii outputs -------------

  if (setup.finalDepthsBandedOutput) {

	  CA::TableReal BANDS(GRID, setup.numberOfBands);
	  BANDS.update(0, setup.numberOfBands, &setup.bands[0], setup.numberOfBands);


	  CA::CellBuffState finalDepthsbandedGrid(GRID);

	  //finalDepthsbandedGrid.fill(realbox, 0.0);

	  
	  // launch the kernel that calculates the time outputs
	  CA::Execute::function(realdomain, 
		  banding, GRID,
		  WD,
		  finalDepthsbandedGrid,
		  BANDS,
		  (*PMASK2), 
		  setup.numberOfBands);
		  

	  CA::AsciiGrid<CA::State> finalDepthsBandedOutputASCII;

	  finalDepthsBandedOutputASCII.ncols = eg.ncols;
	  finalDepthsBandedOutputASCII.nrows = eg.nrows;
	  finalDepthsBandedOutputASCII.xllcorner = eg.xllcorner;
	  finalDepthsBandedOutputASCII.yllcorner = eg.yllcorner;
	  finalDepthsBandedOutputASCII.nodata = eg.nodata;
	  finalDepthsBandedOutputASCII.cellsize = eg.cellsize;
	  //timeGridOutput.data.resize(eg.data.size(), eg.nodata);
	  finalDepthsBandedOutputASCII.data.resize(finalDepthsBandedOutputASCII.ncols*finalDepthsBandedOutputASCII.nrows, eg.nodata);


	  // Retrieve the data
	  finalDepthsbandedGrid.retrieveData(realbox, &finalDepthsBandedOutputASCII.data[0], finalDepthsBandedOutputASCII.ncols, finalDepthsBandedOutputASCII.nrows);


	  std::string filename = ad.output_dir + ad.sdir + setup.short_name + "_" + setup.finalDepthsBandedOutputName + "_MaxDepthBanded.asc";


	  // Write the data.
	  finalDepthsBandedOutputASCII.writeAsciiGrid(filename, setup.rast_places);


  }//------------------------------- end of if finalDepthsBandedOutput


	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Post processing official
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if(preProcessing)
{

	if (setup.output_computation)
	{
		std::cout << "Post-processing : " << setup.sim_name << std::endl;
		std::cout << "------------------------------------------" << std::endl;
	}

	// CREATE ASCII GRID TEMPORARY BUFFERS
	CA::AsciiGrid<CA::Real>& agtmp1 = eg;
	CA::AsciiGrid<CA::Real>  agtmp2 = eg;

	// ----  CELL BUFFERS ----

	// Create two temporary Cell buffer
	CA::CellBuffReal TMP1(GRID);
	CA::CellBuffReal TMP2(GRID);

	// NOTE probably need to clear these temporary grids out from possible revious use


	// ----  SCALAR VALUES ----

	iter = 0;		        // The actual iteration number.
	t = setup.time_start;    // The actual time in seconds
	CA::Real     t_nearest = setup.time_end;      // The next time step
	nodata = agtmp1.nodata; // capture terrain no data value, to display actual no data cells

										 // set output no_data value to zero
	agtmp1.nodata = 0;

	// -- CREATE FULL MASK ---

	

	// ----  INIT RASTER GRID ----

	// List of raster grid data
	//std::vector<RGData> rgdatas(rgs.size());

	// Peak buffers.
	//RGPeak rgpeak;

	// Allocate temporary Ascii file data to ouptut raster.
	if (rgs.size()>0)
	{
		agtmp1.data.resize(agtmp1.ncols * agtmp1.nrows, nodata);
		agtmp2.data.resize(agtmp2.ncols * agtmp2.nrows, nodata);
	}

	//for (size_t i = 0; i<rgs.size(); ++i)
	//{
	//	std::string filename = ad.output_dir + ad.sdir + setup.short_name + "_" + setup.rastergrid_files[i];
	//	initRGData(filename, GRID, nodata, rgs[i], rgdatas[i], rgpeak);
	//}

	// ---- CONTAINER DATA TO REMOVE ---

	// Create the container with the data to remove.
	typedef std::vector< std::pair<std::string, std::string> > RemoveID;
	RemoveID removeIDsCB;
	RemoveID removeIDsEB;

	

	// ----  PEAK RASTER GRID ----


	WD.copy(*rg_manager.getPeak_Depths());

	// Set the water depth to zero if it less than tollerance.
	//CA::State boundary = (setup.rast_boundary) ? 1 : 0;
	//CA::Execute::function(fulldomain, zeroedWD, GRID, WD, MASK, setup.rast_wd_tol, boundary);
    CA::Execute::function(fulldomain, zeroedWD, GRID, WD, (*PMASK1), setup.rast_wd_tol, boundary);




	// Need to realod WD only once for PEAK.
	bool WDloaded = false;

	for (size_t i = 0; i<rgdatas.size(); ++i)
	{
		// Check if the peak values need to be saved.
		if (rgs[i].peak == true)
		{


			


			switch (rgs[i].pv)
			{
			case PV::VEL:
			{
				// Create the name of the file.
				//std::string filenameV = removeExtension(rgdatas[i].filename) + "_V_PEAK.asc";	
				//std::string filenameV = rgdatas[i].filename + "_V_PEAK.asc";
				std::string filenameV = removeExtension(rgdatas[i].filename) + "_V_PEAK";

				// Reset the buffer
				TMP1.fill(fulldomain, nodata);



				TMP1.copy(*rg_manager.getPeak_Velocities());

				// Set the V and A to zero if water depth is less than tollerance.
				//CA::Execute::function(fulldomain, zeroedVA, GRID, TMP1, TMP2, WD, MASK, setup.rast_wd_tol);
                CA::Execute::function(fulldomain, zeroedVA, GRID, TMP1, TMP2, WD, (*PMASK1), setup.rast_wd_tol);

				// Retrieve the data
				TMP1.retrieveData(realbox, &agtmp1.data[0], agtmp1.ncols, agtmp1.nrows);

				if (setup.output_console)
					std::cout << "Write Raster Grid: " << filenameV << std::endl;

				// Write the data.
				//CA::writeAsciiGrid(agtmp1, filenameV, setup.rast_places);
				agtmp1.writeAsciiGrid(filenameV, setup.rast_places);

			}
			break;

			case PV::WL:
			{
				// Create the name of the file.
				//std::string filename = removeExtension(rgdatas[i].filename) + "_PEAK.asc";	
				//std::string filename = rgdatas[i].filename + "_PEAK.asc";
				std::string filename = removeExtension(rgdatas[i].filename) + "_WLPEAK";

				// Reset the buffer
				TMP1.fill(fulldomain, nodata);


				// Make the water depth data and elevation into the water level.
				//CA::Execute::function(fulldomain, makeWL, GRID, TMP1, WD, ELV, MASK);
                CA::Execute::function(fulldomain, makeWL, GRID, TMP1, WD, ELV, (*PMASK1));

				// Retrieve the data
				TMP1.retrieveData(realbox, &agtmp1.data[0], agtmp1.ncols, agtmp1.nrows);

				if (setup.output_console)
					std::cout << "Write Raster Grid: " << filename << std::endl;

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

				if (setup.output_console)
					std::cout << "Write Raster Grid: " << filename << std::endl;

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

	if (setup.remove_data)
	{
		for (size_t i = 0; i<removeIDsCB.size(); i++)
		{
			std::pair <std::string, std::string>& ID = removeIDsCB[i];;
			CA::CellBuffReal::removeData(ad.data_dir, ID.first, ID.second);
		}
		for (size_t i = 0; i<removeIDsEB.size(); i++)
		{
			std::pair <std::string, std::string>& ID = removeIDsEB[i];;
			CA::EdgeBuffReal::removeData(ad.data_dir, ID.first, ID.second);
		}
	}

	if (setup.remove_prec_data)
	{
		// Remove Elevation data.
		CA::CellBuffReal::removeData(ad.data_dir, setup.preproc_name + "_ELV", "0");

		// Remove Grid data.
		CA::Grid::remove(ad.data_dir, setup.preproc_name + "_Grid", "0");

		if (CA::CellBuffReal::existData(ad.data_dir, setup.preproc_name + "_RATIO", "0"))
			CA::CellBuffReal::removeData(ad.data_dir, setup.preproc_name + "_RATIO", "0");

		if (CA::CellBuffState::existData(ad.data_dir, setup.preproc_name + "_ROUGH", "0"))
			CA::CellBuffState::removeData(ad.data_dir, setup.preproc_name + "_ROUGH", "0");

		if (CA::CellBuffState::existData(ad.data_dir, setup.preproc_name + "_INF", "0"))
			CA::CellBuffState::removeData(ad.data_dir, setup.preproc_name + "_INF", "0");
	}

	if (setup.output_console)
		std::cout << "Cleaned data" << std::endl;


	if (setup.output_computation)
	{
		std::cout << "-----------------" << std::endl;
		std::cout << "Total run time taken (s) = " << total_timer.millisecond() / 1000.0 << std::endl;
		std::cout << "-----------------" << std::endl;
	}

}



















	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// end
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  return 0;
}
//---------------------------------------------------------------------------------------------------------------------------
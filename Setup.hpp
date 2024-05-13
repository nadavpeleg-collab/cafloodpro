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

#ifndef _SETUP_HPP_
#define _SETUP_HPP_


//! \file Setup.hpp
//!  Contains the structure with the main setup for the CA2D model.
//! \author Michele Guidolin, University of Exeter, 
//! contact: a.s.chen [at] exeter.ac.uk
//! \date 2024-05

#include"BaseTypes.hpp"
#include<string>
#include<vector>
#include "ArgsData.hpp"


//! Structure with the data that define the main variables used to
//! execute the CA2D model. 
struct Setup
{

	
    bool loadSpatialRainDense;


	//bool hotStart = false;
	bool hotStart;
    std::string startingDepthFileName;

	//CA::AsciiGrid<CA::Real>* startingDepths = NULL;
	CA::AsciiGrid<CA::Real>* startingDepths;

	//bool addHeights = false;
	bool addHeights;
	std::vector<double> addedHeight;
	CA::AsciiGrid<CA::State>* heightAdjustmentMap;

	// array of band values, used to produce the banded output
	double* bands;
	unsigned int numberOfBands;

	//bool finalDepthsBandedOutput = false;
	//bool maxDepthsBandedOutput = false;
	//bool velsBandedOutput = false;
	bool finalDepthsBandedOutput;
	bool maxDepthsBandedOutput;
	bool velsBandedOutput;
	////bool MaxDepthsBandedOutput = false;

	
	//std::string finalDepthsBandedOutputName = "FinalDepthsBanded";
	std::string finalDepthsBandedOutputName;

  //  --- SIMULATION  --- 
  std::vector<CA::Unsigned> sim_version; 
  std::string    sim_name;	//!< Name of the simulation. 
  std::string    short_name;	//!< Short name used as base for output files. 
  std::string    preproc_name;	//!< Name used as base for the pre-processed files.
  MODEL::Type    model_type;	//!< The model used in the simulation

  //bool loadTerrain = true;
  bool loadTerrain;

  //bool loadRoughnesses = true;
  bool loadRoughnesses;

  //bool loadSpatialInfiltrations = true;
  bool loadSpatialInfiltrations;

  //bool loadSpatialRain = true;
  bool loadSpatialRain;

  //std::string timeGridName = "TimeGrid";
  std::string timeGridName;
  CA::AsciiGrid<CA::Real> eg;

  //CA::AsciiGrid<CA::Real> eg;

  //bool timeGridOut = false;
  bool timeGridOut;
  CA::Real thresholdDepth;

//#ifdef CAAPI_GDAL_OPTION

  //std::string gdalDriverName;
  //std::string geoCS;

//#endif



  Setup() :
	  hotStart(false), startingDepths(NULL), addHeights(false), finalDepthsBandedOutput(false), maxDepthsBandedOutput(false), velsBandedOutput(false), finalDepthsBandedOutputName("FinalDepthsBanded"),
	  loadTerrain(true), loadRoughnesses(true), loadSpatialInfiltrations(true), loadSpatialRain(true), timeGridName("TimeGrid"), timeGridOut(false),
	  strough_files(), strain_files(), stinf_files(), rainevent_files(), wlevent_files(), inflowevent_files(), timeplot_files(), rastergrid_files(),
      rainSpatialDense_files(), startingDepthFileName(""), loadSpatialRainDense(true), rainNetCDF(false)
//#ifdef CAAPI_GDAL_OPTION
	  //,gdalDriverName(""), geoCS("")

//#endif
  {

  }



  //  --- TIME VALUES  ---  
  CA::Real      time_start;	//!< Starting time of the simulation (seconds). 
  CA::Real      time_end;	//!< Ending time of the simulation (seconds). 
  CA::Real      time_maxdt;	//!< The initial and maximum time step.
  CA::Real      time_mindt;	//!< The minimum time step.
  CA::Real      time_updatedt;	//!< The time step to update the time step.
  CA::Real      time_alpha;	//!< The proportion of adaptive time step taken. 
  CA::Unsigned  time_maxiters;	//!< The maximum number of iterations. 

  //  --- SIMULATION PARAMETERS ---
  CA::Real      roughness_global; //!< Global default value.  
  CA::Real      infrate_global;   //!< Global default value (mm/hr).  
  CA::Real      tolerance;	  //!< The water difference between cell that can be ignored.
  CA::Real      ignore_wd;	  //!< The water depth that can be ignored.
  CA::Real      tol_slope;	  //!< The slope difference (in %) between cell used to compute dt.
  
  //! The elevation of boundary cell, a high value represent a CLOSED
  //! boundary a low value represents an OPEN boundary.
  CA::Real      boundary_elv;	   

  //  --- PROJECTION PARAMETER ---  
  bool          latlong;	//!< If true various raster file use Latitude/Longitude degrees values.

  //  --- ELEVATION  ---  
  //! ARC/INFO ASCII GRID format file with the specific
  //! elavation value for each cell and the no data value.
  std::string   elevation_ASCII;

  /// --- RAIN SPATIAL DENSE ---    
  //! This list must contains at least two files (in this order) that are used
  //! to create the spatial varying height adjustments:
  //! 1) CSV file with the times value sequence.  
  //! 2) ARC/INFO ASCII GRID map containing the rain values for each time provided in 1).
  std::vector<std::string> rainSpatialDense_files;

  bool rainNetCDF;

  /// --- SPATIAL HEIGH ADJUSTMENT ---    
  //! This list must contains two files (in this order) that are used
  //! to create the spatial varying height adjustments:
  //! 1) ARC/INFO ASCII GRID map containing the indices.
  //! 2) CSV file with the roughness value sequence at various indices.  
  std::vector<std::string> heighAdjustmentMap_files;
 
  /// --- SPATIAL/TEMPORAL ROUGHNESS ---    
  //! This list must contains two files (in this order) that are used
  //! to create the spatial/temporal varying roughness:
  //! 1) ARC/INFO ASCII GRID map containing the indices.
  //! 2) CSV file with the roughness value/time sequence at various indices.  
  //! \attention The model can work with only one single spatial roughness.
  std::vector<std::string> strough_files;    

  /// --- SPATIAL/TEMPORAL RAIN ---    
  //! This list must contains two files (in this order) that are used
  //! to create the spatial/temporal varying rain:
  //! 1) ARC/INFO ASCII GRID map containing the indices.
  //! 2) CSV file with the rain value/time sequence at various indices.  
  std::vector<std::string> strain_files;    

  /// --- SPATIAL/TEMPORAL INFILTRATION ---    
  //! This list must contains two files (in this order) that are used
  //! to create the spatial/temporal varying infiltration:
  //! 1) ARC/INFO ASCII GRID map containing the indices.
  //! 2) CSV file with the infitration value/time sequence at various indices.  
  std::vector<std::string> stinf_files;    

  //  --- RAIN EVENT  ---  
  //! CSV file(s) with the configuration of the rain event(s) to add
  //! into the grid.
  std::vector<std::string> rainevent_files;

  //  --- WATER LEVEL EVENT  ---  
  //! CSV file(s) with the configuration of the water level event(s)
  //! to add into the grid.
  std::vector<std::string> wlevent_files;

  //  --- INFLOW EVENT  ---  
  //! CSV file(s) with the configuration of the inflow event(s)
  //! to add into the grid.
  std::vector<std::string> inflowevent_files;
  
  //  --- TIMEPLOT  ---  
  //! CSV file(s) with the configuration of the time plot(s) to
  //! produce, i.e. the value of a physical variable (WD,WL,VEL) at
  //! specific time and location.
  std::vector<std::string> timeplot_files;

  //  --- RASTERS  ---  
  //! CSV file(s) with the configuration of the raster grid(s) to
  //! produce, i.e. the grid value of a physical variable (WD,WL,VEL) at
  //! specific time.
  std::vector<std::string> rastergrid_files;

  //  --- OUTPUT  ---  
  bool output_console;		//!< If true output info into console.
  CA::Real output_period;	//!< The period to ouput in seconds.
  bool output_computation;	//!< If true ouput computation time.
  bool terrain_info;		//!< If true print terrain info, like slope.
  bool ts_plot;                 //!< If true create a file that plot the time step.
  
  //  --- CHEKS  ---
  bool check_vols;		//!< If true compute the various input/ouput volumes. 

  //  --- REMOVE  ---
  bool remove_data;		//!< If true remove the data created by the model (no pre-proc).
  bool remove_prec_data;	//!< If true remove the data created by the pre-proces

  // ---  RASTER OUTPUT ---
  bool     rast_vel_as_vect;	//!< If true output the velocity raster as a vector field.  
  CA::Real rast_wd_tol;		//!< Ignore cell in raster that are lower than this tolerance.
  bool     rast_boundary;	//!< If true output the value of boundary cell. 
  int      rast_places;		//!< The number of decimal places after the comma.

  // ---  PEAK UPDATE ---
  bool  update_peak_dt;	        //!< If true update the peak at every time step. Default false.

  // --- OPTIMISATIONS ---
  bool     expand_domain;	//!< If true expand the computational domain when needed.
  bool     ignore_upstream;     //!< If true ignore upstream cells.
  CA::Real upstream_reduction;	//!< The amount of elevation to reduce
  bool     decompose_domain;	//!< If true decompose the domain in smaller boxes.
  CA::Real decompose_eff;       //!< The efficiency of the decomposition.

};


//! Initialise the setup structure usign a CSV file. 
//! Each row represents a new "variable" where the 
//! first column is the name of the element 
//! and the following columns have the multiple/single values.
//! \attention The order of elements is not important.
//! \param[in]  filename This is the file where the data is read.
//! \param[out] setup    The structure containing the read data.
//! \return A non zero value if there was an error.
int initSetupFromCSV(const std::string& filename, Setup& setup);


#endif

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

#ifndef _SPATIALTEMPORALDENSE_HPP_
#define _SPATIALTEMPORALDENSE_HPP_


//! \file Spatial.hpp
//!  Contains the structure(s) and classe(s) that are used to manage
//!  spatial/temporal variables.
//! \author Mike Gibson, University of Exeter, 
//! contact: m.j.gibson [at] exeter.ac.uk
//! \date 2024-05


#include"ca2D.hpp"
#include"BaseTypes.hpp"
#include"Box.hpp"
#include<string>
#include<vector>





//! Structure with the configuration value that define a DENSE spatial
//! temporal variable in the CA2D model. 
//! The dense spatial/temporal variables contains a 1) list of grid maps
//! which contain the raw rainfall rates, 2) a list of times at which to
//! apply the grid map rates.
//! \warning For the ASCII/GRID grid map, the nodata value will be
//! considered to be zero. Thus there cannot be a zero index into the CVS.
struct SpatialTemporalDense
{
  std::string                          name;                //!< Name of the spatial variable. 
  std::vector<CA::AsciiGrid<CA::Real>> maps;                //!< The map with the indices. ATTENTION contains only the header, nodata.
  std::vector<CA::Real>                times;               //!< The value/time sequences, one for each index.
  //CA::Unsigned                         l;	                //!< The number of elements in the longest sequence. 
  CA::Real                             end;                 //!< The highest end of a time sequence .

  std::vector<std::string>             map_filenames;       //!< The filename of the ascii grid file.

  std::vector<CA::Real>                rain_maxes;
  std::vector<CA::Real>                rain_totals;


};


//! Initialise the spatial temporal structure.
//! This first arguement list must contain two files or more
//! that are used to create the spatial/temporal-dense variable:
//! 1) CSV file with the roughness value/time sequence at various indices.  
//! 2) at least one; ARC/INFO ASCII GRID map containing the rainfall rates.
//! \warning The ASCI/GRID map dat is not laoded. Only the header.
//! \attention The order of elements is not important.
//! \param[in]  dir       The directory where the files are located.
//! \param[in]  filenames This is the list of two files where the data is read.
//! \param[out] st        The structure containing the read data.
//! \return A non zero value if there was an error.
int initSpatialTemporalDenseFromCSV(const std::string& dir,
			       const std::vector<std::string>& filenames, SpatialTemporalDense& std, bool userNetCDF);


int loadSpatialTemporalDenseData(CA::Grid& GRID, SpatialTemporalDense& std, CA::Box domain, bool loadGrids);


/*
//! Populate the CellBuffer with the indices of the spatial temporal
//! variable. The indices are loaded at specific bits positions and
//! mantain the value of rest of the bits.
//! \param mainid  The main id of the Preprocessed data to load.
//! \param subid   The main id of the Preprocessed data to load.
//! \return A non zero value if there was an error.
int populateSpatialTemporalIndices(CA::Grid&  GRID, const SpatialTemporalDense& st,
				   CA::CellBuffState& DST,  CA::State start, CA::State stop,
				   const std::string& mainid="", const std::string& subid="");


//! Save the spatial temporal map as Preprocessed data
//! \return A non zero value if there was an error.
int saveSpatialTemporalMap(CA::Grid&  GRID, const SpatialTemporalDense& st,
			   const std::string& mainid, const std::string& subid, bool loadMap = true);

//! directly load and populate the statial temporal maps in one go as opposed to saving temporal files
int loadAndPopulateSpatialTemporalMap(
	CA::Grid&  GRID,
	const SpatialTemporalDense& st,
	CA::CellBuffState& DST, CA::State start, CA::State stop,
	bool loadMap = true);
    */



//! Class that manage the spatial/temporal dense rain. This is similar to
//! RainManager. 
class STRainDense
{
private:
  
  //! The structure used during the model computation to store the
  //! spatial/temporal data.
  struct Data
  {

      size_t   index;	        //!< The index of the rain data (rains/times).

      bool firstUploadDone;

      //std::vector<CA::AsciiGrid<CA::Real>> maps;                //!< The map with the indices. ATTENTION contains only the header, nodata.
      //std::vector<CA::Real>                times;               //!< The value/time sequences, one for each index.

      CA::Real currentNoData;



      //std::vector<CA::Real>       volumes;   //!< The total volume of rain of the last period for each index.
      CA::Real volume;
  };

public:

  //! Construct a STD var
  STRainDense(CA::Grid&  GRID, const SpatialTemporalDense& std);

  //! Destroy a ST var.
  ~STRainDense();

  //! Add the computational domain of the rain events into the given domain.
  void addDomain(CA::BoxList& compdomain);

  //! Analyse the area where the various rain event will fall.
  //! \attention This is going to be implemented in future.
  void analyseArea(CA::CellBuffReal& TMP, CA::CellBuffState& MASK, CA::BoxList&  domain);
 
  //! Prepare the rain events for the next update step considering the
  //! simulation time, the lenght of the update step and the next tim
  //! step.
  void prepare(CA::Real t, CA::Real period_time_dt, CA::Real next_dt);

  //! Add the amount of rain 
  void add(CA::CellBuffReal& WD, CA::CellBuffState& MASK, CA::Real t, CA::Real period_time_dt);

  //! Return the volume of rain of the last period_time_dt. 
  //! \attention This is the PERIOD volume.
  CA::Real volume();

  //! Compute the potential velocity that could happen in the next
  //! update/period step.
  //! This is used to limit the time step.
  CA::Real potentialVA(CA::Real t, CA::Real period_time_dt); 

    //! Return the simulation time when the events will not add any
  //! further water.
  CA::Real endTime();

 
private:

  //! Reference to the grid.
  CA::Grid& _grid;

  //! Reference to the Spatial/Temporal rain.
  const SpatialTemporalDense& _std;

  //! Spatial/Temporal rain.
  Data _data;

  CA::CellBuffReal RAIN_STD;

};




#endif

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

#ifndef _SPATIALTEMPORAL_HPP_
#define _SPATIALTEMPORAL_HPP_


//! \file Spatial.hpp
//!  Contains the structure(s) and classe(s) that are used to manage
//!  spatial/temporal variables.
//! \author Michele Guidolin, University of Exeter, 
//! contact: a.s.chen [at] exeter.ac.uk
//! \date 2024-05


#include"ca2D.hpp"
#include"BaseTypes.hpp"
#include"Box.hpp"
#include<string>
#include<vector>


//! Structure that identify a value time sequence. The time represent
//! the end of the value from t=0. For example in the case of the rain
//! a value/time pair 10/3600 represents a 100 mm/hr rain tha fall for
//! an hour. \attention The two list must be of the same size.
struct VT
{
  std::vector<CA::Real>  values;   //!< List of time varying values.
  std::vector<CA::Real>  times;    //!< List of end times.
};


//! Structure with the configuration value that define a spatial
//! temporal variable in the CA2D model. The spatial/temporal
//! variables, like roughness, rain and infiltration are usually
//! composed of two 'items': 1) a grid map which contains various
//! indices, and 2) a table where each index is assigned a value/time
//! sequence. These two items can be build by a spatial map using the
//! CAAPI tool spatialIndexingAG.
//! \warning For the ASCII/GRID grid map, the nodata value will be
//! considered to be zero. Thus there cannot be a zero index into the CVS.
struct SpatialTemporal
{
  std::string    name;		   //!< Name of the spatial variable. 
  CA::AsciiGrid<CA::State> map;	   //!< The map with the indices. ATTENTION contains only the header, nodata.
  std::vector<VT>          seqs;   //!< The value/time sequences, one for each index.
  CA::Unsigned             l;	   //!< The number of elements in the longest sequence. 
  CA::Real                 end;	   //!< The highest end of a time sequence .

  std::string              map_filename;      //!< The filename of the ascii grid file.
};


//! Initialise the spatial temporal structure. 
//! This first argument list must contains two files (in this order)
//! that are used to create the spatial/temporal variable:
//! 1) ARC/INFO ASCII GRID map containing the indices.
//! 2) CSV file with the roughness value/time sequence at various indices.  
//! For the ASCII/GRID grid map, the nodata value will be considered to be 
//! zero. Thus there cannot be a zero index into the CVS.
//! In the CSV file Each row represents a new "element" where the 
//! first column is the name of the element 
//! \warning The ASCI/GRID map dat is not laoded. Only the header.
//! \attention The order of elements is not important.
//! \param[in]  dir       The directory where the files are located.
//! \param[in]  filenames This is the list of two files where the data is read.
//! \param[out] st        The structure containing the read data.
//! \return A non zero value if there was an error.
int initSpatialTemporalFromCSV(const std::string& dir,
			       const std::vector<std::string>& filenames, SpatialTemporal& st);


//! Populate the CellBuffer with the indices of the spatial temporal
//! variable. The indices are loaded at specific bits positions and
//! mantain the value of rest of the bits.
//! \param mainid  The main id of the Preprocessed data to load.
//! \param subid   The main id of the Preprocessed data to load.
//! \return A non zero value if there was an error.
int populateSpatialTemporalIndices(CA::Grid&  GRID, const SpatialTemporal& st, 
				   CA::CellBuffState& DST,  CA::State start, CA::State stop,
				   const std::string& mainid="", const std::string& subid="");


//! Save the spatial temporal map as Preprocessed data
//! \return A non zero value if there was an error.
int saveSpatialTemporalMap(CA::Grid&  GRID, const SpatialTemporal& st, 
			   const std::string& mainid, const std::string& subid, bool loadMap = true);

//! directly load and populate the statial temporal maps in one go as opposed to saving temporal files
int loadAndPopulateSpatialTemporalMap(
	CA::Grid&  GRID,
	const SpatialTemporal& st,
	CA::CellBuffState& DST, CA::State start, CA::State stop,
	bool loadMap = true);

//! Class that manage the spatial/temporal roughness.
class STRough
{
private:
  
  //! The structure used during the model computation to store the
  //! spatial/temporal data.
  struct Data
  {
    std::vector<CA::Real>       values;    //!< The values of the indices for this time period.
    std::vector<size_t>         steps;	   //!< For each sequence, the time step in the sequence.
  };

public:

  //! Construct a ST var
  STRough(CA::Grid&  GRID, const SpatialTemporal& st);

  //! Destroy a ST var.
  ~STRough();
 
  //! Update the given table with the roughness value for each index
  //! of the given update period.
  //! The parameter rough is the default roughness. It used when the
  //! temporal domain is finished.
  //! \attention The values are inverted, i.e. 1/n for all the model
  //! \warning   A part the inertial model where the values are n^2
  void update(CA::TableReal& TAB, CA::Real t, CA::Real period_time_dt, CA::Real next_dt, 
	      CA::Real rough, MODEL::Type model);

 
private:

  //! Reference to the grid.
  CA::Grid& _grid;

  //! Reference to the Spatial/Temporal roughness.
  const SpatialTemporal& _st;

  //! Spatial/Temporal data.
  Data _data;


};


//! Class that manage the spatial/temporal rain. This is similar to
//! RainManager. The difference is that the rain is added unsing a
//! Table and a spatial map.
class STRain
{
private:
  
  //! The structure used during the model computation to store the
  //! spatial/temporal data.
  struct Data
  {
    std::vector<CA::Real>       values;    //!< The values of the indices for this time period.
    std::vector<size_t>         steps;	   //!< For each sequence, the time step in the sequence.

    std::vector<CA::Real>       volumes;   //!< The total volume of rain of the last period for each index.
  };

public:

  //! Construct a ST var
  STRain(CA::Grid&  GRID, const SpatialTemporal& st);

  //! Destroy a ST var.
  ~STRain();

  //! Add the computational domain of the rain events into the given domain.
  void addDomain(CA::BoxList& compdomain);

  //! Analyse the area where the various rain event will fall.
  //! \ttention This is going to be implemented in future.
  void analyseArea(CA::CellBuffReal& TMP, CA::CellBuffState& MASK, CA::BoxList&  domain);
 

  //! Update the given table with the rain value for each index of the
  //! given update period. The equaivalent function in RainManager is
  //! prepare().
  void update(CA::TableReal& TAB, CA::Real t, CA::Real period_time_dt, CA::Real next_dt);

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
  const SpatialTemporal& _st;

  //! Spatial/Temporal rain.
  Data _data;

};


//! Class that manage the spatial/temporal infiltration.
class STInf
{
private:
  
  //! The structure used during the model computation to store the
  //! spatial/temporal data.
  struct Data
  {
    std::vector<CA::Real>       values;    //!< The values of the indices for this time period.
    std::vector<size_t>         steps;	   //!< For each sequence, the time step in the sequence.
  };

public:

  //! Construct a ST var
  STInf(CA::Grid&  GRID, const SpatialTemporal& st);

  //! Destroy a ST var.
  ~STInf();
 
  //! Update the given table with the infiltration value for each index
  //! of the given update period. The values are the infiltration for the update dt.
  void update(CA::TableReal& TAB, CA::Real t, CA::Real period_time_dt, CA::Real next_dt);

  //! Return true if there is any infiltration to perform.
  bool needed() const;
 
private:

  //! Reference to the grid.
  CA::Grid& _grid;

  //! Reference to the Spatial/Temporal roughness.
  const SpatialTemporal& _st;

  //! Spatial/Temporal data.
  Data _data;

};


#endif

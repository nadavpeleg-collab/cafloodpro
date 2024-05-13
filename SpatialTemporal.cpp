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

//! \file SpatialTemporal.cpp
//! contact: a.s.chen [at] exeter.ac.uk
//! \date 2024-05


#include"ca2D.hpp"
#include"ArgsData.hpp"
#include"SpatialTemporal.hpp"
#include"Utilities.hpp"
#include<iostream>
#include<fstream>


// -------------------------//
// Include the CA 2D functions //
// -------------------------//
#include CA_2D_INCLUDE(writeSTIndices)


// Initialise the RainEvents structure usign a CSV file. 
int initSpatialTemporalFromCSV(const std::string& dir,
			       const std::vector<std::string>& filenames, SpatialTemporal& st)
{

  // Default value.
  st.l   = 0;
  st.end = 0;

  // The index zero should always exist and it usually reserved for
  // the nodata.
  st.seqs.resize(1);

  // If the vector is empty we don't do anything.
  if(filenames.size()==0)
    return 0;

  // However if the vector is not empty we nee to check that there is
  // the righ amount of elements.
  if(filenames.size()>0 && filenames.size()!=2)
  {
    std::cerr<<"Spatial Temporal variable setup must contains two files"<<std::endl;
    return 1;
  }

  try
  {
  // ATTENTION. Load only the header here .. not the actual data
  //CA::readAsciiGridHeader(st.map,dir+filenames[0]);
	  st.map.readAsciiGridHeader(dir + filenames[0]);
  }
  catch(std::runtime_error& e)
  {
    std::cerr<<e.what()<<std::endl;
    return 1;
  }

  // Save the Ascii Grid name.
  st.map_filename = dir+filenames[0];


  // Load the CSV file with the indices.
  std::string csvfile (dir+filenames[1]);
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
      
      st.name = CA::trimToken(str," \t\r");
    }

    if(CA::compareCaseInsensitive("Number Seq",tokens[0],true))
    {
      CA::Unsigned num = 0;
      READ_TOKEN(found_tok,num,tokens[1],tokens[0]);

      // Tell the vector of sequence that there is a good change that
      // there are num number of sequences (plus couple extra).
      st.seqs.reserve(num+2);
    }

    if(CA::compareCaseInsensitive("Value",tokens[0],true))
    {
      // Get the first token string without the "Value" substring.
      std::string str = tokens[0].substr(5, tokens[0].size() - 5);

      // Read the number of the index in the first token.
      CA::Unsigned num = 0;
      READ_TOKEN(found_tok,num,str,tokens[0]);

      // The index 0 is reserved.
      if(num == 0)
      {
	std::cerr<<"The index 0 of as spatial temporal variable is reserved on token: "<<tokens[0]<<std::endl;
	return 1;
      }

      // Check that the sequences vector container is big enough.  If
      // it is not the case then we need to increase it size for the
      // new index.
      if(num>=st.seqs.size())
	st.seqs.resize(num+1);

      // Get the sequence object from the vector of sequences.
      VT& vt( st.seqs[num] );

      // Read the sequence of values.
      for (size_t i=1; i<tokens.size(); ++i)
      {
	CA::Real value;
	READ_TOKEN(found_tok,value,tokens[i],tokens[0]);

	vt.values.push_back(value);
      }

      // Check if it is the longest sequence.
      st.l = std::max(vt.values.size(), st.l);
    }

    if(CA::compareCaseInsensitive("Time",tokens[0],true))
    {
      // Get the first token string without the "Time" substring.
      std::string str = tokens[0].substr(4, tokens[0].size() - 4);

      // Read the number of the index in the first token.
      CA::Unsigned num = 0;
      READ_TOKEN(found_tok,num,str,tokens[0]);

      // The index 0 is reserved.
      if(num == 0)
      {
	std::cerr<<"The index 0 of as spatial temporal variable is reserved on token: "<<tokens[0]<<std::endl;
	return 1;
      }

      // Check that the sequences vector container is big enough.  If
      // it is not the case then we need to increase it size for the
      // new index.
      if(num>=st.seqs.size())
	st.seqs.resize(num+1);

      // Get the sequence object from the vector of sequences.
      VT& vt( st.seqs[num] );

      // Read the sequence of values.
      for (size_t i=1; i<tokens.size(); ++i)
      {
	CA::Real value;
	READ_TOKEN(found_tok,value,tokens[i],tokens[0]);

	vt.times.push_back(value);

	// Check if it is the last end
	st.end = std::max(value, st.end);
      }    
    }
    
    // If the token was not identified stop!
    if(!found_tok)
    {
      std::cerr<<"Element '"<<CA::trimToken(tokens[0])<<"' not identified"<<std::endl; \
      return 1;
    }
  }

  // Now we have to check if the time part of each sequence has the
  // same number of items than the value part. If it is not the case
  // expand it with a possible maximum real value to indicate that the
  // value never stop.
  // \attention The index zero is reserved.
  for(size_t i = 1; i< st.seqs.size(); ++i)
  {
    VT& vt( st.seqs[i] );

    if(vt.values.size()>vt.times.size())
    {
      vt.times.resize(vt.values.size(), std::numeric_limits<CA::Real>::max());
      //std::cout<<"WARNING: Expanded the time part of the spatial Temporal sequence "<<i<<std::endl;

      st.end = std::numeric_limits<CA::Real>::max();
    }
  }

  return 0;
}


int populateSpatialTemporalIndices(CA::Grid&  GRID, const SpatialTemporal& st, 
				   CA::CellBuffState& DST, CA::State start, CA::State stop,
				   const std::string& mainid, const std::string& subid)
{
  // Check that there is indices to open.
  if(st.map_filename.empty() || st.map.ncols==0 || st.map.nrows==0 || st.seqs.size()==0)
    return 0;

  // Create a temporary cell buffer.
  CA::CellBuffState TMP(GRID);  

  // Se the default value of the temporary buffer to be nodata.
  TMP.fill(GRID.box(), st.map.nodata);
  
  // Try to load the pre-proc data.
  // If the pre-proc data do not exist use the given map file.
  if(mainid.empty() || subid.empty() || !TMP.loadData(mainid,subid))
  {   
    // Create a temporary ascii grid map
    CA::AsciiGrid<CA::State> map;
    
    // Load the temporary map with the roughness indices.
    map.readAsciiGrid(st.map_filename);
    
    // Find the box of the map. 
    // ATTENTIO This has not been fully tested!!!
    CA::Point     tl( CA::Point::create(GRID, map.xllcorner, map.yllcorner+map.nrows*map.cellsize ) );
    CA::Unsigned  w = map.ncols;
    CA::Unsigned  h = map.nrows;  
    CA::Box mapbox( CA::Box(tl.x(),tl.y(),w,h) );
    
    // limit the size to the maximum Grid box.
    mapbox.limit(GRID.box());
    
    //std::cout<<"MAPBOX ("<<mapbox.x()<<","<<mapbox.y()<<"):("<<mapbox.w()<<","<<mapbox.h()<<")" <<std::endl;
        
    // Insert the map indices into the temporary cell buffer.
    TMP.insertData(mapbox, &map.data[0], map.ncols, map.nrows);        
  }

  // Write the indices from the temporary cell buffer to the specific
  // bits of the MASK.  If the spatial map contains a nodata value,
  // this is replaced to zero since this is the reserved value.  If
  // the spatial map contains an index that is higher than a maximum
  // value, replace it with zero.
  CA::State max  = st.seqs.size()-1;
  CA::Execute::function(GRID.box(), writeSTIndices, GRID, DST, TMP, start, stop,  st.map.nodata, max);     

  return 0;
}

int saveSpatialTemporalMap(CA::Grid&  GRID, const SpatialTemporal& st, 
			   const std::string& mainid, const std::string& subid, bool loadMap)
{
  // Check that there is indices to open.
  if(st.map_filename.empty() || st.map.ncols==0 || st.map.nrows==0 || st.seqs.size()==0)
    return 0;
  
  // Create a temporary cell buffer.
  CA::CellBuffState TMP(GRID);  
  
  // Create a temporary ascii grid map
  CA::AsciiGrid<CA::State> map;
  
  // Load the temporary map with the roughness indices.
  //CA::readAsciiGrid(map, st.map_filename);
  if (loadMap) {
	  map.readAsciiGrid(st.map_filename);
  }else{
	  map.copyData(st.map);
  }
  
  // Find the box of the map. 
  // ATTENTIO This has not been fully tested!!!
  CA::Point     tl( CA::Point::create(GRID, map.xllcorner, map.yllcorner+map.nrows*map.cellsize ) );
  CA::Unsigned  w = map.ncols;
  CA::Unsigned  h = map.nrows;  
  CA::Box mapbox( CA::Box(tl.x(),tl.y(),w,h) );
  
  // limit the size to the maximum Grid box.
  mapbox.limit(GRID.box());
  
  //std::cout<<"MAPBOX ("<<mapbox.x()<<","<<mapbox.y()<<"):("<<mapbox.w()<<","<<mapbox.h()<<")" <<std::endl;
  
  // Se the default value of the temporary buffer to be nodata.
  TMP.fill(GRID.box(), map.nodata);
  
  // Insert the map indices into the temporary cell buffer.
  TMP.insertData(mapbox, &map.data[0], map.ncols, map.nrows);        

#ifdef CA2D_MPI
  // Save the data of the Elevation
  if(!TMP.saveData(mainid,subid))
#else
  // Save the data of the Elevation
  if (!TMP.saveData(mainid, subid))
#endif //CA2D_MPI
  {
    std::cerr<<"Error while saving the Spatial Temporal data"<<std::endl;
    return 1;
  }
  
  return 0;
}



int loadAndPopulateSpatialTemporalMap(
	CA::Grid&  GRID, 
	const SpatialTemporal& st,
	CA::CellBuffState& DST, CA::State start, CA::State stop,
	bool loadMap)
{

	// Check that there is indices to open.
	if (st.map_filename.empty() || st.map.ncols == 0 || st.map.nrows == 0 || st.seqs.size() == 0)
		return 0;

	// Create a temporary cell buffer.
	CA::CellBuffState TMP(GRID);

	// Create a temporary ascii grid map
	CA::AsciiGrid<CA::State> map;

	// Load the temporary map with the roughness indices.
	//CA::readAsciiGrid(map, st.map_filename);
	if (loadMap) {
		map.readAsciiGrid(st.map_filename);
	}
	else {
		map.copyData(st.map);
	}

	// Find the box of the map. 
	// ATTENTIO This has not been fully tested!!!
	CA::Point     tl(CA::Point::create(GRID, map.xllcorner, map.yllcorner + map.nrows*map.cellsize));
	CA::Unsigned  w = map.ncols;
	CA::Unsigned  h = map.nrows;
	CA::Box mapbox(CA::Box(tl.x(), tl.y(), w, h));

	// limit the size to the maximum Grid box.
	mapbox.limit(GRID.box());

	//std::cout<<"MAPBOX ("<<mapbox.x()<<","<<mapbox.y()<<"):("<<mapbox.w()<<","<<mapbox.h()<<")" <<std::endl;

	// Se the default value of the temporary buffer to be nodata.
	TMP.fill(GRID.box(), map.nodata);

	// Insert the map indices into the temporary cell buffer.
	TMP.insertData(mapbox, &map.data[0], map.ncols, map.nrows);


	// Write the indices from the temporary cell buffer to the specific
	// bits of the MASK.  If the spatial map contains a nodata value,
	// this is replaced to zero since this is the reserved value.  If
	// the spatial map contains an index that is higher than a maximum
	// value, replace it with zero.
	CA::State max = st.seqs.size() - 1;
	CA::Execute::function(GRID.box(), writeSTIndices, GRID, DST, TMP, start, stop, st.map.nodata, max);


	return 0;
}// end of loadAndPopulateSpatialTemporalMap




STRough::STRough(CA::Grid&  GRID, const SpatialTemporal& st):
  _grid(GRID),
  _st(st),
  _data()
{
  _data.values.resize(_st.seqs.size());
  _data.steps.resize(_st.seqs.size());
}

STRough::~STRough()
{
}


void STRough::update(CA::TableReal& TAB, CA::Real t, CA::Real period_time_dt, CA::Real next_dt, 
		     CA::Real rough, MODEL::Type model)
{
  switch(model)
  {
  case MODEL::Inertial:    
    // The inertial model use the square of the rough
    std::fill(_data.values.begin(), _data.values.end(), rough*rough);
    break;
  default:
    std::fill(_data.values.begin(), _data.values.end(), 1/rough);
    break;
  }

  // Cycle through each sequence.
  for(size_t i = 0; i<_st.seqs.size(); ++i)
  {    
    // Get the value time structu of this sequence.
    const VT& vt( _st.seqs[i] );

    // If there are no values in this sequence, go to the next
    // sequence.
    if(vt.values.empty() || vt.times.empty())
      continue;

    // Get the the step index of this sequence.
    size_t step = _data.steps[i];

    // Check if the simulation time now is equal or higher than the
    // time when this sequence steps ends. If it is the case,
    // increase the index to the next step.
    if(t >= vt.times[step])
      step++;

    // If the step is larger than the available value/time, do
    // nothing.
    if(step >= vt.values.size() )
      continue;

    // Get the value of the roughness.
    switch(model)
    {
    case MODEL::Inertial:    
      // The inertial model use the square of the rough
      _data.values[i] = vt.values[step]*vt.values[step];      
      break;
    default:
      _data.values[i] = 1/vt.values[step];
      break;
    }
    
    // Update step.
    _data.steps[i] = step;
  }

  // Upload the vector of values into the table.
  TAB.update(0,_data.values.size(),&_data.values[0],_data.values.size());
}



STRain::STRain(CA::Grid&  GRID, const SpatialTemporal& st):
  _grid(GRID),
  _st(st),
  _data()
{
  _data.values.resize(_st.seqs.size());
  _data.steps.resize(_st.seqs.size());
  _data.volumes.resize(_st.seqs.size());
}

STRain::~STRain()
{
}


void STRain::addDomain(CA::BoxList& compdomain)
{


	// NOTE: this needs improving as, currently adding the entire domain for every spatial area.

  // Cycle through each sequence. but ATTENTION we ignore the first one
  // since it is reserved for no rain.
  for(size_t i = 1; i<_st.seqs.size(); ++i)
  {    
    // Get the value time structu of this sequence.
    const VT& vt( _st.seqs[i] );

    // If there are no values in this sequence, go to the next
    // sequence.
    if(vt.values.empty() || vt.times.empty())
      continue;

    // Ok easy here ... we add all the full box domain to be sure.
    compdomain.add(_grid.box());
  }
  
}


void STRain::analyseArea(CA::CellBuffReal& TMP, CA::CellBuffState& MASK, CA::BoxList&  domain)
{

}


void STRain::update(CA::TableReal& TAB, CA::Real t, CA::Real period_time_dt, CA::Real next_dt)
{
  std::fill(_data.values.begin(), _data.values.end(), 0);

  // Cycle through each sequence.
  for(size_t i = 0; i<_st.seqs.size(); ++i)
  {    
    // Get the value time structu of this sequence.
    const VT& vt( _st.seqs[i] );

    // If there are no values in this sequence, go to the next
    // sequence.
    if(vt.values.empty() || vt.times.empty())
      continue;

    // Get the the step index of this sequence.
    size_t step = _data.steps[i];

    // Check if the simulation time now is equal or higher than the
    // time when this sequence steps ends. If it is the case,
    // increase the index to the next step.
    if(t >= vt.times[step])
      step++;

    // If the step is larger than the available value/time, do
    // nothing.
    if(step >= vt.values.size() )      
      continue;

    // Get the value of the rain (transformed in metres from mm) for each dt of the next period.
    _data.values[i] = ( vt.values[step] * 0.001) * (next_dt/3600.0);
    
    // Update step.
    _data.steps[i] = step;
  }

  // Upload the vector of values into the table.
  TAB.update(0,_data.values.size(),&_data.values[0],_data.values.size());
}


CA::Real STRain::volume()
{
  CA::Real rain_volume = 0.0;
  
  return rain_volume;
}


CA::Real STRain::potentialVA(CA::Real t, CA::Real period_time_dt)
{
  CA::Real potential_va = 0.0;

  // Use the maximum amount of rain that is going to fall in the next
  // period to calculate the possible water depth.  Use the critical
  // velocity with this water depth to compute the potential velocity.

  // Loop through the rain sequences.
  for(size_t i = 0; i<_st.seqs.size(); ++i)
  {
    // Get the value time structu of this sequence.
    const VT& vt( _st.seqs[i] );

    // If there are no values in this sequence, go to the next
    // sequence.
    if(vt.values.empty() || vt.times.empty())
      continue;


    // Get the the step index of this sequence.
    size_t step = _data.steps[i];

    // Check if the simulation time now is equal or higher than the
    // time when this sequence steps ends. If it is the case,
    // increase the index to the next step.
    if(t >= vt.times[step])
      step++;

    // If the step is larger than the available value/time, do
    // nothing.
    if(step >= vt.values.size() )
      continue;

    // Get the value of the rain (transformed in metres from mm) for next period.
    CA::Real rain = ( vt.values[step] * 0.001) * (period_time_dt/3600.0);
    
    // Compute the potential velocity.
    potential_va = std::max(potential_va, std::sqrt(static_cast<CA::Real>(9.81)*( rain )) );
  }  
  
  return potential_va;
}


CA::Real STRain::endTime()
{
  CA::Real t_end = 0.0;

  // Loop through the rain sequences.
  for(size_t i = 0; i<_st.seqs.size(); ++i)
  {
    // Get the value time structu of this sequence.
    const VT& vt( _st.seqs[i] );

    // If there are no values in this sequence, go to the next
    // sequence.
    if(vt.values.empty() || vt.times.empty())
      continue;

    // Loop through the time steps.
    for(size_t j = 0; j<vt.times.size(); j++)
    {
      // Check if it is still raining at this time. If yes then
      // updated end_t.
      if(vt.values[j]>0.0)
	t_end = std::max(t_end,vt.times[j]);
    }    
  }
  
  return t_end;
}



STInf::STInf(CA::Grid&  GRID, const SpatialTemporal& st):
  _grid(GRID),
  _st(st),
  _data()
{
  _data.values.resize(_st.seqs.size());
  _data.steps.resize(_st.seqs.size());
}

STInf::~STInf()
{
}


void STInf::update(CA::TableReal& TAB, CA::Real t, CA::Real period_time_dt, CA::Real next_dt)
{
  std::fill(_data.values.begin(), _data.values.end(), 0);

  // Cycle through each sequence.
  for(size_t i = 0; i<_st.seqs.size(); ++i)
  {    
    // Get the value time structu of this sequence.
    const VT& vt( _st.seqs[i] );

    // If there are no values in this sequence, go to the next
    // sequence.
    if(vt.values.empty() || vt.times.empty())
      continue;

    // Get the the step index of this sequence.
    size_t step = _data.steps[i];

    // Check if the simulation time now is equal or higher than the
    // time when this sequence steps ends. If it is the case,
    // increase the index to the next step.
    if(t >= vt.times[step])
      step++;

    // If the step is larger than the available value/time, do
    // nothing.
    if(step >= vt.values.size() )
      continue;

    // Get the value of the infiltration (transformed in metres from mm) for the next period time.
    _data.values[i] = ( vt.values[step] * 0.001) * (period_time_dt/3600.0);
    
    // Update step.
    _data.steps[i] = step;
  }

  // Upload the vector of values into the table.
  TAB.update(0,_data.values.size(),&_data.values[0],_data.values.size());
}


bool STInf::needed() const
{
  bool needed = false;

  // Loop through the rain sequences.
  for(size_t i = 0; i<_st.seqs.size(); ++i)
  {
    // Get the value time structu of this sequence.
    const VT& vt( _st.seqs[i] );

    // If there are no values in this sequence, go to the next
    // sequence.
    if(vt.values.empty() || vt.times.empty())
      continue;

    // Loop through the time steps.
    for(size_t j = 0; j<vt.times.size(); j++)
    {
      // Check if there is a value higher than zero.
      // updated end_t.
      if(vt.values[j]>0.0)
	needed = true;
    }    
  }
  
  return needed;
}

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

//! \file RasterGrid.cpp
//! contact: a.s.chen [at] exeter.ac.uk
//! \date 2024-05


#include"ArgsData.hpp"
#include"RasterGrid.hpp"
#include"Utilities.hpp"
#include<iostream>
#include<fstream>
#include<limits>

#include"Setup.hpp"


// -------------------------//
// Include the CA 2D functions //
// -------------------------//
#include CA_2D_INCLUDE(updatePEAKC)
#include CA_2D_INCLUDE(updatePEAKE)

#include CA_2D_INCLUDE(zeroedWD)
#include CA_2D_INCLUDE(zeroedVA)
#include CA_2D_INCLUDE(makeWL)


// Initialise the RasterGrid structure usign a CSV file. 
int initRasterGridFromCSV(const std::string& filename, RasterGrid& rg)
{
  std::ifstream ifile(filename.c_str());
  
  if(!ifile)
  {
    std::cerr<<"Error opening CSV file: "<<filename<<std::endl;
    return 1;
  }

  rg.pv     = PV::UNKNOWN;
  rg.peak   = false;
  rg.final  = false;
  rg.period = 0;
  
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

    if(CA::compareCaseInsensitive("Raster Grid Name",tokens[0],true))
    {
      std::string str;
      READ_TOKEN(found_tok,str,tokens[1],tokens[0]);
      
      rg.name = CA::trimToken(str);
    }

    if(CA::compareCaseInsensitive("Physical Variable",tokens[0],true))
      READ_TOKEN(found_tok,rg.pv,tokens[1],tokens[0]);

    if(CA::compareCaseInsensitive("Peak",tokens[0],true))
    {
      std::string str( CA::trimToken(tokens[1]) );
      READ_TOKEN(found_tok,rg.peak,str,tokens[0]);
    }

    if(CA::compareCaseInsensitive("Final",tokens[0],true))
    {
      std::string str( CA::trimToken(tokens[1]) );
      READ_TOKEN(found_tok,rg.final,str,tokens[0]);
    }

    if(CA::compareCaseInsensitive("Period",tokens[0],true))
      READ_TOKEN(found_tok,rg.period,tokens[1],tokens[0]);

    // If the token was not identified stop!
    if(!found_tok)
    {
      std::cerr<<"Element '"<<CA::trimToken(tokens[0])<<"' not identified"<<std::endl; \
      return 1;
    }
  }

  return 0;
}


int initRGData(const std::string& filename, CA::Grid& GRID, double nodata, const RasterGrid& rg, 
	       RGData& rgdata, RGPeak& rgpeak)
{

  rgdata.filename = filename;
  
  switch(rg.pv)
  {
  case PV::VEL:	 
    if(rgpeak.V.get() == 0)
    {
      rgpeak.V.reset( new CA::CellBuffReal(GRID) ); 
      rgpeak.V->clear(0.0);
    }
    // ATTENTION! The break is removed since in order to
    // post-process VA we need WD. 
    //break;
  case PV::WD:
  case PV::WL:
    if(rgpeak.WD.get() == 0)
    {
      rgpeak.WD.reset( new CA::CellBuffReal(GRID) ); 
      rgpeak.WD->clear(0.0);
    }
    break;
  default:
    break;
  }

  if(rg.period > 0.0)
    rgdata.time_next = rg.period;
  else
    rgdata.time_next = std::numeric_limits<CA::Real>::max(); 

  return 0;
}


RGManager::RGManager(CA::Grid&  GRID, std::vector<RasterGrid>& rgs, 
		     const std::string& base, std::vector<std::string> names):
  _grid(GRID),
  _rgs(rgs),
  _datas(rgs.size()),
  _peak()
{
  for(size_t i = 0; i<_rgs.size(); ++i)
  {
    std::string filename(base+"_"+names[i]);
    initData(filename, _rgs[i], _datas[i], _peak);
  }
}
  

RGManager::~RGManager()
{

}


cpp11::shared_ptr<CA::CellBuffReal> RGManager::getPeak_Velocities()
{
	return _peak.V;
}// end of getPeak_Velocities


cpp11::shared_ptr<CA::CellBuffReal> RGManager::getPeak_Depths()
{
	return _peak.WD;
}// end of getPeak_Depths


bool RGManager::updatePeak(const CA::BoxList&  domain, 
			   CA::CellBuffReal& WD, CA::CellBuffReal& V, CA::CellBuffState& MASK)
{
  // This variable make sure that the peak are updated only once.
  bool VAPEAKupdated = false;
  bool WDPEAKupdated = false;
  
  for(size_t i = 0; i<_datas.size(); ++i)
  {
    // Check if the peak values need to be updated.
    if(_rgs[i].peak == true)
    {
      switch(_rgs[i].pv)
      {
      case PV::VEL:
	// Update the absolute maximum velocity.
	// ATTENTION need to be tested.
	if(!VAPEAKupdated)
	{
	  CA::Execute::function(domain, updatePEAKC, _grid, (*_peak.V), V, MASK);
	  VAPEAKupdated = true;
	}
	// ATTENTION! The break is removed since in order to
	// post-process VA we need WD. 
	// break;	  
      case PV::WL:
      case PV::WD:	  
	// Update the absolute maximum water depth only once.
	if(!WDPEAKupdated)
	{
	  CA::Execute::function(domain, updatePEAKC, _grid, (*_peak.WD), WD, MASK);
	  WDPEAKupdated = true;
	}
	break;
      default:
	break;
      }
    }
  }

  return (WDPEAKupdated || VAPEAKupdated);
}


bool RGManager::outputPeak(CA::Real t, CA::CellBuffReal& WD, CA::CellBuffReal& V, 
			   const std::string& saveid,  bool output)
{
  // This variables is used to indicates if the output to console
  // happen in the case of time plot.
  bool outputed    = false;
  
  // Since WL variable does not exist. We need to save WD and then use
  // ELV to compute WL. The following variable is used to save WD only
  // once if both WD and WL are requested. The same is the case for
  // VEL, which needs WD.
  bool VAPEAKsaved   = false;
  bool WDPEAKsaved   = false;
  
  for(size_t i = 0; i<_datas.size(); ++i)
  {	
    // Check if the peak values need to be saved.
    if(_rgs[i].peak == true)
    {
      if(!outputed && output)
      {
	std::cout<<"Write Raster Grid (MIN "<<t/60<<"): ";
	outputed = true;
      }

      // Non velocity raster grid.
      switch(_rgs[i].pv)
      {
      case PV::VEL:	  
	if(!VAPEAKsaved)
	{
	  if(output)
	    std::cout<<" VAPEAK";
	  
	  _peak.V->saveData(saveid+"_V","PEAK");
	  VAPEAKsaved = true;
	}
	// ATTENTION! The break is removed since in order to
	// post-process VA we need WD. 
	// break;
      case PV::WL:
      case PV::WD:
	if(!WDPEAKsaved)
	{
	  if(output)
	    std::cout<<" WDPEAK";
	  
	  _peak.WD->saveData(saveid+"_WD","PEAK");
	  WDPEAKsaved = true;
	}
	break;
      default:
	break;
      }
      
    }
    // Update the next time to save a raster grid.
    _datas[i].time_next += _rgs[i].period;
  }

  if(outputed && output)
    std::cout<<std::endl;

  return (WDPEAKsaved || VAPEAKsaved);
}


bool RGManager::output(bool ouputPeak,CA::Real t, CA::Real dt, CA::CellBuffReal& WD,
		       CA::CellBuffReal& V, CA::CellBuffReal& A, CA::EdgeBuffReal &edgeVels,
		       const std::string& saveid,
		       bool output, bool final)
{
  // This variables is used to indicates if the output to console
  // happen in the case of time plot.
  bool outputed    = false;
  
  // Since WL variable does not exist. We need to save WD and then use
  // ELV to compute WL. The following variable is used to save WD only
  // once if both WD and WL are requested. The same is the case for
  // VEL, which needs WD.
  bool VAsaved       = false;
  bool VAPEAKsaved   = false;
  bool WDsaved       = false;
  bool WDPEAKsaved   = false;
  
  if (!ouputPeak) {
	  VAPEAKsaved = true;
	  WDPEAKsaved = true;
  }

  for(size_t i = 0; i<_datas.size(); ++i)
  {
    // Check if it is time to plot or the output is forced!
    if(t >= _datas[i].time_next || (final && _rgs[i].final) || (final && _rgs[i].peak))
    {	
      if(!outputed && output)
      {
	std::cout<<"Write Raster Grid (MIN "<<t/60<<"): ";
	outputed = true;
      }

      // Retrieve the string of the time.
      std::string strtime;
      CA::toString(strtime,std::floor(t+0.5));

      // Save the buffer using direct I/O where the main ID is the
      // buffer name and the subID is the timestep.
      switch(_rgs[i].pv)
      {
      case PV::VEL:	  
	if(!VAsaved)
	{
	  if(output)
	    std::cout<<" VA";	  
	  V.saveData(saveid+"_V",strtime);
	  A.saveData(saveid+"_A",strtime);

	  // save edge data too...
	  //edgeVels.saveData(saveid + "_EV", strtime);

	  // record the time step used, so we can decode to velocities in post processing
	  _rgs[i]._timeStepList[strtime] = dt;

	  VAsaved = true;
	}
	// ATTENTION! The break is removed since in order to
	// post-process VA we need WD. 
	//break;
      case PV::WL:
      case PV::WD:
	if(!WDsaved)
	{
	  if(output)
	    std::cout<<" WD";
	  
	  WD.saveData(saveid+"_WD",strtime);
	  WDsaved = true;
	}
	break;
      default:
	break;
      }
	
      // Check if the peak values need to be saved.
      if(_rgs[i].peak == true)
      {
	// Non velocity raster grid.
	switch(_rgs[i].pv)
	{
	case PV::VEL:	  
	  if(!VAPEAKsaved)
	  {
	    if(output)
	      std::cout<<" VAPEAK";
	  
	    _peak.V->saveData(saveid+"_V","PEAK");


		

	    VAPEAKsaved = true;
	  }
	  // ATTENTION! The break is removed since in order to
	  // post-process VA we need WD. 
	  // break;
	case PV::WL:
	case PV::WD:
	  if(!WDPEAKsaved)
	  {
	    if(output)
	      std::cout<<" WDPEAK";
	  
	    _peak.WD->saveData(saveid+"_WD","PEAK");
	    WDPEAKsaved = true;
	  }
	  break;
	default:
	  break;
	}
	  
      }
      // Update the next time to save a raster grid.
      _datas[i].time_next += _rgs[i].period;
    }
  }

  if(outputed && output)
    std::cout<<std::endl;

  return (WDPEAKsaved || VAPEAKsaved);
}


bool RGManager::writeOutput(bool ouputPeak, CA::Real t, CA::CellBuffReal& WD,
	CA::CellBuffReal& V, CA::CellBuffReal& A,
	bool output, bool final,
	CA::Grid& GRID, Setup setup, CA::BoxList& fulldomain, CA::Box& realbox, CA::AsciiGrid<CA::Real>& eg, CA::CellBuffState& MASK, CA::Real nodata, CA::CellBuffReal& ELV,const ArgsData& ad)
{

	// CREATE ASCII GRID TEMPORARY BUFFERS
	CA::AsciiGrid<CA::Real>& agtmp1 = eg;
	CA::AsciiGrid<CA::Real>  agtmp2 = eg;



	// This variables is used to indicates if the output to console
	// happen in the case of time plot.
	bool outputed = false;

	// Since WL variable does not exist. We need to save WD and then use
	// ELV to compute WL. The following variable is used to save WD only
	// once if both WD and WL are requested. The same is the case for
	// VEL, which needs WD.
	bool VAsaved = false;
	bool VAPEAKsaved = false;
	bool WDsaved = false;
	bool WDPEAKsaved = false;

	if (!ouputPeak) {
		VAPEAKsaved = true;
		WDPEAKsaved = true;
	}

	

	// Create two temporary Cell buffer
	CA::CellBuffReal TMP1(GRID);
	CA::CellBuffReal TMP2(GRID);

	// List of raster grid data
	std::vector<RGData> rgdatas(_rgs.size());

	// Peak buffers.
	RGPeak rgpeak;

	for (size_t i = 0; i<_rgs.size(); ++i)
	{
		std::string filename = ad.output_dir + ad.sdir + setup.short_name + "_" + setup.rastergrid_files[i];
		initRGData(filename, GRID, nodata, _rgs[i], rgdatas[i], rgpeak);
	}

	for (size_t i = 0; i<_datas.size(); ++i)
	{
		// Check if it is time to plot or the output is forced!
		if (t >= _datas[i].time_next || (final && _rgs[i].final))
		{
			if (!outputed && output)
			{
				std::cout << "Write Raster Grid (MIN " << t / 60 << "): ";
				outputed = true;
			}

			// Retrieve the string of the time.
			std::string strtime;
			CA::toString(strtime, std::floor(t + 0.5));

			std::string filename;

			// Save the buffer using direct I/O where the main ID is the
			// buffer name and the subID is the timestep.
			switch (_rgs[i].pv)
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
						filenameV = rgdatas[i].filename + "_V_" + strtime + ".asc";
						filenameA = rgdatas[i].filename + "_A_" + strtime + ".asc";
					}

					TMP1.copy(V);
					TMP2.copy(A);


					// Set the V and A to zero if water depth is less than tollerance.
					CA::Execute::function(fulldomain, zeroedVA, GRID, TMP1, TMP2, WD, MASK, setup.rast_wd_tol);
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
				filename = rgdatas[i].filename + "_" + strtime + ".asc";

				// Reset the buffer
				TMP1.fill(fulldomain, nodata);

				// Make the water depth data and elevation into the water level.
				CA::Execute::function(fulldomain, makeWL, GRID, TMP1, WD, ELV, MASK);

				// Retrieve the data
				TMP1.retrieveData(realbox, &agtmp1.data[0], agtmp1.ncols, agtmp1.nrows);

				if (setup.output_console)
					std::cout << "Write Raster Grid: " << filename << std::endl;

				// Write the data.
				//CA::writeAsciiGrid(agtmp1,filename,setup.rast_places);	      
				agtmp1.writeAsciiGrid(filename, setup.rast_places);


			case PV::WD:



				// Create the name of the file.
				//std::string filename = removeExtension(rgdatas[i].filename) + "_" + strtime + ".asc";	
				filename = rgdatas[i].filename + "_" + strtime + ".asc";

				// Water depth already loaded.

				// Retrieve the data
				WD.retrieveData(realbox, &agtmp1.data[0], agtmp1.ncols, agtmp1.nrows);

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
			_datas[i].time_next += _rgs[i].period;
		}


	}

	

	if (outputed && output)
		std::cout << std::endl;

	return (WDPEAKsaved || VAPEAKsaved);
}



int RGManager::initData(const std::string& filename, const RasterGrid& rg, Data& rgdata, Peak& rgpeak)
{
  rgdata.filename = filename;
  
  switch(rg.pv)
  {
  case PV::VEL:	 
    if(rgpeak.V.get() == 0)
    {
      rgpeak.V.reset( new CA::CellBuffReal(_grid) ); 
      rgpeak.V->clear(0.0);
    }
    // ATTENTION! The break is removed since in order to
    // post-process VA we need WD. 
    //break;
  case PV::WD:
  case PV::WL:
    if(rgpeak.WD.get() == 0)
    {
      rgpeak.WD.reset( new CA::CellBuffReal(_grid) ); 
      rgpeak.WD->clear(0.0);
    }
    break;
  default:
    break;
  }

  if(rg.period > 0.0)
    rgdata.time_next = rg.period;
  else
    rgdata.time_next = std::numeric_limits<CA::Real>::max(); 

  return 0;
}

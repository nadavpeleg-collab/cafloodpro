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

#ifndef _SIMULATIONSTORE_HPP_
#define _SIMULATIONSTORE_HPP_


//! \file SimulationStore.hpp
//!  Contains the structure which contains the setup details and results
//! \author Michael Gibson, University of Exeter, 
//! contact: m.j.gibson [at] exeter.ac.uk
//! \date 2015-06

//#include <boost/thread.hpp>


#include"Setup.hpp"
#include"Arguments.hpp"
#include "AsciiGrid.hpp"
//#include "ESRI_ASCIIGrid.hpp"

#include "SpatialTemporal.hpp"
#include "SpatialTemporalDense.hpp"
#include "Rain.hpp"
#include "WaterLevel.hpp"
#include "Inflow.hpp"
#include "TimePlot.hpp"
#include "RasterGrid.hpp"
#include <thread>



//! Structure with the data that define the main variables used to
//! execute the CA2D model, and the results
struct SimulationStore
{

	bool abortIt;

	//boost::thread* workerPtr;
	std::thread* workerPtr;

	int currentStatus;

	enum STATUS{
		INITIALISED,
		PREPROSSSING,
		RUNNING,
		POSTPROCESSING,
		FINISHED,
		ABORTED,
		FAILED
	};

	enum RESULTCODE{

		SUCCESS,
		FLOOD_IN_NODATA_CELL,
		FLOOD_ADJECT_TO_NODATA, // double check sense ... as water shouldn't be inside a no data cell...
		BURST_ON_NODATA,
		FLOOD_REACH_EDGE_TERRAIN

	};


	double currentTime;

	//------------------------------------------------------------------------------------------
	// constructor
	SimulationStore():
		ad(), setup(), sts(3), res(), wles(), ies(), currentStatus(INITIALISED), currentTime(0.0), workerPtr(NULL), abortIt(false),
		tps(), rgs()//, spatioTemporalDense_rain()
	{
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// ad

		ad.info = true;
		ad.model = "sim";
		ad.pre_proc = true;
		ad.post_proc = true;

		


		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		//setup.



		setup.model_type = MODEL::WCA2Dv2;

		setup.roughness_global = 0.02;

		setup.ignore_wd = 0.001;
		setup.tolerance = 0.001;

		setup.time_updatedt = 60.0;
		
		setup.output_period = 600;
		setup.output_console = true;

		
		
		setup.remove_data = true;
		setup.remove_prec_data = true;
		setup.rast_vel_as_vect = true;

		//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		setup.rast_wd_tol = 0.001;
		//setup.rast_wd_tol = 0.0;

		setup.rast_places = 6;

		setup.upstream_reduction = 1;
		setup.decompose_eff = 0.5;

		//--------------------

		setup.time_alpha = 0.05;
		setup.time_maxiters = 100000000;
		setup.boundary_elv = -1000.0;

		//setup.update_peak_dt = false;
		setup.expand_domain = true;
		setup.check_vols = true;
		setup.output_computation = true;

		//

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

		// default bands values
		setup.bands = new double[1]();
		setup.bands[0] = 1;

		

	}// end of default constructor
	//------------------------------------------------------------------------------------------


	//------------------------------------------------------------------------------------------
	// destructor
	~SimulationStore() 
	{
		// ensure pointer memory is deallocated
		if (timeGrid != NULL){
			delete (timeGrid);
		}
		if (setup.bands != NULL) {
			delete (setup.bands);
		}
		if (setup.startingDepths != NULL) {
			delete (setup.startingDepths);
		}


	}// end of default destructor
	//------------------------------------------------------------------------------------------


	//  --- ARGUMENTS  --- 
	ArgsData ad;
	
	//  --- SETUP  --- 
	Setup setup;

	// --- SETUP DATA  --- 
	
	std::vector<SpatialTemporal> sts;

	SpatialTemporalDense spatioTemporalDense_rain;

	std::vector<RainEvent> res;
	std::vector<WLEvent> wles;
	std::vector<IEvent> ies;
	std::vector<TimePlot> tps;
	std::vector<RasterGrid> rgs;

	// create a pointer to a cell buffer for time output grid
	// in order to minimise memory usuage, using a pointer, which is only declared if appropriated option is selected
	CA::CellBuffReal* timeGrid = NULL;
	double thesholdDepth = 0.0;

	CA::State resultCode = -1;
	
	
	

	


};// end of SimulationStore



#endif

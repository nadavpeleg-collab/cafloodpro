/*
    
Copyright (c) 2015 Centre for Water Systems,
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

//! \file application.h
//! primary applcation function header file of the cafloodpro program that performs fast 2D flood
//! modelling using a cellular automata algorithm and the CAAPI
//! framework.
//! contact: m.j.gibson [at] exeter.ac.uk
//! \date 2015-06

#ifdef CAFUNCSDLL_EXPORTS
#define CAFUNCSDLL_API extern "C" __declspec(dllexport) 
#else
#define CAFUNCSDLL_API extern "C" __declspec(dllimport) 
#endif

// version 0

//CAFUNCSDLL_API int application(int argc, char* argv[]);


//------------------------------------------------------------------------------------------------------------
// version 1
//------------------------------------------------------------------------------------------------------------
// type defs



#ifdef _WIN64
	#if defined __MINGW32__ || defined __MINGW64__
    #include <stdint.h>
        typedef int64_t CA_simulation_mem;
	#else
        typedef _int64 CA_simulation_mem;
    #endif
#else
    #if defined __MINGW32__ || defined __MINGW64__
        typedef int CA_simulation_mem;
	#else
        typedef _int32 CA_simulation_mem;
    #endif
#endif



//------------------------------------------------------------------------------------------------------------
// erros

#define NO_ERRORS			0
#define NOT_ENOUGH_MEMORY	1



#define CANNOT_READ_SETUP_FILE	2

#define SIM_NOT_READY			3

#define NO_LICENSE_FILE			4

#define NO_GDAL_BUILD			5
#define NO_SUCH_DRIVER			6
#define UNKNOWN_EPSG_CODE		7
#define NO_CREATE_COPY			8


//------------------------------------------------------------------------------------------------------------
// functions
// create the simulation data structure and the return the handle to be used by other function to refference
CAFUNCSDLL_API int CADDIES_createSim(CA_simulation_mem& simulation_Handle);


CAFUNCSDLL_API int CADDIES_getGDALDriverCount(int& driverCount);

CAFUNCSDLL_API int CADDIES_getGDALDriver(
	char*& driver, int driverIndex
	);


CAFUNCSDLL_API int CADDIES_setSim_InputOutputDirs(CA_simulation_mem simulation_Handle,
	char* inputDir,
	char* outputDir,
	bool verbose
	);


//CAFUNCSDLL_API int setSim_Inputs(CA_simulation_mem simulation_Handle);
CAFUNCSDLL_API int CADDIES_setSim_Inputs(CA_simulation_mem simulation_Handle,
	char* inputFile,
	bool verbose
	);


//#ifdef CAAPI_GDAL_OPTION


CAFUNCSDLL_API int CADDIES_setGDAL_settings(CA_simulation_mem simulation_Handle,
	char* driverName,
	int epsg
	);


//#endif



CAFUNCSDLL_API int CADDIES_setSim_Names(CA_simulation_mem simulation_Handle,
	char* name,
	char* short_name,
	bool verbose);


CAFUNCSDLL_API int CADDIES_setSim_times(CA_simulation_mem simulation_Handle,
	double start_time,
	double end_time,
	double min_timeStep,
	double max_timeStep,
	double update_dt,
	double update_dt_alpha,
	bool verbose);


// add inflow - via file
CAFUNCSDLL_API int CADDIES_setSim_Inflow_File(CA_simulation_mem simulation_Handle,char* inputFilePath,bool verbose);

// add inflow via paramters
CAFUNCSDLL_API int CADDIES_setSim_Inflow(CA_simulation_mem simulation_Handle,
	char* inflowName,
	double* timeValuesArr,
	double* levelValuesArr,
	unsigned int numOfValues,
	double X,
	double Y,
	bool verbose);


CAFUNCSDLL_API int CADDIES_setSim_RainEvent(CA_simulation_mem simulation_Handle,
	char* eventName,
	double* timeValuesArr,
	double* levelValuesArr,
	unsigned int numOfValues,
	double tlX,
	double tlY,
	double brX,
	double brY,
	bool fullArea,
	bool verbose);

CAFUNCSDLL_API int CADDIES_setSim_WaterLevelEvent(CA_simulation_mem simulation_Handle,
	char* eventName,
	double* timeValuesArr,
	double* levelValuesArr,
	unsigned int numOfValues,
	double tlX,
	double tlY,
	double width,
	double height,
	bool verbose);

CAFUNCSDLL_API int CADDIES_setSim_TimeToThreahold_output(CA_simulation_mem simulation_Handle,
	char* eventName,
	double thresholdDepth,
	bool verbose);



CAFUNCSDLL_API int CADDIES_setSim_Bands(CA_simulation_mem simulation_Handle,
	double* bands,
	unsigned int numberOfBands,
	bool verbose);



//------------------------------------------------------------------------------------------------------------

// set the file to be loaded
CAFUNCSDLL_API int CADDIES_setSim_terrain_ascii(CA_simulation_mem simulation_Handle,
	char* inputFilePath,
	bool verbose);


// set a bunch of files to be loaded, merged (and cut) and used the merged
CAFUNCSDLL_API int CADDIES_setSim_create_terrainExtent_ascii(CA_simulation_mem simulation_Handle,
	char* inputFilePath,
	unsigned int ncols,
	unsigned int nrows,
	double xllcorner,
	double yllcorner,
	double cellsize,
	double nodata,

	unsigned int numberOfTerrains,
	char** inputFilePaths,

	bool verbose);

CAFUNCSDLL_API int CADDIES_setSim_SetTerrainData(CA_simulation_mem simulation_Handle,
	
	unsigned int ncols,
	unsigned int nrows,
	double xllcorner,
	double yllcorner,
	double cellsize,
	double nodata,

	double* data,

	bool verbose);

CAFUNCSDLL_API int CADDIES_setSim_SetInitialWaterDepths(CA_simulation_mem simulation_Handle,

	unsigned int ncols,
	unsigned int nrows,
	double xllcorner,
	double yllcorner,
	double cellsize,
	double nodata,

	double* data,

	bool verbose);



//CAFUNCSDLL_API int CADDIES_setSim_addTerrainData(CA_simulation_mem simulation_Handle,
//	unsigned int numberOfTerrains,
//	char** inputFilePaths,
//------------------------------------------------------------------------------------------------------------
CAFUNCSDLL_API int CADDIES_setSim_SetSpatialRain(CA_simulation_mem simulation_Handle,

	unsigned int ncols,
	unsigned int nrows,
	double xllcorner,
	double yllcorner,
	double cellsize,
	double nodata,

	unsigned int* data,

	double* flowRates,
	double* times,

	unsigned int* numberOfElementsPerSequence,

	unsigned int numberOfSequences,

	bool verbose);

//------------------------------------------------------------------------------------------------------------
CAFUNCSDLL_API int CADDIES_setSim_SetSpatialRainDense(CA_simulation_mem simulation_Handle,

	unsigned int ncols,
	unsigned int nrows,
	double xllcorner,
	double yllcorner,
	double cellsize,
	double nodata,

	double* flowRates,
	double* times,

	unsigned int numberOfSequences,

	bool verbose);


CAFUNCSDLL_API int CADDIES_setSim_SetSpatialInfiltration(CA_simulation_mem simulation_Handle,

	unsigned int ncols,
	unsigned int nrows,
	double xllcorner,
	double yllcorner,
	double cellsize,
	double nodata,

	unsigned int* data,

	double* infiltrationRates,

	double* times,

	unsigned int* numberOfElementsPerSequence,

	unsigned int numberOfSequences,

	bool verbose);

CAFUNCSDLL_API int CADDIES_setSim_SetRoughness(CA_simulation_mem simulation_Handle,

	unsigned int ncols,
	unsigned int nrows,
	double xllcorner,
	double yllcorner,
	double cellsize,
	double nodata,

	double* data,

	double* roughnessFactors,

	double* times,

	unsigned int* numberOfElementsPerSequence,

	unsigned int numberOfSequences,

	bool verbose);

CAFUNCSDLL_API int CADDIES_setSim_SetAdditionalHeights(CA_simulation_mem simulation_Handle,

	unsigned int ncols,
	unsigned int nrows,
	double xllcorner,
	double yllcorner,
	double cellsize,
	double nodata,

	unsigned int* data,

	double* addedHeights,


	unsigned int numberOfElements,

	bool verbose);




//------------------------------------------------------------------------------------------------------------

CAFUNCSDLL_API int CADDIES_setSim_add_raster_output(CA_simulation_mem simulation_Handle,
	char* inputFilePath,
	char* inputFile,
	bool verbose);

CAFUNCSDLL_API int CADDIES_setSim_outputs_Depth_Max(CA_simulation_mem simulation_Handle,
	char* depthsFileName,
	bool verbose);

CAFUNCSDLL_API int CADDIES_setSim_outputs_Level_Max(CA_simulation_mem simulation_Handle,
	char* levelFileName,
	bool verbose);

CAFUNCSDLL_API int CADDIES_setSim_outputs_Vel_Max(CA_simulation_mem simulation_Handle,
	char* velocitiesFileName,
	bool verbose);

//


CAFUNCSDLL_API int CADDIES_setSim_outputs_Depth_Max_banded(CA_simulation_mem simulation_Handle,
	char* depthsMaxbandedFileName,
	bool verbose);




//------------------------------------------------------------------------------------------------------------


CAFUNCSDLL_API int CADDIES_setSim_outputs_PointTimePlot(CA_simulation_mem simulation_Handle,
	char* name,
	char* long_name,
	char* type,
	unsigned int numOfValues,
	double* X,
	double* Y,
	char** point_Names,
	double interval,
	bool verbose);




//------------------------------------------------------------------------------------------------------------

CAFUNCSDLL_API int CADDIES_setSim_outputs_Depth_times(CA_simulation_mem simulation_Handle,
	char* depthsFileName,
	double times,
	bool verbose);

CAFUNCSDLL_API int CADDIES_setSim_outputs_Level_times(CA_simulation_mem simulation_Handle,
	char* depthsFileName,
	double times,
	bool verbose);

CAFUNCSDLL_API int CADDIES_setSim_outputs_Vel_times(CA_simulation_mem simulation_Handle,
	char* depthsFileName,
	double times,
	bool verbose);

//------------------------------------------------------------------------------------------------------------


CAFUNCSDLL_API int CADDIES_setSim_extraParams(CA_simulation_mem simulation_Handle,
	double boundaryElevation,
	double ignore_WD,
	double tolerance,
	double roughnessGlobal,
	int maxIters,
	bool expandedDomain,
	bool verbose);



//------------------------------------------------------------------------------------------------------------

CAFUNCSDLL_API int CADDIES_runSim(CA_simulation_mem simulation_Handle);

CAFUNCSDLL_API int CADDIES_runSim2(CA_simulation_mem simulation_Handle);

CAFUNCSDLL_API int CADDIES_runSim_locking(CA_simulation_mem simulation_Handle);

CAFUNCSDLL_API int CADDIES_runSim_locking2(CA_simulation_mem simulation_Handle);

//------------------------------------------------------------------------------------------------------------


CAFUNCSDLL_API int CADDIES_SimStatus(CA_simulation_mem simulation_Handle,
	int& currentStatus,
	double& currentTime
	);

CAFUNCSDLL_API int CADDIES_getResultCode(CA_simulation_mem simulation_Handle,
	int& code
	);


CAFUNCSDLL_API int CADDIES_Sim_ReIntialise(CA_simulation_mem simulation_Handle);

CAFUNCSDLL_API int CADDIES_SimAbort(CA_simulation_mem simulation_Handle);

//------------------------------------------------------------------------------------------------------------


CAFUNCSDLL_API int CADDIES_freeSim(CA_simulation_mem simulation_Handle);

//------------------------------------------------------------------------------------------------------------

CAFUNCSDLL_API int CADDIES_printSim_detail(CA_simulation_mem simulation_Handle);


//------------------------------------------------------------------------------------------------------------
// test

CAFUNCSDLL_API int CADDIES_Test(bool verbose);

//------------------------------------------------------------------------------------------------------------

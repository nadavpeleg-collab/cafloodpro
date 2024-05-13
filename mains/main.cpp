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

//! \file main.cpp
//! Main file of the cafloodpro program that performs fast 2D flood
//! modelling using a cellular automata algorithm and the CAAPI
//! framework.
//! contact: a.s.chen [at] exeter.ac.uk
//! \date 2024-05


#include"ca2D.hpp"
#include"Arguments.hpp"
#include"Options.hpp"
#include"ArgsData.hpp"
#include"Setup.hpp"
#include"Rain.hpp"
#include"Inflow.hpp"
#include"WaterLevel.hpp"
#include"TimePlot.hpp"
#include"RasterGrid.hpp"
#include"SpatialTemporal.hpp"
#include"SpatialTemporalDense.hpp"
#include<ctime>

// ADD Windows header file. This is meanly used for the API that
// manage the power settings.
#if defined _WIN32 || defined __CYGWIN__   
#include<windows.h>
#endif

#ifdef CAAPI_GDAL_OPTION
//#include "GDAL_Functions.hpp"
#include "gdal_priv.h"

#endif


//! Print the version info to std output.
inline void version()
{
 std::cout<<"Copyright 2024 Centre for Water Systems, University of Exeter"<<std::endl;
 std::cout<<"App                  : "<<CA_QUOTE_MACRO(CAAPI_APP_NAME)<<" ver. "<<CAAPI_APP_VERSION<<std::endl;
 std::cout<<"CA API Version       : "<<caVersion<<std::endl;
 std::cout<<"       Impl Name     : "<<caImplName<<std::endl;
 std::cout<<"       Impl Short    : "<<caImplShortName<<std::endl;
 std::cout<<"       Impl Version  : "<<caImplVersion<<std::endl;
 std::cout<<"       Impl Precision: "<<caImplPrecision<<std::endl;
}


//! Perform the pre-processing using the given parameters (see pre_proc.cpp).
//! \param[in] ad       The arguments data.
//! \param[in] setup    The setup of the simulation.
//! \param[in] ele_file The elevation file.
//! \param[in] sts      The list of spatial temporal imputs.
//! \return A non zero value if there was an error.
int preProc(const ArgsData& ad, const Setup& setup, const std::string& ele_file,
	    const std::vector<SpatialTemporal>& sts,
        SpatialTemporalDense &_std
    );


//! Perform the post processing of the data for a CA 2D model. 
//! \param[in] ad       The arguments data.
//! \param[in] setup    The setup of the simulation.
//! \param[in] eg       The elevation grid.
//! \param[in] tps      The list of time plot outputs.
//! \param[in] rgs      The list of raster grid outputs.
//! \return A non zero value if there was an error.
int postProc(const ArgsData& ad, const Setup& setup, CA::AsciiGrid<CA::Real>& eg,
	const std::vector<TimePlot>& tps, std::vector<RasterGrid>& rgs,
	CA::CellBuffReal* timeGrid);


//! Perform the CADDIES2D flood modelling algorithm using the given
//! parameters (see CADDIES2D.cpp).
//! This function perform the following models:
//! 1) WCA2Dv1
//! 2) WCA2Dv2
//! 3) Inertial
//! 4) Doffusive
//! \param[in] ad       The arguments data.
//! \param[in] setup    The setup of the simulation.
//! \param[in] eg       The elevation grid.
//! \param[in] sts      The list of spatial temporal imputs.
//! \param[in] res      The list of rain event inputs.
//! \param[in] wles     The list of water level event inputs.
//! \param[in] ies      The list of inflow event inputs.
//! \param[in] tps      The list of time plot outputs.
//! \param[in] rgs      The list of raster grid outputs.
//! \return A non zero value if there was an error.
int CADDIES2D(const ArgsData& ad, Setup& setup, CA::AsciiGrid<CA::Real>& eg,
	const std::vector<SpatialTemporal>& sts,
	const std::vector<RainEvent>& res, const std::vector<WLEvent>& wles,
	const std::vector<IEvent>& ies,
	const std::vector<TimePlot>& tps, std::vector<RasterGrid>& rgs,
	CA::CellBuffReal**  timeGrid, double& thresholdDepth,
	double& currentTime, CA::State& resultcode,
	bool& abortIt,
    SpatialTemporalDense& _std,
    bool preProcessing
    );


//! Return the terrrain info of the CA2D model to simulate. 
//! \attention The preProc function should be called before this one.
//! \warning If "Remove Pre-proc data is true, this function remove them."
int terrainInfo(const ArgsData& ad,Setup& setup, const CA::AsciiGrid<CA::Real>& eg);


int main(int argc, char* argv[])
{

#if defined _WIN32 || defined __CYGWIN__   
  // Stop WINDOWS from going to sleep!
  SetThreadExecutionState(ES_CONTINUOUS | ES_SYSTEM_REQUIRED );
#endif

  // Do not sync with stdio
  std::cout.sync_with_stdio(false);
  std::cerr.sync_with_stdio(false);

  // Initialise the 2D caAPI.
  CA::init2D(&argc,&argv);

  // Create the argument structure.
  ArgsData ad;

  // Create the setup structure.
  Setup setup;

  // The number of argument.
  CA::Unsigned na = 0;

  // There is only three compulsory arguments which are the directory
  // where the configuration files and data is found, the initial
  // setup file of the given CA algorithm, and the directory where the
  // output data is saved. The setup file should be in the working
  // directory.
  ad.args.add(na++,"input-dir",   "The input dir which contain the setup file and the data","",  false);
  ad.args.add(na++,"setup-file",  "The setup file of the ca algorithm","",  false);
  ad.args.add(na++,"output-dir",  "The output dir where output files are saved","",  false);

  // Add optional arguments.
  ad.args.add(na++,"test",        "Simple test that check if the executable work.", "",  true, false, true);
  ad.args.add(na++,"version",     "Show the version of the code.", "",  true, false, true);
  ad.args.add(na++,"help",        "Display the help and exit.", "",  true, false, true);
  ad.args.add(na++,"info",        "Print information about the configuration on console",  "",  true, false);
  ad.args.add(na++,"pre-proc",    "Perform the pre-processing of the data (only)",         "",  true, false);
  ad.args.add(na++,"post-proc",   "Perform the post-processing of the data (only)",         "",  true, false);
  ad.args.add(na++,"no-pre-proc", "Do NOT perform the pre-processing of the data",         "",  true, false);
  ad.args.add(na++,"WCA2D",       "Perform the WCA2D flood model (deprecated)",            "",  true, false);
  ad.args.add(na++,"sim",         "Perform the flood model simulation",                    "",  true, false);
  ad.args.add(na++,"terrain-info","Display the terrain info and exit.", "",  true, false, false);
  ad.args.add(na++, "ctime",		"Display the time and date and exit.", "", true, false, false);

#ifdef CAAPI_GDAL_OPTION

  ad.args.add(na++, "driver",	"The GDAL output driver to use", "", true, true, false);
  ad.args.add(na++, "epsg",		"The Geo Co-ordinate epsg system to use", "", true, true, false);
  ad.args.add(na++, "formats",  "List the available GDAL drivers", "", true, false, false);


#endif



  // Add the options from the CA implementation
  ad.args.addList(CA::options());

  // Parse the argument list, if the parse fail (it returns false) then exit!
  if(ad.args.parse(argc,argv,std::cout) == false )
  {
#if defined _WIN32 || defined __CYGWIN__   
    std::cout<<"Use the /help options to show help"<<std::endl;
    //std::cout<<"Press 'Return' to continue"<<std::endl;
    //std::cin.get();
#else
    std::cout<<"Use the -help options to show help"<<std::endl;
    //std::cout<<"Press 'Return' to continue"<<std::endl;
    //std::cin.get();
#endif
    return EXIT_FAILURE;    
  }
  
  // Cycle throught the arguments that were activated
  for(CA::Arguments::List::const_iterator i = ad.args.active().begin(); i != ad.args.active().end();  ++i)
  {
    if((*i)->tag == 0)		// Input dir.
    {
      ad.working_dir = (*i)->value;
    }

    if((*i)->tag == 1)		// setup file.
    {
      ad.setup_file = (*i)->value;
    }

    if((*i)->tag == 2)		// Ouptut dir.
    {
      ad.output_dir = (*i)->value;
    }

    if((*i)->name == "help")
    {
      ad.args.help(std::cout,true);
      return EXIT_SUCCESS;
    }

    if((*i)->name == "version")
    {
      version();
      return EXIT_SUCCESS;
    }
    
    if((*i)->name == "test")
    {
      std::cout<<CAAPI_APP_VERSION<<std::endl;
      return EXIT_SUCCESS;
    }

    if((*i)->name == "info")
      ad.info = true;

    if((*i)->name == "pre-proc")
      ad.pre_proc = true;

    if((*i)->name == "post-proc")
      ad.post_proc = true;

    if((*i)->name == "no-pre-proc")
      ad.no_pre_proc = true;

    if((*i)->name == "WCA2D" || (*i)->name == "sim")
    {
      ad.pre_proc = true;
      ad.model = "sim";
      ad.post_proc = true;
    }

    if((*i)->name == "terrain-info")
    {
      ad.pre_proc = true;
      ad.terrain_info = true;
    }

	if ((*i)->name == "ctime")
	{
		//std::cout << CAAPI_APP_VERSION << std::endl;

		std::time_t t = std::time(0);   // get time now
		std::tm* now = std::localtime(&t);
		std::cout << "DD-MM-YYYY::"
			<< now->tm_mday << '-'
			<< (now->tm_mon + 1) << '-'
			<< (now->tm_year + 1900)
			<< std::endl;

		return EXIT_SUCCESS;
	}


#ifdef CAAPI_GDAL_OPTION

	if ((*i)->name == "driver")
	{
		//setup.gdalDriverName = (*i)->value;
        GDALAllRegister();

        CA::ESRI_ASCIIGrid<CA::Real>::driver = (*i)->value;

        // check driver code exists
        GDALDriver* poDriver;
        poDriver = GetGDALDriverManager()->GetDriverByName(CA::ESRI_ASCIIGrid<CA::Real>::driver.c_str());
        if (poDriver == NULL) {
            std::cout << "Unknown GDAL driver :" << CA::ESRI_ASCIIGrid<CA::Real>::driver << std::endl;
            return EXIT_SUCCESS;
        }

        std::string extension = ".";

        extension += poDriver->GetMetadataItem(GDAL_DMD_EXTENSION);

        std::cout << "GDAL driver["<< CA::ESRI_ASCIIGrid<CA::Real>::driver<<"] extension :" << extension << std::endl;

        //char** metaDomainList = poDriver->GetMetadataDomainList();
        //for (int i = 0; metaDomainList[i] != NULL; ++i) {
        //    std::cout << metaDomainList[i]<<std::endl;
        //}

        const char * createCapable = poDriver->GetMetadataItem(GDAL_DCAP_CREATE);
        

        const char * createCopyCapable = poDriver->GetMetadataItem(GDAL_DCAP_CREATECOPY);

        if (createCapable != NULL) {
            //std::cout << "createCapable[" << CA::ESRI_ASCIIGrid<CA::Real>::driver << "] extension :" << createCapable << std::endl;
            CA::ESRI_ASCIIGrid<CA::Real>::createCapable = true;
        }
        if (createCopyCapable != NULL) {
            //std::cout << "createCopyCapable[" << CA::ESRI_ASCIIGrid<CA::Real>::driver << "] extension :" << createCopyCapable << std::endl;
            CA::ESRI_ASCIIGrid<CA::Real>::createCopyCapable = true;
        }

        if (createCapable == NULL && createCopyCapable == NULL) {
            std::cout << "GDAL driver [" << CA::ESRI_ASCIIGrid<CA::Real>::driver << "] has no create or createCopy capability" << std::endl;
            return EXIT_SUCCESS;
        }

        if (poDriver->GetMetadataItem(GDAL_DCAP_RASTER) == NULL) {
            std::cout << "GDAL driver [" << CA::ESRI_ASCIIGrid<CA::Real>::driver << "] is not a raster driver" << std::endl;
            return EXIT_SUCCESS;
        }


	}

    if ((*i)->name == "formats")
    {

        GDALAllRegister();

        int numOfDrivers = GetGDALDriverManager()->GetDriverCount();


        for (int i = 1,j=1; i < numOfDrivers ; ++i) {
            GDALDriver* poDriver;
            poDriver = GetGDALDriverManager()->GetDriver(i);

            const char* createCapable = poDriver->GetMetadataItem(GDAL_DCAP_CREATE);
            const char* createCopyCapable = poDriver->GetMetadataItem(GDAL_DCAP_CREATECOPY);
            
            std::string shortName = poDriver->GetDescription();

            if (poDriver->GetMetadataItem(GDAL_DCAP_RASTER) != NULL && !(createCapable == NULL && createCopyCapable == NULL) && !shortName.compare("MEM")==0) {



                std::cout << "[" << j << "] -> \"" << shortName << "\" , \"" << poDriver->GetMetadataItem(GDAL_DMD_LONGNAME) << "\"" << std::endl;
                
                ++j;
            }
            //else
            //    std::cout << "NO" << std::endl;

            //std::cout << "-------------------"<< std::endl;


        }

        return EXIT_SUCCESS;

    }

	if ((*i)->name == "epsg")
	{
        GDALAllRegister();
		//setup.geoCS = (*i)->value;
        //CA::ESRI_ASCIIGrid<CA::Real>::driver = (*i)->value;
        std::istringstream stringStream((*i)->value);
        stringStream >> CA::ESRI_ASCIIGrid<CA::Real>::epsgCode;

        OGRSpatialReference oSRS;
        //OSRImportFromEPSG(osr, AsciiGridGeneral<T>::epsgCode);
        char* pszSRS_WKT = NULL;

        //oSRS.SetWellKnownGeogCS("NAD27");
        OGRErr errCode = oSRS.importFromEPSG(CA::ESRI_ASCIIGrid<CA::Real>::epsgCode);
        

        if (errCode != OGRERR_NONE) {
            std::cout << "Unknown EPSG code :" << CA::ESRI_ASCIIGrid<CA::Real>::epsgCode << std::endl;
            return EXIT_SUCCESS;
        }

        //oSRS.exportToPrettyWkt(&pszSRS_WKT);
        
        std::cout << "Geo Co-ordinate system: " << oSRS.GetAttrValue("PROJCS") << std::endl;

        //CPLFree(pszSRS_WKT);

        

	}



#endif



  }// end of cycling through arguements
  //----------------------------------------------------------------------------------------------------

  //----------------------------------------------------------------------------------------------------
  // Set the data directory.
  ad.data_dir = ad.output_dir + ad.sdir;

  
  if(ad.info)
  {
	  version();
	
    std::cout<<"Input Dir            : "<<ad.working_dir<<std::endl;
    std::cout<<"Setup File           : "<<ad.setup_file<<std::endl;
    std::cout<<"Output Dir           : "<<ad.output_dir<<std::endl;
    std::cout<<"Data   Dir           : "<<ad.data_dir<<std::endl;
  }

  if(ad.info)
    std::cout<<std::endl<<"Load simulation setup"<<std::endl;

  


  // use an option to set this.....................................................................
  setup.timeGridOut = true;

  setup.timeGridName = "TemporalOut";

  CA::CellBuffReal* timeGrid = NULL;
  double thesholdDepth = 0.0;
  //..........................
  
  // Load setup information from the CA2D.csv file in the working dir.  
  std::string setup_file = ad.working_dir+ad.sdir+ad.setup_file;

  if(initSetupFromCSV(setup_file, setup)!=0)
  {
    std::cerr<<"Error reading Setup CSV file: "<<setup_file<<std::endl;
    return EXIT_FAILURE;    
  }
    
  if(ad.info)
  {
    std::cout<<"Simulation                : "<<setup.sim_name<<std::endl;
    std::cout<<"Short Name                : "<<setup.short_name<<std::endl;
    std::cout<<"Pre-proc Name             : "<<setup.preproc_name<<std::endl;
    std::cout<<"Version                   : ";
    for(size_t i = 0; i< setup.sim_version.size(); ++i)
      std::cout<<setup.sim_version[i]<<" ";
    std::cout<<std::endl;
    std::cout<<"Model Type                : "<<setup.model_type<<std::endl;
    std::cout<<"Time Start                : "<<setup.time_start<<std::endl;
    std::cout<<"Time End                  : "<<setup.time_end<<std::endl;
    std::cout<<"Max DT                    : "<<setup.time_maxdt<<std::endl;
    std::cout<<"Min DT                    : "<<setup.time_mindt<<std::endl;
    std::cout<<"Update DT                 : "<<setup.time_updatedt<<std::endl;
    std::cout<<"Alpha                     : "<<setup.time_alpha<<std::endl;
    std::cout<<"Max Iterations            : "<<setup.time_maxiters<<std::endl;
    std::cout<<"Roughness Global          : "<<setup.roughness_global<<std::endl;
    std::cout<<"Infiltration Global       : "<<setup.infrate_global<<std::endl;
    std::cout<<"Ignore WD                 : "<<setup.ignore_wd<<std::endl;
    std::cout<<"Tolerance                 : "<<setup.tolerance<<std::endl;
    std::cout<<"Slope Tolerance (%)       : "<<setup.tol_slope<<std::endl;
    std::cout<<"Boundary Ele              : "<<setup.boundary_elv<<std::endl;
    std::cout<<"Geographic Coordinate     : "<<setup.latlong<<std::endl;
    std::cout<<"Elevation ASCII           : "<<setup.elevation_ASCII<<std::endl;
    std::cout<<"Initial Water Depths      : "<<setup.startingDepthFileName << std::endl;

	std::cout<<"Height Adjustment Spatial : ";
	for (size_t i = 0; i< setup.heighAdjustmentMap_files.size(); ++i)
		std::cout << setup.heighAdjustmentMap_files[i] << " ";
	std::cout << std::endl;

    std::cout << "Rain Spatial Dense      : ";
    for (size_t i = 0; i < setup.rainSpatialDense_files.size(); ++i)
        std::cout << setup.rainSpatialDense_files[i] << " ";
    std::cout << std::endl;

    std::cout<<"Roughness Spatial Temp    : ";
    for(size_t i = 0; i< setup.strough_files.size(); ++i)
      std::cout<<setup.strough_files[i]<<" ";
    std::cout<<std::endl;
    std::cout<<"Rain Spatial Temp         : ";
    for(size_t i = 0; i< setup.strain_files.size(); ++i)
      std::cout<<setup.strain_files[i]<<" ";
    std::cout<<std::endl;
    std::cout<<"Infiltration Spatial Temp : ";
    for(size_t i = 0; i< setup.stinf_files.size(); ++i)
      std::cout<<setup.stinf_files[i]<<" ";
    std::cout<<std::endl;
    std::cout<<"Rain Event CSV            : ";
    for(size_t i = 0; i< setup.rainevent_files.size(); ++i)
      std::cout<<setup.rainevent_files[i]<<" ";
    std::cout<<std::endl;
    std::cout<<"Water Level Event CSV     : ";
    for(size_t i = 0; i< setup.wlevent_files.size(); ++i)
      std::cout<<setup.wlevent_files[i]<<" ";
    std::cout<<std::endl;
    std::cout<<"Inflow Event CSV          : ";
    for(size_t i = 0; i< setup.inflowevent_files.size(); ++i)
      std::cout<<setup.inflowevent_files[i]<<" ";
    std::cout<<std::endl;
    std::cout<<"Time Plot CSV             : ";
    for(size_t i = 0; i< setup.timeplot_files.size(); ++i)
      std::cout<<setup.timeplot_files[i]<<" ";
    std::cout<<std::endl;
    std::cout<<"Raster CSV                : ";
    for(size_t i = 0; i< setup.rastergrid_files.size(); ++i)
      std::cout<<setup.rastergrid_files[i]<<" ";
    std::cout<<std::endl;
    std::cout<<"Output Console            : "<<setup.output_console<<std::endl;
    std::cout<<"Output Period             : "<<setup.output_period<<std::endl;
    std::cout<<"Terrain Info              : deprecated"<<std::endl;
    std::cout<<"TS Plot                   : "<<setup.ts_plot<<std::endl;
    std::cout<<"Output Computation Time   : "<<setup.output_computation<<std::endl;
    std::cout<<"Check Volumes             : "<<setup.check_vols<<std::endl;
    std::cout<<"Remove Proc Data          : "<<setup.remove_data<<std::endl;
    std::cout<<"Remove Pre-Proc Data      : "<<setup.remove_prec_data<<std::endl;
    std::cout<<"Raster VEL Vector         : "<<setup.rast_vel_as_vect<<std::endl;
    std::cout<<"Raster WD Tol             : "<<setup.rast_wd_tol<<std::endl;
    std::cout<<"Raster Boundary Output    : "<<setup.rast_boundary<<std::endl;
    std::cout<<"Raster Decimal Places     : "<<setup.rast_places<<std::endl;
    std::cout<<"Update Peak Every DT      : "<<setup.update_peak_dt<<std::endl;
    std::cout<<"Expand Domain             : "<<setup.expand_domain<<std::endl;
    std::cout<<"Ignore Upstream           : "<<setup.ignore_upstream<<std::endl;
    std::cout<<"Upstream Reduction        : "<<setup.upstream_reduction<<std::endl;
    std::cout<<"Domain Decomposition      : "<<setup.decompose_domain<<std::endl;
    std::cout<<"Decomposition Efficiency  : "<<setup.decompose_eff<<std::endl;
  }

  setup.terrain_info = false;

  if(ad.info)
    std::cout<<std::endl<<"Load elevation data "<<std::endl;

  //! Load the eventual elevation.
  std::string ele_file = ad.working_dir+ad.sdir+setup.elevation_ASCII;
  CA::AsciiGrid<CA::Real> eg;

  // ATTENTION. Load only the header here .. not the actual data
  //CA::readAsciiGridHeader(eg,ele_file);
  eg.readAsciiGridHeader(ele_file);


  if(ad.info)
  {
    std::cout<<"Elevation ASCII    : "<<setup.elevation_ASCII<<std::endl;
    std::cout<<"ncols              : "<<eg.ncols<<std::endl;
    std::cout<<"nrows              : "<<eg.nrows<<std::endl;
    std::cout<<"xllcorner          : "<<eg.xllcorner<<std::endl;
    std::cout<<"yllcorner          : "<<eg.yllcorner<<std::endl;
    std::cout<<"cellsize           : "<<eg.cellsize<<std::endl;
    std::cout<<"nodata             : "<<eg.nodata<<std::endl;
  }


  if (!setup.startingDepthFileName.empty())
  {

      setup.hotStart = true;

      // create extent as defined, in memory, will be deleted after use in preproc
      setup.startingDepths = new CA::ESRI_ASCIIGrid<CA::Real>();

      std::string ini_file = ad.working_dir + ad.sdir + setup.startingDepthFileName;

      setup.startingDepths->readAsciiGrid(ini_file);




  }



  // Fixed list of spatial temporal variables. This is expected to be
  // composd of three elements in thi ordeer:
  // 1) Roughness
  // 2) Rain
  // 3) Infiltration.
  // 4) height adjustment
  // These elements must be there also if the simulation does not
  // contain any spatial/temporal element.
  std::vector<SpatialTemporal> sts(4);
  SpatialTemporalDense _std;

  

  if(ad.info)
    std::cout<<std::endl<<"Load spatial temporal roughness configuration "<<std::endl;
  
  if(initSpatialTemporalFromCSV(ad.working_dir+ad.sdir,setup.strough_files, sts[0])!=0)
  {
    std::cerr<<"Error reading Spatial Temporal Rougness from files: ";
    for(size_t i = 0; i<setup.strough_files.size(); ++i)
      std::cerr<<setup.strough_files[i]<<" ";
    std::cerr<<std::endl;

    return EXIT_FAILURE;    
  }
  
  // Make sure that there is at least one sequence.
  if(sts[0].seqs.empty())  sts[0].seqs.resize(1);

  // The index zero of the spatial roughness is reserved for the
  // global roughness value.   
  sts[0].seqs[0].values.resize(1,setup.roughness_global);
  sts[0].seqs[0].times .resize(1,std::numeric_limits<CA::Real>::max());

  if(setup.strough_files.size()==2)
  {
	  
	  
 
    if(ad.info)
    {
      std::cout<<"Name               : "<<sts[0].name<<std::endl;
      std::cout<<"Roughness ASCII    : "<<setup.strough_files[0]<<std::endl;
      std::cout<<"ncols              : "<<sts[0].map.ncols<<std::endl;
      std::cout<<"nrows              : "<<sts[0].map.nrows<<std::endl;
      std::cout<<"xllcorner          : "<<sts[0].map.xllcorner<<std::endl;
      std::cout<<"yllcorner          : "<<sts[0].map.yllcorner<<std::endl;
      std::cout<<"cellsize           : "<<sts[0].map.cellsize<<std::endl;
      std::cout<<"nodata             : "<<sts[0].map.nodata<<std::endl;
      std::cout<<"Roughness Values   : "<<setup.strough_files[1]<<std::endl;
      std::cout<<"Num Indices        : "<<sts[0].seqs.size()-1<<std::endl; // Remove the reservd index 0.
      std::cout<<"Longest Seq.       : "<<sts[0].l<<std::endl;
      std::cout<<"Last Time End      : "<<sts[0].end<<std::endl;
    }
  }

  if(ad.info)
    std::cout<<std::endl<<"Load spatial temporal rain configuration "<<std::endl;
  
  if(initSpatialTemporalFromCSV(ad.working_dir+ad.sdir,setup.strain_files, sts[1])!=0)
  {
    std::cerr<<"Error reading Spatial Temporal Rain from files: ";
    for(size_t i = 0; i<setup.strain_files.size(); ++i)
      std::cerr<<setup.strain_files[i]<<" ";
    std::cerr<<std::endl;

    return EXIT_FAILURE;    
  }
  
  // Make sure that there is at least one sequence.
  if(sts[1].seqs.empty())  sts[1].seqs.resize(1);

  // The index zero of the spatial temporal rain is reserved for the
  // no rain value
  sts[1].seqs[0].values.resize(1,0);
  sts[1].seqs[0].times .resize(1,std::numeric_limits<CA::Real>::max());

  if(setup.strain_files.size()==2)
  {
	  
	  
  
    if(ad.info)
    {
      std::cout<<"Name               : "<<sts[1].name<<std::endl;
      std::cout<<"Rain ASCII         : "<<setup.strain_files[0]<<std::endl;
      std::cout<<"ncols              : "<<sts[1].map.ncols<<std::endl;
      std::cout<<"nrows              : "<<sts[1].map.nrows<<std::endl;
      std::cout<<"xllcorner          : "<<sts[1].map.xllcorner<<std::endl;
      std::cout<<"yllcorner          : "<<sts[1].map.yllcorner<<std::endl;
      std::cout<<"cellsize           : "<<sts[1].map.cellsize<<std::endl;
      std::cout<<"nodata             : "<<sts[1].map.nodata<<std::endl;
      std::cout<<"Rain Values        : "<<setup.strain_files[1]<<std::endl;
      std::cout<<"Num Indices        : "<<sts[1].seqs.size()-1<<std::endl; // Remove the reservd index 0.
      std::cout<<"Longest Seq.       : "<<sts[1].l<<std::endl;
      std::cout<<"Last Time End      : "<<sts[1].end<<std::endl;
    }
  }

  //----------------------

  if(ad.info)
    std::cout<<std::endl<<"Load spatial temporal infiltration configuration "<<std::endl;
  
  if(initSpatialTemporalFromCSV(ad.working_dir+ad.sdir,setup.stinf_files, sts[2])!=0)
  {
    std::cerr<<"Error reading Spatial Temporal Infiltration from files: ";
    for(size_t i = 0; i<setup.stinf_files.size(); ++i)
      std::cerr<<setup.stinf_files[i]<<" ";
    std::cerr<<std::endl;

    return EXIT_FAILURE;    
  }
  
  // Make sure that there is at least one sequence.
  if(sts[2].seqs.empty())  sts[2].seqs.resize(1);

  // The index zero of the spatial temporal infiltration is reserved for the
  // global infiltration value
  sts[2].seqs[0].values.resize(1,setup.infrate_global);
  sts[2].seqs[0].times .resize(1,std::numeric_limits<CA::Real>::max());

  if(setup.stinf_files.size()==2)
  {
	  
  
    if(ad.info)
    {
      std::cout<<"Name               : "<<sts[2].name<<std::endl;
      std::cout<<"Infiltration ASCII : "<<setup.stinf_files[0]<<std::endl;
      std::cout<<"ncols              : "<<sts[2].map.ncols<<std::endl;
      std::cout<<"nrows              : "<<sts[2].map.nrows<<std::endl;
      std::cout<<"xllcorner          : "<<sts[2].map.xllcorner<<std::endl;
      std::cout<<"yllcorner          : "<<sts[2].map.yllcorner<<std::endl;
      std::cout<<"cellsize           : "<<sts[2].map.cellsize<<std::endl;
      std::cout<<"nodata             : "<<sts[2].map.nodata<<std::endl;
      std::cout<<"Infiltration Values: "<<setup.stinf_files[1]<<std::endl;
      std::cout<<"Num Indices        : "<<sts[2].seqs.size()-1<<std::endl; // Remove the reservd index 0.
      std::cout<<"Longest Seq.       : "<<sts[2].l<<std::endl;
      std::cout<<"Last Time End      : "<<sts[2].end<<std::endl;
    }
  }

  //----------------------

  if (ad.info)
	  std::cout << std::endl << "Load spatial height adjustment configuration " << std::endl;

  if (initSpatialTemporalFromCSV(ad.working_dir + ad.sdir, setup.heighAdjustmentMap_files, sts[3]) != 0)
  {
	  std::cerr << "Error reading Spatial Height adjustment from files: ";
	  for (size_t i = 0; i<setup.heighAdjustmentMap_files.size(); ++i)
		  std::cerr << setup.heighAdjustmentMap_files[i] << " ";
	  std::cerr << std::endl;

	  return EXIT_FAILURE;
  }

  // Make sure that there is at least one sequence.
  if (sts[3].seqs.empty())  sts[3].seqs.resize(1);

  // The index zero of the spatial temporal infiltration is reserved for the
  // global infiltration value
  sts[3].seqs[0].values.resize(1, 0);
  sts[3].seqs[0].times.resize(1, std::numeric_limits<CA::Real>::max());

  if (setup.heighAdjustmentMap_files.size() == 2)
  {
	  
	  
	  if (ad.info)
	  {
		  std::cout << "Name               : " << sts[3].name << std::endl;
		  std::cout << "Infiltration ASCII : " << setup.heighAdjustmentMap_files[0] << std::endl;
		  std::cout << "ncols              : " << sts[3].map.ncols << std::endl;
		  std::cout << "nrows              : " << sts[3].map.nrows << std::endl;
		  std::cout << "xllcorner          : " << sts[3].map.xllcorner << std::endl;
		  std::cout << "yllcorner          : " << sts[3].map.yllcorner << std::endl;
		  std::cout << "cellsize           : " << sts[3].map.cellsize << std::endl;
		  std::cout << "nodata             : " << sts[3].map.nodata << std::endl;
		  std::cout << "Infiltration Values: " << setup.heighAdjustmentMap_files[1] << std::endl;
		  std::cout << "Num Indices        : " << sts[3].seqs.size() - 1 << std::endl; // Remove the reservd index 0.
		  std::cout << "Longest Seq.       : " << sts[3].l << std::endl;
		  std::cout << "Last Time End      : " << sts[3].end << std::endl;
	  }

	  setup.addHeights = true;


	  // load the data into the necessary structures
	  //sts[3].map.readAsciiGrid(ad.working_dir + ad.sdir + setup.heighAdjustmentMap_files[1]);
	  setup.heightAdjustmentMap = new CA::AsciiGrid<CA::State>();

	  //std::cout << "Loading height adjustment map: " << ad.working_dir + ad.sdir + setup.heighAdjustmentMap_files[0] << std::endl;

	  setup.heightAdjustmentMap->readAsciiGrid(ad.working_dir + ad.sdir + setup.heighAdjustmentMap_files[0]);

	  //std::cout << "Loaded height adjustment map: " << ad.working_dir + ad.sdir + setup.heighAdjustmentMap_files[0] << std::endl;

	  //setup.addedHeight
	  for (int i = 0; i < sts[3].seqs.size(); ++i) {
		  
		  setup.addedHeight.push_back(sts[3].seqs[i].values[0]);
	  }

	  


  }// end of if heightmap files provided

  //------------------------------------------------------------------------------------

  if (ad.info)
      std::cout << std::endl << "Load Rain spatial dense configuration " << std::endl;

  if (initSpatialTemporalDenseFromCSV(ad.working_dir + ad.sdir, setup.rainSpatialDense_files, _std, setup.rainNetCDF) != 0)
  {
      std::cerr << "Error reading rain Spatial dense from files: ";
      for (size_t i = 0; i < setup.rainSpatialDense_files.size(); ++i)
          std::cerr << setup.rainSpatialDense_files[i] << " ";
      std::cerr << std::endl;

      return EXIT_FAILURE;
  }

  if (setup.rainSpatialDense_files.size() >= 2)
  {
      if (ad.info)
      {
          std::cout << "Name               : " << _std.name << std::endl;
          std::cout << "Rain rates ASCII's : " ;
          for (size_t i = 1; i < _std.map_filenames.size(); ++i)
              std::cout << _std.map_filenames[i] << " ";
          std::cout << std::endl;
          std::cout << "Times              : ";
          for (size_t i = 0; i < _std.times.size(); ++i)
              std::cout << _std.times[i] << " ";
          std::cout << std::endl;

          //std::cout << "cellsize           : " << sts[3].map.cellsize << std::endl;
          
          //std::cout << "rain Times Values  : " << _std.map_filenames[0] << std::endl;
          
          //std::cout << "Longest Seq.       : " << _std.l << std::endl;
          std::cout << "Last Time End      : " << _std.end << std::endl;
      }


  }// end of if rainSpatialDense_files files provided

  //------------------------------------------------------------------------------------

  if(ad.info)
    std::cout<<std::endl<<"Load rain event inputs configuration "<<std::endl;

  // Load any eventual rain event.
  std::vector<RainEvent> res; 
  for(size_t i = 0; i< setup.rainevent_files.size(); ++i)
  {
    std::string file = ad.working_dir+ad.sdir+setup.rainevent_files[i]; 

    RainEvent re;
    if(initRainEventFromCSV(file, re)!=0)
    {
      std::cerr<<"Error reading Rain Event CSV file: "<<file<<std::endl;
      return EXIT_FAILURE;    
    }

    if(ad.info)
    {
      std::cout<<"Rain Event         : "<<re.name<<std::endl;
      std::cout<<"Rain Intensity     : ";
      for(size_t i = 0; i< re.rains.size(); ++i)
	std::cout<<re.rains[i]<<" ";
      std::cout<<std::endl;
      std::cout<<"Time Stop          : ";
      for(size_t i = 0; i< re.times.size(); ++i)
	std::cout<<re.times[i]<<" ";
      std::cout<<std::endl;
      std::cout<<"Area               : ";
      for(size_t i = 0; i< re.area.size(); ++i)
	std::cout<<re.area[i]<<" ";
      std::cout<<std::endl;
    }

    res.push_back(re);    
  } 

  if(ad.info)
    std::cout<<std::endl<<"Load water level event inputs configuration "<<std::endl;

  // Load any eventual water level events
  std::vector<WLEvent> wles; 
  for(size_t i = 0; i< setup.wlevent_files.size(); ++i)
  {
    std::string file = ad.working_dir+ad.sdir+setup.wlevent_files[i]; 

    WLEvent wle;
    if(initWLEventFromCSV(file, wle)!=0)
    {
      std::cerr<<"Error reading Water Level Event CSV file: "<<file<<std::endl;
      return EXIT_FAILURE;    
    }

    if(ad.info)
    {
      std::cout<<"Water Level Event  : "<<wle.name<<std::endl;
      std::cout<<"Water Level        : ";
      for(size_t i = 0; i< wle.wls.size(); ++i)
	std::cout<<wle.wls[i]<<" ";
      std::cout<<std::endl;
      std::cout<<"Time               : ";
      for(size_t i = 0; i< wle.times.size(); ++i)
	std::cout<<wle.times[i]<<" ";
      std::cout<<std::endl;
      std::cout<<"Area               : ";
      for(size_t i = 0; i< wle.area.size(); ++i)
	std::cout<<wle.area[i]<<" ";
      std::cout<<std::endl;
      if(wle.u!=0)
      {
	std::cout<<"U                  : "<<wle.u<<std::endl;
	std::cout<<"N                  : "<<wle.n<<std::endl;
      }
    }

    wles.push_back(wle);    
  } 

  if(ad.info)
    std::cout<<std::endl<<"Load inflow event inputs configuration "<<std::endl;

  // Load any eventual inflow events
  std::vector<IEvent> ies; 
  for(size_t i = 0; i< setup.inflowevent_files.size(); ++i)
  {
    std::string file = ad.working_dir+ad.sdir+setup.inflowevent_files[i]; 

    IEvent ie;

    if(initIEventFromCSV(file, ie)!=0)
    {
      std::cerr<<"Error reading Inflow Event CSV file: "<<file<<std::endl;
      return EXIT_FAILURE;    
    }

    if(ad.info)
    {
      std::cout<<"Inflow Event       : "<<ie.name<<std::endl;
      std::cout<<"Inflow (cumecs)    : ";
      for(size_t i = 0; i< ie.ins.size(); ++i)
	std::cout<<ie.ins[i]<<" ";
      std::cout<<std::endl;
      std::cout<<"Time               : ";
      for(size_t i = 0; i< ie.times.size(); ++i)
	std::cout<<ie.times[i]<<" ";
      std::cout<<std::endl;
      std::cout<<"Area               : ";
      for(size_t i = 0; i< ie.area.size(); ++i)
	std::cout<<ie.area[i]<<" ";
      std::cout<<std::endl;
    }

    ies.push_back(ie);    
  } 

  if(ad.info)
    std::cout<<std::endl<<"Load time plot outputs configuration "<<std::endl;

  // Load any eventual time plots.
  std::vector<TimePlot> tps; 
  for(size_t i = 0; i< setup.timeplot_files.size(); ++i)
  {
    std::string file = ad.working_dir+ad.sdir+setup.timeplot_files[i]; 

    TimePlot tp;
    if(initTimePlotFromCSV(file, tp)!=0)
    {
      std::cerr<<"Error reading Time Plot CSV file: "<<file<<std::endl;
      return EXIT_FAILURE;    
    }

    if(ad.info)
    {
      std::cout<<"Time Plot Name     : "<<tp.name<<std::endl;
      std::cout<<"Physical Variable  : "<<tp.pv<<std::endl;
      std::cout<<"Points X Name      : ";
      for(size_t i = 0; i< tp.pnames.size(); ++i)
	std::cout<<tp.pnames[i]<<" ";
      std::cout<<std::endl;
      std::cout<<"Points X Coo       : ";
      for(size_t i = 0; i< tp.xcoos.size(); ++i)
	std::cout<<tp.xcoos[i]<<" ";
      std::cout<<std::endl;
      std::cout<<"Points Y Coo       : ";
      for(size_t i = 0; i< tp.ycoos.size(); ++i)
	std::cout<<tp.ycoos[i]<<" ";
      std::cout<<std::endl;  
      std::cout<<"Period             : "<<tp.period<<std::endl;
    }

    tps.push_back(tp);    
  } 

  if(ad.info)
    std::cout<<std::endl<<"Load raster grid outputs configuration "<<std::endl;

  // Load any eventual raster grids.
  std::vector<RasterGrid> rgs; 
  for(size_t i = 0; i< setup.rastergrid_files.size(); ++i)
  {
    std::string file = ad.working_dir+ad.sdir+setup.rastergrid_files[i]; 

    RasterGrid rg;
    if(initRasterGridFromCSV(file, rg)!=0)
    {
      std::cerr<<"Error reading Raster Grid CSV file: "<<file<<std::endl;
      return EXIT_FAILURE;    
    }

    if(ad.info)
    {
      std::cout<<"Raster Grid Name   : "<<rg.name<<std::endl;
      std::cout<<"Physical Variable  : "<<rg.pv<<std::endl;
      std::cout<<"Peak               : "<<rg.peak<<std::endl;
      std::cout<<"Final              : "<<rg.final<<std::endl;
      std::cout<<"Period             : "<<rg.period<<std::endl;
    }

    rgs.push_back(rg);    
  } 

  CA::State returnCode = -1;

  // Variable that indicate that something was done.
  bool work_done=false;
  try
  {
    //! Now perform the pre-processing
    if(ad.pre_proc && !ad.no_pre_proc)
    {
      if(ad.info)
	std::cout<<std::endl<<"Starting pre-processing data "<<std::endl;
   
      if(preProc(ad,setup,ele_file,sts,_std)!=0)
      {
	std::cerr<<"Error while performing pre-processing"<<std::endl;
	std::cerr<<"Possible cause is the output directory argument"<<std::endl;
	return EXIT_FAILURE;    
      }
    
      work_done = true;

      if(ad.info)
      {
	std::cout<<std::endl;  
	std::cout<<"Ending pre-processing data "<<std::endl;
      }
    }

    //! Now display the terrain info
    if(ad.terrain_info)
    {
      if(ad.info)
	std::cout<<std::endl<<"Display terrain info "<<std::endl;
  
      if(terrainInfo(ad,setup,eg)!=0)
      {
	std::cerr<<"Error while displaying terrain info"<<std::endl;
	return EXIT_FAILURE;    
      }
      
      work_done = true;
          
      // Make sure the model and post-proc are not performed.
      ad.model.clear();
      ad.post_proc = false;
    }

    //! Now perform the CA2D flood modelling if requested.   
    if(!ad.model.empty())
    {
      if(ad.info)
	std::cout<<std::endl<<"Starting CADDIES2D flood modelling using "<<setup.model_type<<" model"<<std::endl;

	  bool abortIt = false;
	  double localProcessingTime = 0.0;

      //eg.writeAsciiGrid(ad.output_dir + ad.sdir + "elevationTest.tif");

	  if (CADDIES2D(ad, setup, eg, sts, res, wles, ies, tps, rgs,
		  &(timeGrid), thesholdDepth,
		  localProcessingTime, returnCode,
		  abortIt, _std,false) != 0)
      {
	std::cerr<<"Error while performing CADDIES2D flood modelling"<<std::endl;
	return EXIT_FAILURE;    
      }

      work_done = true;
    
      if(ad.info)
      {
	std::cout<<std::endl;  
	std::cout<<"Ending CADDIES2D flood modelling "<<std::endl;
      }
    }

    //! Now perform the post-processing
    if(ad.post_proc)
    {
      if(ad.info)
	std::cout<<std::endl<<"Starting post-processing data "<<std::endl;

	  if (postProc(ad, setup, eg, tps, rgs, timeGrid) != 0)
      {
	std::cerr<<"Error while performing post-processing"<<std::endl;
	return EXIT_FAILURE;    
      }

      work_done = true;

      if(ad.info)
      {
	std::cout<<std::endl;  
	std::cout<<"Ending post-processing data "<<std::endl;
      }
    }
    if(!work_done)
    {
      std::cout<<"\nATTENTION! No work was performed; missing specific options?"<<std::endl;
      std::cout<<"Use the /help options to show help"<<std::endl;
      std::cout<<"Press 'Return' to continue"<<std::endl;
      return EXIT_FAILURE;    
    }

    //eg.writeAsciiGrid(ad.working_dir + ad.sdir + "elevationTest.tif");

  }
  catch(std::exception& e)
  {
    std::cerr<<e.what()<<std::endl;
    std::cout<<"Simulation stopped"<<std::endl;
    return EXIT_FAILURE;    
  }
  // Close the caAPI.
  CA::finalise2D();

#if defined _WIN32 || defined __CYGWIN__   
  // Ok WINDOWS can go back to sleep if he want.
  SetThreadExecutionState(ES_CONTINUOUS);
#endif


  return EXIT_SUCCESS;
}






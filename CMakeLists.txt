########### ADD CAFLOOD application ###############

# Application that performs a 2D flood modelling using cellular
# automata model(s) and it is implemented using the CAAPI.

### Set up the application configuration
cmake_policy(SET CMP0054 NEW)


# The main source of the app.
set(CAAPI_APP_MAIN "mains/main.cpp")

# The source of the app.
set(CAAPI_APP_SOURCES "")
file(GLOB CAAPI_APP_SOURCES "*.cpp")

# The headers of the app.
set(CAAPI_APP_HEADERS "")
file(GLOB CAAPI_APP_HEADERS "*.hpp")

# The data 
set(CAAPI_APP_DATA "NEWS;LICENSE;AUTHORS;Version.txt;CMakeLists.txt")

# Find the CA functions
set(CAAPI_APP_CAFUNCTIONS)
set(CAAPI_IMPL_HEADERS)



####################### Net CDF rainFall ###########################

SET(CAAPI_NETCDF_OPTION OFF CACHE BOOL "Enable GDAL or not?")

if(CAAPI_NETCDF_OPTION)

        
    find_path (NETCDF_INCLUDES_C netcdf.h HINTS NETCDF_DIR NETCDF_DIR)
    message(" inc =  ${NETCDF_INCLUDES_C}  ")
    find_library (NETCDF_LIBRARIES_C  NAMES netcdf)
    message(" lib =  ${NETCDF_LIBRARIES_C}  ")
    
#find_path(NETCDF_INCLUDES_CXX NAMES netcdf
#        HINTS ${NETCDF_CXX_ROOT} PATH_SUFFIXES include)
#    find_library(NETCDF_LIBRARIES_CXX NAMES netcdf_c++4 netcdf-cxx4
#        HINTS ${NETCDF_CXX_ROOT} PATH_SUFFIXES lib)
    
    
    
    # add definition to activate code segments
    if (WIN32)
        add_definitions( /DCAAPI_NETCDF_OPTION=1 )
    else()
        add_definitions( -DCAAPI_NETCDF_OPTION=1 )
    endif()

    # Add the directory with the GDAL includes
    #list(APPEND CAAPI_HEADERS ${GDAL_INCLUDE_DIR})
    # Add the needed directories into the compilation
    #include_directories(${netCDF_INCLUDE_DIR})
    #MESSAGE("[${netCDF_INCLUDE_DIR}]: netCDF_INCLUDE_DIR")
    include_directories(${NETCDF_INCLUDES_C})
    list(APPEND CAAPI_HEADERS ${NETCDF_INCLUDES_C})
    
    #ADD_LIBRARY(netCDF_LIBRARIES STATIC IMPORTED)

    # Add the GDAL libraries 
    #MESSAGE("[${CAAPI_IMPL_LIBRARIES}]: CAAPI_IMPL_LIBRARIES")
    #list(APPEND CAAPI_IMPL_LIBRARIES ${netCDF_LIBRARIES})
    #MESSAGE("[${netCDF_LIBRARIES}]: netCDF_LIBRARIES")
    #MESSAGE("[${CAAPI_IMPL_LIBRARIES}]: CAAPI_IMPL_LIBRARIES")
    
    list(APPEND CAAPI_IMPL_LIBRARIES ${NETCDF_LIBRARIES_C})



ENDIF(CAAPI_NETCDF_OPTION)



########################## GDAL  ##########################

if(CAAPI_NETCDF_OPTION)
    SET(CAAPI_GDAL_OPTION ON BOOL "Enable GDAL as must be used with NetCDF")
else(CAAPI_NETCDF_OPTION)
    SET(CAAPI_GDAL_OPTION OFF CACHE BOOL "Enable GDAL or not?")
ENDIF(CAAPI_NETCDF_OPTION)


if(CAAPI_GDAL_OPTION)

    # Find the required GDAL libraries and headers
    find_package( GDAL REQUIRED )
    
    # add definition to activate code segments
    if (WIN32)
        add_definitions( /DCAAPI_GDAL_OPTION=1 )
    else()
        add_definitions( -DCAAPI_GDAL_OPTION=1 )
    endif()

    # Add the directory with the GDAL includes
    list(APPEND CAAPI_HEADERS ${GDAL_INCLUDE_DIR})
    # Add the needed directories into the compilation
    include_directories(${GDAL_INCLUDE_DIR})

    # Add the GDAL libraries 
    #MESSAGE("[${CAAPI_IMPL_LIBRARIES}]: CAAPI_IMPL_LIBRARIES")
    list(APPEND CAAPI_IMPL_LIBRARIES ${GDAL_LIBRARY})
    #MESSAGE("[${GDAL_LIBRARY}]: GDAL_LIBRARY")
    #MESSAGE("[${CAAPI_IMPL_LIBRARIES}]: CAAPI_IMPL_LIBRARIES")
    


ENDIF(CAAPI_GDAL_OPTION)

###########################################################

# CAFUNCTIONS
CAAPI_add_source_dir(CAfunctions "${CMAKE_CURRENT_SOURCE_DIR}/*.ca" CAAPI_APP_CAFUNCTIONS)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# Execute the module that manage ca for the specific implementation
if(EXISTS "${CAAPI_SPECIFIC_IMPL_DIR}/extra/CAAPI_ca.cmake")

  include("${CAAPI_SPECIFIC_IMPL_DIR}/extra/CAAPI_ca.cmake")

  CAAPI_ca("${CAAPI_SPECIFIC_IMPL_DIR}" CAAPI_IMPL_HEADERS CAAPI_APP_CAFUNCTIONS)
endif()

# check for 32/64bit
if( CMAKE_SIZEOF_VOID_P EQUAL 8 )
    #MESSAGE( "64 bits compiler detected" )
    SET( EX_PLATFORM 64 )
    SET( EX_PLATFORM_NAME "x64" )
	SET( CAAPI_APP_BITNAME "_64" )
else( CMAKE_SIZEOF_VOID_P EQUAL 8 ) 
    #MESSAGE( "32 bits compiler detected" )
    SET( EX_PLATFORM 32 )
    SET( EX_PLATFORM_NAME "x86" )
	SET( CAAPI_APP_BITNAME "" )
endif( CMAKE_SIZEOF_VOID_P EQUAL 8 )

# decide device app name extension
if( CAAPI_DEVICE STREQUAL  "GPU" )
    SET( CAAPI_APP_DEVICENAME "_GPU" )
else( CAAPI_DEVICE STREQUAL  "GPU" ) 
    SET( CAAPI_APP_DEVICENAME "" )
endif( CAAPI_DEVICE STREQUAL  "GPU" )

# The exe name of the app.
set(CAAPI_APP_EXE_SHORT  "cafloodpro")
set(CAAPI_APP_DLL_SHORT "${CAAPI_APP_EXE_SHORT}_API")

set(CAAPI_APP_EXE "${CAAPI_APP_EXE_SHORT}${CAAPI_APP_DEVICENAME}${CAAPI_APP_BITNAME}")
set(CAAPI_APP_DLL "${CAAPI_APP_DLL_SHORT}${CAAPI_APP_DEVICENAME}${CAAPI_APP_BITNAME}")

#MESSAGE("[${CAAPI_APP_EXE}]: CAAPI_APP_EXE")
#MESSAGE("[${CAAPI_APP_DLL}]: CAAPI_APP_DLL")


# Set up the custom target name
set(CAAPI_APP_CUSTOM_TARGET "${CAAPI_APP_EXE}_do_impl_headers")

# Create the target that make the eventual implementation headers
# necessary.
add_custom_target(${CAAPI_APP_CUSTOM_TARGET} DEPENDS 
  ${CAAPI_IMPL_HEADERS} ${CAAPI_HEADERS} ${CAAPI_APP_CAFUNCTIONS})

# Add the executable to the project using the specified source files.
add_executable("${CAAPI_APP_EXE}" "${CAAPI_APP_MAIN}"
  ${CAAPI_APP_SOURCES} ${CAAPI_SOURCES} 
  ${CAAPI_APP_CAFUNCTIONS} ${CAAPI_HEADERS} ${CAAPI_APP_HEADERS})
add_dependencies("${CAAPI_APP_EXE}" ${CAAPI_APP_CUSTOM_TARGET})
target_link_libraries("${CAAPI_APP_EXE}" ${CAAPI_IMPL_LIBRARIES} )
if(WIN32)
  target_link_libraries("${CAAPI_APP_EXE}" wsock32 ws2_32)
endif()
#MESSAGE("[${CAAPI_IMPL_LIBRARIES}]: CAAPI_IMPL_LIBRARIES")

### MPI ###

if(MPI_COMPILE_FLAGS)
	set_target_properties(${CAAPI_APP_EXE} PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
	set_target_properties(${CAAPI_APP_EXE} PROPERTIES
	LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()


########################## VERSION  ##########################

file(STRINGS "Version.txt" CAAPI_APP_VERSION)

# Show the version number
message(STATUS "[${CAAPI_APP_NAME}]: Version ${CAAPI_APP_VERSION}")

if (WIN32)
  add_definitions( /DCAAPI_APP_VERSION=${CAAPI_APP_VERSION} )
  add_definitions( /DCAAPI_APP_NAME=${CAAPI_APP_NAME} )
else()
  add_definitions( -DCAAPI_APP_VERSION=${CAAPI_APP_VERSION} )
  add_definitions( -DCAAPI_APP_NAME=${CAAPI_APP_NAME} )
endif()


########################## INFO  ##########################

# Number of files
list(LENGTH  CAAPI_APP_CAFUNCTIONS   __num_ca)
list(LENGTH  CAAPI_APP_DATA    __num_data)
list(LENGTH  CAAPI_APP_HEADERS __num_headers)
list(LENGTH  CAAPI_APP_SOURCES __num_sources)
MATH(EXPR __num_sources "${__num_sources}+1")
message(STATUS "[${CAAPI_APP_NAME}]: File #: Source (${__num_sources}) Headers (${__num_headers}) CA (${__num_ca}) Data (${__num_data})")



add_definitions(-DCAFUNCSDLL_EXPORTS)

include(GenerateExportHeader)






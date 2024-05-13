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

#ifndef _ARGSDATA_HPP_
#define _ARGSDATA_HPP_

//! \file ArgsData.hpp
//!  Contains the structure with the arguments data.
//! \author Michele Guidolin, University of Exeter, 
//! contact: a.s.chen [at] exeter.ac.uk
//! \date 2024-05

#include"BaseTypes.hpp"
#include"ca2D.hpp"
#include"Utilities.hpp"
#include<string>
#include<vector>

#define READ_TOKEN(_test,_var,_tok,_sect)				\
   {									\
   if(!CA::fromString(_var,_tok))					\
   {									\
     std::cerr<<"Error reading '"<<CA::trimToken(_sect)<<"' element"<<std::endl; \
     return 1;								\
   }									\
   _test = true;							\
   }

//! Remove file extension
inline std::string removeExtension(const std::string& filename) 
{
  size_t lastdot = filename.find_last_of(".");
  if (lastdot == std::string::npos) return filename;
  return filename.substr(0, lastdot); 
}


//! Identifies the model to use in the CADDIES2D flood modelling
struct MODEL
{
  enum Type
  {
    UNKNOWN = 0,
    Diffusive,
    Inertial,
    WCA2Dv1,
    WCA2Dv2
  };
};


//! Transform an input string into a physical variable enum.
std::istream& operator>>(std::istream& in,  MODEL::Type& m);
//! Transform a physical variable enum into an output string.
std::ostream& operator<<(std::ostream& out, MODEL::Type& m);


//! Public data structure that contains all the data set by the
//! arguments.
struct ArgsData
{

  // Create the arguments list and define the prefix which identifies
  // the arguments that are optional. This prefix differs between
  // windows and unix.
  CA::Arguments args;

  //! The working directory where the data and the configuration
  //! files are located.
  std::string working_dir;

  //! The setup file which contain the initial configuration of the ca
  //! algorithm.
  std::string setup_file;

  //! The output directory where the output data files are saved.
  std::string output_dir;

  //! Contain the string used as separator of directory.
  std::string sdir;

  //! The direcotry where temporary GRID/BUFFER data is saved/loaded.
  std::string data_dir;

  //! If true, print any information into console.
  bool info;		

  //! If true, perform preprocessing.
  bool pre_proc;		

  //! If true, do not perform the preprocessing. This take precedence
  //! over pre-proc argument.
  bool no_pre_proc;

  //! If true, perform postprocessing.
  bool post_proc;

  //! Which Model (is a string).
  std::string model;

  //! If true, display the terrrain info and exit.
  bool terrain_info;

  // Constructor
  ArgsData():
#if defined _WIN32 || defined __CYGWIN__   
    args("/"),
#else
    args("-"),
#endif
    working_dir("."),
#if defined _WIN32 || defined __CYGWIN__   
    sdir("\\"),
#else
    sdir("/"),
#endif
    data_dir(),
    info(false),
    pre_proc(false),
    no_pre_proc(false),
    post_proc(false),
    model(),
    terrain_info(false)
  {}
  
  ~ArgsData(){}
  
};


//! Identifies different physical variables, i.e. the buffer of data.
struct PV
{
  enum Type
  {
    UNKNOWN = 0,
    WD,				//!< The water depth. 
    WL,				//!< The water level (depth + elv).
    VEL,			//!< The velocity saved as Vertical, Horizonatl.
  };
};

//! Transform an input string into a physical variable enum.
std::istream& operator>>(std::istream& in,  PV::Type& pv);
//! Transform a physical variable enum into an output string.
std::ostream& operator<<(std::ostream& out, PV::Type& pv);



#endif

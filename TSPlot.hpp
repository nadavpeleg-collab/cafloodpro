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

#ifndef _TSPLOT_HPP_
#define _TSPLOT_HPP_


//! \file TSPlot.hpp
//! Contains the manager that plot the time step(s)
//! \author Michele Guidolin, University of Exeter, 
//! contact: a.s.chen [at] exeter.ac.uk
//! \date 2024-05


#include"ca2D.hpp"
#include"BaseTypes.hpp"
#include"ArgsData.hpp"
#include"PointList.hpp"
#include<string>
#include<vector>
#include<fstream>


//! Class that manage all Time Steps Plots outputs
class TSPlot
{
public:

  //! Construct a TimeSteps Plot manager
  //! \param name This is the name of the file.
  TSPlot(std::string name,bool plot);

  //! Destroy a Time Steps Plot Manager.
  ~TSPlot();
  

  //! Output the time plots.
  //! \params t       The simulation time.
  //! \param  dt      The last dt
  void output(CA::Real t, CA::Real dt);

private:

  cpp11::shared_ptr<std::ofstream> _file;   //!< The file where to output the time plot data.

};


#endif

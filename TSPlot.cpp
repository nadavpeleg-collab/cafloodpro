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

//! \file TSPlot.cpp
//! contact: a.s.chen [at] exeter.ac.uk
//! \date 2024-05


#include"ca2D.hpp"
#include"ArgsData.hpp"
#include"TSPlot.hpp"
#include"Utilities.hpp"
#include<iostream>
#include<fstream>


TSPlot::TSPlot(std::string name, bool plot):
_file()
{
  if(plot)
  {
    // Create file
    _file.reset( new std::ofstream(name.c_str()) );
    
    if(_file->good())
    {
      
      // Set  manipulators 
      _file->setf(std::ios::fixed, std::ios::floatfield);
      _file->precision(6); 
      
      // Write the header
      (*_file)<<"t (s),  dt (s)"<<std::endl;
    }
  }
}
  

TSPlot::~TSPlot()
{

}


void TSPlot::output(CA::Real t, CA::Real dt)
{
  if(_file && _file->good())
  {
    // Write line
    (*_file)<<t<<", "<<dt<<std::endl;
  }
}

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


//! \file ArgsData.cpp
//! contact: a.s.chen [at] exeter.ac.uk
//! \date 2024-05


#include"ArgsData.hpp"
#include"Utilities.hpp"
#include<iostream>
#include<fstream>


// Transform an input string into a model
std::istream& operator>>(std::istream& in, MODEL::Type& m)
{
  CA::String tmp;
  in>>tmp;

  if( CA::compareCaseInsensitive(tmp, "WCA2Dv1" ) )
  {
    m = MODEL::WCA2Dv1; 
    return in;
  }
  if( CA::compareCaseInsensitive(tmp, "WCA2Dv2" ) )
  {
    m = MODEL::WCA2Dv2; 
    return in;
  }
  if( CA::compareCaseInsensitive(tmp, "WCA2D" ) )
  {
    // The default one is WCA2Dv2
    m = MODEL::WCA2Dv2; 
    return in;
  }
  if( CA::compareCaseInsensitive(tmp, "Diffusive" ) )
  {
    m = MODEL::Diffusive; 
    return in;
  }
  if( CA::compareCaseInsensitive(tmp, "Inertial" ) )
  {
    m = MODEL::Inertial; 
    return in;
  }
  m = MODEL::UNKNOWN;
  in.setstate(std::ios::failbit);
  return in;
}


// Transform model variable enum into an output string.
std::ostream& operator<<(std::ostream& out, MODEL::Type& m)
{
  switch(m)
  {
  case MODEL::WCA2Dv1 :
    out<<"WCA2Dv1";
    break;
  case MODEL::WCA2Dv2 :
    out<<"WCA2Dv2";
    break;
  case MODEL::Diffusive :
    out<<"Diffusive";
    break;
  case MODEL::Inertial :
    out<<"Inertial";
    break;
  default:
    out<<"UNKNOWN";
  }

  return out;
}


// Transform an input string into a physical variable enum.
std::istream& operator>>(std::istream& in, PV::Type& pv)
{
  CA::String tmp;
  in>>tmp;

  if( CA::compareCaseInsensitive(tmp, "WD" ) )
  {
    pv = PV::WD; 
    return in;
  }
  if( CA::compareCaseInsensitive(tmp, "WL" ) )
  {
    pv = PV::WL; 
    return in;
  }
  if( CA::compareCaseInsensitive(tmp, "VEL" ) )
  {
    pv = PV::VEL; 
    return in;
  }
  pv = PV::UNKNOWN;
  in.setstate(std::ios::failbit);
  return in;
}


// Transform a physical variable enum into an output string.
std::ostream& operator<<(std::ostream& out, PV::Type& pv)
{
  switch(pv)
  {
  case PV::WD :
    out<<"WD";
    break;
  case PV::WL :
    out<<"WL";
    break;
  case PV::VEL :
    out<<"VEL";
    break;
  default:
    out<<"UNKNOWN";
  }

  return out;
}

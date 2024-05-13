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



//! \file GDAL_Functions.hpp
//!  Contains the structure and function for licensing options
//! \author Mike Gibson, University of Exeter, 
//! contact: m.j.gibson [at] exeter.ac.uk
//! \date 2019-12

#ifdef CAAPI_GDAL_OPTION

#include"ca2D.hpp"
#include"BaseTypes.hpp"
#include"AsciiGrid.hpp"

#include "gdal_priv.h"



int saveRasterAs(CA::AsciiGrid<CA::Real> grid, const std::string driverName, std::string geoCS, std::string fileName);









#endif
//##############################################################################



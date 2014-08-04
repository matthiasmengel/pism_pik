// Copyright (C) 2008-2012 Ed Bueler, Constantine Khroulev, Ricarda Winkelmann,
// Gudfinna Adalgeirsdottir and Andy Aschwanden
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
// version.
//
// PISM is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License
// along with PISM; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

// Implementation of the atmosphere model using constant-in-time precipitation
// and a cosine yearly cycle for near-surface air temperatures.

// This includes the SeaRISE Greenland parameterization.

#include "PATemperaturePIK.hh"
#include "PISMVars.hh"
#include "IceGrid.hh"
#include "pism_options.hh"
#include "PISMTime.hh"

///// PATemperaturePIK

PetscErrorCode PATemperaturePIK::init(PISMVars &vars) {
  PetscErrorCode ierr;

  ierr = verbPrintf(2, grid.com,
     "* Initializing the atmosphere model PATemperaturePIK.\n"
     "  Near-surface air temperature is parameterized as in Huybrechts & De Wolde (1999).\n"
     "  Precipitation is parameterized using a 5 percent increase per degree of warming.\n"); CHKERRQ(ierr);

  ierr = PAYearlyCycle::init(vars); CHKERRQ(ierr);

  // initialize pointers to fields the parameterization depends on:
  usurf = dynamic_cast<IceModelVec2S*>(vars.get("surface_altitude"));
  if (!usurf) SETERRQ(grid.com, 1, "ERROR: surface_altitude is not available");

  lat = dynamic_cast<IceModelVec2S*>(vars.get("latitude"));
  if (!lat) SETERRQ(grid.com, 1, "ERROR: latitude is not available");

  ierr = PISMOptionsIsSet("-paleo_precip", paleo_precipitation_correction); CHKERRQ(ierr);

  if (paleo_precipitation_correction) {
    bool delta_T_set;
    string delta_T_file;

    ierr = PISMOptionsString("-paleo_precip",
                             "Specifies the air temperature offsets file to use with -paleo_precip",
                             delta_T_file, delta_T_set); CHKERRQ(ierr);

    ierr = verbPrintf(2, grid.com, 
                      "  reading delta_T data from forcing file %s for -paleo_precip actions ...\n",
                      delta_T_file.c_str());  CHKERRQ(ierr);

    delta_T = new Timeseries(grid.com, grid.rank, "delta_T",
                             grid.config.get_string("time_dimension_name"));
    ierr = delta_T->set_units("Kelvin", ""); CHKERRQ(ierr);
    ierr = delta_T->set_dimension_units(grid.time->units(), ""); CHKERRQ(ierr);
    ierr = delta_T->set_attr("long_name", "near-surface air temperature offsets");
    CHKERRQ(ierr);
    ierr = delta_T->read(delta_T_file, grid.time->use_reference_date()); CHKERRQ(ierr);
  }

  return 0;
}


// Scale present-day precipitation field as in Pollard & De Conto (2012), Eqn (34b)
PetscErrorCode PATemperaturePIK::mean_precipitation(IceModelVec2S &result) {
  PetscErrorCode ierr;

  ierr = PAYearlyCycle::mean_precipitation(result); CHKERRQ(ierr);

  if ((delta_T != NULL) && paleo_precipitation_correction) {
//    // ... as in Pollard & De Conto (2012):
//    ierr = result.scale(pow (2.0, (0.1* (*delta_T)(t + 0.5 * dt)))); CHKERRQ(ierr); // scale by 2^(0.1*DeltaT)
    ierr = result.scale(pow (1.07, ((*delta_T)(t + 0.5 * dt)))); CHKERRQ(ierr); // scale by 1.05^(DeltaT), i.e., 5% precip increase per degree
  }

  return 0;
}

//! \brief Updates mean annual and mean July near-surface air temperatures.
//! Note that the precipitation rate is time-independent and does not need
//! to be updated.
PetscErrorCode PATemperaturePIK::update(PetscReal my_t, PetscReal my_dt) {
  PetscErrorCode ierr;

  if (lat->has_attr("missing_at_bootstrap")) {
    ierr = PetscPrintf(grid.com, "PISM ERROR: latitude variable was missing at bootstrap;\n"
      "  Atmosphere model depends on latitude and would return nonsense!!\n");
      CHKERRQ(ierr);
    PISMEnd();
  }

  if ((fabs(my_t - t) < 1e-12) &&
      (fabs(my_dt - dt) < 1e-12))
    return 0;

  t  = my_t;
  dt = my_dt;


  PetscScalar **lat_degN, **h;
  ierr = usurf->get_array(h);   CHKERRQ(ierr);
  ierr = lat->get_array(lat_degN); CHKERRQ(ierr);
  ierr = air_temp_mean_annual.begin_access();  CHKERRQ(ierr);
  ierr = air_temp_mean_july.begin_access();  CHKERRQ(ierr); // NOTE: Of course, this is not the mean July, but the mean SUMMER temperature! (Here only denoted as air_temp_mean_july so that the field can be interpreted correctly in the PDD-routines (which were originally meant for the Northern Hemisphere).)

  for (PetscInt i = grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j<grid.ys+grid.ym; ++j) {

//// ALTERNATIVE 1: 
//// annual mean temperature = Martin et al. (2011) parametrization
//// summer mean temperature = anomaly to Huybrechts & DeWolde (1999)
//      air_temp_mean_annual(i,j) = 273.15 + 30 - 0.0075 * h[i][j] - 0.68775 * lat_degN[i][j]*(-1.0);  // surface temperature parameterization as in Martin et al. 2011, Eqn. 2.0.2
//
//      PetscReal gamma_a;
//	if (h[i][j] < 1500.0) {
//	  gamma_a = -0.005102;
//	}else{
//	  gamma_a = -0.014285;
//	}
//
//      PetscReal TMA = 273.15 + 34.46 + gamma_a * h[i][j] - 0.68775 * lat_degN[i][j]*(-1.0); // = TMA, mean annual temperature in Huybrechts & DeWolde (1999)
//      PetscReal TMS = 273.15 + 14.81 - 0.00692 * h[i][j] - 0.27937 * lat_degN[i][j]*(-1.0); // = TMS, mean summer temperature in Huybrechts & DeWolde (1999)
// 
//      air_temp_mean_july(i,j) = air_temp_mean_annual(i,j) + (TMS - TMA);   


//// ALTERNATIVE 2:
//// annual mean temperature = Huybrechts & DeWolde (1	999)
//// summer mean temperature = Huybrechts & DeWolde (1999)
//        PetscReal gamma_a;
//	if (h[i][j] < 1500.0) {
//	  gamma_a = -0.005102;
//	}else{
//	  gamma_a = -0.014285;
//	}
//
//      air_temp_mean_annual(i,j) = 273.15 + 34.46 + gamma_a * h[i][j] - 0.68775 * lat_degN[i][j]*(-1.0); // = TMA, mean annual temperature in Huybrechts & DeWolde (1999)
//      air_temp_mean_july(i,j) = 273.15 + 14.81 - 0.00692 * h[i][j] - 0.27937 * lat_degN[i][j]*(-1.0); // = TMS, mean summer temperature in Huybrechts & DeWolde (1999)  


//// ALTERNATIVE 3:
//// annual mean temperature = Martin et al. (2011)
//// summer mean temperature = Huybrechts & DeWolde (1999)
//      air_temp_mean_annual(i,j) = 273.15 + 30 - 0.0075 * h[i][j] - 0.68775 * lat_degN[i][j]*(-1.0);  // annual mean temperature as in Martin et al. 2011, Eqn. 2.0.2
//      air_temp_mean_july(i,j) = 273.15 + 14.81 - 0.00692 * h[i][j] - 0.27937 * lat_degN[i][j]*(-1.0); // = TMS, mean summer temperature in Huybrechts & DeWolde (1999) 


//// ALTERNATIVE 4:
//// annual mean temperature = new parametrization based on multiple regression analysis of ERA INTERIM data
//// summer mean temperature = new parametrization based on multiple regression analysis of ERA INTERIM data
	air_temp_mean_annual(i,j) = 273.15 + 29.2 - 0.0082 * h[i][j] - 0.576 * lat_degN[i][j]*(-1.0);
	air_temp_mean_july(i,j)   = 273.15 + 16.5 - 0.0068 * h[i][j] - 0.248 * lat_degN[i][j]*(-1.0);


//// ALTERNATIVE 5:
//// annual mean temperature = new parametrization based on multiple regression analysis of ERA INTERIM data with sin(lat)
//// summer mean temperature = new parametrization based on multiple regression analysis of ERA INTERIM data with sin(lat)
//	air_temp_mean_annual(i,j) = 273.15 - 2.0 -0.0082*h[i][j] + 18.4 * (sin(3.1415*lat_degN[i][j]/180)+0.8910)/(1-0.8910);
//	air_temp_mean_july(i,j)   = 273.15 + 3.2 -0.0067*h[i][j] +  8.3 * (sin(3.1415*lat_degN[i][j]/180)+0.8910)/(1-0.8910);

    }
  }

  ierr = usurf->end_access();   CHKERRQ(ierr);
  ierr = lat->end_access(); CHKERRQ(ierr);
  ierr = air_temp_mean_annual.end_access();  CHKERRQ(ierr);
  ierr = air_temp_mean_july.end_access();  CHKERRQ(ierr);

  return 0;
}

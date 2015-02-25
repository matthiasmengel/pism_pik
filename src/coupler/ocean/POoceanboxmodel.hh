// Copyright (C) 2011, 2012 PISM Authors
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

#ifndef _POOCEANBOXMODEL_H_
#define _POOCEANBOXMODEL_H_

#include "PGivenClimate.hh"
#include "POModifier.hh"
#include "Timeseries.hh"

//! A class defining the interface of a PISM ocean model modifier.

//! \brief A class implementing an ocean model.
//! Computes the subshelf melt/refreezing rate based on a simple ocean box model by Olbers & Hellmer (2010).
//class POoceanboxmodel : public PISMOceanModel {
//public:
   //POoceanboxmodel(IceGrid &g, const NCConfigVariable &conf);


class POoceanboxmodel : public PGivenClimate<POModifier,PISMOceanModel>
{
public:
  POoceanboxmodel(IceGrid &g, const NCConfigVariable &conf)
    : PGivenClimate<POModifier,PISMOceanModel>(g, conf, NULL)
  {
    temp_name       = "thetao";
    mass_flux_name  = "salinity"; //NOTE: salinity_name instead of mass_flux_name
    option_prefix   = "-ocean_oceanboxmodel";

    ocean_oceanboxmodel_deltaT_set = false;
    delta_T = NULL;

  }


  virtual ~POoceanboxmodel() {}

  virtual PetscErrorCode init(PISMVars &vars);
  virtual PetscErrorCode update(PetscReal my_t, PetscReal my_dt);
  // virtual PetscErrorCode update(PetscReal my_t, PetscReal my_dt) { t = my_t; dt = my_dt; return 0; } // do nothing //FIXME is this correct?

  virtual PetscErrorCode sea_level_elevation(PetscReal &result) { //FIXME is this obsolete?
    result = sea_level;
    return 0;
  }
  virtual PetscErrorCode identifyICERISES();
  virtual PetscErrorCode extentOfIceShelves();
  virtual PetscErrorCode identifyBOXMODELmask();
  virtual PetscErrorCode extendGLBox();
  virtual PetscErrorCode extendIFBox();
  virtual PetscErrorCode oceanTemperature();
  virtual PetscErrorCode basalMeltRateForGroundingLineBox();
  virtual PetscErrorCode basalMeltRateForIceFrontBox();
  virtual PetscErrorCode basalMeltRateForOtherShelves();

  virtual PetscErrorCode shelf_base_temperature(IceModelVec2S &result);
  virtual PetscErrorCode shelf_base_mass_flux(IceModelVec2S &result);

  virtual PetscErrorCode write_variables(set<string> vars, string filename); // FIXME included by Ronja to write the variables to extra files. Is there a smarter way?

protected:
  IceModelVec2S *ice_thickness, *topg, *lat, *lon, *basins;	// not owned by this class

  IceModelVec2Int *mask;  // not owned by this class

  IceModelVec2S ICERISESmask, BOXMODELmask, DRAINAGEmask, Soc, Soc_base, Toc, Toc_base, Toc_inCelsius, T_star, Toc_anomaly, overturning, heatflux, basalmeltrate_shelf;

  bool ocean_oceanboxmodel_deltaT_set, drainageBasins_set, drainageBasins_OH10_set, exicerises_set;

  Timeseries *delta_T;

  static const int shelf_unidentified, noshelf;
  static const int box_unidentified, box_noshelf, box_GL, box_neighboring, box_IF, maskfloating, maskocean, maskgrounded;

  PetscInt numberOfBasins;

  vector<double>  counter,
                  counter_GLbox,
                  counter_CFbox,
                  k_basins;
  PetscScalar     counter_box_unidentified; 

  vector<double>  mean_salinity_GLbox_vector,
                  mean_meltrate_GLbox_vector,
                  mean_overturning_GLbox_vector; 

  vector<double>  Toc_base_vec,
                  Soc_base_vec,
                  gamma_T_star_vec,
                  C_vec;

};

#endif /* _POOCEANBOXMODEL_H_ */

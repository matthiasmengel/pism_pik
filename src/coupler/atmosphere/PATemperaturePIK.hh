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

#ifndef __PATemperaturePIK_hh
#define __PATemperaturePIK_hh

#include "PAYearlyCycle.hh"
#include "Timeseries.hh"

//! \brief A modification of PAYearlyCycle tailored for the
//! SeaRISE-Greenland assessment. Uses the Fausto [\ref Faustoetal2009]
//! present-day temperature parameterization and stored precipitation data.
//! Adds the precipitation correction for spin-ups.
class PATemperaturePIK : public PAYearlyCycle {
public:
  PATemperaturePIK(IceGrid &g, const NCConfigVariable &conf)
    : PAYearlyCycle(g, conf)
  {
    paleo_precipitation_correction = false;
    delta_T = NULL;
  }

  virtual ~PATemperaturePIK()
  {
    delete delta_T;
  }
  virtual PetscErrorCode init(PISMVars &vars);
  virtual PetscErrorCode update(PetscReal my_t, PetscReal my_dt);
  virtual PetscErrorCode mean_precipitation(IceModelVec2S &result);
protected:
  bool paleo_precipitation_correction;
  Timeseries *delta_T;
  IceModelVec2S *ice_thickness;	//!< current ice thickness produced by IceModel.
  IceModelVec2S *lat, *usurf;
};


#endif	// __PATemperaturePIK_hh

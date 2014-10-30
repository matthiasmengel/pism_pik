// Copyright (C) 2008-2012 Ed Bueler, Constantine Khroulev, Ricarda Winkelmann,
// Gudfinna Adalgeirsdottir, Andy Aschwanden and Torsten Albrecht
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

#include "POoceanboxmodel.hh"
#include "PISMVars.hh"
#include "IceGrid.hh"
#include "pism_options.hh"


/*!
 The aim of these routines is to simulate the ocean circulation underneath the ice shelves and compute the basal melt/refreezing rates according to the ocean box model described in olbers_hellmer10.
*/


// POoceanboxmodel::POoceanboxmodel(IceGrid &g, const NCConfigVariable &conf)
//   : PISMOceanModel(g, conf) {

//   shelfbmassflux.init_2d("shelfbmassflux", g);
//   shelfbmassflux.set_string("pism_intent", "climate_state");
//   shelfbmassflux.set_string("long_name",
//                             "ice mass flux from ice shelf base (positive flux is loss from ice shelf)");
//   shelfbmassflux.set_units("m s-1");
//   shelfbmassflux.set_glaciological_units("m year-1");

//   shelfbtemp.init_2d("shelfbtemp", g);
//   shelfbtemp.set_string("pism_intent", "climate_state");
//   shelfbtemp.set_string("long_name",
//                         "absolute temperature at ice shelf base");
//   shelfbtemp.set_units("Kelvin");
// }

PetscErrorCode POoceanboxmodel::init(PISMVars &vars) {
  PetscErrorCode ierr;

  ierr = verbPrintf(2, grid.com,
                    "* Initializing the ocean box model (based on Olbers & Hellmer (2010)...\n"); CHKERRQ(ierr);

  ierr = process_options(); CHKERRQ(ierr);

  ierr = set_vec_parameters("", ""); CHKERRQ(ierr);

  ierr = temp.create(grid, temp_name, false); CHKERRQ(ierr);
  ierr = mass_flux.create(grid, mass_flux_name, false); CHKERRQ(ierr); //NOTE: salinity instead of mass_flux

  ierr = temp.set_attrs("climate_forcing", "absolute temperature at ice shelf base",
                        "Kelvin", ""); CHKERRQ(ierr); // temp_boundlayer is not a standard name
  ierr = mass_flux.set_attrs("climate_forcing", "ocean salinity at ice shelf",
                             "g/kg",
                             ""); CHKERRQ(ierr); // salinity_boundlayer is not a standard name

  ierr = temp.init(filename); CHKERRQ(ierr);
  ierr = mass_flux.init(filename); CHKERRQ(ierr);

  // mask to identify different basins
  ierr = SHELFmask.create(grid, "SHELFmask", true); CHKERRQ(ierr);
  ierr = SHELFmask.set_attrs("model_state", "mask displaying shelf regions","", ""); CHKERRQ(ierr);

  // mask to identify the ocean boxes
  ierr = BOXMODELmask.create(grid, "BOXMODELmask", true); CHKERRQ(ierr);
  ierr = BOXMODELmask.set_attrs("model_state", "mask displaying ocean box model grid","", ""); CHKERRQ(ierr);

  // salinity
  ierr = Soc.create(grid, "Soc", true); CHKERRQ(ierr);
  ierr = Soc.set_attrs("model_state", "ocean salinity field","", "ocean salinity field"); CHKERRQ(ierr);  //NOTE unit=psu
  ierr = Soc_base.create(grid, "Soc_base", true); CHKERRQ(ierr);
  ierr = Soc_base.set_attrs("model_state", "ocean base salinity field","", "ocean base salinity field"); CHKERRQ(ierr);  //NOTE unit=psu

  // temperature
  ierr = Toc.create(grid, "Toc", true); CHKERRQ(ierr);
  ierr = Toc.set_attrs("model_state", "ocean temperature field","K", "ocean temperature field"); CHKERRQ(ierr);
  ierr = Toc_base.create(grid, "Toc_base", true); CHKERRQ(ierr);
  ierr = Toc_base.set_attrs("model_state", "ocean base temperature","K", "ocean base temperature"); CHKERRQ(ierr);

  ierr = Toc_inCelsius.create(grid, "Toc_inCelsius", true); CHKERRQ(ierr);
  ierr = Toc_inCelsius.set_attrs("model_state", "ocean box model temperature field","degree C", "ocean box model temperature field"); CHKERRQ(ierr);

  ierr = T_star.create(grid, "T_star", true); CHKERRQ(ierr);
  ierr = T_star.set_attrs("model_state", "T_star field","degree C", "T_star field"); CHKERRQ(ierr);

  ierr = Toc_anomaly.create(grid, "Toc_anomaly", true); CHKERRQ(ierr);
  ierr = Toc_anomaly.set_attrs("model_state", "ocean temperature anomaly","K", "ocean temperature anomaly"); CHKERRQ(ierr);


  // overturning rate
  ierr = overturning.create(grid, "overturning", true); CHKERRQ(ierr);
  ierr = overturning.set_attrs("model_state", "cavity overturning","m^3 s-1", "cavity overturning"); CHKERRQ(ierr); // no CF standard_name ??

  // heat flux
  ierr = heatflux.create(grid, "ocean heat flux", true); CHKERRQ(ierr);
  ierr = heatflux.set_attrs("climate_state", "ocean heat flux", "W/m^2", ""); CHKERRQ(ierr);

  // basal melt rate
  ierr = basalmeltrate_shelf.create(grid, "basal melt rate from ocean box model", true); CHKERRQ(ierr);
  ierr = basalmeltrate_shelf.set_attrs("climate_state", "basal melt rate from ocean box model", "m/s", ""); CHKERRQ(ierr);


  ice_thickness = dynamic_cast<IceModelVec2S*>(vars.get("land_ice_thickness"));
  if (!ice_thickness) { SETERRQ(grid.com, 1, "ERROR: ice thickness is not available"); }

  topg = dynamic_cast<IceModelVec2S*>(vars.get("bedrock_altitude"));
  if (topg == NULL) SETERRQ(grid.com, 1, "bedrock topography is not available");

  lat = dynamic_cast<IceModelVec2S*>(vars.get("latitude"));
  if (lat == NULL) SETERRQ(grid.com, 1, "latitude is not available");

  lon = dynamic_cast<IceModelVec2S*>(vars.get("longitude"));
  if (!lon) { SETERRQ(grid.com, 1, "ERROR: longitude is not available"); }

  ierr = PISMOptionsIsSet("-drainageBasins", drainageBasins_set); CHKERRQ(ierr);

  if (drainageBasins_set){
    basins = dynamic_cast<IceModelVec2S*>(vars.get("drainage_basins"));
    if (!basins) { SETERRQ(grid.com, 1, "ERROR: drainage basins is not available"); }
  }
  ierr = PISMOptionsIsSet("-obm_deltaT", ocean_oceanboxmodel_deltaT_set); CHKERRQ(ierr);

  if (ocean_oceanboxmodel_deltaT_set) {
    bool delta_T_set;
    string delta_T_file;

    ierr = PISMOptionsString("-obm_deltaT",
                             "Specifies the ocean temperature offsets file to use with -obm_deltaT",
                             delta_T_file, delta_T_set); CHKERRQ(ierr);

    ierr = verbPrintf(2, grid.com,
                      "  reading delta_T data from forcing file %s for -obm_deltaT actions ...\n",
                      delta_T_file.c_str());  CHKERRQ(ierr);

    delta_T = new Timeseries(grid.com, grid.rank, "delta_T",
                             grid.config.get_string("time_dimension_name"));
    ierr = delta_T->set_units("Kelvin", ""); CHKERRQ(ierr);
    ierr = delta_T->set_dimension_units(grid.time->units(), ""); CHKERRQ(ierr);
    ierr = delta_T->set_attr("long_name", "ocean temperature offsets");
    CHKERRQ(ierr);
    ierr = delta_T->read(delta_T_file, grid.time->use_reference_date()); CHKERRQ(ierr);
  }

  return 0;
}

const int POoceanboxmodel::shelf_unidentified = -99.0; // This should never show up in the .nc-files.
// const int POoceanboxmodel::basin_RossSea = -1.0; // This should never show up in the .nc-files.
// const int POoceanboxmodel::basin_WeddellSea = -2.0; // This should never show up in the .nc-files.
// const int POoceanboxmodel::basin_EastAntarctica = -3.0; // This should never show up in the .nc-files.
// const int POoceanboxmodel::basin_AmundsenSea = -4.0;  // This should never show up in the .nc-files.
const int POoceanboxmodel::noshelf = 0.0;
const int POoceanboxmodel::shelf_RossSea = 1.0;
const int POoceanboxmodel::shelf_WeddellSea = 2.0;
const int POoceanboxmodel::shelf_EastAntarctica = 3.0;
const int POoceanboxmodel::shelf_AmundsenSea = 4.0;

const int POoceanboxmodel::box_unidentified = -99.0;     // This should never show up in the .nc-files.
const int POoceanboxmodel::box_near_GL = -1.0; // This should never show up in the .nc-files.
const int POoceanboxmodel::box_noshelf = 0.0;
const int POoceanboxmodel::box_GL = 1.0;  // ocean box covering the grounding line region
const int POoceanboxmodel::box_IF = 2.0;  // ocean box covering the rest of the ice shelf

// TODO read shelf mask from file
// TODO why is there no ice front box?
// TODO Toc_base and Soc_base need to be read from file instead! Jippie!

//! Find all ice shelves in four pre-defined regions: Ross Sea, Weddell Sea, Amundsen Sea, East Antarctica
PetscErrorCode POoceanboxmodel::AntarcticBasins() {
  PetscErrorCode ierr;
  ierr = verbPrintf(2, grid.com,
                    "NOW in AntarcticBasins rountine\n"); CHKERRQ(ierr);
  const PetscScalar rhoi = config.get("ice_density");
  const PetscScalar rhow = config.get("sea_water_density");

  ierr = ice_thickness->begin_access();   CHKERRQ(ierr);
  ierr = topg->begin_access();   CHKERRQ(ierr);
  ierr = lat->begin_access();   CHKERRQ(ierr);
  ierr = lon->begin_access();   CHKERRQ(ierr);
  ierr = SHELFmask.begin_access();   CHKERRQ(ierr);

  // STEP 1: Declare floating cells which should certainly belong to a specific shelf
  // shelf regions are given as lon, lat - ranges; mark each floating box accordingly
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      const PetscScalar hgrounded = (*topg)(i,j) + (*ice_thickness)(i,j);
      const PetscScalar hfloating = (1.0 - rhoi/rhow) * (*ice_thickness)(i,j); // FIXME This assumes currentSeaLevel=0 !
      bool herefloating=false;

      if ((*ice_thickness)(i,j) > 0.0 && hgrounded < hfloating) {herefloating = true;}


      if (herefloating){ // FIXME adjust regions!
        if(( (((*lon)(i,j)+180.0)>=0.0 && ((*lon)(i,j)+180.0)<50.0) || (((*lon)(i,j)+180.0)>=335.0 && ((*lon)(i,j)+180.0)<=360.0)) &&
	       (*lat)(i,j)>=-85.0 && (*lat)(i,j)<-80.0 ){
			SHELFmask(i,j) = shelf_RossSea;
        }
        else if(((*lon)(i,j)+180.0)>=90.0 && ((*lon)(i,j)+180.0)<145.0 && (*lat)(i,j)>=-85.0 && (*lat)(i,j)<-78.0){
			SHELFmask(i,j) = shelf_WeddellSea;
        }
        else if(((*lon)(i,j)+180.0)>=245.0 && ((*lon)(i,j)+180.0)<255.0 && (*lat)(i,j)>=-75.0 && (*lat)(i,j)<-70.0){
			SHELFmask(i,j) = shelf_EastAntarctica;
        }
        else if(((*lon)(i,j)+180.0)>=76.0 && ((*lon)(i,j)+180.0)<82.0 && (*lat)(i,j)>=-76.0 && (*lat)(i,j)<-74.0){
			SHELFmask(i,j) = shelf_AmundsenSea;
	    }
	    else{
			SHELFmask(i,j) = shelf_RossSea;//shelf_unidentified;
		}

      }
      if(!herefloating){
		  SHELFmask(i,j) = noshelf;
	  }
    } // end j
  } // end i

  ierr = ice_thickness->end_access();   CHKERRQ(ierr);
  ierr = topg->end_access();   CHKERRQ(ierr);
  ierr = lat->end_access();   CHKERRQ(ierr);
  ierr = lon->end_access();   CHKERRQ(ierr);
  ierr = SHELFmask.end_access();   CHKERRQ(ierr);
  ierr = SHELFmask.beginGhostComm(); CHKERRQ(ierr);
  ierr = SHELFmask.endGhostComm(); CHKERRQ(ierr);


  //// STEP 2: Find the rest of each shelf recursively (this is in order to find the WHOLE ice shelf whose extent we cannot know beforehand).
  //// NOTE This excludes the case that a (BOXMODELmask=box_IF)-part belongs to one shelf and the (BOXMODELmask=box_GL)-part to another which MUST NOT occur. If it did occur, the following would happen: Since the basal melt rate of the (BOXMODELmask=box_GL)-part is needed to compute the basal melt rate of the (BOXMODELmask=box_IF)-part, this case would lead to basalmeltrate_shelf=NaN which would lead to ice_thickness=NaN which would lead to ks=NaN in temperatureStep() which would lead to a Segmentation Violation!
  //bool done = false;
  //PetscInt loopcount = 0;
  //while(! done){
    //done = true;

    //ierr = SHELFmask.begin_access();   CHKERRQ(ierr);
    //for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
      //for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {

	//if ((SHELFmask(i,j) == shelf_unidentified) &&
	    //(SHELFmask(i-1,j) == shelf_RossSea || SHELFmask(i,j-1) == shelf_RossSea || SHELFmask(i,j+1) == shelf_RossSea || SHELFmask(i+1,j) == shelf_RossSea)){
		//SHELFmask(i,j) = shelf_RossSea;
		//done = false;
	//}
	//if ((SHELFmask(i,j) == shelf_unidentified) &&
	    //(SHELFmask(i-1,j) == shelf_WeddellSea || SHELFmask(i,j-1) == shelf_WeddellSea || SHELFmask(i,j+1) == shelf_WeddellSea || SHELFmask(i+1,j) == shelf_WeddellSea)){
		//SHELFmask(i,j) = shelf_WeddellSea;
		//done = false;
	//}
	//if ((SHELFmask(i,j) == shelf_unidentified) &&
	    //(SHELFmask(i-1,j) == shelf_EastAntarctica || SHELFmask(i,j-1) == shelf_EastAntarctica || SHELFmask(i,j+1) == shelf_EastAntarctica || SHELFmask(i+1,j) == shelf_EastAntarctica)){
		//SHELFmask(i,j) = shelf_EastAntarctica;
		//done = false;
	//}
	//if ((SHELFmask(i,j) == shelf_unidentified) &&
	    //(SHELFmask(i-1,j) == shelf_AmundsenSea || SHELFmask(i,j-1) == shelf_AmundsenSea || SHELFmask(i,j+1) == shelf_AmundsenSea || SHELFmask(i+1,j) == shelf_AmundsenSea)){
		//SHELFmask(i,j) = shelf_AmundsenSea;
		//done = false;
	//} // NOTE The ordering above puts a certain emphasis on RossSea > WeddellSea > AmundsenSea because it could also be that, e.g., one cell has BOTH a RossSea- and a WeddellSea-neighbour. Then it is declared "near Ross" in the first if-loop and is off-limits for Weddell. This case will rarely happen.
      //}
    //}

    //ierr = SHELFmask.end_access(); CHKERRQ(ierr); // FIXME I deleted a second i,j-loop here because I think it is obsolete - CHECK!

    //ierr = SHELFmask.beginGhostComm(); CHKERRQ(ierr);
    //ierr = SHELFmask.endGhostComm(); CHKERRQ(ierr);

    //// We're "done" only if we are done on *all* processor sub-domains:
    //int flag = done;
    //MPI_Allreduce(MPI_IN_PLACE, &flag, 1, MPI_INT, MPI_LAND, grid.com);
    //done = flag;

    //loopcount += 1;
  //}

  ierr = extentOfIceShelves(); CHKERRQ(ierr);
  ierr = identifyGroundingLineBox(); CHKERRQ(ierr);
  ierr = identifyIceFrontBox(); CHKERRQ(ierr);

  return 0;
}


//! Compute the extent of each ice shelf
PetscErrorCode POoceanboxmodel::extentOfIceShelves() {
  PetscErrorCode ierr;
  ierr = verbPrintf(2, grid.com,
                    "NOW in extent of ice shelves rountine\n"); CHKERRQ(ierr);
  const PetscScalar rhoi = config.get("ice_density");
  const PetscScalar rhow = config.get("sea_water_density");

  // see olbers_hellmer10, table 1 for observed lengths of the ice shelves // FIXME config!
  const PetscScalar RossLengthInitial=ceil( 800 / (grid.dx* 1e-3)); // length is given in km   // NOTE: This only works when initializing from a present-day configuration!
  const PetscScalar WeddellLengthInitial=ceil( 620 / (grid.dx* 1e-3));
  const PetscScalar EastAntarcticaLengthInitial=ceil( 540 / (grid.dx* 1e-3));
  const PetscScalar AmundsenLengthInitial=ceil( 70 / (grid.dx* 1e-3));

  PetscScalar lcounterRoss=0.0; PetscScalar lcounterWeddell=0.0; PetscScalar lcounterEastAntarctica=0.0; PetscScalar lcounterAmundsen=0.0;

  ierr = ice_thickness->begin_access();   CHKERRQ(ierr);
  ierr = topg->begin_access();   CHKERRQ(ierr);
  ierr = SHELFmask.begin_access(); CHKERRQ(ierr);
  ierr = BOXMODELmask.begin_access(); CHKERRQ(ierr);

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {

      const PetscScalar hgrounded = (*topg)(i,j) + (*ice_thickness)(i,j);
      const PetscScalar hfloating = (1.0 - rhoi/rhow) * (*ice_thickness)(i,j); // FIXME This assumes currentSeaLevel=0 !
      bool herefloating=false;
      if ((*ice_thickness)(i,j) > 0.0 && hgrounded < hfloating) {herefloating = true;}

      bool grounded_N=false; bool grounded_NE=false; bool grounded_NW=false;
      bool grounded_S=false; bool grounded_SE=false; bool grounded_SW=false;
      bool grounded_E=false; bool grounded_W=false;

      if (herefloating){
	  //BOXMODELmask(i,j) = box_unidentified;
	  // Count cells for specific shelves
	  if (SHELFmask(i,j) == shelf_RossSea)        lcounterRoss++;
	  if (SHELFmask(i,j) == shelf_WeddellSea)     lcounterWeddell++;
	  if (SHELFmask(i,j) == shelf_EastAntarctica) lcounterEastAntarctica++;
	  if (SHELFmask(i,j) == shelf_AmundsenSea)    lcounterAmundsen++;

	  PetscScalar hg_N = (*topg)(i,j+1) + (*ice_thickness)(i,j+1),   hf_N = (1.0 - rhoi/rhow) * (*ice_thickness)(i,j+1);
	  if((*ice_thickness)(i,j+1) > 0.0 && hg_N > hf_N){grounded_N = true;}
	  PetscScalar hg_NE = (*topg)(i+1,j+1) + (*ice_thickness)(i+1,j+1),   hf_NE = (1.0 - rhoi/rhow) * (*ice_thickness)(i+1,j+1);
	  if((*ice_thickness)(i+1,j+1) > 0.0 && hg_NE > hf_NE){grounded_NE = true;}
	  PetscScalar hg_NW = (*topg)(i-1,j+1) + (*ice_thickness)(i-1,j+1),   hf_NW = (1.0 - rhoi/rhow) * (*ice_thickness)(i-1,j+1);
	  if((*ice_thickness)(i-1,j+1) > 0.0 && hg_NW > hf_NW){grounded_NW = true;}
	  PetscScalar hg_S = (*topg)(i,j-1) + (*ice_thickness)(i,j-1),   hf_S = (1.0 - rhoi/rhow) * (*ice_thickness)(i,j-1);
	  if((*ice_thickness)(i,j-1) > 0.0 && hg_S > hf_S){grounded_S = true;}
	  PetscScalar hg_SE = (*topg)(i+1,j-1) + (*ice_thickness)(i+1,j-1),   hf_SE = (1.0 - rhoi/rhow) * (*ice_thickness)(i+1,j-1);
	  if((*ice_thickness)(i+1,j-1) > 0.0 && hg_SE > hf_SE){grounded_SE = true;}
	  PetscScalar hg_SW = (*topg)(i-1,j-1) + (*ice_thickness)(i-1,j-1),   hf_SW = (1.0 - rhoi/rhow) * (*ice_thickness)(i-1,j-1);
	  if((*ice_thickness)(i-1,j-1) > 0.0 && hg_SW > hf_SW){grounded_SW = true;}
	  PetscScalar hg_E = (*topg)(i+1,j) + (*ice_thickness)(i+1,j),   hf_E = (1.0 - rhoi/rhow) * (*ice_thickness)(i+1,j);
	  if((*ice_thickness)(i+1,j) > 0.0 && hg_E > hf_E){grounded_E = true;}
	  PetscScalar hg_W = (*topg)(i-1,j) + (*ice_thickness)(i-1,j),   hf_W = (1.0 - rhoi/rhow) * (*ice_thickness)(i-1,j);
	  if((*ice_thickness)(i-1,j) > 0.0 && hg_W > hf_W){grounded_W = true;}

	  if (grounded_N || grounded_NE || grounded_NW || grounded_S || grounded_SE || grounded_SW || grounded_E || grounded_W) { // i.e., a floating box next to grounded box
		BOXMODELmask(i,j) = box_GL;

	  }else{ // i.e., all other floating boxes
		//ierr = verbPrintf(2, grid.com, "setting box model mask to box_unidentified \n"); CHKERRQ(ierr);
	    BOXMODELmask(i,j) = box_unidentified;
	  }

	}else{ // i.e., not floating
	  BOXMODELmask(i,j) = box_noshelf;
	}
	//if (BOXMODELmask(i,j)== box_GL){ ierr = verbPrintf(2, grid.com, "BOXMODELmask of i=%d, j=%d is %e \n", i, j, BOXMODELmask(i,j)); CHKERRQ(ierr); }
    }
  }

  ierr = ice_thickness->end_access();   CHKERRQ(ierr);
  ierr = topg->end_access();   CHKERRQ(ierr);
  ierr = SHELFmask.end_access(); CHKERRQ(ierr);
  ierr = BOXMODELmask.end_access(); CHKERRQ(ierr);

  ierr = PISMGlobalSum(&lcounterRoss, &counterRoss, grid.com); CHKERRQ(ierr);
  ierr = PISMGlobalSum(&lcounterWeddell, &counterWeddell, grid.com); CHKERRQ(ierr);
  ierr = PISMGlobalSum(&lcounterEastAntarctica, &counterEastAntarctica, grid.com); CHKERRQ(ierr);
  ierr = PISMGlobalSum(&lcounterAmundsen, &counterAmundsen, grid.com); CHKERRQ(ierr);

  // This gives the initial area in terms of number of boxes
  // FIXME the counterXXX_init is 0, so k_XXX is infinity, at the moment this is not used,
  // TOD find a new routine to set the amount of iterations for the groundling line box definition
  //if (firstOceanBoxModelStep==true){
    counterRoss_init=counterRoss;
    counterWeddell_init=counterWeddell;
    counterEastAntarctica_init=counterEastAntarctica;
    counterAmundsen_init=counterAmundsen;
  //}

  // The extent of the GL- and CF boxes depends on these parameters:  // FIXME config!
  k_Ross=2;//ceil((counterRoss/counterRoss_init) * 0.4 * 0.33 * RossLengthInitial);
  k_Weddell=2;//ceil((counterWeddell/counterWeddell_init) * 0.2 * 0.33 * WeddellLengthInitial);
  k_EastAntarctica=2;//ceil((counterEastAntarctica/counterEastAntarctica_init) * 0.2 * 0.33 * EastAntarcticaLengthInitial);
  k_Amundsen=2;//ceil((counterAmundsen/counterAmundsen_init) * 0.1 * 0.33 * AmundsenLengthInitial);

  //firstOceanBoxModelStep = false;
  return 0;
}


//! Identify all cells which belong to the grounding line box
PetscErrorCode POoceanboxmodel::identifyGroundingLineBox() {
  PetscErrorCode ierr;
  ierr = verbPrintf(2, grid.com,"NOW in identifying gl box rountine\n"); CHKERRQ(ierr);
  bool done=false;
  PetscScalar kcounter=0.0;

  while(done == false && kcounter < PetscMax(PetscMax(PetscMax(k_Ross,k_Weddell),k_EastAntarctica),k_Amundsen)){
    done = true;
    kcounter++;

    ierr = SHELFmask.begin_access(); CHKERRQ(ierr);
    ierr = BOXMODELmask.begin_access(); CHKERRQ(ierr);

    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
	if (kcounter < k_Ross && SHELFmask(i,j) == shelf_RossSea){ // this implies that i,j is floating
	  if (BOXMODELmask(i,j) == box_unidentified &&
	    (BOXMODELmask(i-1,j) == box_GL || BOXMODELmask(i,j-1) == box_GL || BOXMODELmask(i,j+1) == box_GL || BOXMODELmask(i+1,j) == box_GL || BOXMODELmask(i+1,j+1) == box_GL ||  BOXMODELmask(i+1,j-1) == box_GL ||  BOXMODELmask(i-1,j+1) == box_GL || BOXMODELmask(i-1,j-1) == box_GL )){
		BOXMODELmask(i,j) = box_near_GL;
		done = false;
	  }
	}
	if (kcounter < k_Weddell && SHELFmask(i,j) == shelf_WeddellSea){ // this implies that i,j is floating
	  if (BOXMODELmask(i,j) == box_unidentified &&
	    (BOXMODELmask(i-1,j) == box_GL || BOXMODELmask(i,j-1) == box_GL || BOXMODELmask(i,j+1) == box_GL || BOXMODELmask(i+1,j) == box_GL || BOXMODELmask(i+1,j+1) == box_GL ||  BOXMODELmask(i+1,j-1) == box_GL ||  BOXMODELmask(i-1,j+1) == box_GL || BOXMODELmask(i-1,j-1) == box_GL )){
		BOXMODELmask(i,j) = box_near_GL;
		done = false;
	  }
	}
	if (kcounter < k_EastAntarctica && SHELFmask(i,j) == shelf_EastAntarctica){ // this implies that i,j is floating
	  if (BOXMODELmask(i,j) == box_unidentified &&
      	    (BOXMODELmask(i-1,j) == box_GL || BOXMODELmask(i,j-1) == box_GL || BOXMODELmask(i,j+1) == box_GL || BOXMODELmask(i+1,j) == box_GL || BOXMODELmask(i+1,j+1) == box_GL ||  BOXMODELmask(i+1,j-1) == box_GL ||  BOXMODELmask(i-1,j+1) == box_GL || BOXMODELmask(i-1,j-1) == box_GL )){
		BOXMODELmask(i,j) = box_near_GL;
		done = false;
	  }
	}
	if (kcounter < k_Amundsen && SHELFmask(i,j) == shelf_AmundsenSea){ // this implies that i,j is floating
	  if (BOXMODELmask(i,j) == box_unidentified &&
	    (BOXMODELmask(i-1,j) == box_GL || BOXMODELmask(i,j-1) == box_GL || BOXMODELmask(i,j+1) == box_GL || BOXMODELmask(i+1,j) == box_GL || BOXMODELmask(i+1,j+1) == box_GL ||  BOXMODELmask(i+1,j-1) == box_GL ||  BOXMODELmask(i-1,j+1) == box_GL || BOXMODELmask(i-1,j-1) == box_GL )){
		BOXMODELmask(i,j) = box_near_GL;
		done = false;
	  }
	}

      }
    }

    ierr = SHELFmask.end_access(); CHKERRQ(ierr);
    ierr = BOXMODELmask.end_access(); CHKERRQ(ierr);
    ierr = BOXMODELmask.beginGhostComm(); CHKERRQ(ierr);
    ierr = BOXMODELmask.endGhostComm(); CHKERRQ(ierr);

    ierr = BOXMODELmask.begin_access(); CHKERRQ(ierr);
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
	if (BOXMODELmask(i,j) == box_near_GL){
	  BOXMODELmask(i,j) = box_GL;
	}
      }
    }
    ierr = BOXMODELmask.end_access(); CHKERRQ(ierr);
    ierr = BOXMODELmask.beginGhostComm(); CHKERRQ(ierr);
    ierr = BOXMODELmask.endGhostComm(); CHKERRQ(ierr);

    // We're "done" only if we are done on *all* processor sub-domains:
    int flag = done;
    MPI_Allreduce(MPI_IN_PLACE, &flag, 1, MPI_INT, MPI_LAND, grid.com);
    done = flag;
  }

  return 0;
}


/*!
Define the rest of each ice shelf as the respective near-ice front-box i.
*/
PetscErrorCode POoceanboxmodel::identifyIceFrontBox() {
  PetscErrorCode ierr;
    ierr = verbPrintf(2, grid.com,
                    "NOW in identifying ice front box rountine\n"); CHKERRQ(ierr);
  // counter to determine the size of the grounding line and ice front boxes for each shelf
  PetscScalar lcounterRoss_CFbox=0.0;
  PetscScalar lcounterWeddell_CFbox=0.0;
  PetscScalar lcounterEastAntarctica_CFbox=0.0;
  PetscScalar lcounterAmundsen_CFbox=0.0;

  ierr = SHELFmask.begin_access(); CHKERRQ(ierr);
  ierr = BOXMODELmask.begin_access(); CHKERRQ(ierr);

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (BOXMODELmask(i,j) == box_unidentified){
	BOXMODELmask(i,j) = box_IF;
	if (SHELFmask(i,j) == shelf_RossSea) lcounterRoss_CFbox++;
	if (SHELFmask(i,j) == shelf_WeddellSea) lcounterWeddell_CFbox++;
	if (SHELFmask(i,j) == shelf_EastAntarctica) lcounterEastAntarctica_CFbox++;
	if (SHELFmask(i,j) == shelf_AmundsenSea) lcounterAmundsen_CFbox++;
      }
    }
  }

  ierr = SHELFmask.end_access(); CHKERRQ(ierr);
  ierr = BOXMODELmask.end_access(); CHKERRQ(ierr);

  ierr = PISMGlobalSum(&lcounterRoss_CFbox, &counterRoss_CFbox, grid.com); CHKERRQ(ierr);
  ierr = PISMGlobalSum(&lcounterWeddell_CFbox, &counterWeddell_CFbox, grid.com); CHKERRQ(ierr);
  ierr = PISMGlobalSum(&lcounterEastAntarctica_CFbox, &counterEastAntarctica_CFbox, grid.com); CHKERRQ(ierr);
  ierr = PISMGlobalSum(&lcounterAmundsen_CFbox, &counterAmundsen_CFbox, grid.com); CHKERRQ(ierr);

  counterRoss_GLbox=counterRoss-counterRoss_CFbox;
  counterWeddell_GLbox=counterWeddell-counterWeddell_CFbox;
  counterEastAntarctica_GLbox=counterEastAntarctica-counterEastAntarctica_CFbox;
  counterAmundsen_GLbox=counterAmundsen-counterAmundsen_CFbox;

  return 0;
}



/*!
Compute ocean temperature outside of the ice shelf cavities.
*/
PetscErrorCode POoceanboxmodel::oceanTemperature() { // FIXME unnecessary when everything is read from files?
  PetscErrorCode ierr;
  ierr = verbPrintf(2, grid.com,
                    "NOW in ocean temp rountine\n"); CHKERRQ(ierr);
  const PetscScalar rhoi = config.get("ice_density");
  const PetscScalar rhow = config.get("sea_water_density");

  ierr = ice_thickness->begin_access();   CHKERRQ(ierr);
  ierr = topg->begin_access();   CHKERRQ(ierr);
  ierr = lat->begin_access();   CHKERRQ(ierr);
  ierr = lon->begin_access();   CHKERRQ(ierr);
  ierr = SHELFmask.begin_access();   CHKERRQ(ierr);
  ierr = Soc_base.begin_access();   CHKERRQ(ierr);
  ierr = Toc_base.begin_access();   CHKERRQ(ierr);
  ierr = Toc_anomaly.begin_access();   CHKERRQ(ierr);
  ierr = Toc.begin_access();   CHKERRQ(ierr);

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {

      // make sure all temperatures are zero at the beginning of each timestep
      Toc(i,j) = 273.15; // in K
      Toc_base(i,j) = 273.15;  // in K
      Toc_anomaly(i,j) = 0.0;  // in K or °C
      Soc_base(i,j) = 0.0; // in psu

      const PetscScalar shelfbaseelev = - (rhoi / rhow) * (*ice_thickness)(i,j);
      const PetscScalar hgrounded = (*topg)(i,j) + (*ice_thickness)(i,j);
      const PetscScalar hfloating = (1.0 - rhoi/rhow) * (*ice_thickness)(i,j); // FIXME This assumes currentSeaLevel=0 !
      bool herefloating=false;
      if ((*ice_thickness)(i,j) > 0.0 && hgrounded < hfloating) {herefloating = true;}

      if (herefloating){
	  // olbers_hellmer10, table 4: observed temperatures of the northern reservoirs = T0 := (T2min+T2max)
	  // from these, subtract the temperature difference between year 2010 and the pre-industrial temperature (data from RCP scenario)
	  if (SHELFmask(i,j)==shelf_RossSea){
// 	    Toc_base(i,j) = 273.15 -1.7-1.0457091e-01;
	    Toc_base(i,j) = 273.15 -1.8475;  // FIXME config! NO! READ FROM FILE!
	    Soc_base(i,j) = 34.83;
	  }else if (SHELFmask(i,j)==shelf_WeddellSea){
	    // Here, weighted sum, i.e. T0 = (area(Filchner)/area(Weddell))*((-2-1.95)/2) + (area(Ronne)/area(Weddell))*((-1.9-1.5)/2)
// 	    Toc_base(i,j) = 273.15 -1.76-1.7764463e-01;
	    Toc_base(i,j) = 273.15 -1.8464;
	    Soc_base(i,j) = 34.74;
	  }else if (SHELFmask(i,j)==shelf_EastAntarctica){
// 	    Toc_base(i,j) = 273.15 -1.8-1.5475473e-01;
	    Toc_base(i,j) = 273.15 -1.8344;
	    Soc_base(i,j) = 34.55;
	  }else if (SHELFmask(i,j)==shelf_AmundsenSea){  // NOTE Pine Island northern reservoir has positive temperature!
// 	    Toc_base(i,j) = 273.15 +1.1-1.1295637e-01;
	    Toc_base(i,j) = 273.15 +0.8427;
	    Soc_base(i,j) = 34.67;
	  }else{
	    Toc_base(i,j) = 273.15 -1.9;
	    Soc_base(i,j) = 34.67;
	  }

	// Add temperature anomalies from given nc-file  // FIXME different nc-files for each basin!
	if ((delta_T != NULL) && ocean_oceanboxmodel_deltaT_set) {
	    Toc_anomaly(i,j) = (*delta_T)(t + 0.5 * dt);
	}else{
	  Toc_anomaly(i,j) = 0.0;
	}

	Toc(i,j) = Toc_base(i,j) + Toc_anomaly(i,j); // in K

      } // end if herefloating
    } // end j
  } // end i

  ierr = ice_thickness->end_access();   CHKERRQ(ierr);
  ierr = topg->end_access();   CHKERRQ(ierr);
  ierr = lat->end_access();   CHKERRQ(ierr);
  ierr = lon->end_access();   CHKERRQ(ierr);
  ierr = SHELFmask.end_access();   CHKERRQ(ierr);
  ierr = Soc_base.end_access();   CHKERRQ(ierr);
  ierr = Toc_base.end_access();   CHKERRQ(ierr);
  ierr = Toc_anomaly.end_access();   CHKERRQ(ierr);
  ierr = Toc.end_access();   CHKERRQ(ierr);

  return 0;
}



// NOTE This will compute the mean over ALL shelves in one region! That's why I need to reintroduce the "compute melt rate in other shelves!
// FIXME exclude ice rises!

//! Compute the basal melt / refreezing rates for each shelf cell bordering the grounding line box
PetscErrorCode POoceanboxmodel::basalMeltRateForGroundingLineBox() {
  PetscErrorCode ierr;
  ierr = verbPrintf(2, grid.com,
                    "NOW in basal melt rate gl rountine\n"); CHKERRQ(ierr);
  //! constants // FIXME config
  const PetscScalar earth_grav = config.get("standard_gravity");
  const PetscScalar rhoi = config.get("ice_density");
  const PetscScalar latentHeat = config.get("water_latent_heat_fusion");
  const PetscScalar c_p_ocean = 3974.0;       // J/(K*kg), specific heat capacity of ocean mixed layer
  const PetscScalar rho_star  = 1033;         // kg/m^3
  const PetscScalar a         = -0.057;       // °C/psu
  const PetscScalar b         = 0.0832;       // °C
  const PetscScalar c         = 7.64e-4;      // °C/dbar
  const PetscScalar alpha     = 7.5e-5;       // 1/°C, NOTE K vs °C
  const PetscScalar beta      = 7.7e-4;       // 1/psu
  const PetscScalar nu        = rhoi / rho_star;       // no unit
  const PetscScalar lambda    = latentHeat / c_p_ocean;   // °C, NOTE K vs °C

  //! parameters // FIXME config
//   const PetscScalar gamma_T_star=0.7*3.5e-5;
  const PetscScalar gamma_T_star_Ross=4.5080e-07;  // FIXME use one standard value for all shelves!
  const PetscScalar gamma_T_star_Weddell=2.4*5.8555e-07;
  const PetscScalar gamma_T_star_EastAntarctica=1.1*1.0265e-06;
  const PetscScalar gamma_T_star_Amundsen=1.6*1.1134e-05;
  const PetscScalar C_Ross = 10e6;           // m^6/(kg*s), advection parameter
  const PetscScalar C_Weddell = 2e6;   // m^6/(kg*s), advection parameter
  const PetscScalar C_EastAntarctica = 8e6;           // m^6/(kg*s), advection parameter
  const PetscScalar C_Amundsen = 3e6;      // m^6/(kg*s), advection parameter

  PetscScalar lcounter_edge_of_GLbox_Ross=0.0;   // only a local counter
  PetscScalar lcounter_edge_of_GLbox_Weddell=0.0;
  PetscScalar lcounter_edge_of_GLbox_EastAntarctica=0.0;
  PetscScalar lcounter_edge_of_GLbox_Amundsen=0.0;

  PetscScalar lmean_salinity_Ross_GLbox=0.0, lmean_salinity_Weddell_GLbox=0.0, lmean_salinity_EastAntarctica_GLbox=0.0, lmean_salinity_Amundsen_GLbox=0.0;
  PetscScalar lmean_meltrate_Ross_GLbox=0.0, lmean_meltrate_Weddell_GLbox=0.0, lmean_meltrate_EastAntarctica_GLbox=0.0, lmean_meltrate_Amundsen_GLbox=0.0;
  PetscScalar lmean_overturning_Ross_GLbox=0.0, lmean_overturning_Weddell_GLbox=0.0, lmean_overturning_EastAntarctica_GLbox=0.0, lmean_overturning_Amundsen_GLbox=0.0;

  // TODO: numberOfDrainageBasins = 17 (header) ; fehlermeldung, wenn max(nc-basins) != numberOfDrainageBasins; abbrechen
  // vektoren der länge 17 für mean salinity,...

  ierr = ice_thickness->begin_access(); CHKERRQ(ierr);
  ierr = SHELFmask.begin_access(); CHKERRQ(ierr);
  ierr = BOXMODELmask.begin_access(); CHKERRQ(ierr);
  ierr = T_star.begin_access(); CHKERRQ(ierr);
  ierr = Toc_base.begin_access(); CHKERRQ(ierr);
  ierr = Toc_anomaly.begin_access(); CHKERRQ(ierr);
  ierr = Toc_inCelsius.begin_access(); CHKERRQ(ierr);
  ierr = Toc.begin_access(); CHKERRQ(ierr);
  ierr = Soc_base.begin_access(); CHKERRQ(ierr);
  ierr = Soc.begin_access(); CHKERRQ(ierr);
  ierr = overturning.begin_access(); CHKERRQ(ierr);
  ierr = basalmeltrate_shelf.begin_access(); CHKERRQ(ierr);

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {

      // Make sure everything is at default values at the beginning of each timestep
      // Toc_base,Toc_anomaly,Toc,Soc_base are set to default in computeShelfBaseTemperature()
      T_star(i,j) = 0.0; // in °C
      Toc_inCelsius(i,j) = 0.0; // in °C
      Soc(i,j) = 0.0; // in psu

      if (BOXMODELmask(i,j) == box_GL &&
	  (SHELFmask(i,j) == shelf_RossSea || SHELFmask(i,j) == shelf_WeddellSea ||
	   SHELFmask(i,j) == shelf_EastAntarctica || SHELFmask(i,j) == shelf_AmundsenSea)){ // this implies floating
	  const PetscScalar pressure = rhoi * earth_grav * (*ice_thickness)(i,j) * 1e-4; // MUST be in dbar  // NOTE 1dbar = 10000 Pa = 1e4 kg m-1 s-2
									       // FIXME need to include atmospheric pressure?
	  T_star(i,j) = a*Soc_base(i,j) + b - c*pressure - (Toc_base(i,j)-273.15+Toc_anomaly(i,j)); // in °C

	  PetscScalar gamma_T_star,C1,g1;
	  if (SHELFmask(i,j) == shelf_RossSea){
	    gamma_T_star=gamma_T_star_Ross;
	    C1 = C_Ross;
	    g1 = (counterRoss_GLbox * grid.dx * grid.dy)*gamma_T_star / (C_Ross*rho_star);   // in (m^2*m/s) / ((m^6/(kg*s)) * kg/m^3) = 1
	  } else if (SHELFmask(i,j) == shelf_WeddellSea){
	    gamma_T_star=gamma_T_star_Weddell;
	    C1 = C_Weddell;
	    g1 = (counterWeddell_GLbox * grid.dx * grid.dy)*gamma_T_star / (C_Weddell*rho_star); // in (m^2*m/s) / ((m^6/(kg*s)) * kg/m^3) = 1
	  } else if (SHELFmask(i,j) == shelf_EastAntarctica){
	    gamma_T_star=gamma_T_star_EastAntarctica;
	    C1 = C_EastAntarctica;
	    g1 = (counterEastAntarctica_GLbox * grid.dx * grid.dy)*gamma_T_star / (C_EastAntarctica*rho_star);   // in (m^2*m/s) / ((m^6/(kg*s)) * kg/m^3) = 1
	  } else if (SHELFmask(i,j) == shelf_AmundsenSea){
	    gamma_T_star=gamma_T_star_Amundsen;
	    C1 = C_Amundsen;
	    g1 = (counterAmundsen_GLbox * grid.dx * grid.dy)*gamma_T_star / (C_Amundsen*rho_star);   // in (m^2*m/s) / ((m^6/(kg*s)) * kg/m^3) = 1
	  } else { // This must not happen, since SHELFmask needs to be either RIS, FRIS, AIS or PIG at this point
	    ierr = verbPrintf(2, grid.com,"PISM_ERROR: [rank %d] at %d, %d  -- SHELFmask(i,j)=%f causes problems.\n   Aborting... \n",grid.rank, i, j, SHELFmask(i,j)); CHKERRQ(ierr);
	    PISMEnd();
	  }

          //! temperature for grounding line box
          PetscScalar helpterm1, helpterm2;
	  helpterm1 = g1/(beta*(Soc_base(i,j) / (nu*lambda)) - alpha);                  // in 1 / (1/°C) = °C
	  helpterm2 = (g1*T_star(i,j)) / (beta*(Soc_base(i,j) / (nu*lambda)) - alpha); // in °C / (1/°C) = °C^2
	  if ((0.25*PetscSqr(helpterm1) -helpterm2) < 0.0) {
	    ierr = verbPrintf(2, grid.com,"PISM_ERROR: square-root is negative! %f at %d, %d\n...with 0.25*helpterm^2=%f,helpterm2=%f,g1=%f,(beta*(Soc_base(i,j) / (nu*lambda)) - alpha)=%f,Tstar=%f\n   Aborting... \n", 0.25*PetscSqr(helpterm1) -helpterm2, i, j, 0.25*PetscSqr(helpterm1),helpterm2,g1,(beta*(Soc_base(i,j) / (nu*lambda)) - alpha),T_star(i,j)); CHKERRQ(ierr);
	    PISMEnd();
	  }

	  // NOTE Careful, Toc_base(i,j) is in K, Toc_inCelsius(i,j) NEEDS to be in °C!
	  Toc_inCelsius(i,j) = (Toc_base(i,j)-273.15+Toc_anomaly(i,j)) - ( -0.5*helpterm1 + sqrt(0.25*PetscSqr(helpterm1) -helpterm2) );

	  //! salinity for grounding line box
	  Soc(i,j) = Soc_base(i,j) - (Soc_base(i,j) / (nu*lambda)) * ((Toc_base(i,j)-273.15+Toc_anomaly(i,j)) - Toc_inCelsius(i,j));  // in psu

	  //! basal melt rate for grounding line box
	  basalmeltrate_shelf(i,j) = (-gamma_T_star/(nu*lambda)) * (a*Soc(i,j) + b - c*pressure - Toc_inCelsius(i,j));  // in m/s

	  //! overturning
	  // FIXME NOTE Actually, there is of course no overturning-FIELD, it is only a scalar for each shelf.
	  // Here, I compute overturning as 	MEAN[C1*rho_star* (beta*(Soc_base(i,j)-Soc(i,j)) - alpha*((Toc_base(i,j)-273.15+Toc_anomaly(i,j))-Toc_inCelsius(i,j)))]
	  // while in fact it should be 	C1*rho_star* (beta*(Soc_base-MEAN[Soc(i,j)]) - alpha*((Toc_base-273.15+Toc_anomaly)-MEAN[Toc_inCelsius(i,j)]))
	  // which is the SAME since Soc_base, Toc_base and Toc_anomaly are the same FOR ALL i,j CONSIDERED, so this is just nomenclature!
	  overturning(i,j) = C1*rho_star* (beta*(Soc_base(i,j)-Soc(i,j)) - alpha*((Toc_base(i,j)-273.15+Toc_anomaly(i,j))-Toc_inCelsius(i,j))); // in m^3/s


	  // TODO loop über alle basins
	  if (BOXMODELmask(i-1,j)==box_IF ||
	    BOXMODELmask(i+1,j)==box_IF ||
	    BOXMODELmask(i,j-1)==box_IF ||
	    BOXMODELmask(i,j+1)==box_IF){ // i.e., if this cell is from the GL box and one of the neighbours is from the CF box - It is important to only take the border of the grounding line box to the calving front box into account, because the following mean value will be used to compute the value for the calving front box. I.e., this helps avoiding discontinuities!
	      if (SHELFmask(i,j) == shelf_RossSea){
		lcounter_edge_of_GLbox_Ross++;
		lmean_salinity_Ross_GLbox += Soc(i,j);
		lmean_meltrate_Ross_GLbox += basalmeltrate_shelf(i,j);
		lmean_overturning_Ross_GLbox += overturning(i,j);
	      } else if (SHELFmask(i,j) == shelf_WeddellSea){
		lcounter_edge_of_GLbox_Weddell++;
		lmean_salinity_Weddell_GLbox += Soc(i,j);
		lmean_meltrate_Weddell_GLbox += basalmeltrate_shelf(i,j);
		lmean_overturning_Weddell_GLbox += overturning(i,j);
	      } else if (SHELFmask(i,j) == shelf_EastAntarctica){
		lcounter_edge_of_GLbox_EastAntarctica++;
		lmean_salinity_EastAntarctica_GLbox += Soc(i,j);
		lmean_meltrate_EastAntarctica_GLbox += basalmeltrate_shelf(i,j);
		lmean_overturning_EastAntarctica_GLbox += overturning(i,j);
	      } else if (SHELFmask(i,j) == shelf_AmundsenSea){
		lcounter_edge_of_GLbox_Amundsen++;
		lmean_salinity_Amundsen_GLbox += Soc(i,j);
		lmean_meltrate_Amundsen_GLbox += basalmeltrate_shelf(i,j);
		lmean_overturning_Amundsen_GLbox += overturning(i,j);
	      } else { // This must not happen, since SHELFmask needs to be either RIS, FRIS, AIS or PIG at this point
		ierr = verbPrintf(2, grid.com,"PISM_ERROR: [rank %d] at %d, %d  -- SHELFmask(i,j)=%f causes problems.\n   Aborting... \n",grid.rank, i, j, SHELFmask(i,j)); CHKERRQ(ierr);
		PISMEnd();
	      }
	  } // no else-case necessary since all variables are set to zero at the beginning of this routine

	  // 	ierr = verbPrintf(5, grid.com, "ERROR DETECTION GLBOX: SHELFmask=%f, BOXMODELmask=%f, H=%f, T_star=%f, C1=%e, g1=%e, helpterm1=%e, helpterm2=%e, Toc_base=%f, Toc_anomaly=%f, Toc_inCelsius=%f, Toc=%f, Soc_base=%f, Soc=%f, overturning=%e, basalmeltrate_shelf=%e\n", SHELFmask(i,j), BOXMODELmask(i,j), (*ice_thickness)(i,j), T_star(i,j), C1,g1, helpterm1, helpterm2, Toc_base(i,j), Toc_anomaly(i,j), Toc_inCelsius(i,j), Toc(i,j), Soc_base(i,j), Soc(i,j), overturning(i,j), basalmeltrate_shelf(i,j)); CHKERRQ(ierr);

// 	  if (grid.year >= 1916.00) {
// 	  ierr = verbPrintf(2, grid.com,"PISM_INFO: [rank %d] at %d, %d  -- T_star(i,j)=%e,C1=%e,g1=%e,helpterm=%e1,helpterm2=%e,Toc_inCelsius(i,j)=%e,Soc(i,j)=%e,basalmeltrate_shelf(i,j)=%e,overturning(i,j)=%e\n",grid.rank, i, j, T_star(i,j),C1,g1,helpterm1, helpterm2,Toc_inCelsius(i,j),Soc(i,j),basalmeltrate_shelf(i,j),overturning(i,j)); CHKERRQ(ierr);
// 	  }
      }else { // i.e., neither Ross nor Weddell nor EastAntarctica nor Pine Island
	basalmeltrate_shelf(i,j) = 0.0;
      }
   } // end j
  } // end i

  ierr = ice_thickness->end_access(); CHKERRQ(ierr);
  ierr = SHELFmask.end_access(); CHKERRQ(ierr);
  ierr = BOXMODELmask.end_access(); CHKERRQ(ierr);
  ierr = T_star.end_access(); CHKERRQ(ierr);
  ierr = Toc_base.end_access(); CHKERRQ(ierr);
  ierr = Toc_anomaly.end_access(); CHKERRQ(ierr);
  ierr = Toc_inCelsius.end_access(); CHKERRQ(ierr);
  ierr = Toc.end_access(); CHKERRQ(ierr);
  ierr = Soc_base.end_access(); CHKERRQ(ierr);
  ierr = Soc.end_access(); CHKERRQ(ierr);
  ierr = overturning.end_access(); CHKERRQ(ierr);
  ierr = basalmeltrate_shelf.end_access(); CHKERRQ(ierr);


  PetscScalar counter_edge_of_GLbox_Ross, counter_edge_of_GLbox_Weddell, counter_edge_of_GLbox_EastAntarctica, counter_edge_of_GLbox_Amundsen;
  ierr = PISMGlobalSum(&lcounter_edge_of_GLbox_Ross, &counter_edge_of_GLbox_Ross, grid.com); CHKERRQ(ierr);
  ierr = PISMGlobalSum(&lcounter_edge_of_GLbox_Weddell, &counter_edge_of_GLbox_Weddell, grid.com); CHKERRQ(ierr);
  ierr = PISMGlobalSum(&lcounter_edge_of_GLbox_EastAntarctica, &counter_edge_of_GLbox_EastAntarctica, grid.com); CHKERRQ(ierr);
  ierr = PISMGlobalSum(&lcounter_edge_of_GLbox_Amundsen, &counter_edge_of_GLbox_Amundsen, grid.com); CHKERRQ(ierr);

  ierr = PISMGlobalSum(&lmean_meltrate_Ross_GLbox, &mean_meltrate_Ross_GLbox, grid.com); CHKERRQ(ierr);
  ierr = PISMGlobalSum(&lmean_meltrate_Weddell_GLbox, &mean_meltrate_Weddell_GLbox, grid.com); CHKERRQ(ierr);
  ierr = PISMGlobalSum(&lmean_meltrate_EastAntarctica_GLbox, &mean_meltrate_EastAntarctica_GLbox, grid.com); CHKERRQ(ierr);
  ierr = PISMGlobalSum(&lmean_meltrate_Amundsen_GLbox, &mean_meltrate_Amundsen_GLbox, grid.com); CHKERRQ(ierr);

  ierr = PISMGlobalSum(&lmean_salinity_Ross_GLbox, &mean_salinity_Ross_GLbox, grid.com); CHKERRQ(ierr);
  ierr = PISMGlobalSum(&lmean_salinity_Weddell_GLbox, &mean_salinity_Weddell_GLbox, grid.com); CHKERRQ(ierr);
  ierr = PISMGlobalSum(&lmean_salinity_EastAntarctica_GLbox, &mean_salinity_EastAntarctica_GLbox, grid.com); CHKERRQ(ierr);
  ierr = PISMGlobalSum(&lmean_salinity_Amundsen_GLbox, &mean_salinity_Amundsen_GLbox, grid.com); CHKERRQ(ierr);

  ierr = PISMGlobalSum(&lmean_overturning_Ross_GLbox, &mean_overturning_Ross_GLbox, grid.com); CHKERRQ(ierr);
  ierr = PISMGlobalSum(&lmean_overturning_Weddell_GLbox, &mean_overturning_Weddell_GLbox, grid.com); CHKERRQ(ierr);
  ierr = PISMGlobalSum(&lmean_overturning_EastAntarctica_GLbox, &mean_overturning_EastAntarctica_GLbox, grid.com); CHKERRQ(ierr);
  ierr = PISMGlobalSum(&lmean_overturning_Amundsen_GLbox, &mean_overturning_Amundsen_GLbox, grid.com); CHKERRQ(ierr);

  if (counter_edge_of_GLbox_Ross>0.0){
    mean_salinity_Ross_GLbox=mean_salinity_Ross_GLbox/counter_edge_of_GLbox_Ross;
    mean_meltrate_Ross_GLbox=mean_meltrate_Ross_GLbox/counter_edge_of_GLbox_Ross;
    mean_overturning_Ross_GLbox=mean_overturning_Ross_GLbox/counter_edge_of_GLbox_Ross;
  } else { // This means that there is no [cell from the GLbox neighboring a cell from the CFbox], NOT necessarily that there is no GLbox!
    mean_salinity_Ross_GLbox=0.0; mean_meltrate_Ross_GLbox=0.0; mean_overturning_Ross_GLbox=0.0;
  }
  if (counter_edge_of_GLbox_Weddell>0.0){
    mean_salinity_Weddell_GLbox=mean_salinity_Weddell_GLbox/counter_edge_of_GLbox_Weddell;
    mean_meltrate_Weddell_GLbox=mean_meltrate_Weddell_GLbox/counter_edge_of_GLbox_Weddell;
    mean_overturning_Weddell_GLbox=mean_overturning_Weddell_GLbox/counter_edge_of_GLbox_Weddell;
  } else { // This means that there is no [cell from the GLbox neighboring a cell from the CFbox], NOT necessarily that there is no GLbox!
    mean_salinity_Weddell_GLbox=0.0; mean_meltrate_Weddell_GLbox=0.0; mean_overturning_Weddell_GLbox=0.0;
  }
  if (counter_edge_of_GLbox_EastAntarctica>0.0){
    mean_salinity_EastAntarctica_GLbox=mean_salinity_EastAntarctica_GLbox/counter_edge_of_GLbox_EastAntarctica;
    mean_meltrate_EastAntarctica_GLbox=mean_meltrate_EastAntarctica_GLbox/counter_edge_of_GLbox_EastAntarctica;
    mean_overturning_EastAntarctica_GLbox=mean_overturning_EastAntarctica_GLbox/counter_edge_of_GLbox_EastAntarctica;
  } else { // This means that there is no [cell from the GLbox neighboring a cell from the CFbox], NOT necessarily that there is no GLbox!
    mean_salinity_EastAntarctica_GLbox=0.0; mean_meltrate_EastAntarctica_GLbox=0.0; mean_overturning_EastAntarctica_GLbox=0.0;
  }
  if (counter_edge_of_GLbox_Amundsen>0.0){
    mean_salinity_Amundsen_GLbox=mean_salinity_Amundsen_GLbox/counter_edge_of_GLbox_Amundsen;
    mean_meltrate_Amundsen_GLbox=mean_meltrate_Amundsen_GLbox/counter_edge_of_GLbox_Amundsen;
    mean_overturning_Amundsen_GLbox=mean_overturning_Amundsen_GLbox/counter_edge_of_GLbox_Amundsen;
  } else { // This means that there is no [cell from the GLbox neighboring a cell from the CFbox], NOT necessarily that there is no GLbox!
    mean_salinity_Amundsen_GLbox=0.0; mean_meltrate_Amundsen_GLbox=0.0; mean_overturning_Amundsen_GLbox=0.0;
  }

//   ierr = verbPrintf(5, grid.com, "          GLbox to CFbox: mean melt rate for Ross = %f m/a, Weddell = %f m/a, EastAntarctica = %f m/a, PIG = %f m/a \n                 overturning for Ross = %f Sv, Weddell = %f Sv, EastAntarctica = %f Sv, PIG = %f Sv \n",
//   mean_meltrate_Ross_GLbox*secpera,mean_meltrate_Weddell_GLbox*secpera,mean_meltrate_EastAntarctica_GLbox*secpera,mean_meltrate_Amundsen_GLbox*secpera,
//   mean_overturning_Ross_GLbox*1e-6,mean_overturning_Weddell_GLbox*1e-6,mean_overturning_EastAntarctica_GLbox*1e-6,mean_overturning_Amundsen_GLbox*1e-6); CHKERRQ(ierr);
  return 0;
}

//! Compute the basal melt / refreezing rates for each shelf cell bordering the ice front box
PetscErrorCode POoceanboxmodel::basalMeltRateForIceFrontBox() {
  PetscErrorCode ierr;  // FIXME redo all verbprintfs!
  ierr = verbPrintf(2, grid.com,
                    "NOW in bm ice front rountine\n"); CHKERRQ(ierr);
  //! constants
  const PetscScalar earth_grav = config.get("standard_gravity");
  const PetscScalar rhoi = config.get("ice_density");
  const PetscScalar latentHeat = config.get("water_latent_heat_fusion");
  const PetscScalar c_p_ocean	 = 3974.0;       // J/(K*kg), specific heat capacity of ocean mixed layer
  const PetscScalar rho_star = 1033;         // kg/m^3
  const PetscScalar a        = -0.057;       // °C/psu
  const PetscScalar b        = 0.0832;       // °C
  const PetscScalar c        = 7.64e-4;      // °C/dbar
  const PetscScalar nu       = rhoi / rho_star;       // no unit
  const PetscScalar lambda   = latentHeat / c_p_ocean;   // °C, NOTE K vs °C

  //! parameters
//   const PetscScalar gamma_T_star=0.7*3.5e-5;
//   const PetscScalar gamma_T_star_Ross=4.5080e-07;
//   const PetscScalar gamma_T_star_Weddell=5.8555e-07;
//   const PetscScalar gamma_T_star_EastAntarctica=1.0265e-06;
//   const PetscScalar gamma_T_star_Amundsen=1.1134e-05;
  const PetscScalar gamma_T_star_Ross=4.5080e-07;
  const PetscScalar gamma_T_star_Weddell=2.4*5.8555e-07;
  const PetscScalar gamma_T_star_EastAntarctica=1.1*1.0265e-06;
  const PetscScalar gamma_T_star_Amundsen=1.6*1.1134e-05;

  ierr = SHELFmask.begin_access(); CHKERRQ(ierr);
  ierr = BOXMODELmask.begin_access(); CHKERRQ(ierr);
  ierr = ice_thickness->begin_access(); CHKERRQ(ierr);
  ierr = T_star.begin_access(); CHKERRQ(ierr);
  ierr = Toc_base.begin_access(); CHKERRQ(ierr);
  ierr = Toc_anomaly.begin_access(); CHKERRQ(ierr);
  ierr = Toc_inCelsius.begin_access(); CHKERRQ(ierr);
  ierr = Toc.begin_access(); CHKERRQ(ierr);
  ierr = Soc_base.begin_access(); CHKERRQ(ierr);
  ierr = Soc.begin_access(); CHKERRQ(ierr);
  ierr = overturning.begin_access(); CHKERRQ(ierr);
  ierr = basalmeltrate_shelf.begin_access(); CHKERRQ(ierr);


  //! The ice front box = BOX I
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {  // FIXME REPAIR
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {

      if (BOXMODELmask(i,j)==box_IF &&
	  (SHELFmask(i,j) == shelf_RossSea || SHELFmask(i,j) == shelf_WeddellSea ||
	   SHELFmask(i,j) == shelf_EastAntarctica || SHELFmask(i,j) == shelf_AmundsenSea)){ // this implies floating
	  const PetscScalar pressure = rhoi * earth_grav * (*ice_thickness)(i,j) * 1e-4; // MUST be in dbar  // NOTE 1dbar = 10000 Pa = 1e4 kg m-1 s-2
									       // FIXME need to include atmospheric pressure?
	  T_star(i,j) = a*Soc_base(i,j) + b - c*pressure - (Toc_base(i,j)-273.15+Toc_anomaly(i,j)); // in °C


	  PetscScalar  gamma_T_star,area_GLbox,area_CFbox,mean_salinity_in_GLbox,mean_meltrate_in_GLbox,mean_overturning_in_GLbox;
	  if (SHELFmask(i,j) == shelf_RossSea){
	    gamma_T_star=gamma_T_star_Ross;
	    area_CFbox = (counterRoss_CFbox * grid.dx * grid.dy);
	    area_GLbox = (counterRoss_GLbox * grid.dx * grid.dy);
	    mean_salinity_in_GLbox = mean_salinity_Ross_GLbox;
	    mean_meltrate_in_GLbox = mean_meltrate_Ross_GLbox;
	    mean_overturning_in_GLbox = mean_overturning_Ross_GLbox;
	  } else if (SHELFmask(i,j) == shelf_WeddellSea){
	    gamma_T_star=gamma_T_star_Weddell;
	    area_CFbox = (counterWeddell_CFbox * grid.dx * grid.dy);
	    area_GLbox = (counterWeddell_GLbox * grid.dx * grid.dy);
	    mean_salinity_in_GLbox = mean_salinity_Weddell_GLbox;
	    mean_meltrate_in_GLbox = mean_meltrate_Weddell_GLbox;
	    mean_overturning_in_GLbox = mean_overturning_Weddell_GLbox;
	  } else if (SHELFmask(i,j) == shelf_EastAntarctica){
	    gamma_T_star=gamma_T_star_EastAntarctica;
	    area_CFbox = (counterEastAntarctica_CFbox * grid.dx * grid.dy);
	    area_GLbox = (counterEastAntarctica_GLbox * grid.dx * grid.dy);
	    mean_salinity_in_GLbox = mean_salinity_EastAntarctica_GLbox;
	    mean_meltrate_in_GLbox = mean_meltrate_EastAntarctica_GLbox;
	    mean_overturning_in_GLbox = mean_overturning_EastAntarctica_GLbox;
	  } else if (SHELFmask(i,j) == shelf_AmundsenSea){
	    gamma_T_star=gamma_T_star_Amundsen;
	    area_CFbox = (counterAmundsen_CFbox * grid.dx * grid.dy);
	    area_GLbox = (counterAmundsen_GLbox * grid.dx * grid.dy);
	    mean_salinity_in_GLbox = mean_salinity_Amundsen_GLbox;
	    mean_meltrate_in_GLbox = mean_meltrate_Amundsen_GLbox;
	    mean_overturning_in_GLbox = mean_overturning_Amundsen_GLbox;
	  } else { // This must not happen, since SHELFmask needs to be either RIS, FRIS, AIS or PIG at this point
		ierr = verbPrintf(2, grid.com,"PISM_ERROR: SHELFmask(i,j)=%f is wrong! [rank %d] at %d, %d \n   Aborting... \n",grid.rank, i, j, SHELFmask(i,j)); CHKERRQ(ierr);
		PISMEnd();
	  }

          PetscScalar k1,k2,k3,k4,k5,k6;
	  k1 = (area_CFbox*gamma_T_star)/(nu*lambda);                                                 // in (m^2*m/s)/(°C) = m^3/(s*°C)
	  k2 = 1/(mean_overturning_in_GLbox + area_CFbox*gamma_T_star);                               // in s/m^3
	  k3 = (area_CFbox*gamma_T_star*T_star(i,j) - nu*lambda*area_GLbox*mean_meltrate_in_GLbox);  // in m^2 * m/s * °C = m^3 * °C / s
          k4 = (-k1*k2*area_CFbox*gamma_T_star*a + k1*a);              // in m^3/(s*°C) * s/m^3 * m^2 * m/s * °C/psu = m^3/(s*psu)
          k5 = (mean_overturning_in_GLbox + Soc_base(i,j)*k1*k2*area_CFbox*gamma_T_star*a  - k1*Soc_base(i,j)*a - k1*T_star(i,j) + k1*k2*k3); // in m^3/s
          k6 = (k1*Soc_base(i,j)*T_star(i,j) - k1*k2*Soc_base(i,j)*k3  - area_GLbox*mean_meltrate_in_GLbox*mean_salinity_in_GLbox);  // in psu * m^3/s

          //! salinity for calving front box
	  if (k4 == 0.0) {
	    ierr = verbPrintf(2, grid.com,"PISM_ERROR: Division by zero! k4=%f at %d, %d\n   Aborting... \n", k4, i, j); CHKERRQ(ierr);
	    PISMEnd();
	  }
	  if ((0.25*PetscSqr(k5/k4) - (k6/k4)) < 0.0) {
	    ierr = verbPrintf(2, grid.com,"PISM_ERROR: Square-root is negative! %f at %d, %d\n...with 0.25*PetscSqr((k5/k4))=%f,(k6/k4)=%f,k5=%f,k6=%f,k4=%f\n   Aborting... \n", 0.25*PetscSqr(k5/k4) - (k6/k4), i, j,0.25*PetscSqr((k5/k4)), (k6/k4),k5,k6,k4) ; CHKERRQ(ierr);
	    PISMEnd();
	  }
          Soc(i,j) = Soc_base(i,j) - ( -0.5*(k5/k4) - sqrt(0.25*PetscSqr(k5/k4) - (k6/k4)) ); // in psu

	  //! temperature for calving front box
	  // NOTE Careful, Toc_base(i,j) is in K, Toc_inCelsius(i,j) NEEDS to be in °C!
	  Toc_inCelsius(i,j) = (Toc_base(i,j)-273.15+Toc_anomaly(i,j)) - ( k2*area_CFbox*gamma_T_star*a*(Soc_base(i,j)-Soc(i,j)) - k2*k3 );

	  //! basal melt rate for calving front box
	  basalmeltrate_shelf(i,j) = (-gamma_T_star/(nu*lambda)) * (a*Soc(i,j) + b - c*pressure - Toc_inCelsius(i,j)); // in m/s

	  if (mean_salinity_in_GLbox == 0.0 || mean_meltrate_in_GLbox == 0.0 || mean_overturning_in_GLbox == 0.0){
	    // Basically, this must not occur in a cell which belongs to the CFbox of RIS, FRIS, AIS or PIG - so the run is terminated.
	    ierr = verbPrintf(2, grid.com, "PISM_ERROR: DETECTION CFBOX: There is no neighbouring grounding line box for this calving front box at %d,%d! \nThis will lead to a zero k4 and in turn to NaN in Soc, Toc_inCelsius and basalmeltrate_shelf. After the next massContExplicitStep(), H will be NaN, too! This will cause ks in temperatureStep() to be NaN and lead to a Segmentation Violation! \nIn particular: SHELFmask=%f, BOXMODELmask=%f, H=%f, T_star=%f, \narea_GLbox=%e, area_CFbox=%e, mean_salinity_in_GLbox=%f, mean_meltrate_in_GLbox=%e, mean_overturning_in_GLbox=%e, \nk1=%e,k2=%e,k3=%e,k4=%e,k5=%e,k6=%e, \nToc_base=%f, Toc_anomaly=%f, Toc_inCelsius=%f, Toc=%f, Soc_base=%f, Soc=%f, basalmeltrate_shelf=%e \n   Aborting... \n", i,j, SHELFmask(i,j), BOXMODELmask(i,j), (*ice_thickness)(i,j), T_star(i,j), area_GLbox,area_CFbox,mean_salinity_in_GLbox,mean_meltrate_in_GLbox,mean_overturning_in_GLbox,k1,k2,k3,k4,k5,k6, Toc_base(i,j), Toc_anomaly(i,j), Toc_inCelsius(i,j), Toc(i,j), Soc_base(i,j), Soc(i,j), basalmeltrate_shelf(i,j)); CHKERRQ(ierr);
	    PISMEnd();
	  }

// 	    if (grid.year >= 1916.00) {
// 	    ierr = verbPrintf(2, grid.com, "PISM_INFO: [rank %d] at %d, %d  -- T_star(i,j)=%e, area_GLbox=%e, area_CFbox=%e, mean_salinity_in_GLbox=%e, mean_meltrate_in_GLbox=%e, mean_overturning_in_GLbox=%e, k1=%e,k2=%e,k3=%e,k4=%e,k5=%e,k6=%e,Soc(i,j)=%e,Toc_inCelsius(i,j)=%e,basalmeltrate_shelf(i,j)=%e\n",grid.rank, i, j, T_star(i,j), area_GLbox,area_CFbox,mean_salinity_in_GLbox,mean_meltrate_in_GLbox,mean_overturning_in_GLbox,k1,k2,k3,k4,k5,k6,Soc(i,j),Toc_inCelsius(i,j),basalmeltrate_shelf(i,j)); CHKERRQ(ierr);
// 	    }

    } // NOTE NO else-case, since  basalMeltRateForGroundingLineBox() and basalMeltRateForOtherShelves() cover all other cases and we would overwrite those results here.
   } // end j
  } // end i

  ierr = ice_thickness->end_access(); CHKERRQ(ierr);
  ierr = SHELFmask.end_access(); CHKERRQ(ierr);
  ierr = BOXMODELmask.end_access(); CHKERRQ(ierr);
  ierr = T_star.end_access(); CHKERRQ(ierr);
  ierr = Toc_base.end_access(); CHKERRQ(ierr);
  ierr = Toc_anomaly.end_access(); CHKERRQ(ierr);
  ierr = Toc_inCelsius.end_access(); CHKERRQ(ierr);
  ierr = Toc.end_access(); CHKERRQ(ierr);
  ierr = Soc_base.end_access(); CHKERRQ(ierr);
  ierr = Soc.end_access(); CHKERRQ(ierr);
  ierr = overturning.end_access(); CHKERRQ(ierr);
  ierr = basalmeltrate_shelf.end_access(); CHKERRQ(ierr);

  return 0;
}



//! Convert Toc_inCelsius from °C to K and write into Toc for the .nc-file; NOTE It is crucial, that Toc_inCelsius is in °C for the computation of the basal melt rate
//! Compute the melt rate for all other ice shelves.
PetscErrorCode POoceanboxmodel::basalMeltRateForOtherShelves() {
  PetscErrorCode ierr;
  ierr = verbPrintf(2, grid.com,
                    "NOW in bm others rountine\n"); CHKERRQ(ierr);
  const PetscScalar rhoi = config.get("ice_density");
  const PetscScalar rhow = config.get("sea_water_density");
  const PetscScalar latentHeat = config.get("water_latent_heat_fusion");
  const PetscScalar c_p_ocean	 = 3974.0;       // J/(K*kg), specific heat capacity of ocean mixed layer
  const PetscScalar gamma_T = 0.178567865873;  // FIXME!!!! (Wrong value!) FIXME config
  const PetscScalar meltFactor = 0.002; // FIXME!!!! (Wrong value!) FIXME config

  ierr = SHELFmask.begin_access(); CHKERRQ(ierr);
  ierr = ice_thickness->begin_access(); CHKERRQ(ierr);
  ierr = Toc_base.begin_access(); CHKERRQ(ierr);
  ierr = Toc_anomaly.begin_access(); CHKERRQ(ierr);
  ierr = Toc_inCelsius.begin_access(); CHKERRQ(ierr);
  ierr = Toc.begin_access(); CHKERRQ(ierr);
  ierr = overturning.begin_access(); CHKERRQ(ierr);
  ierr = basalmeltrate_shelf.begin_access(); CHKERRQ(ierr); // NOTE meltrate has units:   J m-2 s-1 / (J kg-1 * kg m-3) = m s-1


  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {

      if (SHELFmask(i,j) == shelf_RossSea || SHELFmask(i,j) == shelf_WeddellSea || SHELFmask(i,j) == shelf_EastAntarctica || SHELFmask(i,j) == shelf_AmundsenSea){
          Toc(i,j) = 273.15 + Toc_inCelsius(i,j) + Toc_anomaly(i,j); // in K  // FIXME I think Toc should not occur in any of the routines before!
      } else if (SHELFmask(i,j) == shelf_unidentified){
	Toc(i,j) = Toc_base(i,j) + Toc_anomaly(i,j); // in K, NOTE: Toc_base is already in K, so no (+273.15)
	// default: compute the melt rate from the temperature field according to beckmann_goosse03 (see below)
	const PetscScalar shelfbaseelev = - (rhoi / rhow) * (*ice_thickness)(i,j);
	const PetscScalar c_p_ocean = 3974.0;       // J/(K*kg), specific heat capacity of ocean mixed layer
        const PetscScalar gamma_T   = 1e-4;     // m/s, thermal exchange velocity
        PetscScalar T_f = 273.15 + (0.0939 -0.057*35.0 + 7.64e-4* shelfbaseelev); // add 273.15 to get it in Kelvin... 35 is the salinity
	heatflux(i,j) = meltFactor * rhow * c_p_ocean * gamma_T * (Toc(i,j) - T_f);  // in W/m^2
	basalmeltrate_shelf(i,j) = heatflux(i,j) / (latentHeat * rhoi); // in m s-1
      } else if (SHELFmask(i,j) == noshelf) {
	basalmeltrate_shelf(i,j) = 0.0;
      } else { // This must not happen, since SHELFmask needs to be one of the abovenoshelf
	ierr = verbPrintf(2, grid.com,"PISM_ERROR: [rank %d] at %d, %d  -- SHELFmask(i,j)=%f causes problems.\n   Aborting... \n",grid.rank, i, j, SHELFmask(i,j)); CHKERRQ(ierr);
	PISMEnd();
      }
    } // end j
  } // end i

  ierr = SHELFmask.end_access(); CHKERRQ(ierr);
  ierr = ice_thickness->end_access(); CHKERRQ(ierr);
  ierr = Toc_base.end_access(); CHKERRQ(ierr);
  ierr = Toc_anomaly.end_access(); CHKERRQ(ierr);
  ierr = Toc_inCelsius.end_access(); CHKERRQ(ierr);
  ierr = Toc.end_access(); CHKERRQ(ierr);
  ierr = overturning.end_access(); CHKERRQ(ierr);
  ierr = basalmeltrate_shelf.end_access(); CHKERRQ(ierr);

  return 0;
}



PetscErrorCode POoceanboxmodel::shelf_base_temperature(IceModelVec2S &result) {
  PetscErrorCode ierr = AntarcticBasins(); CHKERRQ(ierr);
  ierr = oceanTemperature(); CHKERRQ(ierr);
  ierr = Toc.copy_to(result); CHKERRQ(ierr);
  ierr = verbPrintf(2, grid.com,
                    "NOW in shelf_base_temperature rountine\n"); CHKERRQ(ierr);
  return 0;
}

//! \brief Computes mass flux in ice-equivalent m s-1.
/*!
 * Assumes that mass flux is proportional to the shelf-base heat flux.
 */
PetscErrorCode POoceanboxmodel::shelf_base_mass_flux(IceModelVec2S &result) {
  PetscErrorCode ierr = basalMeltRateForGroundingLineBox(); CHKERRQ(ierr);
  ierr = basalMeltRateForIceFrontBox(); CHKERRQ(ierr); // TODO Diese Routinen woanders aufrufen (um Dopplung zu vermeiden)
  ierr = basalMeltRateForOtherShelves(); CHKERRQ(ierr);
  ierr = basalmeltrate_shelf.copy_to(result); CHKERRQ(ierr);
  ierr = verbPrintf(2, grid.com,
                    "NOW in shelf_base_mass_flux rountine\n"); CHKERRQ(ierr);

  return 0;
}



//  void POoceanboxmodel::add_vars_to_output(string keyword, map<string,NCSpatialVariable> &result) {
//    if (keyword == "medium" || keyword == "big") {
//      result["shelfbtemp"] = shelfbtemp;
//      result["shelfbmassflux"] = shelfbmassflux;
//    }
//  }

//  PetscErrorCode POoceanboxmodel::define_variables(set<string> vars, const PIO &nc,
//                                                 PISM_IO_Type nctype) {
//    PetscErrorCode ierr;

//    if (set_contains(vars, "shelfbtemp")) {
//      ierr = shelfbtemp.define(nc, nctype, true); CHKERRQ(ierr);
//    }

//    if (set_contains(vars, "shelfbmassflux")) {
//      ierr = shelfbmassflux.define(nc, nctype, true); CHKERRQ(ierr);
//    }

//    return 0;
//  }

// PetscErrorCode POoceanboxmodel::write_variables(set<string> vars, string filename) {
//   PetscErrorCode ierr;
//   IceModelVec2S tmp;

//   if (set_contains(vars, "shelfbtemp")) {
//     if (!tmp.was_created()) {
//       ierr = tmp.create(grid, "tmp", false); CHKERRQ(ierr);
//     }

//     ierr = tmp.set_metadata(shelfbtemp, 0); CHKERRQ(ierr);
//     ierr = shelf_base_temperature(tmp); CHKERRQ(ierr);
//     ierr = tmp.write(filename.c_str()); CHKERRQ(ierr);
//   }

//   if (set_contains(vars, "shelfbmassflux")) {
//     if (!tmp.was_created()) {
//       ierr = tmp.create(grid, "tmp", false); CHKERRQ(ierr);
//     }

//     ierr = tmp.set_metadata(shelfbmassflux, 0); CHKERRQ(ierr);
//     tmp.write_in_glaciological_units = true;
//     ierr = shelf_base_mass_flux(tmp); CHKERRQ(ierr);
//     ierr = tmp.write(filename.c_str()); CHKERRQ(ierr);
//   }

//   if (set_contains(vars, "BOXMODELmask")) {  ierr = BOXMODELmask.write(filename.c_str()); CHKERRQ(ierr);  }
//   if (set_contains(vars, "SHELFmask")) {  ierr = SHELFmask.write(filename.c_str()); CHKERRQ(ierr);  }
//   if (set_contains(vars, "Soc")) {  ierr = Soc.write(filename.c_str()); CHKERRQ(ierr); }
//   if (set_contains(vars, "Soc_base")) { ierr = Soc_base.write(filename.c_str()); CHKERRQ(ierr); }
//   if (set_contains(vars, "Toc")) { ierr = Toc.write(filename.c_str()); CHKERRQ(ierr); }
//   if (set_contains(vars, "Toc_base")) {  ierr = Toc_base.write(filename.c_str()); CHKERRQ(ierr);  }
//   if (set_contains(vars, "Toc_inCelsius")) { ierr = Toc_inCelsius.write(filename.c_str()); CHKERRQ(ierr);  }
//   if (set_contains(vars, "T_star")) {  ierr = T_star.write(filename.c_str()); CHKERRQ(ierr);  }
//   if (set_contains(vars, "Toc_anomaly")) {  ierr = Toc_anomaly.write(filename.c_str()); CHKERRQ(ierr); }
//   if (set_contains(vars, "overturning")) {  ierr = overturning.write(filename.c_str()); CHKERRQ(ierr); }
//   if (set_contains(vars, "heatflux")) {  ierr = heatflux.write(filename.c_str()); CHKERRQ(ierr);  }
//   if (set_contains(vars, "basalmeltrate_shelf")) { ierr = basalmeltrate_shelf.write(filename.c_str()); CHKERRQ(ierr); }
//   return 0;
// }

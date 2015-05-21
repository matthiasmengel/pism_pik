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


PetscErrorCode POoceanboxmodel::init(PISMVars &vars) {
  PetscErrorCode ierr;

  ierr = verbPrintf(2, grid.com,
                    "* Initializing the ocean box model (based on Olbers & Hellmer (2010)...\n"); CHKERRQ(ierr);

  ierr = process_options(); CHKERRQ(ierr);

  ierr = set_vec_parameters("", ""); CHKERRQ(ierr);

  ierr = temp.create(grid, temp_name, false); CHKERRQ(ierr);
  ierr = mass_flux.create(grid, mass_flux_name, false); CHKERRQ(ierr); //NOTE: salinity instead of mass_flux

  ierr = temp.set_attrs("climate_forcing", "absolute temperature at ice shelf base", "Kelvin", ""); CHKERRQ(ierr); 
  // temp_boundlayer is not a standard name
  ierr = mass_flux.set_attrs("climate_forcing", "ocean salinity at ice shelf", "g/kg", ""); CHKERRQ(ierr); 
  // salinity_boundlayer is not a standard name

  ierr = temp.init(filename); CHKERRQ(ierr);
  ierr = mass_flux.init(filename); CHKERRQ(ierr);

  // mask to identify the ocean boxes
  ierr = BOXMODELmask.create(grid, "BOXMODELmask", true); CHKERRQ(ierr);
  ierr = BOXMODELmask.set_attrs("model_state", "mask displaying ocean box model grid","", ""); CHKERRQ(ierr);

  // mask to identify the grounded ice rises
  ierr = ICERISESmask.create(grid, "ICERISESmask", true); CHKERRQ(ierr);
  ierr = ICERISESmask.set_attrs("model_state", "mask displaying ice rises","", ""); CHKERRQ(ierr);
  ierr = PISMOptionsIsSet("-exclude_icerises", exicerises_set); CHKERRQ(ierr);

  // mask to calculate the mean salinity and ocean temperature in front of each basin 
  ierr = OCEANMEANmask.create(grid, "OCEANMEANmask", true); CHKERRQ(ierr);
  ierr = OCEANMEANmask.set_attrs("model_state", "mask displaying ocean region for parameter input","", ""); CHKERRQ(ierr);

  // FIXME delete CHECKmask 
  ierr = CHECKmask.create(grid, "CHECKmask", true); CHKERRQ(ierr);
  ierr = CHECKmask.set_attrs("model_state", "mask displaying rounded basins","", ""); CHKERRQ(ierr);


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

  mask = dynamic_cast<IceModelVec2Int*>(vars.get("mask"));
  if (!mask) { SETERRQ(grid.com, 1, "ERROR: mask is not available"); }

  basins = dynamic_cast<IceModelVec2S*>(vars.get("drainage_basins"));
  if (!basins) { SETERRQ(grid.com, 1, "ERROR: drainage basins is not available"); }

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


  bool gamma_T_set;
  PetscReal gamma_T = 1e-6;
  ierr = PISMOptionsReal("-gamma_T","-gamma_T",gamma_T, gamma_T_set); CHKERRQ(ierr);

  bool value_C_set;
  PetscReal value_C = 5e6;
  ierr = PISMOptionsReal("-value_C","-value_C",value_C, value_C_set); CHKERRQ(ierr);

  //set number of basins per option
  bool number_of_basins_set;
  numberOfBasins = 18;
  ierr = PISMOptionsInt("-number_of_basins","Number of Drainage Basins",numberOfBasins, number_of_basins_set); CHKERRQ(ierr);

  Toc_base_vec.resize(numberOfBasins);
  Soc_base_vec.resize(numberOfBasins);
  gamma_T_star_vec.resize(numberOfBasins);
  C_vec.resize(numberOfBasins);

  counter.resize(numberOfBasins);
  counter_GLbox.resize(numberOfBasins);
  counter_CFbox.resize(numberOfBasins);
  k_basins.resize(numberOfBasins);;

  mean_salinity_GLbox_vector.resize(numberOfBasins);
  mean_meltrate_GLbox_vector.resize(numberOfBasins);
  mean_overturning_GLbox_vector.resize(numberOfBasins);


  //FIXME calculate the number of basins 
  
  
  for(int i=0;i<numberOfBasins;i++) {
      Toc_base_vec[i] = -1.5; //dummy, FIXME why these values? 
      Soc_base_vec[i] = 34.5; //dummy
      gamma_T_star_vec[i]= gamma_T; 
      C_vec[i]           = value_C;
      if(i==1){
        ierr = verbPrintf(2, grid.com,"Using %d drainage basins and default values: gamma_T_star= %f, C = %f, calculate Soc and Toc from thetao and salinity...\n ", numberOfBasins, gamma_T_star_vec[i], C_vec[i]  ); CHKERRQ(ierr);       
      }
  }


  return 0;
}

const int POoceanboxmodel::shelf_unidentified = -99.0; // FIXME clean up and delete This should never show up in the .nc-files.
const int POoceanboxmodel::noshelf = 0.0; // FIXME clean up and delete 

const int POoceanboxmodel::box_unidentified = -99.0;     // This should never show up in the .nc-files.
const int POoceanboxmodel::box_neighboring = -1.0; // This should never show up in the .nc-files.
const int POoceanboxmodel::box_noshelf = 0.0;
const int POoceanboxmodel::box_GL = 1.0;  // ocean box covering the grounding line region
const int POoceanboxmodel::box_IF = 2.0;  // ocean box covering the rest of the ice shelf

const int POoceanboxmodel::maskfloating = MASK_FLOATING;  
const int POoceanboxmodel::maskocean = MASK_ICE_FREE_OCEAN;  
const int POoceanboxmodel::maskgrounded = MASK_GROUNDED;  


PetscErrorCode POoceanboxmodel::update(PetscReal my_t, PetscReal my_dt) {
  PetscErrorCode ierr = update_internal(my_t, my_dt); CHKERRQ(ierr);

  ierr = temp.at_time(t); CHKERRQ(ierr);
  ierr = mass_flux.at_time(t); CHKERRQ(ierr);

  ierr = verbPrintf(4, grid.com,"######### POoceanboxmodel::update() start\n"); CHKERRQ(ierr);
  ierr = verbPrintf(4, grid.com,"0  : calculating mean salinity and temperatures\n"); CHKERRQ(ierr);
  ierr = computeOCEANMEANS(); CHKERRQ(ierr);     

  ierr = verbPrintf(4, grid.com,"A  : calculating shelf_base_temperature\n"); CHKERRQ(ierr);
  if (exicerises_set) {
    ierr = identifyICERISES(); CHKERRQ(ierr);}
  ierr = extentOfIceShelves(); CHKERRQ(ierr); 
  ierr = identifyBOXMODELmask(); CHKERRQ(ierr); 
  ierr = oceanTemperature(); CHKERRQ(ierr);
  ierr = Toc.copy_to(temp); CHKERRQ(ierr);

  ierr = verbPrintf(4, grid.com,"B  : calculating shelf_base_mass_flux\n"); CHKERRQ(ierr);

  //Assumes that mass flux is proportional to the shelf-base heat flux.
  ierr = basalMeltRateForGroundingLineBox(); CHKERRQ(ierr);
  ierr = basalMeltRateForIceFrontBox(); CHKERRQ(ierr); // TODO Diese Routinen woanders aufrufen (um Dopplung zu vermeiden)
  ierr = basalMeltRateForOtherShelves(); CHKERRQ(ierr);
  ierr = basalmeltrate_shelf.copy_to(mass_flux); CHKERRQ(ierr);

  return 0;
}



//! Round basin mask non integer values to an integral value of the next neigbor
PetscErrorCode POoceanboxmodel::roundBasins(PetscInt i, PetscInt j) {
  PetscErrorCode ierr;

  ierr = basins->begin_access();   CHKERRQ(ierr);

  PetscInt  id_rounded = static_cast<int>(round((*basins)(i,j)));
  PetscReal id_fractional = (*basins)(i,j);
  PetscInt  id = -1;
  
  if(id_fractional-static_cast<double>(id_rounded) != 0.0){ //if id_fractional differns from integer value  
    //FIXME needs ghost comunication?
    if ((i-1 >= grid.xs) && (j-1 >= grid.ys)&& ((*basins)(i-1,j-1)-static_cast<double>(static_cast<int>(round((*basins)(i-1,j-1)))) == 0.0)){
      id = static_cast<int>(round((*basins)(i-1,j-1)));
    } else if((i+1 <= grid.xs+grid.xm) && (j-1 >= grid.ys)&& ((*basins)(i+1,j-1)-static_cast<double>(static_cast<int>(round((*basins)(i+1,j-1)))) == 0.0)){
      id = static_cast<int>(round((*basins)(i+1,j-1)));
    } else if ((i-1 >= grid.xs) && (j+1 <= grid.ys+grid.ym)&& ((*basins)(i-1,j+1)-static_cast<double>(static_cast<int>(round((*basins)(i-1,j+1)))) == 0.0)){
      id = static_cast<int>(round((*basins)(i-1,j+1)));
    } else if ((i+1 <= grid.xs+grid.xm) && (j+1 <= grid.ys+grid.ym)&& ((*basins)(i+1,j+1)-static_cast<double>(static_cast<int>(round((*basins)(i+1,j+1)))) == 0.0)){
      id = static_cast<int>(round((*basins)(i+1,j+1)));
    } else { //if no neigbour has an integer id
      id = id_rounded;
      //ierr = verbPrintf(2, grid.com,"no neighbour has an integer id\n"); CHKERRQ(ierr);
    }
  } else { // if id_rounded is id_fractional
    id = id_rounded;
  }
  id = id_rounded;

  //ierr = DRAINAGEmask.end_access();   CHKERRQ(ierr);
  ierr = basins->end_access();   CHKERRQ(ierr); 

  return id;
}



//! When ocean_given is set compute mean salinity and temperature in each basin.
PetscErrorCode POoceanboxmodel::computeOCEANMEANS() {
  /*NOTE: 
  -basin 0 (which is not existing) is the margin of the domain-it has currently seven boxes which lie in the OCEANMEANmask
    -> shall we exclude this case? Or do we continue with it, e.g. for the case of just one basin?
  -we call this rountine before we call extentOfIceShelves, where the Drainage mask is set to basins' values. 
    -> should we use basins instead? Or define the drainage mask earlier? This could be the reason why
  -during the initialisation we have wrong values for the salinity and temperature
  -currently, by iterating, we may include ocean cells which neighbor directly bedrock or grounded ice
    -> do we want to exclude them?
  -how big do we want to have the region over which we calculate the ocean means? 
    -> search in literature? currently 200km
  -currently we enter this routine for all calculations
    -> do we want to have an if(there is an ocean filed) {call routine?}
  -check the values for basin 14 (there are just 6 cells in the mask): Calculation is correct!
  */

  PetscErrorCode ierr;
  ierr = verbPrintf(4, grid.com,"0a : in computeOCEANMEANS routine \n"); CHKERRQ(ierr);

  PetscInt number_of_iterations = 0; // to get circa 200km before the ice shelf, FIXME do we want this to be chosable?
  number_of_iterations = round(200e3 / ((grid.dx +grid.dy)/2.0)); // meter 
  //ierr = verbPrintf(2, grid.com,"0a  : number of Iterations = %d , %f\n", number_of_iterations, 200e3 / ((grid.dx +grid.dy)/2.0)); CHKERRQ(ierr);

  PetscInt ocean_mean_region = 2,
           ocean_mean_region_candidate = 1, 
           unidentified = -1;

  PetscScalar lm_count[numberOfBasins]; //count cells to take mean over for each basin
  PetscScalar m_count[numberOfBasins];
  PetscScalar lm_Sval[numberOfBasins]; //add salinity for each basin
  PetscScalar lm_Tval[numberOfBasins]; //add temperature for each basin
  PetscScalar m_Tval[numberOfBasins]; //FIXME delete?
  PetscScalar m_Sval[numberOfBasins]; //FIXME delete?

  for(int k=0;k<numberOfBasins;k++){
    m_count[k]=0.0;
    lm_count[k]=0.0;
    lm_Sval[k]=0.0;
    lm_Tval[k]=0.0;
    m_Tval[k]=0.0; //FIXME delete?
    m_Sval[k]=0.0; //FIXME delete?
  }

  ierr = OCEANMEANmask.begin_access();   CHKERRQ(ierr); 
  ierr = mask->begin_access();   CHKERRQ(ierr); 
  ierr = temp.begin_access();   CHKERRQ(ierr); 
  ierr = mass_flux.begin_access();   CHKERRQ(ierr); //salinitiy

  //label ocean cells neighboring ice shelf cells as ocean_mean_region
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
        if ((*mask)(i,j)==maskocean && 
          ((*mask)(i,j+1)==maskfloating || (*mask)(i,j-1)==maskfloating || 
           (*mask)(i+1,j)==maskfloating || (*mask)(i-1,j)==maskfloating)){ // if the cell is ocean directly at the ice shelf front
            
            PetscInt shelf_id = roundBasins(i,j);
            OCEANMEANmask(i,j)=ocean_mean_region;
            lm_count[shelf_id]+=1;
            lm_Sval[shelf_id]+=mass_flux(i,j);
            lm_Tval[shelf_id]+=temp(i,j);

        }
        else { // not ocean or not neighboring the calving front
          OCEANMEANmask(i,j)=unidentified;
        }
      }
    }


  ierr = OCEANMEANmask.end_access();   CHKERRQ(ierr); //FIXME necessary? 
  ierr = mask->end_access();   CHKERRQ(ierr);  
  ierr = temp.end_access();   CHKERRQ(ierr);  
  ierr = mass_flux.end_access();   CHKERRQ(ierr);  //salinity
  //FIXME GhostComm?
  // FIXME do we need to do something with the lcounters for parallel computing reasons?

  ierr = OCEANMEANmask.begin_access();   CHKERRQ(ierr); 
  ierr = mask->begin_access();   CHKERRQ(ierr); 

  for(PetscInt k=0;k<number_of_iterations;k++){
    //NOTE at the moment we do not only extend the mask towards the ocean, but also along the non-ice-shelf coasts. Do we care? 

    // Find candidates for ocean_mean_region:
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
        if ((*mask)(i,j)==maskocean && OCEANMEANmask(i,j)==unidentified &&
            (OCEANMEANmask(i,j+1)==ocean_mean_region || OCEANMEANmask(i,j-1)==ocean_mean_region || 
             OCEANMEANmask(i+1,j)==ocean_mean_region || OCEANMEANmask(i-1,j)==ocean_mean_region)){ // if the cell is ocean directly at the ice shelf front
             OCEANMEANmask(i,j)=ocean_mean_region_candidate;
        } //if
      } //j
    } //i

    //Label candidads for ocean_mean_region:
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
        if ( OCEANMEANmask(i,j)==ocean_mean_region_candidate){

            PetscInt shelf_id = roundBasins(i,j);
            OCEANMEANmask(i,j)=ocean_mean_region;
            lm_count[shelf_id]+=1;
            lm_Sval[shelf_id]+=mass_flux(i,j);
            lm_Tval[shelf_id]+=temp(i,j);
        } //if
      } //j
    } //i

  } // number of iteration

  ierr = OCEANMEANmask.end_access();   CHKERRQ(ierr); 
  ierr = mask->end_access();   CHKERRQ(ierr);  

  for(int i=0;i<numberOfBasins;i++) {
    ierr = PISMGlobalSum(&lm_count[i], &m_count[i], grid.com); CHKERRQ(ierr);
    ierr = PISMGlobalSum(&lm_Sval[i], &m_Sval[i], grid.com); CHKERRQ(ierr);
    ierr = PISMGlobalSum(&lm_Tval[i], &m_Tval[i], grid.com); CHKERRQ(ierr);
  } 

  for(int i=0;i<numberOfBasins;i++) {
        
    if(i>0 && m_count[i]==0){ //if basin is not dummy basin 0 and there are no cells of this basin in the OCEANMEANmask 
      ierr = verbPrintf(2, grid.com,"PISM_WARNING: basin %d contains no cells in OCEANMEANmask, no mean salinity od temperature values are computed! \n ", i); CHKERRQ(ierr);   
    }

    m_Sval[i] = m_Sval[i] / m_count[i];
    m_Tval[i] = m_Tval[i] / m_count[i]; 
    
    Toc_base_vec[i]=m_Tval[i] - 273.15;
    Soc_base_vec[i]=m_Sval[i];
    ierr = verbPrintf(5, grid.com,"0: basin= %d, temp =%.3f, salinity=%.3f\n", i, Toc_base_vec[i], Soc_base_vec[i]); CHKERRQ(ierr);
  } 

  return 0;
}



//! Identify if grounded area is continent or ice rise

PetscErrorCode POoceanboxmodel::identifyICERISES() {
  
  PetscErrorCode ierr;
  ierr = verbPrintf(4, grid.com,"A1a: in identifyICERISES rountine\n"); CHKERRQ(ierr);

  PetscInt seed_x = (grid.Mx - 1)/2,
           seed_y = (grid.My - 1)/2;

  PetscScalar lcontinent_identified = 0.0,
              all_continent_identified = 1.0,
              previous_step_identified = 0.0;

  PetscInt continent = 2,
           ocean = 0,
           icerise = 1,
           unidentified = -1;

  PetscInt iteration_round = 0;

  ierr = ICERISESmask.begin_access();   CHKERRQ(ierr);
  ierr = mask->begin_access();   CHKERRQ(ierr); 

  ICERISESmask.set(unidentified);
  
  if ((seed_x >= grid.xs) && (seed_x < grid.xs+grid.xm) && (seed_y >= grid.ys)&& (seed_y <= grid.ys+grid.ym)){
    ICERISESmask(seed_x,seed_y)=continent;
  }

  ierr = ICERISESmask.end_access();   CHKERRQ(ierr);
  ierr = mask->end_access();   CHKERRQ(ierr);  


  // find continent first
  while(all_continent_identified > previous_step_identified){

    iteration_round+=1;
    previous_step_identified = all_continent_identified;

    ierr = ICERISESmask.begin_access();   CHKERRQ(ierr);
    ierr = mask->begin_access();   CHKERRQ(ierr); 

    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
        if ((*mask)(i,j)==maskgrounded && ICERISESmask(i,j)==unidentified &&
          (ICERISESmask(i,j+1)==continent || ICERISESmask(i,j-1)==continent || 
           ICERISESmask(i+1,j)==continent || ICERISESmask(i-1,j)==continent)){
           ICERISESmask(i,j)=continent;
           lcontinent_identified+=1;
        }
        else if ((*mask)(i,j)!=maskgrounded){
          ICERISESmask(i,j)=ocean;
        }
      }
    }

    ierr = ICERISESmask.end_access();   CHKERRQ(ierr);
    ierr = mask->end_access();   CHKERRQ(ierr);   
    ierr = ICERISESmask.beginGhostComm(); CHKERRQ(ierr);
    ierr = ICERISESmask.endGhostComm(); CHKERRQ(ierr);

    ierr = PISMGlobalSum(&lcontinent_identified, &all_continent_identified, grid.com); CHKERRQ(ierr);
  
  }

   //set ice rises value
   ierr = ICERISESmask.begin_access();   CHKERRQ(ierr);
    ierr = mask->begin_access();   CHKERRQ(ierr); 

    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
        if (ICERISESmask(i,j)==unidentified){
          ICERISESmask(i,j)=icerise;
        }
      }
    }
    ierr = ICERISESmask.end_access();   CHKERRQ(ierr);
    ierr = mask->end_access();   CHKERRQ(ierr);   


  return 0;
}



//! Compute the extent of the ice shelves of each basin/region (i.e. counter)
/*  Start to fill in the BOXMODELmask: -identify the boxes directly at the groundingline
                                       -identify the boxes directly at the calving front 
                                       -Set all other shelf_boxes to shelf_unidentified*/
PetscErrorCode POoceanboxmodel::extentOfIceShelves() {
  
  PetscErrorCode ierr;
  ierr = verbPrintf(4, grid.com,"A1b: in extent of ice shelves rountine\n"); CHKERRQ(ierr);

  PetscScalar lcounter_box_unidentified = 0; //count the total amount of unidentified shelf boxes.
  PetscScalar lcounter[numberOfBasins];
  PetscScalar lcounter_CFbox[numberOfBasins];
  PetscScalar lcounter_GLbox[numberOfBasins];

  for (int i=0;i<numberOfBasins;i++){ 
    lcounter[i]=0.0; 
    lcounter_CFbox[i]=0.0;
    lcounter_GLbox[i]=0.0;
  }

  PetscInt continent = 2;
  
  ierr = mask->begin_access();   CHKERRQ(ierr); 
  ierr = BOXMODELmask.begin_access(); CHKERRQ(ierr);
  if (exicerises_set) { ierr = ICERISESmask.begin_access(); CHKERRQ(ierr);}


  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if ((*mask)(i,j)==maskfloating){ //if this is a shelf cell
        lcounter[roundBasins(i,j)]++;

        bool neighbor_to_land;
        if (exicerises_set) {
          neighbor_to_land = (  ICERISESmask(i,j+1)==continent || ICERISESmask(i,j-1)==continent || 
                                ICERISESmask(i+1,j)==continent || ICERISESmask(i-1,j)==continent || 
                                ICERISESmask(i+1,j+1)==continent || ICERISESmask(i+1,j-1)==continent ||
                                ICERISESmask(i-1,j+1)==continent || ICERISESmask(i-1,j-1)==continent );
        } else {
          neighbor_to_land = (  (*mask)(i,j+1)<maskfloating || (*mask)(i,j-1)<maskfloating || 
                                (*mask)(i+1,j)<maskfloating || (*mask)(i-1,j)<maskfloating || 
                                (*mask)(i+1,j+1)<maskfloating || (*mask)(i+1,j-1)<maskfloating ||
                                (*mask)(i-1,j+1)<maskfloating || (*mask)(i-1,j-1)<maskfloating );
        }

        if (neighbor_to_land){
             // i.e. there is a grounded neighboring cell (which is not ice rise!)
            BOXMODELmask(i,j) = box_GL;
            lcounter_GLbox[roundBasins(i,j)]++;
        } else if ((*mask)(i,j+1)==maskocean || (*mask)(i,j-1)==maskocean || 
            (*mask)(i+1,j)==maskocean || (*mask)(i-1,j)==maskocean || 
            (*mask)(i+1,j+1)==maskocean || (*mask)(i+1,j-1)==maskocean ||
            (*mask)(i-1,j+1)==maskocean || (*mask)(i-1,j-1)==maskocean ){ // i.e. there is an ocean neighboring cell 
            BOXMODELmask(i,j) = box_IF;
            lcounter_CFbox[roundBasins(i,j)]++;
		    } else { // i.e., all other floating boxes
            BOXMODELmask(i,j) = box_unidentified;
            lcounter_box_unidentified++;
		    }

      }else{ // i.e., not floating
        BOXMODELmask(i,j) = box_noshelf;
      }
    }
  }

if (exicerises_set) { ierr = ICERISESmask.end_access();   CHKERRQ(ierr);}
  ierr = mask->end_access();   CHKERRQ(ierr); 
  ierr = BOXMODELmask.end_access(); CHKERRQ(ierr);

  ierr = BOXMODELmask.beginGhostComm(); CHKERRQ(ierr);
  ierr = BOXMODELmask.endGhostComm(); CHKERRQ(ierr);
  

  ierr = PISMGlobalSum(&lcounter_box_unidentified, &counter_box_unidentified, grid.com); CHKERRQ(ierr);
  for(int i=0;i<numberOfBasins;i++) {
  	ierr = PISMGlobalSum(&lcounter[i], &counter[i], grid.com); CHKERRQ(ierr);
    ierr = PISMGlobalSum(&lcounter_CFbox[i], &counter_CFbox[i], grid.com); CHKERRQ(ierr);
    ierr = PISMGlobalSum(&lcounter_GLbox[i], &counter_GLbox[i], grid.com); CHKERRQ(ierr);
  } 
  for(int i=0;i<numberOfBasins;i++){ ierr = verbPrintf(5, grid.com,"AfterExtentOfIceShelves: basin= %d, counter[i] = %.0f, counter_CFbox= %.0f, counter_GLbox = %.0f\n", i, counter[i], counter_CFbox[i], counter_GLbox[i]); CHKERRQ(ierr);}
 
  return 0;
}



//! extent the grouningline box with neighboring unidentified shelf cells
PetscErrorCode POoceanboxmodel::extendGLBox() {
  PetscErrorCode ierr;

  PetscScalar lcounter_box_unidentified = 0.0; 
  PetscScalar all_counter_box_unidentified = 0.0;
  PetscScalar lcounter_GLbox[numberOfBasins];
  PetscScalar all_counter_GLbox[numberOfBasins];
  for (int i=0;i<numberOfBasins;i++){ lcounter_GLbox[i]=0.0; all_counter_GLbox[i]=0.0; }; 

  ierr = mask->begin_access();   CHKERRQ(ierr); 
  ierr = BOXMODELmask.begin_access(); CHKERRQ(ierr);

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (BOXMODELmask(i,j)==box_unidentified && 
      	  (BOXMODELmask(i,j+1)==box_GL || BOXMODELmask(i,j-1)==box_GL || 
           BOXMODELmask(i+1,j)==box_GL || BOXMODELmask(i-1,j)==box_GL // || BOXMODELmask(i+1,j+1)==box_GL || BOXMODELmask(i+1,j-1)==box_GL || BOXMODELmask(i-1,j+1)==box_GL || BOXMODELmask(i-1,j-1)==box_GL
          ) ){ // i.e. this is an unidentified shelf cell with a neighbor that is in the GLbox
      	  BOXMODELmask(i,j)=box_neighboring; 
      }
    }
  }

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (BOXMODELmask(i,j)==box_neighboring){ 
      	  BOXMODELmask(i,j)=box_GL;
      	  lcounter_box_unidentified++; 
          lcounter_GLbox[roundBasins(i,j)]++;
      }
    }
  }

  ierr = mask->end_access();   CHKERRQ(ierr); 
  ierr = BOXMODELmask.end_access(); CHKERRQ(ierr);
  ierr = BOXMODELmask.beginGhostComm(); CHKERRQ(ierr);
  ierr = BOXMODELmask.endGhostComm(); CHKERRQ(ierr);
  
  ierr = PISMGlobalSum(&lcounter_box_unidentified, &all_counter_box_unidentified, grid.com); CHKERRQ(ierr);
  counter_box_unidentified -= all_counter_box_unidentified;

  for(int i=0;i<numberOfBasins;i++) { 
    ierr = PISMGlobalSum(&lcounter_GLbox[i], &all_counter_GLbox[i], grid.com); CHKERRQ(ierr);
    counter_GLbox[i] += all_counter_GLbox[i];
  } 

  return 0;
}



//! extend the ice_front box with neighboring unidentified shelf cells.
PetscErrorCode POoceanboxmodel::extendIFBox() {
  
  PetscErrorCode ierr;

  PetscScalar lcounter_box_unidentified = 0.0; 
  PetscScalar all_counter_box_unidentified = 0.0;
  PetscScalar lcounter_CFbox[numberOfBasins];
  PetscScalar all_counter_CFbox[numberOfBasins];
  for (int i=0;i<numberOfBasins;i++){ lcounter_CFbox[i] = 0.0; all_counter_CFbox[i] = 0.0;}

  ierr = mask->begin_access();   CHKERRQ(ierr); 
  ierr = BOXMODELmask.begin_access(); CHKERRQ(ierr);
  
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (BOXMODELmask(i,j)==box_unidentified && 
      	  (BOXMODELmask(i,j+1)==box_IF || BOXMODELmask(i,j-1)==box_IF || 
           BOXMODELmask(i+1,j)==box_IF || BOXMODELmask(i-1,j)==box_IF  // || BOXMODELmask(i+1,j+1)==box_IF || BOXMODELmask(i+1,j-1)==box_IF || //FIXME eight neigbors?BOXMODELmask(i-1,j+1)==box_IF || BOXMODELmask(i-1,j-1)==box_IF
          ) ){ // i.e. this is an unidentified shelf cell with a neighbor that is in the IFbox
      	  BOXMODELmask(i,j)=box_neighboring;  
      }
    }
  }

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (BOXMODELmask(i,j)==box_neighboring){ 
      	  BOXMODELmask(i,j)=box_IF; 
      	  lcounter_box_unidentified++; 
          lcounter_CFbox[roundBasins(i,j)]++;
      }
    }
  }

  ierr = mask->end_access();   CHKERRQ(ierr); 
  ierr = BOXMODELmask.end_access(); CHKERRQ(ierr);
  ierr = BOXMODELmask.beginGhostComm(); CHKERRQ(ierr);
  ierr = BOXMODELmask.endGhostComm(); CHKERRQ(ierr);

  ierr = PISMGlobalSum(&lcounter_box_unidentified, &all_counter_box_unidentified, grid.com); CHKERRQ(ierr);
  counter_box_unidentified -= all_counter_box_unidentified;
  
  for(int i=0;i<numberOfBasins;i++) { 
    ierr = PISMGlobalSum(&lcounter_CFbox[i], &all_counter_CFbox[i], grid.com); CHKERRQ(ierr);
    counter_CFbox[i] += all_counter_CFbox[i];
  } 

  return 0;
}


//! Compute the boxmodelmask, calculate the extent of each box in each region
PetscErrorCode POoceanboxmodel::identifyBOXMODELmask() {

  PetscErrorCode ierr;
  ierr = verbPrintf(4, grid.com,"A1c: in identify boxmodel mask rountine\n"); CHKERRQ(ierr);
  
  PetscScalar lcounter_box_unidentified = counter_box_unidentified+1.0; 
  
  while((counter_box_unidentified > 0.0)&&(lcounter_box_unidentified != counter_box_unidentified)){

      lcounter_box_unidentified = counter_box_unidentified;
      ierr = verbPrintf(5, grid.com,"A1b: counter_box_unidentified=%.0f, lcounter_box_unidentified=%.0f\n", counter_box_unidentified, lcounter_box_unidentified); CHKERRQ(ierr);

      //ierr = extendIFBox(); CHKERRQ(ierr); // FIXME size depends on how often this routine is called
      //ierr = extendIFBox(); CHKERRQ(ierr);
      //ierr = extendIFBox(); CHKERRQ(ierr);
      ierr = extendIFBox(); CHKERRQ(ierr);
      ierr = extendIFBox(); CHKERRQ(ierr);
      ierr = extendIFBox(); CHKERRQ(ierr);
      ierr = extendGLBox(); CHKERRQ(ierr); 

  }

  if(counter_box_unidentified>0.0){ //if there are still floating cells which were not labels before, i.e. a shelf connected to an island (in the ex_icerises case)
  	ierr = BOXMODELmask.begin_access(); CHKERRQ(ierr);
  	for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
        if (BOXMODELmask(i,j)==box_unidentified ){ // i.e. this is an unidentified shelf cell with a neighbor that is in the IFbox
      	  	ierr = verbPrintf(4, grid.com,"A1b: left over i=%d, j=%d \n", i,j); CHKERRQ(ierr);
      	  	BOXMODELmask(i,j)=box_GL;   
            //FIXME counter für Anzahl der Groundingline boxen hochzählen!
      	}
      }
  	}
  	ierr = BOXMODELmask.end_access(); CHKERRQ(ierr);
  }

  //FIXME delete!
  ierr = CHECKmask.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
        CHECKmask(i,j) = roundBasins(i,j);
      }
  }
  ierr = CHECKmask.end_access(); CHKERRQ(ierr);
  //FIXME end of delete!

  for(int i=0;i<numberOfBasins;i++){ ierr = verbPrintf(5, grid.com,"A1b: basin= %d, counter[i] = %.0f, counter_CFbox= %.0f, counter_GLbox = %.0f, ratio_CF_box= %.3f, ratio_GL_box= %.3f\n", i, counter[i], counter_CFbox[i], counter_GLbox[i], counter_CFbox[i]/counter[i], counter_GLbox[i]/counter[i]); CHKERRQ(ierr);}
  return 0;
}





/*!
Compute ocean temperature outside of the ice shelf cavities.
*/
PetscErrorCode POoceanboxmodel::oceanTemperature() { 
	/* IN THIS ROUTINE
	Before: Set the ocean temperature in front of the ice shelf for each region (This is T_0 in the paper). Note that here, we have a field.
	Want: Set the ocean temperature in front of the ice shelf for each basin (This is T_0 in the paper). Note that here, we have a field.
	TODO : calculate the mean over the ocean temperature in front of the ice shelf for each basin*/
  PetscErrorCode ierr;
  ierr = verbPrintf(4, grid.com,"A2 : in ocean temp rountine\n"); CHKERRQ(ierr);

  ierr = mask->begin_access();   CHKERRQ(ierr);
  ierr = Soc_base.begin_access();   CHKERRQ(ierr);
  ierr = Toc_base.begin_access();   CHKERRQ(ierr);
  ierr = Toc_anomaly.begin_access();   CHKERRQ(ierr);
  ierr = Toc.begin_access();   CHKERRQ(ierr);
  ierr = temp.begin_access();   CHKERRQ(ierr);
  ierr = mass_flux.begin_access();   CHKERRQ(ierr);


  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {

      // make sure all temperatures are zero at the beginning of each timestep
      Toc(i,j) = 273.15; // in K
      Toc_base(i,j) = 273.15;  // in K
      Toc_anomaly(i,j) = 0.0;  // in K or °C
      Soc_base(i,j) = 0.0; // in psu


      if ((*mask)(i,j)==maskfloating){
        PetscInt shelf_id = roundBasins(i,j);
        Toc_base(i,j) = 273.15 + Toc_base_vec[shelf_id ];
      	Soc_base(i,j) =  Soc_base_vec[shelf_id ];

        //! salinity and temperature for grounding line box
        if ( Soc_base(i,j) == 0.0 || Toc_base_vec[shelf_id] == 0.0 ) {
          ierr = verbPrintf(2, grid.com,"PISM_ERROR: Missing Soc_base and Toc_base for %d, %d, basin %d \n   Aborting... \n", i, j, shelf_id); CHKERRQ(ierr);
          PISMEnd();
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

  ierr = mask->end_access();   CHKERRQ(ierr);
  ierr = Soc_base.end_access();   CHKERRQ(ierr);
  ierr = Toc_base.end_access();   CHKERRQ(ierr);
  ierr = Toc_anomaly.end_access();   CHKERRQ(ierr);
  ierr = Toc.end_access();   CHKERRQ(ierr);
  ierr = temp.end_access();   CHKERRQ(ierr);
  ierr = mass_flux.end_access();   CHKERRQ(ierr);

  return 0;
}



// NOTE Mean Gl_box meltrate is needed for basalMeltRateForIceFrontBox(). Here, mean is taken over all shelves in one drainage basin!

//! Compute the basal melt / refreezing rates for each shelf cell bordering the grounding line box
PetscErrorCode POoceanboxmodel::basalMeltRateForGroundingLineBox() {
  PetscErrorCode ierr;
  ierr = verbPrintf(4, grid.com,"B1 : in basal melt rate gl rountine\n"); CHKERRQ(ierr);

  //! constants 
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

  PetscScalar lcounter_edge_of_GLbox_vector[numberOfBasins], 
              lmean_salinity_GLbox_vector[numberOfBasins], 
              lmean_meltrate_GLbox_vector[numberOfBasins], 
              lmean_overturning_GLbox_vector[numberOfBasins];
  for (int i= 0;i<numberOfBasins;i++){ 
    lcounter_edge_of_GLbox_vector[i]=0.0;
    lmean_salinity_GLbox_vector[i]=0.0;
    lmean_meltrate_GLbox_vector[i]=0.0;
    lmean_overturning_GLbox_vector[i]=0.0;
  }
  

  ierr = ice_thickness->begin_access(); CHKERRQ(ierr);
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

      PetscInt shelf_id = roundBasins(i,j);

      // Make sure everything is at default values at the beginning of each timestep
      T_star(i,j) = 0.0; // in °C
      Toc_inCelsius(i,j) = 0.0; // in °C
      Soc(i,j) = 0.0; // in psu

      if (BOXMODELmask(i,j) == box_GL && shelf_id > 0.0){

	  	  const PetscScalar pressure = rhoi * earth_grav * (*ice_thickness)(i,j) * 1e-4; // MUST be in dbar  // NOTE 1dbar = 10000 Pa = 1e4 kg m-1 s-2
		    // FIXME need to include atmospheric pressure?
	  	  T_star(i,j) = a*Soc_base(i,j) + b - c*pressure - (Toc_base(i,j)-273.15+Toc_anomaly(i,j)); // in °C

	  	  PetscScalar gamma_T_star,C1,g1;
        PetscInt shelf_id = roundBasins(i,j);
	  	  gamma_T_star = gamma_T_star_vec[shelf_id];
	  	  C1 = C_vec[shelf_id];
	  	  g1 = (counter_GLbox[shelf_id] * grid.dx * grid.dy) * gamma_T_star / (C1*rho_star); 

        //! temperature for grounding line box
	  	  PetscScalar helpterm1 = g1/(beta*(Soc_base(i,j) / (nu*lambda)) - alpha);                  // in 1 / (1/°C) = °C
	  	  PetscScalar helpterm2 = (g1*T_star(i,j)) / (beta*(Soc_base(i,j) / (nu*lambda)) - alpha); // in °C / (1/°C) = °C^2

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
	  	  // NOTE Actually, there is of course no overturning-FIELD, it is only a scalar for each shelf.
	  	  // Here, I compute overturning as 	MEAN[C1*rho_star* (beta*(Soc_base(i,j)-Soc(i,j)) - alpha*((Toc_base(i,j)-273.15+Toc_anomaly(i,j))-Toc_inCelsius(i,j)))]
	  	  // while in fact it should be 	C1*rho_star* (beta*(Soc_base-MEAN[Soc(i,j)]) - alpha*((Toc_base-273.15+Toc_anomaly)-MEAN[Toc_inCelsius(i,j)]))
	  	  // which is the SAME since Soc_base, Toc_base and Toc_anomaly are the same FOR ALL i,j CONSIDERED, so this is just nomenclature!
	  	  overturning(i,j) = C1*rho_star* (beta*(Soc_base(i,j)-Soc(i,j)) - alpha*((Toc_base(i,j)-273.15+Toc_anomaly(i,j))-Toc_inCelsius(i,j))); // in m^3/s

        if (BOXMODELmask(i-1,j)==box_IF || BOXMODELmask(i+1,j)==box_IF || BOXMODELmask(i,j-1)==box_IF || BOXMODELmask(i,j+1)==box_IF){ 
        // i.e., if this cell is from the GL box and one of the neighbours is from the CF box - It is important to only take the border of the grounding line box 
        // to the calving front box into account, because the following mean value will be used to compute the value for the calving front box. I.e., this helps avoiding discontinuities!
          lcounter_edge_of_GLbox_vector[shelf_id]++;
          lmean_salinity_GLbox_vector[shelf_id] += Soc(i,j);
          lmean_meltrate_GLbox_vector[shelf_id] += basalmeltrate_shelf(i,j);
          lmean_overturning_GLbox_vector[shelf_id] += overturning(i,j);

	  	  } // no else-case necessary since all variables are set to zero at the beginning of this routine

      }else { // i.e., not GL_box
		      basalmeltrate_shelf(i,j) = 0.0;
      }

    } // end j
  } // end i

  ierr = ice_thickness->end_access(); CHKERRQ(ierr);
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
     


  for(int i=0;i<numberOfBasins;i++) {
    PetscScalar counter_edge_of_GLbox_vector=0.0;
    ierr = PISMGlobalSum(&lcounter_edge_of_GLbox_vector[i], &counter_edge_of_GLbox_vector, grid.com); CHKERRQ(ierr);
    ierr = PISMGlobalSum(&lmean_meltrate_GLbox_vector[i], &mean_meltrate_GLbox_vector[i], grid.com); CHKERRQ(ierr);
    ierr = PISMGlobalSum(&lmean_salinity_GLbox_vector[i], &mean_salinity_GLbox_vector[i], grid.com); CHKERRQ(ierr);
    ierr = PISMGlobalSum(&lmean_overturning_GLbox_vector[i], &mean_overturning_GLbox_vector[i], grid.com); CHKERRQ(ierr);

    if (counter_edge_of_GLbox_vector>0.0){
      mean_salinity_GLbox_vector[i] = mean_salinity_GLbox_vector[i]/counter_edge_of_GLbox_vector;
      mean_meltrate_GLbox_vector[i] = mean_meltrate_GLbox_vector[i]/counter_edge_of_GLbox_vector;
      mean_overturning_GLbox_vector[i] = mean_overturning_GLbox_vector[i]/counter_edge_of_GLbox_vector;
    } else { // This means that there is no [cell from the GLbox neighboring a cell from the CFbox], NOT necessarily that there is no GLbox!
      mean_salinity_GLbox_vector[i]=0.0; mean_meltrate_GLbox_vector[i]=0.0; mean_overturning_GLbox_vector[i]=0.0;
    }
      ierr = verbPrintf(5, grid.com,"   %d: cnt=%f, sal=%.3f, melt=%.3f, over=%.3f \n", i,counter_edge_of_GLbox_vector,mean_salinity_GLbox_vector[i],mean_meltrate_GLbox_vector[i],mean_overturning_GLbox_vector[i]) ; CHKERRQ(ierr); 
  }
  return 0;
}




//! Compute the basal melt / refreezing rates for each shelf cell bordering the ice front box
PetscErrorCode POoceanboxmodel::basalMeltRateForIceFrontBox() {
  PetscErrorCode ierr;  // FIXME redo all verbprintfs!
  ierr = verbPrintf(4, grid.com,"B2 : in bm ice front rountine\n"); CHKERRQ(ierr);
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

      PetscInt shelf_id = roundBasins(i,j);

      if (BOXMODELmask(i,j)==box_IF && shelf_id > 0.0){

	  	  const PetscScalar pressure = rhoi * earth_grav * (*ice_thickness)(i,j) * 1e-4; // MUST be in dbar  // NOTE 1dbar = 10000 Pa = 1e4 kg m-1 s-2
									       // FIXME need to include atmospheric pressure?
	  	  T_star(i,j) = a*Soc_base(i,j) + b - c*pressure - (Toc_base(i,j) - 273.15 + Toc_anomaly(i,j)); // in °C

	  	  PetscScalar  gamma_T_star,area_GLbox,area_CFbox,mean_salinity_in_GLbox,mean_meltrate_in_GLbox,mean_overturning_in_GLbox;

        gamma_T_star = gamma_T_star_vec[shelf_id]; 
        area_CFbox = (counter_CFbox[shelf_id] * grid.dx * grid.dy); 
        area_GLbox = (counter_GLbox[shelf_id] * grid.dx * grid.dy);  
        mean_salinity_in_GLbox = mean_salinity_GLbox_vector[shelf_id];
        mean_meltrate_in_GLbox = mean_meltrate_GLbox_vector[shelf_id];
        mean_overturning_in_GLbox = mean_overturning_GLbox_vector[shelf_id];


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
          ierr = verbPrintf(2, grid.com,"PISM_ERROR: Probably mean_overturning_in_GLbox = %f is zero, check if there is a grounding line box in basin %d , \n   ", mean_overturning_in_GLbox, shelf_id); CHKERRQ(ierr);
          // FIXME rewrite this warning? Do not stop but calculate melt rates according to Beckmann-Gose?
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
	    	  // this must not occur since there must always be a GL_box neighbor 
	    	  ierr = verbPrintf(2, grid.com, "PISM_ERROR: DETECTION CFBOX: There is no neighbouring grounding line box for this calving front box at %d,%d! \nThis will lead to a zero k4 and in turn to NaN in Soc, Toc_inCelsius and basalmeltrate_shelf. After the next massContExplicitStep(), H will be NaN, too! This will cause ks in temperatureStep() to be NaN and lead to a Segmentation Violation! \nIn particular: basin_id=%f, BOXMODELmask=%f, H=%f, T_star=%f, \narea_GLbox=%e, area_CFbox=%e, mean_salinity_in_GLbox=%f, mean_meltrate_in_GLbox=%e, mean_overturning_in_GLbox=%e, \nk1=%e,k2=%e,k3=%e,k4=%e,k5=%e,k6=%e, \nToc_base=%f, Toc_anomaly=%f, Toc_inCelsius=%f, Toc=%f, Soc_base=%f, Soc=%f, basalmeltrate_shelf=%e \n   Aborting... \n", i,j, roundBasins(i,j), BOXMODELmask(i,j), (*ice_thickness)(i,j), T_star(i,j), area_GLbox,area_CFbox,mean_salinity_in_GLbox,mean_meltrate_in_GLbox,mean_overturning_in_GLbox,k1,k2,k3,k4,k5,k6, Toc_base(i,j), Toc_anomaly(i,j), Toc_inCelsius(i,j), Toc(i,j), Soc_base(i,j), Soc(i,j), basalmeltrate_shelf(i,j)); CHKERRQ(ierr);
	    	  PISMEnd();
	  	  }
      } // NOTE NO else-case, since  basalMeltRateForGroundingLineBox() and basalMeltRateForOtherShelves() cover all other cases and we would overwrite those results here.
    } // end j
  } // end i

  ierr = ice_thickness->end_access(); CHKERRQ(ierr);
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
  ierr = verbPrintf(4, grid.com,"B3 : in bm others rountine\n"); CHKERRQ(ierr);

  const PetscScalar rhoi = config.get("ice_density");
  const PetscScalar rhow = config.get("sea_water_density");
  const PetscScalar latentHeat = config.get("water_latent_heat_fusion");
  const PetscScalar c_p_ocean	 = 3974.0;       // J/(K*kg), specific heat capacity of ocean mixed layer
  const PetscScalar gamma_T = 0.178567865873;  // FIXME!!!! (Wrong value!) FIXME config
  const PetscScalar meltFactor = 0.002; // FIXME!!!! (Wrong value!) FIXME config

  ierr = ice_thickness->begin_access(); CHKERRQ(ierr);
  ierr = Toc_base.begin_access(); CHKERRQ(ierr);
  ierr = Toc_anomaly.begin_access(); CHKERRQ(ierr);
  ierr = Toc_inCelsius.begin_access(); CHKERRQ(ierr);
  ierr = Toc.begin_access(); CHKERRQ(ierr);
  ierr = overturning.begin_access(); CHKERRQ(ierr);
  ierr = basalmeltrate_shelf.begin_access(); CHKERRQ(ierr); // NOTE meltrate has units:   J m-2 s-1 / (J kg-1 * kg m-3) = m s-1


  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {

      PetscInt shelf_id = roundBasins(i,j);

      if (shelf_id > 0.0){
        Toc(i,j) = 273.15 + Toc_inCelsius(i,j) + Toc_anomaly(i,j);} // in K  // FIXME I think Toc should not occur in any of the routines before!

      else if (shelf_id == shelf_unidentified){
		    Toc(i,j) = Toc_base(i,j) + Toc_anomaly(i,j); // in K, NOTE: Toc_base is already in K, so no (+273.15)
		    // default: compute the melt rate from the temperature field according to beckmann_goosse03 (see below)
		    const PetscScalar shelfbaseelev = - (rhoi / rhow) * (*ice_thickness)(i,j);
	      const PetscScalar c_p_ocean = 3974.0;       // J/(K*kg), specific heat capacity of ocean mixed layer
        const PetscScalar gamma_T   = 1e-4;     // m/s, thermal exchange velocity
        PetscScalar T_f = 273.15 + (0.0939 -0.057*35.0 + 7.64e-4* shelfbaseelev); // add 273.15 to get it in Kelvin... 35 is the salinity
		    heatflux(i,j) = meltFactor * rhow * c_p_ocean * gamma_T * (Toc(i,j) - T_f);  // in W/m^2
		    basalmeltrate_shelf(i,j) = heatflux(i,j) / (latentHeat * rhoi); // in m s-1

      } else if (shelf_id == noshelf) {
		    basalmeltrate_shelf(i,j) = 0.0;

      } else { // This must not happen, since SHELFmask needs to be one of the abovenoshelf
		    ierr = verbPrintf(2, grid.com,"PISM_ERROR: [rank %d] at %d, %d  -- basins(i,j)=%f causes problems.\n   Aborting... \n",grid.rank, i, j, roundBasins(i,j)); CHKERRQ(ierr);
		    PISMEnd();
      }
    } // end j
  } // end i


  //ierr = DRAINAGEmask.end_access(); CHKERRQ(ierr);
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
  PetscErrorCode ierr = temp.copy_to(result); CHKERRQ(ierr);
  return 0;
}

//! \brief Computes mass flux in ice-equivalent m s-1.
/*!
 * Assumes that mass flux is proportional to the shelf-base heat flux.
 */
PetscErrorCode POoceanboxmodel::shelf_base_mass_flux(IceModelVec2S &result) {
  PetscErrorCode ierr = mass_flux.copy_to(result); CHKERRQ(ierr);
  return 0;
}


//FIXME Included again by Ronja to write the variables to extra file. Is there a smarter way?? 
PetscErrorCode POoceanboxmodel::write_variables(set<string> vars, string filename) {
  PetscErrorCode ierr;

  if (set_contains(vars, "BOXMODELmask")) {  ierr = BOXMODELmask.write(filename.c_str()); CHKERRQ(ierr);  }
  if (set_contains(vars, "OCEANMEANmask")) {  ierr = OCEANMEANmask.write(filename.c_str()); CHKERRQ(ierr);  }
  if (set_contains(vars, "Soc")) {  ierr = Soc.write(filename.c_str()); CHKERRQ(ierr); }
  if (set_contains(vars, "Soc_base")) { ierr = Soc_base.write(filename.c_str()); CHKERRQ(ierr); }
  if (set_contains(vars, "Toc")) { ierr = Toc.write(filename.c_str()); CHKERRQ(ierr); }
  if (set_contains(vars, "Toc_base")) {  ierr = Toc_base.write(filename.c_str()); CHKERRQ(ierr);  }
  if (set_contains(vars, "Toc_inCelsius")) { ierr = Toc_inCelsius.write(filename.c_str()); CHKERRQ(ierr);  }
  if (set_contains(vars, "T_star")) {  ierr = T_star.write(filename.c_str()); CHKERRQ(ierr);  }
  if (set_contains(vars, "Toc_anomaly")) {  ierr = Toc_anomaly.write(filename.c_str()); CHKERRQ(ierr); }
  if (set_contains(vars, "overturning")) {  ierr = overturning.write(filename.c_str()); CHKERRQ(ierr); }
  if (set_contains(vars, "heatflux")) {  ierr = heatflux.write(filename.c_str()); CHKERRQ(ierr);  }
  if (set_contains(vars, "basalmeltrate_shelf")) { ierr = basalmeltrate_shelf.write(filename.c_str()); CHKERRQ(ierr); }
  if (set_contains(vars, "CHECKmask")) { ierr = CHECKmask.write(filename.c_str()); CHKERRQ(ierr); } //FIXME delete
  if (exicerises_set) {
    if (set_contains(vars, "ICERISESmask")) {  ierr = ICERISESmask.write(filename.c_str()); CHKERRQ(ierr);  }
  }

  return 0;
}

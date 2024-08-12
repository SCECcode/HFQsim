#include <vector>//Vector hearder, defines the vector container class
#include <iostream>//Header that defines the standard input/output stream objects, e.g. cin, cout, cerr
#include <fstream>
#include <string.h>//functions to manipulate C strings and arrays
#include <stdio.h>//file related input/output
#include <algorithm>//vector copy function, binary search
#include <stdlib.h>// srand, rand
#include <numeric> // std::iota
#include <iterator> // std::distance
#include <omp.h> // openmp

#include <math.h> 
#include "physical_constants.h"

#include "input_parameters.h" //Input value 

#include "fileio.h"
#include "file_array.h"
#include "temperature_profile.h"
#include "strength_profile.h"
#include "creep_law.h"
#include "side_creep.h"
#include "stiffness_matrix.h"
#include "fault.h"
#include "earthquake.h"

#include "romberg_integration.h"
#include "fault_cooling.h"
#include "activ_energy.h"

int Romberg::JMAX = 20;
int Romberg::K = 5;
double Romberg::EPS = 1.0e-10;

/* PREM density of crustal rock kg/km^3*/
double FaultCoolingModel::density = 2.6*Units::g/(Units::cm*Units::cm*Units::cm);

/* Fault Width */
double FaultCoolingModel::faultWidth = 5.0*Units::cm;   // 10 cm was initially used

/* Diffusivity of rock, can't remember origin */
double FaultCoolingModel::diffusivity = 1.0E-2*Units::cm*Units::cm/Units::second;    // 1.0e-2 is the physical value. 

/* Specific heat capacity for granite*/
double FaultCoolingModel::specHeat = 790*Units::J/(Units::kg*Units::K);

/* Distance away from fault plane, where temperature remains the background temperature */
double FaultCoolingModel::cooling_distance = 100.0*Units::m;

/** Make operator for multiplying a vector by constant scalar */
std::vector<double> operator*( const std::vector<double>& v, double alpha)
{
  std::vector<double> v2 = v;
  for( unsigned int i = 0; i<v2.size(); i++ ) v2[i] *= alpha; 
  return v2;
}

//-------------------------------------------------------------------------
int main(int argc, char * argv[]) //argc -> argument count, char *argv[] -> argument vector
{

  std::vector<int> saveTime(int(In::maxTime/Units::year));
  std::iota(saveTime.begin(),saveTime.end(),1);
  for( unsigned int i = 0; i<int(In::maxTime/Units::year); i++){
    saveTime[i] = int(saveTime[i] * 2000);
  }

  std::cerr << "Program: " << argv[0] << " No. Args: " << argc << std::endl;

  /*Background Temperatre*/
  BackgroundTemperatureProfile T( In::Tsurface, In::Tgradient);

  /*Static strength profile*/
  int nCells = In::nL*In::nD; //128*32
  std::vector<double> vector_fs(nCells,In::fs);
  //if( !txtfile2faultvec( vector_fs, In::ifile_frictionCoef, In::nL, In::nD)){
  //  return -1;
  //}
  StaticStrengthProfile StatStrength(In::tau0, vector_fs, In::dsigmadz, In::nL, In::nD, In::cellHeight);
  std::vector<double> taus(nCells,0.0);
  StatStrength.getStaticStrength(taus); 

  /*Dynamic stress drop*/
  std::vector<double> strengthDrops (nCells,In::dynStrengthDrop); // τs-τa
  // The original files have 128*32 points. Need to be modified if higher resolution is applied.
  if( !txtfile2faultvec( strengthDrops, In::ifile_strengthdrops, In::nL, In::nD ) ){   
    return -1;
  }

  /* Initial slip deficit */
  std::vector<double> slipDeficit(nCells, 0.0);
  if( !txtfile2faultvec( slipDeficit, In::ifile_initSlipDef, In::nL, In::nD) ){
    return -1;
  }

  /*Stiffness matrix*/
  StiffnessMatrix K(In::cellLength, In::nL, In::cellHeight, In::nD, In::rigidity); //cellLength = 70 km/nL, cellHeight = 17.5 km/nD, rigidity of the half space 90 GPa

  FaultActivationEnergy FAE;

  /* set up creep parameters */
  std::vector<double> faultE(nCells, In::activationEnergy);
  std::vector<double> faultA(nCells, In::arrhAmplitude);
  std::vector<double> faultn(nCells, In::stressExponent );
  std::vector<double> faultdE(nCells); // activation energy after EQ

  /* setup a linear activation energy gradient with depth until zBD*/
  FAE.setAE_linear( faultE, In::surfaceE, In::activationEnergy, In::zBD,In::nL, In::nD, In::cellHeight);

  std::vector<double> strengthProfile(In::nD);
  std::vector<double> dosMatrix(nCells);
  for( int j=0; j<In::nD; j++){
    double depth = ((double)j+0.5)*In::cellHeight;
    for( int i=0; i<In::nL; i++){
      int iCell = i*In::nD+j;
      dosMatrix[iCell] = In::dosCoef;
      CreepLaw C0( faultA[iCell], faultE[iCell], faultn[iCell] );
      strengthProfile[j] = std::min( taus[i*In::nD+j],
                                     C0.stress( In::loadingStrainRate, T(depth) ) );
    }
  }

  /* Get the  background temperature */
  std::vector<double> faultBkgdT = T.faultgrid_co( In::nL, In::nD, In::cellHeight );
  /* Side creep*/
  sidecreepmask( faultA, faultE, faultn, faultBkgdT,
		                      In::nL, In::nD, In::cellLength, In::xBD);
  /* add random fluctuation to activation energy */
  for( int j = 0; j<In::nD; j++){
    for( int i = 0; i<In::nL; i++){
      int fluc = rand() % 5;
      int iCell = i*In::nD+j;
      if( fluc == 0 ){
	    faultE[iCell] = 0.9 * faultE[iCell];
      }
      else if( fluc == 1 ){
	    faultE[iCell] = 1.1 * faultE[iCell];
      }
      else{
        faultE[iCell] = faultE[iCell];
      }
      faultdE[iCell] = 0.8 * faultE[iCell];
    }
  }  

  std::vector<double> bkgdT(nCells), initStress(nCells), taua(nCells), faultE0(nCells);
    
  /* container for the time of last EQ at each cell */
  std::vector<double> fault_EQtime(nCells, 0.0);

  /* initial strength profile */
  for( int j=0; j<In::nD; j++ ){
    double depth = ((double)j + 0.5)*In::cellHeight;
    for( int i=0; i<In::nL; i++ ){
      int iCell = i*In::nD+j;
	  taua[iCell] = taus[iCell] - strengthDrops[iCell]; // arrest stress
      faultE0[iCell] = faultE[iCell];
      CreepLaw C( faultA[iCell], faultE[iCell], faultn[iCell] );
	  bkgdT[iCell] = T(depth);
	  initStress[iCell] = std::min( taus[iCell],C.stress( In::loadingStrainRate, bkgdT[iCell] ));
	}
  }

  /* make healing mechanism stronger with depth, proportional to normal stress, therefore depth. */
  std::vector<double> tauHeal(nCells, 0.0);
  // std::vector<double> tauMin(nCells, In::tau_min);
  for( unsigned int j=0; j<In::nD; j++){
    //double depth = ((double)j + 0.5)*In::cellHeight;
    for(unsigned int i=0; i<In::nL; i++){
      tauHeal[i*In::nD+j] = In::tau_max - (In::tau_max-In::tau_min)*(double)j/(double)In::nD; // depth has unit of km
    }
  }

  /* Save the initial temperature */
  const std::string ofilename_faultBkgdT(argv[0] + In::ofilesuffix_faultBkgdT );
  faultvec2txtfile( ofilename_faultBkgdT, faultBkgdT*(1.0/Units::K), In::nL, In::nD );

  /* Save the Arrhenius amplitude and activation energy */
  const std::string ofilename_Arrheamp( argv[0] + In::ofilesuffix_Arrheamp );
  faultvec2txtfile( ofilename_Arrheamp, faultA*(1.0/( pow(Units::Pa, 3) * Units::second)), In::nL, In::nD );

  const std::string ofilename_activE( argv[0] + In::ofilesuffix_activE );
  faultvec2txtfile( ofilename_activE, faultE*(1.0/Units::kJ/Units::mole), In::nL, In::nD );

  /* Save the dosCoef distribution */
  const std::string ofilename_DOS( argv[0] + In::ofilesuffix_DOS );
  faultvec2txtfile( ofilename_DOS, dosMatrix, In::nL, In::nD);

  /* save tau_max distribution */
  const std::string ofilename_tau_heal( argv[0] + In::ofilesuffix_tau_heal );
  faultvec2txtfile( ofilename_tau_heal, tauHeal*(1.0/Units::year), In::nL, In::nD );

  /* Save the initial stress */
  const std::string ofilename_initStress( argv[0] + In::ofilesuffix_initStress );
  faultvec2txtfile( ofilename_initStress, initStress*(1.0/Units::MPa), In::nL, In::nD );

  /* dynamic strength */
  std::vector<double> taud(nCells,0.0);
  for( int k=0; k<nCells; k++){
    taud[k] = taus[k] - ( taus[k] - taua[k] ) / dosMatrix[k];
  }
  
  /* zero initial stress, calculate stress only from slip deficit */
  for( unsigned int k=0; k<nCells; k++){
    initStress[k] = 0.0 * In::tau0;
  }

  /** Initialize the fault */
  Fault F( In::nL, In::nD, In::cellLength, In::cellHeight, In::faultWidth,
           taus, taua, taud, dosMatrix, faultA, faultn, K, initStress );
  
  F.setSlipDeficit(slipDeficit);

  /* Set containers for storing the fault info */
  std::vector<double> stress(nCells), creepVel(nCells);
  F.getStress( stress ); 

  /* Containers for recording the slip */
  std::vector<double> totalCreepSlip( nCells, 0.0), totalEqkSlip( nCells, 0.0);

  /* Containers for recording the temperature on the fault. Each loop T_start = bkgdT + fault_dT */
  std::vector<double> T_start(nCells);
  
  /*starting time of each event*/
  std::vector<double> tHist;  

  /* slip deficit on the fault at each step  */
  std::vector<double> slipdefHist;

  /*duration of each event*/
  std::vector<double> tStep;
          
  //std::vector<double> qHist(2*In::nTimeMax, 0.0);
  std::vector<double> qHist;

  /*container for recording dT on the fault*/
  std::vector<double> fault_dT(nCells, 0.0);

  /*dT on the fault from each time step*/
  std::vector<double> hist_dT;

  /*stress on the fault at the end of each step*/
  std::vector<double> hist_tau;

  /*creep velocity on the fault at each step*/
  std::vector<double> hist_vc;

  /*heat rate at each cell*/
  std::vector<double> dqdt(nCells, 0.0);

  /*record each step as either EQ or creep slip*/
  std::vector<int> slip_type; //1 for creep, 2 for EQ

  /*record slip on the fault at each time step*/
  std::vector<double> faultSlip;

  /* activation energy evolution on the fault at each step */
  std::vector<double> histActE;
  
  /* Container of stress on the fault before loading */
  std::vector<double> stressPre(nCells,0.0);
  
  /* Calculate the stress */
  F.getStress( stress );

  /* Initialize the cooling model */
  FaultCoolingModel FCM;

  /* Print header for the earthquake catalog */
  std::cout << "Time_yr, x_km, z_km, Mag_P, Mag_W, Area_km2, StressDrop_MPa\n";

  /* Run the algorithm */
  int i = 0;
  double time = 0.0;
  /* index for saving files */
  int iSave = 0;
  while( i<In::nTimeMax && time < In::maxTime )
    {
    /* 1. apply log activation energy healing
       2. update temperature on the fault
    */
    for(unsigned int k=0; k<nCells; k++){
      T_start[k] = bkgdT[k] + fault_dT[k];
    }
    FAE.logHeal(faultE,faultdE,faultE0,fault_EQtime,tauHeal,time,In::nL,In::nD);
    histActE.insert(histActE.end(),faultE.begin(),faultE.end());
    
    /* update slip deficit */
    F.getSlipDeficit( slipDeficit );
    slipdefHist.insert(slipdefHist.end(),slipDeficit.begin(),slipDeficit.end());

    /* Get the creep velocity on the fault */
    F.getCreepVelocity( creepVel, stress, T_start, faultE);
    hist_vc.insert(hist_vc.end(),creepVel.begin(),creepVel.end());
 
    /* Compute the time to failure */
    double timeStep = F.estimateTimeToFailure( stress, creepVel, In::plateVelocity );

    /* Adjust the time step */
    if( timeStep < 0 ){
      /* Negative implies something went wrong */
        std::cerr << "ERROR: negative time step\n";
        return -1;
      } else if( timeStep > In::maxTimeStep ){
        /* Don't let it get too big or creep rates will be inaccurate */
        timeStep = In::maxTimeStep;
	std::cerr << "Max time step. Model time: " << time/Units::year  << std::endl;
      } else if( timeStep < In::minTimeStep ){
        /* Don't let it get too small or we will be waiting forever */
        timeStep = In::minTimeStep;
	std::cerr << "Min time step. " << std::endl;
      } else{
        std::cerr << "Time step: " << timeStep << std::endl;
      }

    /* Load the fault */
    F.loadFault( In::plateVelocity, creepVel, timeStep );
    /* update slip deficit */
    F.getSlipDeficit( slipDeficit );

    /*step count +1, creep slip type = 1 */
    slip_type.push_back(1);

    /* record creep slip */
    for( unsigned int k=0; k<nCells; k++ ){
      double c_slip = creepVel[k]*timeStep;
      totalCreepSlip[k] += c_slip;
      faultSlip.push_back(c_slip);
    }
    
    /* record time, time step, stress */
    tHist.push_back(time);
    tStep.push_back(timeStep);
    stressPre = stress;

    /* Get the new stress */
    F.getStress( stress );
    hist_tau.insert(hist_tau.end(),stress.begin(),stress.end());

    /* set heat generated by creep to zero */
    for( unsigned int k=0; k<nCells; k++) {
	  // dqdt[k] = 0.0*(stressPre[k] + stress[k]) * creepVel[k] * timeStep;
      // qHist.push_back(dqdt[k]);
      dqdt[k] = 0.0;
      qHist.push_back(0.0);
    }

    /* update temperature on the fault */
    FCM.finite_cooling(fault_dT, timeStep, dqdt, nCells);
    hist_dT.insert(hist_dT.end(),fault_dT.begin(),fault_dT.end());
        
    /* update model time*/
    time += timeStep;

    /* Calculate whether there are any hypocenters */
    unsigned int iHypo, jHypo;
    unsigned int nCrit = F.nCriticalCells(iHypo, jHypo, stress );
    if( nCrit > 1 ){
      std::cerr << "\nWARNING: Multiple hypocenters (" << nCrit << ") ";
    }
    if( nCrit > 0 ){
      /* Output a dot for progress to terminal */
      std::cerr <<"time "<< time/Units::year<< "  **EQ** "<<std::endl;

      F.getSlipDeficit( slipDeficit );
      slipdefHist.insert(slipdefHist.end(),slipDeficit.begin(),slipDeficit.end());

      /* Compute the slip deficit before the earthquake */
      std::vector<double> slipDeficitPre;
      F.getSlipDeficit( slipDeficitPre );

      /* Store the prior stress */
      std::vector<double> stressPre = stress;

      /* record activation energy for plotting purpose */
      histActE.insert(histActE.end(),faultE.begin(),faultE.end());

      /* Store the prior time */
      double timePre = time;
      tHist.push_back(timePre);

      /* Compute the earthquake */
      std::vector<double> eqkSlip( nCells, 0.0 );
      F.computeEarthquake( stress, time, iHypo, jHypo, In::slipVelocity, eqkSlip, fault_EQtime, In::recordTime );

      /*record the duration*/
      tStep.push_back(time-timePre);

      /* Compute the new slip deficit */
      std::vector<double> slipDeficitPost;
      F.getSlipDeficit( slipDeficitPost );

      /* Compute the earthquake properties */
      Earthquake thisEQ( iHypo, jHypo, timePre, In::cellLength*In::cellHeight,
                    slipDeficitPre, slipDeficitPost, stressPre, stress );
        
      /* record EQ slip */
      for( unsigned int k=0; k<nCells; k++ ){
        double e_slip = slipDeficitPre[k] - slipDeficitPost[k];
	    if( e_slip > 0 ){
	      totalEqkSlip[k] += e_slip;
          dqdt[k] = stress[k]*e_slip;
        }
	    else{
	      dqdt[k] = 0;
        }
        qHist.push_back(dqdt[k]);
        faultSlip.push_back(e_slip);
      }

      /* record type of slip, EQ slip type = 2 */
      slip_type.push_back(2);
      /* update temperature on the fault */
      FCM.finite_cooling(fault_dT, time-timePre, dqdt, nCells);

      /* record temperature and stress*/
        hist_dT.insert(hist_dT.end(),fault_dT.begin(),fault_dT.end());
        hist_tau.insert(hist_tau.end(),stress.begin(),stress.end());

      /* check real time max on-fault temperature */
      std::cerr <<"Current max dT on the fault, index= "
        <<std::max_element(fault_dT.begin(),fault_dT.end())-fault_dT.begin()
        <<" , dT = "<<*std::max_element(fault_dT.begin(),fault_dT.end())<<std::endl;
        
      /* record post-seismic creep velocity */
      hist_vc.insert(hist_vc.end(),creepVel.begin(),creepVel.end());
        
      /* Record the total earthquake slip, heat, and temperature increase*/
      if( i > 0 ){
        /* Output the earthquake catalog: time, x, z, magP, magW, area, stressdrop*/
        printf( "%9.6E, %6.2f, %5.2f, %4.2f, %4.2f, %6.2f, %5.2f\n",
                  timePre/Units::year,
                  (iHypo+0.5)*In::cellLength/Units::km ,
                  (jHypo+0.5)*In::cellHeight/Units::km ,
                  thisEQ.potMagnitude() ,
                  thisEQ.momMagnitude( In::rigidity ) ,
                  thisEQ.ruptureArea() / (Units::km*Units::km) ,
                  thisEQ.staticStressDrop()/Units::MPa );
        } /* end if */
        
      } /*END earthquake */
      
      /* save output files and clear memory */
      if(i > 0 && i % 30000 == 0){

	// produce some input slip deficit
	F.getSlipDeficit( slipDeficit );
        const std::string ofilename_slipdef( argv[0] + std::to_string(iSave) + In::ofilesuffix_slipdef );
        faultvec2txtfile( ofilename_slipdef, slipDeficit*(1.0/Units::m), In::nL, In::nD ); 

	if ( time/Units::year >= 300 && time/Units::year < 500){          
        /* model time */
        const std::string ofilename_timeStep( argv[0] + std::to_string(iSave) + In::ofilesuffix_timeStep );
        x_txtfile( ofilename_timeStep, tHist*(1.0/Units::year));

	/* slip deficit */        
	const std::string ofilename_slipdefHist( argv[0] + std::to_string(iSave) + In::ofilesuffix_slipdefHist );
        x_txtfile( ofilename_slipdefHist, slipdefHist*(1.0/Units::m));

	/* slip */
        const std::string ofilename_faultSlip( argv[0] + std::to_string(iSave) + In::ofilesuffix_faultSlip );
        x_txtfile( ofilename_faultSlip, faultSlip*(1.0/Units::m));

	/* stress */
        const std::string ofilename_hist_tau( argv[0] + std::to_string(iSave) + In::ofilesuffix_hist_tau );
        x_txtfile( ofilename_hist_tau, hist_tau*(1.0/Units::MPa));

	/* temperature change */
        const std::string ofilename_hist_dT( argv[0] + std::to_string(iSave) + In::ofilesuffix_hist_dT );
        x_txtfile( ofilename_hist_dT, hist_dT*(1.0/Units::K));

	/* activation energy */
        const std::string ofilename_hist_actE( argv[0] + std::to_string(iSave) + In::ofilesuffix_hist_actE );
        x_txtfile( ofilename_hist_actE, histActE*(1.0/Units::kJ/Units::mole));

	/* creep velocity */
        const std::string ofilename_hist_vc( argv[0] + std::to_string(iSave) + In::ofilesuffix_hist_vc );
        x_txtfile( ofilename_hist_vc, hist_vc*(1.0*Units::second/Units::m));
	}

	tHist.clear();
        slipdefHist.clear();
        faultSlip.clear();
        hist_tau.clear();
        hist_dT.clear();
        histActE.clear();
        hist_vc.clear();

	/* heat */
        qHist.clear();

        /* model time */
        tHist.clear();

        /* slip type */
	// const std::string ofilename_slip_type(argv[0] + std::to_string(iSave) + In::ofilesuffix_slip_type);
	// x_txtfile( ofilename_slip_type, slip_type);
        slip_type.clear();

        /* tStep */
        tStep.clear();
        
	    iSave++;
      }

      i++;

    } /*END fault algorithm */

  /* Break the line */
  std::cerr << std::endl;

  /* Report why we finished the loop */
  if( i == In::nTimeMax ) std::cerr << "Max number of iterations reached\n";
  else std::cerr << "Max time reached\n";

  /* Output files */
  /* total creep slip */
  const std::string ofilename_creepSlip( argv[0] + In::ofilesuffix_creepSlip );
  faultvec2txtfile( ofilename_creepSlip, totalCreepSlip*(1.0/Units::m), In::nL, In::nD );

  /* total EQ slip */
  const std::string ofilename_eqSlip( argv[0] + In::ofilesuffix_eqSlip );
  faultvec2txtfile( ofilename_eqSlip, totalEqkSlip*(1.0/Units::m), In::nL, In::nD );

  /* final slip deficit */
  F.getSlipDeficit( slipDeficit );
  const std::string ofilename_slipdef( argv[0] + In::ofilesuffix_slipdef );
  faultvec2txtfile( ofilename_slipdef, slipDeficit*(1.0/Units::m), In::nL, In::nD );
    
  /* slip deficit at each step */
  const std::string ofilename_slipdefHist( argv[0] + std::to_string(iSave) + In::ofilesuffix_slipdefHist );
  x_txtfile( ofilename_slipdefHist, slipdefHist*(1.0/Units::m));

  /* final stress */
  const std::string ofilename_finalStress( argv[0] + In::ofilesuffix_finalStress );
  faultvec2txtfile( ofilename_finalStress, stress*(1.0/Units::MPa), In::nL, In::nD );

  /* slip */
  const std::string ofilename_faultSlip( argv[0] + std::to_string(iSave) + In::ofilesuffix_faultSlip );
  x_txtfile( ofilename_faultSlip, faultSlip*(1.0/Units::m));
    
  /* temperature change */
  const std::string ofilename_hist_dT( argv[0] + std::to_string(iSave) + In::ofilesuffix_hist_dT );
  x_txtfile( ofilename_hist_dT, hist_dT*(1.0/Units::K));

  /* type of slip at each step */
  const std::string ofilename_slip_type( argv[0] + In::ofilesuffix_slip_type);
  x_txtfile(ofilename_slip_type, slip_type);
    
  /* starting time of each step */
  const std::string ofilename_timeStep( argv[0] + std::to_string(iSave) + In::ofilesuffix_timeStep );
  x_txtfile( ofilename_timeStep, tHist*(1.0/Units::year));
    
  /* stress at the end of each step */
  const std::string ofilename_hist_tau( argv[0] + std::to_string(iSave) + In::ofilesuffix_hist_tau );
  x_txtfile( ofilename_hist_tau, hist_tau*(1.0/Units::MPa));
    
  /* creep velocity */
  const std::string ofilename_hist_vc( argv[0] + std::to_string(iSave) + In::ofilesuffix_hist_vc );
  x_txtfile( ofilename_hist_vc, hist_vc*(1.0*Units::second/Units::m));
 
  /* activation energy */
  const std::string ofilename_hist_actE( argv[0] + std::to_string(iSave) + In::ofilesuffix_hist_actE );
  x_txtfile( ofilename_hist_actE, histActE*(1.0/Units::kJ/Units::mole));
}

/* END */         
                 












  

#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#include <vector>
#include <string.h>
#include <math.h>
#include "units.h"


//----------------------------------------
/**  Set parameters used in main function
*/
namespace In
{
  /* Fault properties */
  const int nL(256), nD(64); // number of cells along strike and down-dip

  const double cellLength( 70.0*Units::km/(double)nL ), // slip cell length
    cellHeight( 17.5*Units::km/(double)nD ), // slip cell height
    faultWidth( 1.0*Units::km ), // width of the fault zone for strain calc, 10 cm is the initial value.
    rigidity( 30*Units::GPa ); // rigidity of the halfspace, value for granite is 24 GPa. 

  /* Loading and creep parameters */
  const double plateVelocity( 35*Units::mm/Units::year ), // plate velocity
    loadingStrainRate( plateVelocity/faultWidth ) , // strain rate of fault loading
    // du_max( 100.0*Units::cm ), // maximum creep distance
    arrhAmplitude( 6.31E-24 / ( pow(Units::Pa, 3) * Units::second) ), // Amplitude in the Arrhenius relation
    activationEnergy( 130.0*Units::kJ/Units::mole), // Activation energy at fault surface = 130 initially            Hirth et al(2001): 135+/-15 for quartz. Check West.   
    surfaceE( 130.0*Units::kJ/Units::mole), // Activation energy at z = 0   
    // srufaceEdamaged( 50*Units::kJ/Units::mole), // not using#################   
    AEasp_ratio( 0.9), // scaling factor of the asperity activation energy    
    stressExponent( 3.0 ), // exponent to stress in arrhenius equation
    xBD( 7.5*Units::km ), // creeping boundary width
    zBD( 10*Units::km), // BD boundary depth 
    zShallow( 3*Units::km), //shallow depth 
    tau_ratio( 4.0 ),
    Bterm_length( 7.5*Units::km); // parameters used in BZ96

  /* healing mechanism constants */
  const double tau_max( 1.0*Units::year), // max healing time of activation energy
    tau_min( 1.0*Units::day); // min healing time
    // tau_maxBD( 0.3*Units::day), //at zBD
    // tau_minBD( 0.1*Units::day),
    // scale_factor(-0.15/log(tau_max/tau_min) ),  // scale factor for activation energy healing function, not using ################
    // tau_tHeal( 500.0 * Units::year); // full recovery time of static strength

  /* Background temperature profile */
  const double Tsurface( 278.15*Units::K ),  // Temperature at surface, 273.15 is the original value
    Tgradient( 20.0 * Units::K/Units::km ); // Change in temperature with depth

  /*  Thermal parameters */
  const double density = 2.6*Units::g/(Units::cm*Units::cm*Units::cm), // density of crustal rock
    diffusivity = 1.0E-2*Units::cm*Units::cm/Units::second, // diffusivity of granite
    specHeat = 790*Units::J/(Units::kg*Units::K), // specific heat capacity for granite
    SlipZoneWidth = 5.0*Units::cm, // fault slip zone width
    Cooling_distance = 100.0*Units::m; //normal distance from the fault plane where temperature field is constant

  /* Strength and initial stress parameters */
  const double
    tau0( 10.0 * Units::MPa ), // cohesion 6MPa is the initial value
    fs(0.75), // coefficient of friction
    dsigmadz( 18.0 * Units::MPa/Units::km ), // effective normal stress gradient
    dsigmadz_dry( 28.0 * Units::MPa/Units::km ), // contact zone, zero pore pressure 
    dynStrengthDrop( 2.0*Units::MPa ), // tau_static - tau_dynamic
    dosCoef(1.25);

  /* Distribution of friction coef on the fault. */
  // const std::string ifile_frictionCoef("../../inputs/frictionCoef75.txt"); // 0.75 +/- 0.05
  // const std::string ifile_frictionCoef("../../inputs/frictionCoef75_4plus.txt"); // 0.75 +/- 0.05, 4 times resolution
  // const std::string ifile_frictionCoef("../../inputs/frictionCoef75_16plus.txt");
  // const std::string ifile_frictionCoef("../../inputs/frictionCoef70.txt"); // 0.7 +/- 0.1 
  // const std::string ifile_frictionCoef("../../inputs/frictionCoef75_512_128.txt"); // uniform distribution [0.7,0.8], 16 times resolution
  const std::string ifile_frictionCoef("../../inputs/frictionCoef75_256_64.txt");
  // const std::string ifile_frictionCoef("../../inputs/frictionCoef75_128_32.txt");

  /* Distribution of tau_s - tau_d, should be in Pa */
  // const std::string ifile_strengthdrops("../../inputs/stressdrops_rwalk.txt");
  // const std::string ifile_strengthdrops("../../inputs/stressdrops_unif.txt");
  // const std::string ifile_strengthdrops("../../inputs/stressdrops_frac.txt");
  // const std::string ifile_strengthdrops("../../inputs/arrestUnif.txt");
  // const std::string ifile_strengthdrops("../../inputs/arrestUnif_4plus.txt");
  // const std::string ifile_strengthdrops("../../inputs/stressdrops_2to6MPa_128_32.txt");
  const std::string ifile_strengthdrops("../../inputs/stressdrops_2to6MPa_256_64.txt"); // 2 to 6 MPa, uniform distribution, 256x64 cells
  // const std::string ifile_strengthdrops("../../inputs/arrestUnif_16plus.txt");
  // const std::string ifile_strengthdrops("../../inputs/stressdrops_2to6MPa_512_128.txt"); // 2 to 6 MPa, uniform distribution, 512*128 cells

  /* Distribution of initial slip deficit */
  // const std::string ifile_initSlipDef("../../inputs/initSlipDef_2000.txt"); // slip deficit of 2000 yrs loading
  // const std::string ifile_initSlipDef("../../inputs/initSlipDef_5000_128_32.txt"); //slip deficit of 5000 yrs loading, E=130 kJ/mol
  // const std::string ifile_initSlipDef("../../inputs/initSlipDef_4plus.txt");
  // const std::string ifile_initSlipDef("../../inputs/initSlipDef_zero_512_128.txt"); 
  // const std::string ifile_initSlipDef("../../inputs/initSlipDef_zero_128_32.txt");
  const std::string ifile_initSlipDef("../../inputs/init_slipdef_2500.txt");

  /* Time taken for 1 m of slip */
  const double slipVelocity(3.0 * Units::km/Units::second);   

  /* Truncate threshold */
  const double threshold = 1.0*Units::year;

  /* Algorithm time constraints */
  const double maxTimeStep = 3.0*Units::day;
  const double minTimeStep = 0.001*Units::second;        
  const int nTimeMax = 1e9; // maximum number of time steps
  const double maxTime = 80.0*Units::year;
  const double recordTime = 0.1*Units::year; // time to start applying heating and healing
  // const double coolTime = 0.5*Units::second; // quench temperature with fluid cooling, reset dT to zero

  /* Output file names */
  const std::string ofilesuffix_initStress = "_initstress.txt";
  const std::string ofilesuffix_finalStress = "_finalstress.txt";
  const std::string ofilesuffix_creepSlip = "_creepslip.txt";
  const std::string ofilesuffix_eqSlip = "_eqkslip.txt";
  const std::string ofilesuffix_slipdef = "_finalslipdef.txt";
  const std::string ofilesuffix_faultBkgdT = "_faultBkgdT.txt";
  const std::string ofilesuffix_faultSlip = "_faultSlip.txt";
  const std::string ofilesuffix_hist_dT = "_hist_dT.txt";
  const std::string ofilesuffix_slip_type = "_slip_type.txt";
  const std::string ofilesuffix_timeStep = "_timeStep.txt";
  const std::string ofilesuffix_hist_tau = "_hist_tau.txt";
  const std::string ofilesuffix_hist_vc = "_hist_vc.txt";
  //const std::string ofilesuffix_stressavg = "_stressavg.txt";
  //const std::string ofilesuffix_stressdev = "_stressdev.txt";
  const std::string ofilesuffix_activE = "_activE.txt";
  const std::string ofilesuffix_Arrheamp = "_Arrheamp.txt";
  const std::string ofilesuffix_DOS = "_DOS.txt";
  const std::string ofilesuffix_hist_actE = "_hist_actE.txt";
  const std::string ofilesuffix_tau_heal = "_tau_heal.txt";
  const std::string ofilesuffix_slipdefHist = "_hist_slipdef.txt";
  const std::string ofilesuffix_hist_taus = "_hist_taus.txt";

  const std::string ofilesuffix_MHPStress = "_MHPstress.txt";
  const std::string ofilesuffix_MHPTemp = "_MHPtemp.txt";
  const std::string ofilesuffix_MHPVcreep = "_MHPvcreep.txt";
}

#endif /* PARAMETERS_H_ */

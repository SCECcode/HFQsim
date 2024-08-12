#ifndef FAULT_COOLING_H_
#define FAULT_COOLING_H_

#include <math.h>
#include <vector>
#include "romberg_integration.h"
//#include "units.h"

//------------------------------------------------------------------------------
/**
    A halfspace with constant heat decay properties surrounds a fault
    with constant finite width.

    Cardwell et al (1978)
 */
class FaultCoolingModel{

 public:
  /* All members are static and must be defined at in the initial part
     of the main source file.  Use static because should have only one
     decay function, and makes passing a function for integration
     easier. */
  static double faultWidth; // default 10 cm
  static double density; // default is 2.6 g/cm^3
  static double diffusivity; // default is 0.01 cm^2/s
  static double specHeat; // default is 790 J/(kg K)
  static double cooling_distance; // distance away from fault plane where temperature remains constant. 

  /* Temperature generated immediately from heat rate spike*/
  double temperature(
		     const double &heatRate); // heat rate in J/s

  /* Temperature change at time t from constant heat rate between t1 and t2*/
  double temperature(
		     const double &heatRate, // heat rate in J/s
		     const double &t1, // start time for heat rate
		     const double &t2, // end time for heat rate
		     const double &tNow ); // time to calculate temperature at

  /* Temperature decay function for generated heat  */
  double integTempDecay(
			const double &tSlipStart, // start time of heat generation
			const double &tSlipEnd, // end time of heat generation
			const double &tNow ); // measurement time

  //double tempDecay( const double &tSinceSlip); //error function 

  /*update dT on the fault*/ 
  void hist_cooling(
      const std::vector<double> &heatHist, //dqdt from all the previous events
      const std::vector<double> &timeHist, //tStart from all the previous events
      const std::vector<double> &tStep,    //duration of each event (time of faulting in Cardwell 1977)
      std::vector<double> &faultdT,  //temperature profile on the fault 
      const double &tNow, //time to calculate temperature at
      const int &nEvent,
      const int &nSkip,   //take out early events to reduce calculation time 
      const int &nCells);

  /*truncate data by calculating all the previous dT at one single point. If dT < threshold, throw out the event.
  A cell with higher heat generation should be used. 
*/

  /*overloaded heat recording funtion, 1) creep heat. 2) earthquake heat*/
 /* void heat_cal(
      std::vector<double> &creepVel, //creep velocity at each cell
      std::vector<double> &stress,   //stress at each cell at the beginning of the heating process
      std:;vector<double> &qHist,    //container for recording dqdt
      const double &nCells,
      const double &nEvent,
      const double &qSkip);

  void heat_cal(
      std::vector<double> &slipDeficitPre,
      std::vector<double> &slipDeficitPost,
      std::vector<double> &totalEqkSlip,
      std::vector<double> &stress,
      std::vector<double> &qHist,
      const double &tStep,
      const double &nCells,
      const double &nEvent,
      const double &qSkip);
*/
  int truncate(
      const std::vector<double> &timeHist,
      const double &tNow,
      int &nSkip, //total number of event skipped
      const double &threshold);
  /*
  int truncate_q(
      std::vector<double> &timeHist,
      const double &tNow,
      const int &nSkip,
      int &qSkip);
  */
  void finite_cooling( std::vector<double> &fault_dT,
                       const double &timeStep,
                       const std::vector<double> &dq,
                       const int &nCells
                       );

  double finite_cooling_cell( const double &dT,
                              const double &dq,
                              const double &duration 
                              );  
   
   

private:
  /* Romberg object for performing the integration */
  Romberg R;



  /* Temperature decay function for a unit pulse of heat, at
     t=timeSinceSlip ago */
  static double tempDecay( double tSinceSlip );

};
//------------------------------------------------------------------------------
#endif /*FAULT_COOLING_H_*/

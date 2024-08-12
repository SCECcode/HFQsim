/**
    Source file for the member functions of FaultCoolingModel
*/
#include "fault_cooling.h"

//------------------------------------------------------------------------------
/* Imediate temperature from heat rate over a certain area */         //there is another area in the volume term
double FaultCoolingModel::temperature( const double &heatRate)
{
  return 0.5*heatRate/( faultWidth*density*specHeat ); 
}

//------------------------------------------------------------------------------
/* Temperature from heat rate */
double FaultCoolingModel::temperature( const double &heatRate,
				       const double &t1, const double &t2,
				       const double &tNow )
{
  return temperature(heatRate)*integTempDecay( t1, t2, tNow );
}

//------------------------------------------------------------------------------
/* Temperature decay function for generated heat  */
double FaultCoolingModel::integTempDecay(
		      const double &tSlipStart, // start time of heat generation
		      const double &tSlipEnd, // end time of heat generation
		      const double &tNow ) // measurement time
{
  /* Check what time to return heat for */
  if( tNow <= tSlipStart){
    /*Before t=0*/
    return 0.0;
  }
  else if ( tNow < tSlipEnd){
    double tSinceStart = tNow - tSlipStart;
    return R.integrate( &tempDecay, tNow-tSinceStart, tNow);
  }
  else{
    /* After t=0 */
    //double tSinceStart = tNow - tSlipStart;
    //double tSinceEnd = fmax(0.0, tNow - tSlipEnd);
    double tSinceStart = tSlipEnd - tSlipStart;
  
    //if( fabs( tSinceEnd-tSinceStart ) < Romberg::EPS*tSinceEnd){
      /* Slip pulse */
    //  return (tSinceStart-tSinceEnd)*tempDecay( 0.5*(tSinceStart+tSinceEnd) );
    //}
    //else{
      /* Period of slip */
      //return R.integrate( &tempDecay, tSinceEnd, tSinceStart );
    return R.integrate( &tempDecay, tNow-tSinceStart, tNow);
    //}
  }
}

//------------------------------------------------------------------------------
/* Temperature decay function */
double FaultCoolingModel::tempDecay( double tSinceSlip )
{
  /* Check greater than 0.0 */
  if( fabs(tSinceSlip) >= Romberg::EPS ){
    double tmp = 0.25*faultWidth/sqrt( diffusivity*tSinceSlip );
    return erf(tmp) - erf(-tmp);
  }
  else{
    /* Case where tSinceSlip = 0.0: erf(Inf) - erf(Inf) = 1 - -1 = 2 */
    return 2.0;
  }
}
//------------------------------------------------------------------------------
// add on 10/23/2018

//Use current time to update dT on the fault. 
void FaultCoolingModel::hist_cooling(  const std::vector<double> &heatHist,    //all the previous heat pulse
                                       const std::vector<double> &timeHist,    //the starting time of each heat pulse
                                       const std::vector<double> &tStep,       //duration of each step, therefore tEnd[i] = timeHist[i]+tStep[i]
                                       std::vector<double> &faultdT,     //temperature profile on the fault
                                       const double &tNow, //current time 
                                       const int &nEvent, //number of events in the calculation 
                                       const int &nSkip, //number of steps skipped from beginning
                                       const int &nCells //total number of cells
                                       )
{
  for(unsigned int k=0; k<nCells; k++){
    double dT = 0.0; //record the total temperature change
    for(int i =nSkip; i<=nEvent-1; i++){  //loop till the current event
      //double tEnd = timeHist[i] + tStep[i];
      //double dTstep = temperature(heatHist[nCells*(i-qSkip)+k],timeHist[i],timeHist[i]+tStep[i],tNow);
      double dTstep = temperature(heatHist[nCells*i+k],0.0,tStep[i],tNow-timeHist[i]);
      dT += dTstep;
      }
    faultdT[k] = dT;
  }
  return;
}

/*truncate the data        
  1. dT goes to zero roughly after a year. Check tNow-tStart below threshold = 10,000 hr. 
*/
int FaultCoolingModel::truncate(const std::vector<double> &timeHist,
                                const double &tNow,
                                int &nSkip,  //total number of event skipped
                                const double &threshold )
{
  int low_events = 0;
  for(int i = nSkip; i<(nSkip+50); i++) //loop through the first 50 events
  {
    if((tNow-timeHist[i]) > threshold){
      low_events += 1;
    }
  }
  nSkip += low_events;
  return nSkip;
}

/*
int FaultCoolingModel::truncate_q(std::vector<double> &timeHist,
                                  const double &tNow,
                                  const int &nSkip,
                                  int &qSkip)
{
  double threshold = 36000000;
  int low_events = 0;
  for(int i = nSkip; i<(nSkip+50); i++) //loop through the first 50 events
  {
    if((tNow-timeHist[i]) > threshold){
      low_events += 1;
    }
  }
  qSkip += low_events;
  return qSkip;
}
*/
                                  
void FaultCoolingModel::finite_cooling( std::vector<double> &fault_dT,
                                        const double &timeStep,
                                        const std::vector<double> &dq,
                                        const int &nCells)
{
  #pragma omp parallel
  {
    #pragma omp for
    for( unsigned int k=0; k<nCells; k++){
      fault_dT[k] = finite_cooling_cell( fault_dT[k], dq[k], timeStep);
    }
  }
  return;
}

double FaultCoolingModel::finite_cooling_cell( const double &dT,
                                               const double &dq,
                                               const double &duration ) 
{
  return (dT + temperature(dq) ) * exp(-diffusivity / (faultWidth*cooling_distance) * duration );
}

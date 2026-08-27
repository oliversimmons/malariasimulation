/*
 * adult_mosquito_eqs_fe.h
 *
 *  Created on: 27 Aug 2026
 */

#ifndef SRC_ADULT_MOSQUITO_EQS_FE_H_
#define SRC_ADULT_MOSQUITO_EQS_FE_H_

#include <queue>
#include "timeseries.h"
#include "adult_mosquito_eqs.h"
#include "aquatic_mosquito_eqs.h"

/*
 * AdultMosquitoModelFE
 *
 * A data structure for storing equation parameters for each timestep
 * 
 * foim, mu and lagged_incubating are updated each timestep
 */
struct AdultMosquitoModelFE {
    Rcpp::XPtr<Timeseries> emergence_timeseries;
    std::deque<double> lagged_incubating; //last tau values for incubating mosquitoes
    double mu; //death rate for adult female mosquitoes
    const double tau; //extrinsic incubation period
    double foim; //force of infection towards mosquitoes
    size_t total_M;
    AdultMosquitoModelFE(Rcpp::XPtr<Timeseries> emergence_timeseries, double, double, double, double);
};

// create a system of equations for the solver
integration_function_t create_eqs(AdultMosquitoModelFE& model);

#endif /* SRC_ADULT_MOSQUITO_EQS_FE_H_ */
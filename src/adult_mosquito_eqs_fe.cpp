/*
 * adult_mosquito_eqs_fe.cpp
 *
 *  Created on: 27 Aug 2026
 */

#include <Rcpp.h>
#include "adult_mosquito_eqs_fe.h"

AdultMosquitoModelFE::AdultMosquitoModelFE(
  Rcpp::XPtr<Timeseries> emergence_timeseries,
  double mu,
  double tau,
  double incubating,
  double foim
) : emergence_timeseries(emergence_timeseries), mu(mu), tau(tau), foim(foim)
{
  for (auto i = 0u; i < tau; ++i) {
    lagged_incubating.push_back(incubating);
  }
}

integration_function_t create_eqs(AdultMosquitoModelFE& model) {
  return [&model](const state_t& x, state_t& dxdt, double t) {
    auto emergence = model.emergence_timeseries->at(t, false); 
    //set submodel total_M
    model.total_M =
      x[get_idx(AdultState::S)] +
      x[get_idx(AdultState::E)] +
      x[get_idx(AdultState::I)];
    
    //run the adult ode
    auto incubation_survival = exp(-model.mu * model.tau);
    
    dxdt[get_idx(AdultState::S)] =
      emergence //growth to adult female
    - x[get_idx(AdultState::S)] * model.foim //infections
    - x[get_idx(AdultState::S)] * model.mu; //deaths   
    
    dxdt[get_idx(AdultState::E)] =
    x[get_idx(AdultState::S)] * model.foim  //infections
    - model.lagged_incubating.front() * incubation_survival //survived incubation period
    - x[get_idx(AdultState::E)] * model.mu; // deaths
    
    dxdt[get_idx(AdultState::I)] = model.lagged_incubating.front() * incubation_survival //survived incubation period
    - x[get_idx(AdultState::I)] * model.mu; // deaths
  };
}

//[[Rcpp::export]]
Rcpp::XPtr<AdultMosquitoModelFE> create_adult_mosquito_model_fe(
    Rcpp::XPtr<Timeseries> emergence_timeseries,
    double mu,
    double tau,
    double susceptible,
    double foim
) {
  auto model = new AdultMosquitoModelFE(
    emergence_timeseries,
    mu,
    tau,
    susceptible,
    foim
  );
  return Rcpp::XPtr<AdultMosquitoModelFE>(model, true);
}

//[[Rcpp::export]]
void adult_mosquito_model_update_fe(
    Rcpp::XPtr<AdultMosquitoModelFE> model,
    double mu,
    double foim,
    double susceptible
) {
  model->mu = mu;
  model->foim = foim;
  model->lagged_incubating.push_back(susceptible * foim);
  if (model->lagged_incubating.size() > 0) {
    model->lagged_incubating.pop_front();
  }
}

//[[Rcpp::export]]
std::vector<double> adult_mosquito_model_fe_save_state(
    Rcpp::XPtr<AdultMosquitoModelFE> model
) {
  // Only the lagged_incubating needs to be saved. The rest of the model
  // state is reset at each time step by a call to update before it gets
  // used.
  return {model->lagged_incubating.begin(), model->lagged_incubating.end()};
}

//[[Rcpp::export]]
void adult_mosquito_model_fe_restore_state(
    Rcpp::XPtr<AdultMosquitoModelFE> model,
    std::vector<double> state
) {
  model->lagged_incubating.assign(state.begin(), state.end());
}

//[[Rcpp::export]]
Rcpp::XPtr<Solver> create_adult_fe_solver(
    Rcpp::XPtr<AdultMosquitoModelFE> model,
    std::vector<double> init,
    double r_tol,
    double a_tol,
    size_t max_steps
) {
  return Rcpp::XPtr<Solver>(
    new Solver(
        init,
        create_eqs(*model),
        r_tol,
        a_tol,
        max_steps
    ),
    true
  );
}
// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
using namespace Rcpp;
// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// simple example of creating two matrices and
// returning the result of an operatioon on them
//
// via the exports attribute we tell Rcpp to make this function
// available from R
//
// [[Rcpp::export]]
arma::umat get_locations(arma::sp_mat& B)
{

  // Make const iterator
  arma::sp_mat::const_iterator start = B.begin();
  arma::sp_mat::const_iterator end   = B.end();

  // Calculate number of points
  int n = std::distance(start, end);

  // Kill process if no values are found (very sparse matrix)
  if (n <= 0) { Rcpp::stop("No values found!"); }

  // Build a location storage matrix
  arma::umat locs(2, n);

  // Create a vector to store each row information in. (Row, Col)
  arma::uvec temp(2);

  // Start collecting locations
  arma::sp_mat::const_iterator it = start;
  for(int i = 0; i < n; ++i)
  {
    temp(0) = it.row();
    temp(1) = it.col();
    locs.col(i) = temp;
    ++it; // increment
  }

  return locs;
}


//This is the basic loop.
// This function outputs an updated net tension
// [[Rcpp::export]]
arma::mat all_dims_one_loc(const arma::sp_mat& ten_mat,const arma::vec &kvect,const arma::vec &Hvect,const arma::vec &dvect, 
                            const arma::mat &dzvect, const arma::mat &net_tension)
{
  const int dzvect_size = dzvect.n_cols;
  //std::cout << dzvect_size;
  arma::sp_mat temp = ten_mat;
  arma::mat result = net_tension;
  //calculate locs outside the loop
  const arma::umat locs = get_locations(temp);
  //Cycle through each dimension and and update the tension
  for (int i = 0; i < dzvect_size; ++i) {

    // the '%' sign means element wise multiplication
    temp = std::move( arma::sp_mat(locs, kvect % (Hvect-dvect) % dzvect.col(i) / Hvect) );
    //std::cout << temp <<"\n";
    result.col(i) = arma::vec(arma::sum(temp, 1));
  }

  return result;
}


// It should be noted that the algo can either be in timeshift OR noisy termination mode.
// in the case both variables are true noisy termination takes precedence
//This function does not modify the inputs
// [[Rcpp::export]]
List core_while_loop_sparse_cpp(const double &max_iter, const double &sample, const arma::mat &network_dynamics,
                                const arma::mat &elevation, const arma::uvec &non_empty_t_vect, const arma::uvec &non_empty_vect, const arma::vec &dvect,
                                const arma::mat &velocity, const arma::mat &acceleration, const arma::mat &static_force,const arma::mat &force, 
                                const arma::mat &net_tension, const double tstep, const arma::sp_mat& ten_mat, const arma::vec &kvect, const double mass,
                                const bool dynamic_reset, const double tstep_change, const arma::mat &net_force, 
                                const double coef_drag, const double static_limit, const double tol,
                                const bool timeshift, const bool noisy_termination, const bool verbose_reporting){


  //
  //Create new variables that can be changed from constant inputs
  //
  
  arma::mat acceleration_temp = acceleration;
  arma::mat elevation_temp = elevation;
  arma::mat network_dynamics_temp = network_dynamics;
  arma::mat net_force_temp;
  arma::mat net_tension_temp = net_tension;
  arma::mat static_force_temp = static_force;
  arma::sp_mat ten_mat_temp  = ten_mat;
  double tstep_temp = tstep;
  arma::mat velocity_temp = velocity;

  //
  //Create variables for internal function use
  // All temporary vectors/matrices even internal are described as such for clarity
  //

  arma::mat dzvect_temp = elevation.rows(non_empty_t_vect)-elevation.rows(non_empty_vect);
  arma::mat Hvect_temp = sqrt(pow(dvect,2) + arma::vec(arma::sum(pow(dzvect_temp, 2), 1)));
  arma::mat friction_temp =  coef_drag*velocity;
  double current_time = 0;
  double is_noisy = 0;
  bool system_stable = false;
  int sample_counter  = 0;
  double static_force_scaler = 0;

  //
  // Create the saved values in case of noisy convergence
  //

  arma::mat elevation_saved = elevation;
  arma::mat net_tension_saved = net_tension;
  arma::mat velocity_saved = velocity;
  arma::mat friction_saved =  friction_temp;
  arma::mat static_force_saved = static_force;
  arma::mat net_force_saved = net_force;
  arma::mat acceleration_saved = acceleration;

  //
  //Inside the loop variables use modify in place where possible
  //
  int i = 0;
  while((i < max_iter) && !system_stable){

    dzvect_temp = elevation_temp.rows(non_empty_t_vect)-elevation_temp.rows(non_empty_vect);

    //Hvect <- sqrt(rowSums(dzvect^2) + dvect^2)
    Hvect_temp = sqrt(pow(dvect,2) + arma::vec(arma::sum(pow(dzvect_temp, 2), 1)));


    elevation_temp = (velocity_temp * tstep_temp) + (0.5 * acceleration_temp * tstep_temp * tstep_temp) + elevation_temp; //Distance s1 = ut+0.5at^2+s0
    velocity_temp = velocity_temp + (acceleration_temp * tstep_temp); //velocity v1 = v0 +at
    static_force_temp = force + net_tension_temp;


    net_tension_temp = all_dims_one_loc(ten_mat_temp, kvect, Hvect_temp, dvect, dzvect_temp, net_tension_temp );

    friction_temp = coef_drag * velocity_temp; //friction of an object in a viscous fluid under laminar flow
    net_force_temp = static_force_temp - friction_temp; //net force
    acceleration_temp = net_force_temp/mass; //acceleration

    sample_counter = floor((i+1)/sample)*sample;
    if(i+1 == sample_counter ){

      //The value of the sum of the staic force is used three times this reduces re-calcs
      static_force_scaler = arma::accu(arma::abs(static_force_temp));
      if(verbose_reporting){

        std::cout << "Iteration "<< i +1 << " total static force " << static_force_scaler << "\n";

      }
      
      
      double dynamics_row = (i+1)/sample -1;

      if(dynamics_row > 0 && (timeshift || noisy_termination)){//
        //The convergence is noisy if the static force at time t is greater than the static force at t-1
        is_noisy = static_force_scaler  > network_dynamics_temp(dynamics_row-1,2);

      }
 
      if(timeshift && !noisy_termination){ //only activates if in timeshift mode and noisy termination is off
        if(is_noisy){
          //change the tstep
          tstep_temp = tstep_temp * tstep_change;
          std::cout << "Iteration "<< i +1 << " total static force " << static_force_scaler << 
          " greater than previous static force of "<<  network_dynamics_temp(dynamics_row-1,2) << " reducing timestep to " << tstep_temp<<"\n";

          //overwrite data with previously saved info
          //If dynamic reset is zero then all the dynamics are set to zero
          elevation_temp = elevation_saved;
          net_tension_temp = net_tension_saved;
          velocity_temp = velocity_saved * (!dynamic_reset);
          friction_temp = friction_saved * (!dynamic_reset);
          static_force_temp = static_force_saved;
          //some net force and accerlation is required. if the system uses dynamic reset then
          //only the static components are used. As a result the the below expressions have two components only
          //one of which will be used depending on the dynamic reset
          net_force_temp = net_force_saved * (!dynamic_reset) + static_force_saved * (dynamic_reset);
          acceleration_temp = acceleration_saved * (!dynamic_reset) + (net_force_temp/mass) * (dynamic_reset);

        } else {

          //as the time can change, it needs to be tracked
          //current time is not updated if the system reverts to the previous saved point
          current_time = current_time + sample * tstep_temp;


          // Should these saved variables be declared at the beggining?
          //someone who knows more about c++ needs to tell me this
          elevation_saved = elevation_temp;
          net_tension_saved = net_tension_temp;
          velocity_saved = velocity_temp;
          friction_saved = friction_temp;
          static_force_saved = static_force_temp;
          net_force_saved = net_force_temp;
          acceleration_saved = acceleration_temp;
        }

      }

      if(!timeshift){
          //as the time can change, it needs to be tracked
          //current time is not updated if the system reverts to the previous saved point
          current_time = current_time + sample * tstep_temp;

      }


      //These will be the same as the previous row if the reset has been tripped.
      network_dynamics_temp(dynamics_row, 0) = i+1; //Iteration
      network_dynamics_temp(dynamics_row, 1) = current_time; //time in seconds
      network_dynamics_temp(dynamics_row, 2) = static_force_scaler;  //static force. The force exerted on the node
      network_dynamics_temp(dynamics_row, 3) = arma::accu(arma::abs(0.5 * mass * velocity_temp / tstep_temp)); //kinetic_force #I am not sure how I justify this value
      network_dynamics_temp(dynamics_row, 4) = arma::accu( 0.5 * kvect % pow((Hvect_temp-dvect),2));     //spring potential_energy
      network_dynamics_temp(dynamics_row, 5) = arma::accu(0.5 * mass * pow(velocity_temp,2));    //kinetic_energy
      
      
      if(timeshift && !noisy_termination){
        
        // The finite check has been removed!!!!!
        //This needs to be debugged and put back in!
        ///
        system_stable = //!(network_dynamics(dynamics_row,2).is_finite())||
          (network_dynamics(dynamics_row,2) > static_limit) || //if bigger than the static limit the system is diverging and the process should be stopped
          (network_dynamics(dynamics_row,2) < tol);

      } else {

        //Checks for early termination conditions. There are three or conditions
        //1 If the static force is not a finite value, this covers NA, NaN and infinite.
        //2 The static force exceeds the static limit
        //3 The system is in the noisy zone
        //4 If the static force is less than the required tolerance then the system is stable and the process can terminate
        system_stable =// !is.finite(network_dynamics[dynamics_row,3])| 
          (network_dynamics(dynamics_row,2)>static_limit) ||
          is_noisy ||
          (network_dynamics(dynamics_row,2) < tol);
      
      }


    }
    i++;
  }

  return List::create(_["dzvect"] = dzvect_temp, _["Hvect"] = Hvect_temp,
                      _["elevation"] = elevation_temp, _["velocity"] = velocity_temp, _["static_force"] = static_force_temp,
                        _["net_tension"] = net_tension_temp, _["friction"] = friction_temp, _["net_force"] = net_force_temp, 
                        _["acceleration"] = acceleration_temp, _["network_dynamics"] = network_dynamics_temp, _["Iter"] = i+1, _["tstep"] = tstep_temp);
}



//This is the basic loop for the dense function on each dimension
// This function outputs an updated net tension
// [[Rcpp::export]]
arma::mat all_dims_one_loc_dense(const arma::mat& ten_mat,const arma::vec &kvect,const arma::vec &Hvect,
const arma::vec &dvect, const arma::mat &dzvect, const arma::mat &net_tension, const arma::uvec &non_empty_index_vec)
{
    const int dzvect_size = dzvect.n_cols;
    //std::cout << dzvect_size;
    arma::mat temp = ten_mat;
    arma::mat result = net_tension;
    //calculate locs outside the loop
    //const arma::umat locs = get_locations_dense(temp);
    //Cycle through each dimension and and update the tension
    for (int i = 0; i < dzvect_size; ++i) {

        // the '%' sign means element wise multiplication
        temp.elem(non_empty_index_vec) = kvect % (Hvect-dvect) % dzvect.col(i) / Hvect;

        //std::cout << temp <<"\n";
        result.col(i) = arma::vec(arma::sum(temp, 1));
    }

    return result;
}

//         net_tension_temp = all_dims_one_loc_dense(ten_mat_temp, kvect, Hvect_temp, dvect, dzvect_temp, net_tension_temp, non_empty_index_vec );


//This function does not modify the inputs
// [[Rcpp::export]]
    List core_while_loop_dense_cpp(const double &max_iter, const double &sample, const arma::mat &network_dynamics,
    const arma::mat &elevation, const arma::uvec &non_empty_t_vect, const arma::uvec &non_empty_vect, const arma::vec &dvect,
    const arma::mat &velocity, const arma::mat &acceleration, const arma::mat &static_force,const arma::mat &force, const arma::mat 
    &net_tension, const double tstep, const arma::mat& ten_mat, const arma::vec &kvect, const double mass, const bool dynamic_reset, 
    const double tstep_change, const arma::mat &net_force, const double coef_drag, const double static_limit, const double tol, 
    const arma::uvec &non_empty_index_vec, const bool timeshift, const bool noisy_termination, const bool verbose_reporting ){


    //
    //Create new variables that can be changed from constant inputs
    //

    arma::mat acceleration_temp = acceleration;
    arma::mat elevation_temp = elevation;
    arma::mat network_dynamics_temp = network_dynamics;
    arma::mat net_force_temp;
    arma::mat net_tension_temp = net_tension;
    arma::mat static_force_temp = static_force;
    arma::mat ten_mat_temp  = ten_mat;
    double tstep_temp = tstep;
    arma::mat velocity_temp = velocity;

    //
    //Create variables for internal function use
    // All temporary vectors/matrices even internal are described as such for clarity
    //

    arma::mat dzvect_temp = elevation.rows(non_empty_t_vect)-elevation.rows(non_empty_vect);
    arma::mat Hvect_temp = sqrt(pow(dvect,2) + arma::vec(arma::sum(pow(dzvect_temp, 2), 1)));
    arma::mat friction_temp =  coef_drag*velocity;
    double current_time = 0;
    double is_noisy = 0;
    bool system_stable = false;
    int sample_counter  = 0;
    double static_force_scaler = 0;

    //
    // Create the saved values in case of noisy convergence
    //

    arma::mat elevation_saved = elevation;
    arma::mat net_tension_saved = net_tension;
    arma::mat velocity_saved = velocity;
    arma::mat friction_saved =  friction_temp;
    arma::mat static_force_saved = static_force;
    arma::mat net_force_saved = net_force;
    arma::mat acceleration_saved = acceleration;

    //
    //Inside the loop variables use modify in place where possible
    //
    int i = 0;
    while((i < max_iter) && !system_stable){
    //std::cout << "start loop " << i <<"\n";
    dzvect_temp = elevation_temp.rows(non_empty_t_vect)-elevation_temp.rows(non_empty_vect);

    //Hvect <- sqrt(rowSums(dzvect^2) + dvect^2)
    Hvect_temp = sqrt(pow(dvect,2) + arma::vec(arma::sum(pow(dzvect_temp, 2), 1)));


    elevation_temp = (velocity_temp * tstep_temp) + (0.5 * acceleration_temp * tstep_temp * tstep_temp) + elevation_temp; //Distance s1 = ut+0.5at^2+s0
    velocity_temp = velocity_temp + (acceleration_temp * tstep_temp); //velocity v1 = v0 +at
    static_force_temp = force + net_tension_temp;


    net_tension_temp = all_dims_one_loc_dense(ten_mat_temp, kvect, Hvect_temp, dvect, dzvect_temp, net_tension_temp, non_empty_index_vec );

    friction_temp = coef_drag * velocity_temp; //friction of an object in a viscous fluid under laminar flow
    net_force_temp = static_force_temp - friction_temp; //net force
    acceleration_temp = net_force_temp/mass; //acceleration

    sample_counter = floor((i+1)/sample)*sample;
    if(i+1 == sample_counter ){

    //The value of the sum of the staic force is used three times this reduces re-calcs
    static_force_scaler = arma::accu(arma::abs(static_force_temp));
    if(verbose_reporting){

        std::cout << "Iteration "<< i +1 << " total static force " << static_force_scaler << "\n";

    }


    double dynamics_row = (i+1)/sample -1;

    if(dynamics_row > 0 && (timeshift || noisy_termination)){//
        //The convergence is noisy if the static force at time t is greater than the static force at t-1
        is_noisy = static_force_scaler  > network_dynamics_temp(dynamics_row-1,2);

    }

    if(timeshift && !noisy_termination){ //only activates if in timeshift mode and noisy termination is off
        if(is_noisy){
        //change the tstep
        tstep_temp = tstep_temp * tstep_change;
        std::cout << "Iteration "<< i +1 << " total static force " << static_force_scaler << 
        " greater than previous static force of "<<  network_dynamics_temp(dynamics_row-1,2) << " reducing timestep to " << tstep_temp<<"\n";

        //overwrite data with previously saved info
        //If dynamic reset is zero then all the dynamics are set to zero
        elevation_temp = elevation_saved;
        net_tension_temp = net_tension_saved;
        velocity_temp = velocity_saved * (!dynamic_reset);
        friction_temp = friction_saved * (!dynamic_reset);
        static_force_temp = static_force_saved;
        //some net force and accerlation is required. if the system uses dynamic reset then
        //only the static components are used. As a result the the below expressions have two components only
        //one of which will be used depending on the dynamic reset
        net_force_temp = net_force_saved * (!dynamic_reset) + static_force_saved * (dynamic_reset);
        acceleration_temp = acceleration_saved * (!dynamic_reset) + (net_force_temp/mass) * (dynamic_reset);

        } else {

        //as the time can change, it needs to be tracked
        //current time is not updated if the system reverts to the previous saved point
        current_time = current_time + sample * tstep_temp;


        // Should these saved variables be declared at the beggining?
        //someone who knows more about c++ needs to tell me this
        elevation_saved = elevation_temp;
        net_tension_saved = net_tension_temp;
        velocity_saved = velocity_temp;
        friction_saved = friction_temp;
        static_force_saved = static_force_temp;
        net_force_saved = net_force_temp;
        acceleration_saved = acceleration_temp;
        }

    }

    if(!timeshift){
        //as the time can change, it needs to be tracked
        //current time is not updated if the system reverts to the previous saved point
        current_time = current_time + sample * tstep_temp;

    }


    //These will be the same as the previous row if the reset has been tripped.
    network_dynamics_temp(dynamics_row, 0) = i+1; //Iteration
    network_dynamics_temp(dynamics_row, 1) = current_time; //time in seconds
    network_dynamics_temp(dynamics_row, 2) = static_force_scaler;  //static force. The force exerted on the node
    network_dynamics_temp(dynamics_row, 3) = arma::accu(arma::abs(0.5 * mass * velocity_temp / tstep_temp)); //kinetic_force #I am not sure how I justify this value
    network_dynamics_temp(dynamics_row, 4) = arma::accu( 0.5 * kvect % pow((Hvect_temp-dvect),2));     //spring potential_energy
    network_dynamics_temp(dynamics_row, 5) = arma::accu(0.5 * mass * pow(velocity_temp,2));    //kinetic_energy


    if(timeshift && !noisy_termination){
        
        // The finite check has been removed!!!!!
        //This needs to be debugged and put back in!
        ///
        system_stable = //!(network_dynamics(dynamics_row,2).is_finite())||
        (network_dynamics(dynamics_row,2) > static_limit) || //if bigger than the static limit the system is diverging and the process should be stopped
        (network_dynamics(dynamics_row,2) < tol);

    } else {

        //Checks for early termination conditions. There are three or conditions
        //1 If the static force is not a finite value, this covers NA, NaN and infinite.
        //2 The static force exceeds the static limit
        //3 The system is in the noisy zone
        //4 If the static force is less than the required tolerance then the system is stable and the process can terminate
        system_stable =// !is.finite(network_dynamics[dynamics_row,3])| 
        (network_dynamics(dynamics_row,2)>static_limit) ||
        is_noisy ||
        (network_dynamics(dynamics_row,2) < tol);

    }


    }
    i++;
    }

    return List::create(_["dzvect"] = dzvect_temp, _["Hvect"] = Hvect_temp,
                    _["elevation"] = elevation_temp, _["velocity"] = velocity_temp, _["static_force"] = static_force_temp,
                        _["net_tension"] = net_tension_temp, _["friction"] = friction_temp, _["net_force"] = net_force_temp, 
                        _["acceleration"] = acceleration_temp, _["network_dynamics"] = network_dynamics_temp, _["Iter"] = i+1, _["tstep"] = tstep_temp);
}

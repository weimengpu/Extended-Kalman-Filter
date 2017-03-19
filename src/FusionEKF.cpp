#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>
#include <cmath>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  /**
  TODO:
    * Finish initializing the FusionEKF.
  */
  //const float c = 0.0225;

  R_laser_ << 0.0225, 0,
        0, 0.0225;

  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  H_laser_ << 1.0, 0, 0, 0,
              0, 1.0, 0, 0;
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    auto Q = MatrixXd(4, 4);

    auto P = MatrixXd(4, 4);
	  P << 1.0, 0, 0, 0,
			    0, 1.0, 0, 0,
			    0, 0, 1000.0, 0,
			    0, 0, 0, 1000.0;

    auto F = MatrixXd(4, 4);
    F << 1.0, 0, 1.0, 0,
    	   0, 1.0, 0, 1.0,
      	 0, 0, 1.0, 0,
      	 0, 0, 0, 1.0;

    previous_timestamp_ = measurement_pack.timestamp_;

    // first measurement
    auto x_ = VectorXd(4);
    //x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      double rho = measurement_pack.raw_measurements_(0,0);
      double phi = measurement_pack.raw_measurements_(1,0);
      double rho_dot = measurement_pack.raw_measurements_(2,0);
      auto pos = PolarToCartesian(rho, phi);
      auto vel = PolarToCartesian(rho_dot, phi);
      x_ << pos(0,0), pos(1,0), vel(0,0), vel(1,0);
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      x_ << measurement_pack.raw_measurements_(0,0), measurement_pack.raw_measurements_(1,0), 0, 0;
    }

    ekf_.Init(x_, P, F, Q);

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
   */
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds

  previous_timestamp_ = measurement_pack.timestamp_;

  float noise_ax = 8.5;
	float noise_ay = 8.5;

  float dt_2 = dt * dt;
	float dt_3 = dt_2 * dt;
	float dt_4 = dt_3 * dt;

	ekf_.Q_ <<  dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
			         0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,
			         dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
			         0, dt_3/2*noise_ay, 0, dt_2*noise_ay;

  ekf_.F_(0,2) = dt;
  ekf_.F_(1,3) = dt;

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    float px = ekf_.x_(0,0);
    float py = ekf_.x_(1,0);
    float vx = ekf_.x_(2,0);
    float vy = ekf_.x_(3,0);
    float eps = 1e-5;
    float rho = sqrt(px * px + py * py);
    float phi = atan2(py, px);
    float rho_dot = (px * vx + py * vy) / (eps + rho);
    VectorXd result(3);
    result << rho, phi, rho_dot;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_, result, R_radar_);
  } else {
    // Laser updates
    ekf_.Update(measurement_pack.raw_measurements_, H_laser_, R_laser_);
  }
}

/**
* A helper method to calculate Cartesian.
*/
Vector2f FusionEKF::PolarToCartesian(double rho, double phi){
  auto cart = Vector2f();
  cart << rho * cos(phi), rho * sin(phi);
  return cart;
}

#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;
#include <math.h>

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

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  H_laser_ << 1,0,0,0,
		      0,1,0,0;

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */
  ekf_.F_ = MatrixXd(4,4);
  ekf_.F_ <<  1, 0, 1, 0,
  			  0, 1, 0, 1,
  			  0, 0, 1, 0,
  			  0, 0, 0, 1;


  	//set the acceleration noise components
  	noise_ax = 9;
  	noise_ay = 9;

  	//create a 4D state vector, we don't know yet the values of the x state
	ekf_.x_ = VectorXd(4);

	//state covariance matrix P
	ekf_.P_ = MatrixXd(4, 4);
	ekf_.P_ << 1, 0, 0, 0,
			   0, 1, 0, 0,
			   0, 0, 1000, 0,
			   0, 0, 0, 1000;
	ekf_.Q_ = MatrixXd(4, 4);
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
    // first measurement
    cout << "EKF: " << endl;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
       float rho = measurement_pack.raw_measurements_[0];
       float phi = measurement_pack.raw_measurements_[1];
       ekf_.x_ << rho*cos(phi), rho*sin(phi), 0, 0;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
    }
    previous_timestamp_ = measurement_pack.timestamp_;
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
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */

    float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
	previous_timestamp_ = measurement_pack.timestamp_;
	if (dt == 0) {
	      dt = 1e-5;
	  }
	float dt_2 = dt * dt;
	float dt_3 = dt_2 * dt;
	float dt_4 = dt_3 * dt;

	//Modify the F matrix so that the time is integrated
	ekf_.F_(0, 2) = dt;
	ekf_.F_(1, 3) = dt;

	//set the process covariance matrix Q
	ekf_.Q_ <<  dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
			   0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,
			   dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
			   0, dt_3/2*noise_ay, 0, dt_2*noise_ay;

	ekf_.x_ = ekf_.F_ * ekf_.x_;
	MatrixXd Ft = ekf_.F_.transpose();
	ekf_.P_ = ekf_.F_ * ekf_.P_ * Ft + ekf_.Q_;

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */
  float pi = 3.14159;
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
	    VectorXd z = VectorXd(3);
		z << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], measurement_pack.raw_measurements_[2];
		ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
		ekf_.R_ = R_radar_;
		//VectorXd z_pred = ekf_.H_ * ekf_.x_;
		float ro = sqrt(ekf_.x_(0)*ekf_.x_(0) + ekf_.x_(1)*ekf_.x_(1));
		float phi = atan2(ekf_.x_(1), ekf_.x_(0));
		float ro_dot = (ekf_.x_(0)*ekf_.x_(2)+ekf_.x_(1)*ekf_.x_(3))/ro;
		VectorXd z_pred(3);
	    z_pred << ro, phi, ro_dot;
		VectorXd y = z - z_pred;

		if( abs(y(1)) > pi){
			if (y(1)>0)
				y(1) = y(1) - ceil((y(1)-pi)/(2*pi))*2*pi;
			else
				y(1) = y(1) + floor((y(1)+pi)/2*pi)*2*pi;
		}
		MatrixXd Ht = ekf_.H_.transpose();
		MatrixXd S = ekf_.H_ * ekf_.P_ * Ht + ekf_.R_;
		MatrixXd Si = S.inverse();
		MatrixXd PHt = ekf_.P_ * Ht;
		MatrixXd K = PHt * Si;

		//new estimate
		ekf_.x_ = ekf_.x_ + (K * y);
		long x_size = ekf_.x_.size();
		MatrixXd I = MatrixXd::Identity(x_size, x_size);
		ekf_.P_ = (I - K * ekf_.H_) * ekf_.P_;

  } else {
    // Laser updates
        VectorXd z = VectorXd(2);
        z << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1];
        ekf_.H_ = H_laser_;
        ekf_.R_ = R_laser_;
	    VectorXd z_pred = ekf_.H_ * ekf_.x_;
	  	VectorXd y = z - z_pred;
	  	MatrixXd Ht = ekf_.H_.transpose();
	  	MatrixXd S = ekf_.H_ * ekf_.P_ * Ht + ekf_.R_;
	  	MatrixXd Si = S.inverse();
	  	MatrixXd PHt = ekf_.P_ * Ht;
	  	MatrixXd K = PHt * Si;

	  	//new estimate
	  	ekf_.x_ = ekf_.x_ + (K * y);
	  	long x_size = ekf_.x_.size();
	  	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	  	ekf_.P_ = (I - K * ekf_.H_) * ekf_.P_;

  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}

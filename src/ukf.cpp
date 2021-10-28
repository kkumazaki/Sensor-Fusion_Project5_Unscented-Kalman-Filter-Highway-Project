#include "ukf.h"
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 3;
  //std_a_ = 30; // 30 is too big. Should be less than 3.

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 2*M_PI;
  //std_yawdd_ = 30; // 30 is too big. Should be less than 2Pi.
  
  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  
  /**
   * End DO NOT MODIFY section for measurement noise values 
   */
  
  /**
   * TODO: Complete the initialization. See ukf.h for other member properties.
   * Hint: one or more values initialized above might be wildly off...
   */

  is_initialized_ = false;
  time_us_ = 0;
  n_x_ = 5;
  n_aug_ = 7;
  lambda_ = 3 - n_aug_;
  Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_+1);

  // Initialize x_
  x_ << 0, 0, 0, 0, 0; // *** There may be a room to modify it.

  // Initialize P_
  P_  = MatrixXd::Identity(5, 5);

  // Initialize weights_
  weights_ = VectorXd(2*n_aug_+1);
  for (int i = 0; i < 2*n_aug_+1; i++){
    if (i == 0){
        weights_(i) = lambda_ / (lambda_+n_aug_);
    } else {
        weights_(i) = 1.0 / (2.0*(lambda_+n_aug_));
    }
  }

  // Initialize the Counter of ProcessMeasurement
  processCnt_ = 0;
  processCntHalf_= 0;
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */
  //std::cout << "-----ProcessMeasurement---------------" << std::endl;
  //std::cout << "count of process measurement : " << processCnt_ << std::endl;
  //std::cout << "meas_package.timestamp_ : " << meas_package.timestamp_ << std::endl;

  processCntHalf_ = processCnt_/2;

  /* 0. Initialize */
  if (processCntHalf_ == 0){
  //if (!is_initialized_){
    std::cout << "0. Initialize Start-----" << std::endl;
    std::cout << "timestamp [sec] " << time_us_ / 1000000 << std::endl;

    // Initialize the state x_ with the first measurement.
    // Convert radar from polar to cartesian coordinates
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      float rho = meas_package.raw_measurements_(0);
      float phi = meas_package.raw_measurements_(1);
      float rhodot = meas_package.raw_measurements_(2);

      x_(0) = rho * cos(phi);
      x_(1) = rho * sin(phi);
      x_(2) = sqrt(rhodot * rhodot);
      x_(3) = phi;
    }

    // Initialize the state x_ with the first measurement.
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {    
      x_(0) = meas_package.raw_measurements_(0);
      x_(1) = meas_package.raw_measurements_(1);
     }

    // update timestamp
    time_us_ = meas_package.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    std::cout << "0. Initialize End-----" << std::endl;
  }
  else {
    double delta_t;
    // (1)Lidar
    if ((use_laser_ == true)&&(processCntHalf_%2==1)&&(meas_package.sensor_type_ == MeasurementPackage::LASER)){
      std::cout << "----- count of process measurement : " << processCntHalf_ << std::endl;
      delta_t = (double)(meas_package.timestamp_ - time_us_) / 1000000.;
      time_us_ = meas_package.timestamp_;
      std::cout << "timestamp [sec] " << (double)time_us_ / 1000000. << std::endl;
      std::cout << "delta_t [sec] " << delta_t << std::endl;
      Prediction(delta_t);
      UpdateLidar(meas_package);
    }
    // (2)Ridar
    if ((use_radar_ == true)&&(processCntHalf_%2==0)&&(meas_package.sensor_type_ == MeasurementPackage::RADAR)){
      std::cout << "----- count of process measurement : " << processCntHalf_ << std::endl;
      delta_t = (double)(meas_package.timestamp_ - time_us_) / 1000000.;
      time_us_ = meas_package.timestamp_;
      std::cout << "timestamp [sec] " << (double)time_us_ / 1000000. << std::endl;
      std::cout << "delta_t [sec] " << delta_t << std::endl;
      Prediction(delta_t);
      UpdateRadar(meas_package);
    }

  }

  processCnt_ += 1;

}

void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */

  std::cout << "1. Predict Start-----" << std::endl;

  /* 1. Generating Sigma Points with augmentation*/
  // create augmented mean vector
  VectorXd x_aug = VectorXd(7);

  // create augmented state covariance
  MatrixXd P_aug = MatrixXd(7, 7);

  // create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
   
  // create Process Noise Covariant Matrix
  MatrixXd Q = MatrixXd(2, 2);
  Q << std_a_*std_a_,   0,
       0,       std_yawdd_*std_yawdd_;
 
  // create augmented mean state
  x_aug.head(n_x_) = x_;
  
  // create augmented covariance matrix
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug.bottomRightCorner(2, 2) = Q;

  // create square root matrix
  MatrixXd A = P_aug.llt().matrixL();

  // create augmented sigma points
  Xsig_aug.col(0) = x_aug;
  
  for (int i = 0; i < n_aug_; ++i){
      Xsig_aug.col(i+1) = x_aug + sqrt(lambda_ + n_aug_)*A.col(i);
      Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_ + n_aug_)*A.col(i);
  }

  // print result
  //std::cout << "Xsig_aug = " << std::endl << Xsig_aug << std::endl;

  /* 2. Sigma Point Prediction */
  // create matrix with predicted sigma points as columns
  //MatrixXd Xsig_pred = MatrixXd(n_x, 2 * n_aug + 1);

  //double delta_t = 0.1; // time diff in sec

  // predict sigma points
  for (int i = 0; i< 2*n_aug_+1; ++i) {
    // extract values for better readability
    double p_x = Xsig_aug(0,i);
    double p_y = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);

    // predicted state values
    double px_p, py_p;

    // avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    } else {
        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    // add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    // write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }

  // print result
  //std::cout << "Xsig_pred_ = " << std::endl << Xsig_pred_ << std::endl;

  /* 3. Predicted Mean and Covariance */
   // create vector for weights
  //VectorXd weights = VectorXd(2*n_aug+1);
  
  // create vector for predicted state
  //VectorXd x = VectorXd(n_x);

  // create covariance matrix for prediction
  //MatrixXd P = MatrixXd(n_x, n_x); 

  // set weights
  // predict state mean
  for (int i = 0; i < 2*n_aug_+1; i++){
      x_ += weights_(i) * Xsig_pred_.col(i);
  }

  // predict state covariance matrix
  for (int i = 0; i < 2*n_aug_+1; i++){
      P_ += weights_(i) * (Xsig_pred_.col(i) - x_)*(Xsig_pred_.col(i) - x_).transpose();
  }

  // print result
  //std::cout << "Predicted state" << std::endl;
  //std::cout << x_ << std::endl;
  //std::cout << "Predicted covariance matrix" << std::endl;
  //std::cout << P_ << std::endl; 

  std::cout << "2. Predict End-----" << std::endl;

}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */

  std::cout << "3a. UpdateLidar Start-----" << std::endl;

  /* 1. Predict Radar Measurement */
  // set measurement dimension, lidar can measure px and py
  int n_z = 2;

  // create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  // mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  
  // measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);

  // transform sigma points into measurement space
  for (int i = 0; i < 2*n_aug_+1; i++){
      double px = Xsig_pred_.col(i)(0);
      double py = Xsig_pred_.col(i)(1);
      
      Zsig.col(i)(0) = px;
      Zsig.col(i)(1) = py;
  }
  
  // calculate mean predicted measurement
  for (int i = 0; i < 2*n_aug_+1; i++){
      z_pred += weights_(i) * Zsig.col(i);
  }
  // calculate innovation covariance matrix S
  for (int i = 0; i < 2*n_aug_+1; i++){
      S += weights_(i) * (Zsig.col(i) - z_pred)*(Zsig.col(i) - z_pred).transpose();
  }
  S(0,0) += std_laspx_*std_laspx_;
  S(1,1) += std_laspy_*std_laspy_;

  // print result
  //std::cout << "z_pred: " << std::endl << z_pred << std::endl;
  //std::cout << "S: " << std::endl << S << std::endl;

  /* 2. UKF Update */
  // create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  // calculate cross correlation matrix
  for (int i = 0; i < 2 * n_aug_ + 1; i++){
      Tc += weights_(i) * (Xsig_pred_.col(i) - x_) * (Zsig.col(i) - z_pred).transpose();
  }

  // calculate Kalman gain K;
  MatrixXd K = MatrixXd(n_x_, n_z);
  
  K = Tc * S.inverse();

  // create vector for incoming radar measurement
  VectorXd z = VectorXd(n_z);
  z = meas_package.raw_measurements_;

  // update state mean and covariance matrix
  x_ = x_+ K*(z - z_pred);
  P_= P_ - K * S * K.transpose();

  // print result
  //std::cout << "Updated state x: " << std::endl << x_ << std::endl;
  //std::cout << "Updated state covariance P: " << std::endl << P_ << std::endl;

  std::cout << "3a. UpdateLidar End-----" << std::endl;

}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */

  std::cout << "3b. UpdateRadar Start-----" << std::endl;

  /* 1. Predict Radar Measurement */
  // set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;

  // create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  // mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  
  // measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);

  // transform sigma points into measurement space
  for (int i = 0; i < 2*n_aug_+1; i++){
      double px = Xsig_pred_.col(i)(0);
      double py = Xsig_pred_.col(i)(1);
      double v = Xsig_pred_.col(i)(2);
      double psi = Xsig_pred_.col(i)(3);
      
      Zsig.col(i)(0) = sqrt(px*px + py*py);
      Zsig.col(i)(1) = atan(py/px);
      Zsig.col(i)(2) = (px*cos(psi)+py*sin(psi))*v/sqrt(px*px + py*py);
  }
  
  // calculate mean predicted measurement
  for (int i = 0; i < 2*n_aug_+1; i++){
      z_pred += weights_(i) * Zsig.col(i);
  }
  // calculate innovation covariance matrix S
  for (int i = 0; i < 2*n_aug_+1; i++){
      S += weights_(i) * (Zsig.col(i) - z_pred)*(Zsig.col(i) - z_pred).transpose();
  }
  S(0,0) += std_radr_*std_radr_;
  S(1,1) += std_radphi_*std_radphi_;
  S(2,2) += std_radrd_*std_radrd_;

  // print result
  //std::cout << "z_pred: " << std::endl << z_pred << std::endl;
  //std::cout << "S: " << std::endl << S << std::endl;

  /* 2. UKF Update */
  // create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  // calculate cross correlation matrix
  for (int i = 0; i < 2 * n_aug_ + 1; i++){
      Tc += weights_(i) * (Xsig_pred_.col(i) - x_) * (Zsig.col(i) - z_pred).transpose();
  }

  // calculate Kalman gain K;
  MatrixXd K = MatrixXd(n_x_, n_z);
  
  K = Tc * S.inverse();

  // create vector for incoming radar measurement
  VectorXd z = VectorXd(n_z);
  z = meas_package.raw_measurements_;

  // update state mean and covariance matrix
  x_ = x_+ K*(z - z_pred);
  P_= P_ - K * S * K.transpose();

  // print result
  //std::cout << "Updated state x: " << std::endl << x_ << std::endl;
  //std::cout << "Updated state covariance P: " << std::endl << P_ << std::endl;
  std::cout << "3b. UpdateRadar End-----" << std::endl;

}
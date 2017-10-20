#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = false;

  ///* State dimension
  n_x_ = 5;

  ///* Augmented state dimension
  n_aug_ = 7;

  ///* Sigma point spreading parameter
  lambda_= 3 - n_x_;

  // initial state vector
  x_ = VectorXd(n_x_);
  x_.setZero();

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);
  P_.setIdentity();

  // calculate the number of sigma points
  n_sig_ = 2 * n_aug_ + 1;

  ///* predicted sigma points matrix
  Xsig_pred_ = MatrixXd(n_x_, n_sig_);
  Xsig_pred_.setZero();

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 10; // TODO: experiment with this value

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 20; // TODO: experiment with this value

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
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

  //create augmented mean vector
  x_aug_ = VectorXd(n_aug_);
  x_aug_.setZero();

  //create augmented state covariance
  P_aug_ = MatrixXd(n_aug_, n_aug_);
  P_aug_.setZero();

  //create augmented sigma point matrix
  Xsig_aug_ = MatrixXd(n_aug_, n_sig_);
  Xsig_aug_.setZero();

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */

  /* INITIALIZATION */
  /* Initialization structure similar to EKF Project */

  if (!is_initialized_) {
    // we can now finish initialization that is dependent on first measurement

    if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      /**
       Initialize state.
       */
      // TODO: set ekf_.x_[0] to x or raw_measurements[0]
      x_[0] = meas_package.raw_measurements_[0];
      // TODO: set ekf_.x_[1] to y or raw_measurements[1]
      x_[1] = meas_package.raw_measurements_[1];
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      /**
       Convert radar from polar to cartesian coordinates and initialize state.
       */
      // TODO: ro is raw_measurements[0], theta is raw_measurements[1]
      // TODO: set x_[0] to ro*cos(theta)
      x_[0] =
      meas_package.raw_measurements_[0] *
      cos(meas_package.raw_measurements_[1]);
      // TODO: set x_[1] to ro*sin(theta)
      x_[1] =
      meas_package.raw_measurements_[0] *
      sin(meas_package.raw_measurements_[1]);
    }

    // done initializing, no need to predict or update
    ///* time when the state is true, in us
    time_us_ = meas_package.timestamp_;
    is_initialized_ = true;
    return;
  }

  /* PREDICTION */
  /* Control structure similar to EKF Project */

  float dt = (meas_package.timestamp_ - time_us_)/1000000.0;
  if (dt > 0.001) {
    time_us_ = meas_package.timestamp_;
    Prediction(dt);
  }

  /* UPDATE */

  if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
    UpdateLidar(meas_package);
  }
  else if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    UpdateRadar(meas_package);
  }

}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  /* CREATE AUGMENTED SIGMA POINTS */
  /* Lesson 7, Section 17 */

  //create augmented mean state
  x_aug_.setZero();
  x_aug_.head(n_x_) = x_;

  //create augmented covariance matrix
  P_aug_.setZero();
  P_aug_.topLeftCorner(P_.rows(), P_.cols()) = P_;
  P_aug_(P_aug_.rows() - 2, P_aug_.cols() - 2) = std_a_ * std_a_;
  P_aug_(P_aug_.rows() - 1, P_aug_.cols() - 1) = std_yawdd_ * std_yawdd_;

  //create square root matrix
  MatrixXd A = P_aug_.llt().matrixL();

  //create augmented sigma points

  //set first column of sigma point matrix
  Xsig_aug_.setZero();
  Xsig_aug_.col(0) = x_aug_;

  // calculate square root of (lambda + n_x)
  double factor = sqrt(lambda_ + n_aug_);

  // calculate remaining sigma points
  for (int i = 0; i < n_aug_; ++i) {
    Xsig_aug_.col(i + 1) = x_aug_ + factor * A.col(i);
    Xsig_aug_.col(i + 1 + n_aug_) = x_aug_ - factor * A.col(i);
  }

  /* PREDICT SIGMA POINTS */
  /* Lesson 7, Section 20 */

  //predict sigma points
  //avoid division by zero
  //write predicted sigma points into right column
  
  for (int i = 0; i < n_sig_; ++i) {
    
    double vel = Xsig_aug_(2, i); // magnituded of velocity
    double psi = Xsig_aug_(3, i); // yaw angle
    double psi_dot = Xsig_aug_(4, i); // rate of change of yaw angle
    double nu_a = Xsig_aug_(5, i); // from Lesson 7 Section 8, longitudinal acceleration noise
    double nu_psi_dot_dot = Xsig_aug_(6, i); // from Lesson 7 Section 8, yaw acceleration noise
    double psi_dot_mult_delta_t = psi_dot * delta_t;
    double psi_plus_psi_dot_delta_t = psi + psi_dot_mult_delta_t;
    double one_half_delta_t_squared = 0.5 * (delta_t * delta_t);
    double cos_psi = cos(psi);
    double sin_psi = sin(psi);
    
    if (fabs(psi_dot) > 0.001) {
      double vel_div_psi_dot = vel / psi_dot;
      // yaw rate psi dot is not zero; use first formula
      Xsig_pred_(0, i) =
        Xsig_aug_(0, i) +
        vel_div_psi_dot * (sin(psi_plus_psi_dot_delta_t) - sin_psi) +
        one_half_delta_t_squared * cos_psi * nu_a;
      Xsig_pred_(1, i) =
        Xsig_aug_(1, i) +
        vel_div_psi_dot * (-cos(psi_plus_psi_dot_delta_t) + cos_psi) +
        one_half_delta_t_squared * sin_psi * nu_a;
      Xsig_pred_(2, i) =
        Xsig_aug_(2, i) +
        0 +
        delta_t * nu_a;
      Xsig_pred_(3, i) =
        Xsig_aug_(3, i) +
        psi_dot_mult_delta_t +
        one_half_delta_t_squared * nu_psi_dot_dot;
      Xsig_pred_(4, i) =
        Xsig_aug_(4, i) +
        0 +
        delta_t * nu_psi_dot_dot;
    } else {
      // yaw rate psi dot is zero; use second formula
      Xsig_pred_(0, i) =
        Xsig_aug_(0, i) +
        vel * cos_psi * delta_t +
        one_half_delta_t_squared * cos_psi * nu_a;
      Xsig_pred_(1, i) =
        Xsig_aug_(1, i) +
        vel * sin_psi * delta_t +
        one_half_delta_t_squared * sin_psi * nu_a;
      Xsig_pred_(2, i) =
        Xsig_aug_(2, i) +
        0 +
        delta_t * nu_a;
      Xsig_pred_(3, i) =
        Xsig_aug_(3, i) +
        psi_dot_mult_delta_t +
        one_half_delta_t_squared * nu_psi_dot_dot;
      Xsig_pred_(4, i) =
        Xsig_aug_(4, i) +
        0 +
        delta_t * nu_psi_dot_dot;
    }

    /* PREDICT MEAN AND COVARIANCE */
    /* Lesson 7, Section 23 */

    // weights are already initialized
    ///* Weights of sigma points
    weights_ = VectorXd(n_sig_);
    weights_(0) = lambda_ / (lambda_ + n_aug_);
    for (int i = 1; i < n_sig_; ++i) {
      weights_(i) = 1.0 / (2.0 * (lambda_ + n_aug_));
    }

    // predicted state mean
    x_.setZero();
    for (int i = 0; i < n_sig_; ++i) {
      x_ += weights_(i) * Xsig_pred_.col(i);
    }

    // predicted state covariance matrix

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3) > M_PI) x_diff(3) -= 2. * M_PI;
    while (x_diff(3) < -M_PI) x_diff(3) += 2. * M_PI;

    P_.setZero();
    for (int i = 0; i < n_sig_; ++i) {
      P_ += weights_(i) * x_diff * x_diff.transpose();
    }
  }
  std::cout << "x_" << endl;
  std::cout << x_ << endl;
  std::cout << "P_" << endl;
  std::cout << P_ << endl;
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */

  //measurement noise covariance
  MatrixXd R_ = MatrixXd(2, 2);
  R_ <<
    std_laspx_, 0,
    0, std_laspx_;

  //measurement matrix
  MatrixXd H_ = MatrixXd(2, n_x_);
  H_ <<
    1, 0, 0, 0, 0,
    0, 1, 0, 0, 0;

  int n_z = 2;
  VectorXd z = VectorXd(n_z);
  z = meas_package.raw_measurements_;
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */

  /* PREDICT RADAR SIGMA POINTS */

  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, n_sig_);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);

  //transform sigma points into measurement space
  for (int i = 0; i < n_sig_; ++i) {
    double px = Xsig_pred_(0, i);
    double py = Xsig_pred_(1, i);
    double vel = Xsig_pred_(2, i);
    double rho = Xsig_pred_(3, i);
    double sqrt_px2_py2 = sqrt(px * px + py * py);
    Zsig(0, i) = sqrt_px2_py2;
    Zsig(1, i) = atan(py / px);
    Zsig(2, i) = (px * cos(rho) * vel + py * sin(rho) * vel) / sqrt_px2_py2;
  }

  //calculate mean predicted measurement
  for (int i = 0; i < n_sig_; ++i) {
    z_pred += weights_(i) * Zsig.col(i);
  }

  //calculate measurement covariance matrix S
  for (int i = 0; i < n_sig_; ++i) {
    S += weights_(i) * (Zsig.col(i) - z_pred) * (Zsig.col(i) - z_pred).transpose();
  }

  MatrixXd R = MatrixXd(n_z, n_z);
  R <<
    std_radr_ * std_radr_, 0, 0,
    0, std_radphi_ * std_radphi_, 0,
    0, 0, std_radrd_ * std_radrd_;

  S = S + R;

  /* UPDATE RADAR */
  /* Lesson 7, Section 29 */

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  //calculate cross correlation matrix
  for (int i = 0; i < n_sig_; ++i) {
    Tc += weights_(i) * (Xsig_pred_.col(i) - x_) * (Zsig.col(i) - z_pred).transpose();
  }

  //calculate Kalman gain K;
  MatrixXd K = MatrixXd(n_z, n_x_);
  K = Tc * S.inverse();

  //residual
  VectorXd z = VectorXd(n_z);
  z = meas_package.raw_measurements_;
  VectorXd z_diff = z - z_pred;

  //angle normalization
  while (z_diff(1) > M_PI) z_diff(1) -= 2. * M_PI;
  while (z_diff(1) < -M_PI) z_diff(1) += 2. * M_PI;

  //update state mean and covariance matrix
  x_ += K * z_diff;
  P_ -= K * S * K.transpose();
}

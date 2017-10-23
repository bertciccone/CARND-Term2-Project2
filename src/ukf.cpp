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
  // true after initializing state vector with first available measurement
  is_initialized_ = false;

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  ///* State dimension
  n_x_ = 5;

  ///* Augmented state dimension
  n_aug_ = 7;

  // calculate the number of sigma points
  n_sig_ = 2 * n_aug_ + 1;

  // initial state vector
  x_ = VectorXd(n_x_).setZero();

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_).setIdentity();

  ///* predicted sigma points matrix
  Xsig_pred_ = MatrixXd(n_x_, n_sig_).setZero();

  ///* Sigma point spreading parameter
  lambda_= 3 - n_aug_;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 2.35; // TODO: experiment with this value

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.75; // TODO: experiment with this value

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // lidar measurement noise covariance
  R_laser_ = MatrixXd(n_z_laser_, n_z_laser_);
  R_laser_ <<
    std_laspx_ * std_laspx_, 0,
    0, std_laspy_ * std_laspy_;

  // lidar measurement matrix
  H_ = MatrixXd(n_z_laser_, n_x_);
  H_ <<
    1, 0, 0, 0, 0,
    0, 1, 0, 0, 0;

  // lidar measurement matrix transpose
  Ht = H_.transpose();

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  // radar measurement noise covariance
  R_radar_ = MatrixXd(n_z_radar_, n_z_radar_).setZero();
  R_radar_ <<
    std_radr_ * std_radr_, 0, 0,
    0, std_radphi_ * std_radphi_, 0,
    0, 0, std_radrd_ * std_radrd_;

  //create augmented mean vector
  x_aug_ = VectorXd(n_aug_).setZero();

  //create augmented state covariance
  P_aug_ = MatrixXd(n_aug_, n_aug_).setZero();

  //create augmented sigma point matrix
  Xsig_aug_ = MatrixXd(n_aug_, n_sig_).setZero();

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
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

    ///* Weights of sigma points
    weights_ = VectorXd(n_sig_);
    weights_(0) = lambda_ / (lambda_ + n_aug_);
    for (int i = 1; i < n_sig_; ++i) {
      weights_(i) = 1.0 / (2.0 * (lambda_ + n_aug_));
    }

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

  // print NIS values with spacing and caps for values above 95% threshold
  std::cout <<
    ((NIS_laser_ > 9.991) ? "        NIS LIDAR " : " nis lidar ") << NIS_laser_ <<
    ((NIS_radar_ > 7.815) ? "        NIS RADAR " : " nis radar ") << NIS_radar_ <<
    endl;

}

/**
 * Normalizes the angle phi.
 * @param {double} phi the angle that needs normalization
 */
void NormalizeAngle(double& phi)
{
  phi = atan2(sin(phi), cos(phi));
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
  const double sqrt_lambda_n_aug = sqrt(lambda_ + n_aug_);

  // calculate remaining sigma points
  for (int i = 0; i < n_aug_; ++i) {
    Xsig_aug_.col(i + 1) = x_aug_ + sqrt_lambda_n_aug * A.col(i);
    Xsig_aug_.col(i + 1 + n_aug_) = x_aug_ - sqrt_lambda_n_aug * A.col(i);
  }

  /* PREDICT SIGMA POINTS */
  /* Lesson 7, Section 20 */

  //predict sigma points
  //avoid division by zero
  //write predicted sigma points into right column
  
  for (int i = 0; i < n_sig_; ++i) {
    
    double p_x = Xsig_aug_(0,i);
    double p_y = Xsig_aug_(1,i);
    double vel = Xsig_aug_(2, i); // magnituded of velocity
    double yaw = Xsig_aug_(3, i); // yaw angle
    double yawd = Xsig_aug_(4, i); // rate of change of yaw angle
    double nu_a = Xsig_aug_(5, i); // from Lesson 7 Section 8, longitudinal acceleration noise
    double nu_yawdd = Xsig_aug_(6, i); // from Lesson 7 Section 8, yaw acceleration noise
    double yawd_dt = yawd * delta_t;
    double yaw_yawd_dt = yaw + yawd_dt;
    double half_delta_t_squared = 0.5 * (delta_t * delta_t);
    double cos_yaw = cos(yaw);
    double sin_yaw = sin(yaw);
    
    if (fabs(yawd) > 0.001) {
      double vel_div_psi_dot = vel / yawd;
      // yaw rate psi dot is not zero; use first formula
      Xsig_pred_(0, i) =
        p_x +vel_div_psi_dot * (sin(yaw_yawd_dt) - sin_yaw) +
        half_delta_t_squared * cos_yaw * nu_a;
      Xsig_pred_(1, i) =
        p_y +
        vel_div_psi_dot * (-cos(yaw_yawd_dt) + cos_yaw) +
        half_delta_t_squared * sin_yaw * nu_a;
      Xsig_pred_(2, i) = vel + 0 + delta_t * nu_a;
      Xsig_pred_(3, i) = yaw + yawd_dt + half_delta_t_squared * nu_yawdd;
      Xsig_pred_(4, i) = yawd + 0 + delta_t * nu_yawdd;
    } else {
      // yaw rate psi dot is zero; use second formula
      Xsig_pred_(0, i) =
        p_x +
        vel * cos_yaw * delta_t +
        half_delta_t_squared * cos_yaw * nu_a;
      Xsig_pred_(1, i) =
        p_y +
        vel * sin_yaw * delta_t +
        half_delta_t_squared * sin_yaw * nu_a;
      Xsig_pred_(2, i) = vel + 0 + delta_t * nu_a;
      Xsig_pred_(3, i) = yaw + yawd_dt + half_delta_t_squared * nu_yawdd;
      Xsig_pred_(4, i) = yawd + 0 + delta_t * nu_yawdd;
    }
  }

    /* PREDICT MEAN AND COVARIANCE */
    /* Lesson 7, Section 23 */

    // predicted state mean
    x_.setZero();
    for (int i = 0; i < n_sig_; ++i) {
      x_ += weights_(i) * Xsig_pred_.col(i);
    }

    // predicted state covariance matrix

    P_.setZero();
    for (int i = 0; i < n_sig_; ++i) {
      // state difference
      VectorXd x_diff = Xsig_pred_.col(i) - x_;
      //angle normalization
      NormalizeAngle(x_diff(3));
      P_ += weights_(i) * x_diff * x_diff.transpose();
    }

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

  /* Lesson 5, Section 12 */

  VectorXd z = VectorXd(n_z_laser_) = meas_package.raw_measurements_;
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd S = H_ * P_ * Ht + R_laser_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ += (K * y);
  MatrixXd I = MatrixXd(x_.size(), x_.size()).setIdentity();
  P_ -= K * H_ * P_;

  // lidar NIS calculation
  NIS_laser_ = (z - z_pred).transpose() * Si * (z - z_pred);

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
  /* Lesson 7, Section 26 */

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z_radar_, n_sig_).setZero();

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z_radar_).setZero();

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z_radar_, n_z_radar_).setZero();

  //transform sigma points into measurement space
  for (int i = 0; i < n_sig_; ++i) {
    double px = Xsig_pred_(0, i);
    double py = Xsig_pred_(1, i);
    if (px == 0 and py == 0) return; // atan2 would be undefined; skip update
    double vel = Xsig_pred_(2, i);
    double rho = Xsig_pred_(3, i);
    double sqrt_px2_py2 = sqrt(px * px + py * py);
    double v1 = cos(rho) * vel;
    double v2 = sin(rho) * vel;
    Zsig(0, i) = sqrt_px2_py2;
    Zsig(1, i) = atan2(py, px);
    Zsig(2, i) = (px * v1 + py * v2) / sqrt_px2_py2;
  }

  //calculate mean predicted measurement
  for (int i = 0; i < n_sig_; ++i) {
    z_pred += weights_(i) * Zsig.col(i);
  }

  //calculate measurement covariance matrix S
  for (int i = 0; i < n_sig_; ++i) {
    S +=
      weights_(i) * (Zsig.col(i) - z_pred) * (Zsig.col(i) - z_pred).transpose();
  }

  S += R_radar_;

  /* UPDATE RADAR */
  /* Lesson 7, Section 29 */

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z_radar_).setZero();

  //calculate cross correlation matrix
  for (int i = 0; i < n_sig_; ++i) {
    Tc +=
      weights_(i) * (Xsig_pred_.col(i) - x_) * (Zsig.col(i) - z_pred).transpose();
  }

  //calculate Kalman gain K;
  MatrixXd K = MatrixXd(n_z_radar_, n_x_).setZero();
  K = Tc * S.inverse();

  //residual
  VectorXd z = VectorXd(n_z_radar_).setZero();
  z = meas_package.raw_measurements_;
  VectorXd z_diff = z - z_pred;

  //angle normalization
  NormalizeAngle(z_diff(1));

  //update state mean and covariance matrix
  x_ += K * z_diff;
  P_ -= K * S * K.transpose();

  // radar NIS calculation
  MatrixXd Si = S.inverse();
  NIS_radar_ = (z - z_pred).transpose() * Si * (z - z_pred);

}


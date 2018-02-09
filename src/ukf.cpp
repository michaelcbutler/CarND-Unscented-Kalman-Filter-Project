#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */

// normalize angle to range [-pi, pi]
inline double normalize_range(double phi)
{
	while (phi > M_PI)
		phi -= 2*M_PI;
	while (phi < -M_PI)
		phi += 2*M_PI;
	return phi;
}


const int N_X = 5; // state dimension: <px, py, v, yaw, yawd>

UKF::UKF() {
  // need first measurement to initialize
	is_initialized_ = false;

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(N_X); // intialize at first measurement

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 6.0; 

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1.0; 
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
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
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  
  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  n_x_ = N_X;
  n_aug_ = N_X + 2; // + long and yaw accel noise
  lambda_ = 3 - n_aug_;
  const int n_sigma = 2*n_aug_ + 1;

  // covariance matrix
  P_ = MatrixXd(N_X, N_X);

  // sigma points matrix
  Xsig_pred_ = MatrixXd(n_x_, n_sigma);

  // sigma point weights
  weights_ = VectorXd(n_sigma);
  weights_(0) = lambda_/(lambda_ + n_aug_);
  const double weight = 0.5/(lambda_ + n_aug_);
  for (int i = 1; i < n_sigma; ++i) {
    weights_(i) = weight;
  }

  // measurement noise covariance matrices
  R_laser_ = MatrixXd(2, 2);
  R_laser_ << std_laspx_*std_laspx_, 0,
              0, std_laspy_*std_laspy_;
  R_radar_ = MatrixXd(3, 3);
  R_radar_ << std_radr_*std_radr_, 0, 0,
              0, std_radphi_*std_radphi_, 0,
              0, 0,std_radrd_*std_radrd_;

  lidar_nis_.open("../lidar_nis.txt");
  radar_nis_.open("../radar_nis.txt");

  if (radar_nis_.is_open() && lidar_nis_.is_open())
      cout << "nis files opened" << endl;
}

UKF::~UKF() {
    lidar_nis_.close();
    radar_nis_.close();
    cout << "nis files closed" << endl;
}

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
  if (!is_initialized_) {
    time_us_ = meas_package.timestamp_;
    elapsed_time_ = 0.0;
    radar_nis_ << "std_a_ =" << std_a_ << "\tstd_yawdd_ =" << std_yawdd_ << endl;
    lidar_nis_ << "std_a_ =" << std_a_ << "\tstd_yawdd_ =" << std_yawdd_ << endl;

    if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      double px = meas_package.raw_measurements_[0];
      double py = meas_package.raw_measurements_[1];
      double v = 0.0;
      double yaw = 0.0;
      double yawd = 0.0;
      x_ << px, py, v, yaw, yawd;

      // covariance matrix - initial value estimates
      P_.fill(0.0);
      P_(0, 0) = std_laspx_*std_laspx_;
      P_(1, 1) = std_laspy_*std_laspx_;
      P_(2, 2) = 0.5;               // 0.7 m/s std dev = 1.56 mph
      P_(3, 3) = 0.016;             // 0.126 rad std dev - 7.25 degrees
      P_(4, 4) = 0.005;             // 0.0.7 rad/s std dev - 4 degrees/s
    }
    else {
      assert(meas_package.sensor_type_ == MeasurementPackage::RADAR);
      double rho = meas_package.raw_measurements_[0];
      double phi = meas_package.raw_measurements_[1];
      double rho_dot = meas_package.raw_measurements_[2];

      phi = normalize_range(phi); // normalize phi to range [-pi, pi]

      double px = rho*cos(phi);
      double py = rho*sin(phi);

      double v = 0.0;    // TODO: improve this estimate
      double yaw = 0.0;
      double yawd = 0.0;
      x_ << px, py, v, yaw, yawd;

      // covariance matrix - initial value estimates
      P_.fill(0.0);
      P_(0, 0) = P_(1, 1) = std_radr_*std_radr_;
      P_(2, 2) = std_radrd_*std_radrd_;
      P_(3, 3) = 0.016;             // 0.126 rad std dev - 7.25 degrees
      P_(4, 4) = 0.005;             // 0.0.7 rad/s std dev - 4 degrees/s
    }

    is_initialized_ = true;

    return;
  }

  //compute the time elapsed between the current and previous measurements
	double dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
  elapsed_time_ += dt;
	time_us_ = meas_package.timestamp_;

  Prediction(dt);

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      UpdateRadar(meas_package);
    }
    else {
      assert(meas_package.sensor_type_ == MeasurementPackage::LASER);
      UpdateLidar(meas_package);
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
   const int n_sigma = 2*n_aug_ + 1;

  // create augmented mean state
  VectorXd x_aug = VectorXd(n_aug_);
  x_aug.head(5) = x_;
  x_aug(5) = 0.0;
  x_aug(6) = 0.0;

  // create augmented covariance matrix
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;

  // create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  // create augmented sigma points
  MatrixXd Xsig_aug = MatrixXd(n_aug_, n_sigma);
  Xsig_aug.col(0)  = x_aug;
  const double coef = sqrt(lambda_ + n_aug_);
  for (int i = 0; i < n_aug_; ++i)
  {
    Xsig_aug.col(i + 1)          = x_aug + coef*L.col(i);
    Xsig_aug.col(i + 1 + n_aug_) = x_aug - coef*L.col(i);
  }

  // predict sigma points
  for (int i = 0; i < n_sigma; ++i)
  {
    // extract values for better readability
    double p_x =      Xsig_aug(0, i);
    double p_y =      Xsig_aug(1, i);
    double v   =      Xsig_aug(2, i);
    double yaw =      Xsig_aug(3, i);
    double yawd =     Xsig_aug(4, i);
    double nu_a =     Xsig_aug(5, i);
    double nu_yawdd = Xsig_aug(6, i);

    // predicted state values
    double px_p, py_p;

    // avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = p_x + v/yawd*(sin(yaw + yawd*delta_t) - sin(yaw));
        py_p = p_y + v/yawd*(cos(yaw) - cos(yaw + yawd*delta_t));
    }
    else {
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
    Xsig_pred_(0, i) = px_p;
    Xsig_pred_(1, i) = py_p;
    Xsig_pred_(2, i) = v_p;
    Xsig_pred_(3, i) = yaw_p;
    Xsig_pred_(4, i) = yawd_p;
  }
  
  // predicted state mean
  x_ = weights_(0)*Xsig_pred_.col(0);
  for (int i = 1; i < n_sigma; ++i) {
    x_ = x_ + weights_(i)*Xsig_pred_.col(i);
  }

  // predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < n_sigma; ++i) {

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    x_diff(3) = normalize_range(x_diff(3));

    P_ = P_ + weights_(i)*x_diff*x_diff.transpose();
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
    
  const int n_sigma = 2*n_aug_ + 1;
  const int n_z = 2; // px, py

  VectorXd z = VectorXd(n_z);
  double px = meas_package.raw_measurements_(0);
  double py = meas_package.raw_measurements_(1);
  z << px, py;

  MatrixXd Zsig = MatrixXd(n_z, n_sigma);

  // transform sigma points into measurement space
  for (int i = 0; i < n_sigma; ++i) {

    // extract values for better readibility
    double p_x = Xsig_pred_(0, i);
    double p_y = Xsig_pred_(1, i);

    // measurement model
    Zsig(0,i) = p_x; 
    Zsig(1,i) = p_y;
  }

  // mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i=0; i < n_sigma; ++i) {
      z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  // innovation covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);
  for (int i = 0; i < n_sigma; ++i) {
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  // add measurement noise covariance matrix
  S = S + R_laser_;

  // calculate cross correlation matrix
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);
  for (int i = 0; i < n_sigma; ++i) {

    // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
 
    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  // Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  // residual
  VectorXd z_diff = z - z_pred;

  // update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();

  // compute NIS
  double error = z_diff.transpose()*S.inverse()*z_diff;
  lidar_nis_ << elapsed_time_ << "\t" << error << "\n";
  lidar_nis_.flush();
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
  
  const int n_sigma = 2*n_aug_ + 1;
  const int n_z = 3; // rho, phi, rho_dot

  VectorXd z = VectorXd(n_z);
  double rho = meas_package.raw_measurements_(0);
  double phi = meas_package.raw_measurements_(1);
  double rho_dot = meas_package.raw_measurements_(2);
  z << rho, phi, rho_dot;

  MatrixXd Zsig = MatrixXd(n_z, n_sigma);

  // transform sigma points into measurement space
  for (int i = 0; i < n_sigma; ++i) {

    // extract values for better readibility
    double p_x = Xsig_pred_(0, i);
    double p_y = Xsig_pred_(1, i);
    double v   = Xsig_pred_(2, i);
    double yaw = Xsig_pred_(3 ,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    double r = sqrt(p_x*p_x + p_y*p_y);

    // measurement model
    Zsig(0,i) = r;                                  // rho
    Zsig(1,i) = atan2(p_y,p_x);                     // phi
    Zsig(2,i) = (p_x*v1 + p_y*v2 ) / max(0.1, r);   // rho_dot // TODO: better zero-divide?
  }

  // mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i=0; i < n_sigma; ++i) {
      z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  // innovation covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);
  for (int i = 0; i < n_sigma; ++i) {
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    // angle normalization
    z_diff(1) = normalize_range(z_diff(1));

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  // add measurement noise covariance matrix
  S = S + R_radar_;

  // calculate cross correlation matrix
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);
  for (int i = 0; i < n_sigma; ++i) {

    // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    // angle normalization
    z_diff(1) = normalize_range(z_diff(1));

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // angle normalization
    x_diff(3) = normalize_range(x_diff(3));

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  // Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  // residual
  VectorXd z_diff = z - z_pred;

  // angle normalization
  z_diff(1) = normalize_range(z_diff(1));

  // update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();

  // compute NIS
  double error = z_diff.transpose()*S.inverse()*z_diff;
  radar_nis_ << elapsed_time_ << "\t" << error << "\n";
  radar_nis_.flush();

}

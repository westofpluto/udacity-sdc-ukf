#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif
#define M_2PI            (2.0*M_PI)

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  nis_radar_file.open("../NIS_radar.txt", ios::out);
  nis_laser_file.open("../NIS_laser.txt", ios::out);

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

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;
  
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

  ///* State dimension
  n_x_ = 5;

  ///* Augmented state dimension
  n_aug_ = 7;

  ///* Sigma point spreading parameter
  lambda_ = 3 - n_aug_;

  ///* augmented sigma points matrix
  Xsig_aug_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  ///* predicted sigma points matrix
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  Zsig_radar_ = MatrixXd(3, 2 * n_aug_ + 1);
  Zsig_laser_ = MatrixXd(2, 2 * n_aug_ + 1);

  // weights
  weights_ = VectorXd(2*n_aug_+1);
  //set vector for weights
  weights_(0) = lambda_/(lambda_+n_aug_);
  for (int i=1;i<2*n_aug_+1;i++) {
      weights_(i)=0.5/(lambda_+n_aug_);
  }

  //mean predicted measurement radar
  z_pred_radar_ = VectorXd(3);
  
  //measurement covariance matrix S radar
  S_radar_ = MatrixXd(3,3);

  //mean predicted measurement laser
  z_pred_laser_ = VectorXd(2);
  
  //measurement covariance matrix S laser
  S_laser_ = MatrixXd(2,2);

  R_radar_ = MatrixXd(3,3);
  R_radar_ <<    std_radr_*std_radr_, 0, 0,
          0, std_radphi_*std_radphi_, 0,
          0, 0,std_radrd_*std_radrd_;

  R_laser_ = MatrixXd(2,2);
  R_laser_ <<    std_laspx_*std_laspx_, 0, 
                 0, std_laspy_*std_laspy_;
}

UKF::~UKF() {
  std::cout << "Closing NIS files" << std::endl;
  nis_radar_file.close();
  nis_laser_file.close();
  std::cout << "Finished closing NIS files" << std::endl;
}

void UKF::AugmentedSigmaPoints() {

  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  //create augmented mean state
  x_aug.head(n_x_)=x_;
  x_aug(5)=0.0;
  x_aug(6)=0.0;
  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_, n_x_)=P_;
  P_aug(5,5)=std_a_*std_a_;
  P_aug(5,6)=0.0;
  P_aug(6,5)=0.0;
  P_aug(6,6)=std_yawdd_*std_yawdd_;
  //create square root matrix
  MatrixXd sqrtP_aug=P_aug.llt().matrixL();
  //std::cout << "sqrtP_aug = " << sqrtP_aug << std::endl;

  //create augmented sigma points
  double f=sqrt(lambda_+n_aug_);
  Xsig_aug_.col(0)=x_aug;
  for (int i = 0; i< n_aug_; i++) {
    Xsig_aug_.col(i+1)       = x_aug + f * sqrtP_aug.col(i);
    Xsig_aug_.col(i+1+n_aug_) = x_aug - f * sqrtP_aug.col(i);
  }

  //print result
  //std::cout << "Xsig_aug_ = " << std::endl << Xsig_aug_ << std::endl;

  // result is in Xsig_aug_
}

void UKF::SigmaPointPrediction(double delta_t) {

  //predict sigma points
  //avoid division by zero
  //write predicted sigma points into right column
  
  //predict sigma points
  for (int i = 0; i< 2*n_aug_+1; i++) {
    //
    // extract ith column of Xk into local vars
    //
    double p_x = Xsig_aug_(0,i);
    double p_y = Xsig_aug_(1,i);
    double v = Xsig_aug_(2,i);
    double yaw = Xsig_aug_(3,i);
    double yawd = Xsig_aug_(4,i);
    double nu_a = Xsig_aug_(5,i);
    double nu_yawdd = Xsig_aug_(6,i);

    std::cout << "p_x, p_y are " << p_x << ", " << p_y << std::endl;
    std::cout << "v is " << v << std::endl;
    std::cout << "yaw, yawd are " << yaw << ", " << yawd << std::endl;
    std::cout << "nu_a, nu_yawdd are " << nu_a << ", " << nu_yawdd << std::endl;
    std::cout << std::endl;

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > 0.00001) {
        px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    } else {
        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    //while (yaw_p > M_PI) yaw_p-=2.*M_PI;
    //while (yaw_p<-M_PI) yaw_p+=2.*M_PI;

    //write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }

  //print result
  //std::cout << "Xsig_pred_ = " << std::endl << Xsig_pred_ << std::endl;

  // result stored in Xsig_pred_
}

void UKF::PredictMeanAndCovariance() {
  
  //create vector for predicted state
  VectorXd x = VectorXd(n_x_);

  //create covariance matrix for prediction
  MatrixXd P = MatrixXd(n_x_, n_x_);

  //predict state mean
  //predict state covariance matrix

  //predicted state mean
  x.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    x = x+ weights_(i) * Xsig_pred_.col(i);
  }
  //while (x(3)> M_PI) x(3)-=2.*M_PI;
  //while (x(3)<-M_PI) x(3)+=2.*M_PI;

  //predicted state covariance matrix
  P.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

    VectorXd x_diff = Xsig_pred_.col(i) - x;
    //angle normalization
    //std::cout << "x_diff(3) is " << x_diff(3) << std::endl;
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P = P + weights_(i) * x_diff * x_diff.transpose() ;
  }

  //print result
  std::cout << "Predicted state" << std::endl;
  std::cout << x << std::endl;
  std::cout << "Predicted covariance matrix" << std::endl;
  std::cout << P << std::endl;

  //write result
  x_ = x;
  P_ = P;
}

void UKF::PredictRadarMeasurement() {

  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;

  //transform sigma points into measurement space
  Zsig_radar_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model
    Zsig_radar_(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
    Zsig_radar_(1,i) = atan2(p_y,p_x);                                 //phi
    Zsig_radar_(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
  }

  //mean predicted measurement
  z_pred_radar_.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
      z_pred_radar_ = z_pred_radar_ + weights_(i) * Zsig_radar_.col(i);
  }

  while (z_pred_radar_(1)> M_PI) z_pred_radar_(1)-=2.*M_PI;
  while (z_pred_radar_(1)<-M_PI) z_pred_radar_(1)+=2.*M_PI;

  //measurement covariance matrix S
  S_radar_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig_radar_.col(i) - z_pred_radar_;

    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S_radar_ = S_radar_ + weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  S_radar_ = S_radar_ + R_radar_;

  // result stored in z_pred_radar_ and S_radar
}

void UKF::PredictLaserMeasurement() {

  //set measurement dimension, laser measures px and py directly   
  int n_z = 2;

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);

    // measurement model
    Zsig_laser_(0,i) = p_x;
    Zsig_laser_(1,i) = p_y;
  }

  //mean predicted measurement
  z_pred_laser_.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
      z_pred_laser_ = z_pred_laser_ + weights_(i) * Zsig_laser_.col(i);
  }

  //measurement covariance matrix S
  S_laser_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig_laser_.col(i) - z_pred_laser_;

    S_laser_ = S_laser_ + weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  S_laser_ = S_laser_ + R_laser_;

  //print result
  //std::cout << "z_pred: " << std::endl << z_pred_laser_ << std::endl;
  //std::cout << "S_laser_: " << std::endl << S_laser_ << std::endl;

  // result stored in z_pred_laser_ and S_laser_
}


/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
    std::cout << "ProcessMeasurement called " << std::endl;

    if (!is_initialized_) {
        if (meas_package.sensor_type_ == MeasurementPackage::RADAR && !use_radar_) {
            return;
        } else if (meas_package.sensor_type_ == MeasurementPackage::LASER && !use_laser_) {
            return;
        }
        /**
          * Initialize the state x_ with the first measurement.
          * Create the covariance matrix.
        */
        // first measurement
        cout << "First measurement UKF: " << endl;

        if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
            /**
            Convert radar from polar to cartesian coordinates and initialize state.
            */
            float range = meas_package.raw_measurements_[0];
            float phi=meas_package.raw_measurements_[1];
            float drdt=meas_package.raw_measurements_[2];
            float px=range*cos(phi);
            float py=range*sin(phi);

            float vx = drdt*cos(phi);
            float vy = drdt*sin(phi);
            float v = sqrt(vx*vx + vy*vy);

            x_ << px, py, v, 0, 0;

        } else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
            x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
        }
        //
        // we also need to create the initial state covariance matrix.
        // This represents the initial state uncertainty and usually requires some physical intuition regarding
        // the problem at hand. In our case, the measurement error covariance components are pretty small
        // so we should have a pretty good idea of the initial positions. We can use 1.0 for the initial error
        // uncertainty there. We really don't know what the velocity is, but we can be conservative and use 50*50=2500
        // for the velocity error uncertainties.
        //
        P_ = MatrixXd(n_x_, n_x_);
        double sigsq=1.0;
        P_ << sigsq, 0, 0, 0, 0,
             0, sigsq, 0, 0, 0,
             0, 0, sigsq, 0, 0,
             0, 0, 0, sigsq, 0,
             0, 0, 0, 0, sigsq;

        previous_timestamp_ = meas_package.timestamp_;
        is_initialized_ = true;
        std::cout << "End of first measurement UKF: " << std::endl;
        return;
    }

    /************************************************************************
    ** If we are here then filter is already initialized, so do filtering
    ************************************************************************/

    //
    // PREDICTION STEP
    //
    std::cout << "start of prediction" << std::endl;
    double delta_t = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;   //delta_t - expressed in seconds
    previous_timestamp_ = meas_package.timestamp_;

    std::cout << "just before Prediction, delta_t is " << delta_t << std::endl;
    Prediction(delta_t);
    std::cout << "just after Prediction" << std::endl;

    //
    // At this point we have the predicted state x_ and the predicted state error covariance matrix P_
    // so now we can do updates from measurements
    //

    //
    // UPDATE STEP
    //
    std::cout << "start of update" << std::endl;
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
        std::cout << "RADAR" << std::endl;
        PredictRadarMeasurement();
        UpdateStateRadar(meas_package);
    } else if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) {
        std::cout << "LASER" << std::endl;
        PredictLaserMeasurement();
        UpdateStateLaser(meas_package);
    }
    std::cout << "after update" << std::endl;

    // print the output
    std::cout << "x_ = " << x_ << std::endl;
    std::cout << "P_ = " << P_ << std::endl;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
    std::cout << "Prediction0" << std::endl;
    AugmentedSigmaPoints();
    std::cout << "Prediction1" << std::endl;
    SigmaPointPrediction(delta_t);
    std::cout << "Prediction2" << std::endl;
    PredictMeanAndCovariance();
    std::cout << "Prediction3" << std::endl;
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateStateLaser(MeasurementPackage meas_package) {
    UpdateStateFromMeasurement(meas_package);
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateStateRadar(MeasurementPackage meas_package) {
    UpdateStateFromMeasurement(meas_package);
}

void UKF::UpdateStateFromMeasurement(MeasurementPackage meas_package) {
  VectorXd z_diff;
  MatrixXd S;
  MatrixXd K;

  int n_z = 0;
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      n_z = 3;
  } else {
      n_z = 2;
  }
  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    //residual
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
        z_diff = Zsig_radar_.col(i) - z_pred_radar_;
        while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
        while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
    } else {
        z_diff = Zsig_laser_.col(i) - z_pred_laser_;
    }

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  VectorXd z = meas_package.raw_measurements_;
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      S = S_radar_;
      K = Tc * S.inverse();
      z_diff = z - z_pred_radar_;
      //angle normalization
      while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
      while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
  } else {
      S = S_laser_;
      K = Tc * S.inverse();
      z_diff = z - z_pred_laser_;
  }

  std::cout << "## x_ is " << x_ << std::endl;
  std::cout << "## P_ is " << P_ << std::endl;
  std::cout << "## K is " << K << std::endl;
  std::cout << "## z_diff is " << z_diff << std::endl;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();

  while (x_(3)> M_PI) x_(3)-=2.*M_PI;
  while (x_(3)<-M_PI) x_(3)+=2.*M_PI;
  //
  // Finally compute appropriate NIS value
  //
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      nis_radar_ = z_diff.transpose()*S.inverse()*z_diff;
      nis_radar_file << nis_radar_ << std::endl;
  } else {
      nis_laser_ = z_diff.transpose()*S.inverse()*z_diff;
      nis_laser_file << nis_laser_ << std::endl;
  }
}

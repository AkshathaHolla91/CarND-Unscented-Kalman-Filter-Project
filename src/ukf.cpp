#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#include "tools.h"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  is_initialized_ = false;
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.3;
  
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
  //set state dimension
  n_x_= 5;
  // Set augmented dimension
  n_aug_ = 7;
  lambda_ = 3 - n_aug_;
  
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
  weights_ = VectorXd(2*n_aug_+1);


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
  if(!is_initialized_){
    if(meas_package.sensor_type_==MeasurementPackage::RADAR){
      float rho=meas_package.raw_measurements_(0);
      float phi=meas_package.raw_measurements_(1);
      float rho_dot=meas_package.raw_measurements_(2);

      float px=rho*cos(phi);
      float py=rho*sin(phi);

      x_<<px, py,rho_dot, 0, 0;
      P_<< 1, 0, 0, 0, 0,
          0, 1, 0, 0, 0,
          0, 0, 1, 0, 0,
          0, 0, 0, 1, 0,
          0, 0, 0, 0, 1;


    }
    else if(meas_package.sensor_type_==MeasurementPackage::LASER){
      float px=meas_package.raw_measurements_(0);
      float py=meas_package.raw_measurements_(1);
      x_<<px, py, 0 , 0, 0;
      P_<< 1, 0, 0, 0, 0,
          0, 1, 0, 0, 0,
          0, 0, 1, 0, 0,
          0, 0, 0, 1, 0,
          0, 0, 0, 0, 1;

    }
    time_us_=meas_package.timestamp_;
    is_initialized_=true;
    return;
  }
  float delta_t=(meas_package.timestamp_- time_us_)/1000000.0;
  time_us_=meas_package.timestamp_;


  Prediction(delta_t);

  if(use_radar_ && meas_package.sensor_type_==MeasurementPackage::RADAR){
    UpdateRadar(meas_package);
  }
  else if(use_laser_ && meas_package.sensor_type_==MeasurementPackage::LASER){
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
  VectorXd x_aug=VectorXd(7);
  MatrixXd P_aug=MatrixXd(7, 7);
  x_aug.head(5)=x_;
  x_aug(5)=0;
  x_aug(6)=0;

  P_aug.fill(0);
  P_aug.topLeftCorner(5, 5)=P_;
  P_aug(5,5)= std_a_*std_a_;
  P_aug(6,6)= std_yawdd_* std_yawdd_;
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  MatrixXd p_root=P_aug.llt().matrixL();
  // Augumented sigma values are generated
  Xsig_aug.col(0)=x_aug;
  for(int i=0;i<n_aug_;i++){
    Xsig_aug.col(i+1)=x_aug+sqrt(lambda_+n_aug_) * p_root.col(i);
    Xsig_aug.col(i+1+n_aug_)=x_aug-sqrt(lambda_+n_aug_) * p_root.col(i);
  }

  //predict sigma points
  for(int i=0;i<2*n_aug_+1;i++){
    double px=Xsig_aug(0,i);
    double py=Xsig_aug(1,i);
    double v=Xsig_aug(2,i);
    double yaw=Xsig_aug(3,i);
    double yawd=Xsig_aug(4,i);
    double nu_a=Xsig_aug(5,i);
    double nu_yawdd=Xsig_aug(6,i);

    double px_pred, py_pred;
    if(fabs(yawd)>0.001){
      px_pred=px+v/yawd* (sin (yaw+yawd*delta_t)- sin(yaw));
      py_pred = py + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    }
    else{
      px_pred=px+v*delta_t*cos(yaw);
      py_pred=py+v*delta_t*sin(yaw);
    }
    double v_pred=v;
    double yaw_pred=yaw+yawd*delta_t;
    double yawd_pred=yawd;

    px_pred+=0.5*nu_a*delta_t*delta_t*cos(yaw);
    py_pred+=0.5*nu_a*delta_t*delta_t*sin(yaw);
    v_pred+=nu_a*delta_t;
    yaw_pred+=0.5*nu_yawdd*delta_t*delta_t;
    yawd_pred+=nu_yawdd*delta_t;

    Xsig_pred_(0,i)=px_pred;
    Xsig_pred_(1,i)=py_pred;
    Xsig_pred_(2,i)=v_pred;
    Xsig_pred_(3,i)=yaw_pred;
    Xsig_pred_(4,i)=yawd_pred;

  }
  double weights_n=0.5/(lambda_+n_aug_);
  weights_.fill(weights_n);
  double weights_0=lambda_/(lambda_+n_aug_);
  weights_(0)=weights_0;

  x_.fill(0.0);
  for(int i=0;i<2*n_aug_+1;i++){
    x_+=weights_(i)*Xsig_pred_.col(i);
  }

  P_.fill(0.0);
  Tools tools;
  for(int i=0;i<2*n_aug_+1;i++){
    VectorXd x_diff=Xsig_pred_.col(i)-x_;
    x_diff(3)=tools.NormalizeAngle(x_diff(3));
    P_+=weights_(i)*x_diff*x_diff.transpose();
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
   Tools tools;
  int n_z=2;
  MatrixXd Zsig=MatrixXd(n_z, 2*n_aug_+1);

  for(int i=0;i<2*n_aug_+1;i++){
    double px=Xsig_pred_(0,i);
    double py=Xsig_pred_(1,i);
    
    Zsig(0,i)=px;
    Zsig(1,i)=py;
    
  }
 
  //mean predicted measurement
  VectorXd z_pred=VectorXd(n_z);
  z_pred.fill(0);
  for(int i=0;i<2*n_aug_+1;i++){
    z_pred+=weights_(i)*Zsig.col(i);
  }
 //innovation covariance matrix S
 
  MatrixXd S=MatrixXd(n_z,n_z);
  S.fill(0);
  for(int i=0;i<2*n_aug_+1;i++){
    VectorXd z_diff=Zsig.col(i)-z_pred;
    S+=weights_(i)*z_diff*z_diff.transpose();
  }
  MatrixXd R=MatrixXd(n_z,n_z);
  R << std_laspx_*std_laspx_, 0,
      0, std_laspy_*std_laspy_;
  S+=R;
//create vector for incoming radar measurement
VectorXd z=VectorXd(n_z);
z(0)=meas_package.raw_measurements_(0);
z(1)=meas_package.raw_measurements_(1);


//create matrix for cross correlation Tc
MatrixXd tc=MatrixXd(n_x_,n_z);
tc.fill(0);
for(int i=0;i<2*n_aug_+1;i++){
  VectorXd z_diff=Zsig.col(i)-z_pred;
  

  VectorXd x_diff=Xsig_pred_.col(i)-x_;
  

  tc+=weights_(i)*x_diff*z_diff.transpose();
}

MatrixXd K=tc*S.inverse();

VectorXd y=z-z_pred;


x_+=K*y;
P_-=K*S*K.transpose();
  
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
  Tools tools;
  int n_z=3;
  MatrixXd Zsig=MatrixXd(n_z, 2*n_aug_+1);

  for(int i=0;i<2*n_aug_+1;i++){
    double px=Xsig_pred_(0,i);
    double py=Xsig_pred_(1,i);
    double v=Xsig_pred_(2,i);
    double yaw=Xsig_pred_(3,i);

    double vx=v*cos(yaw);
    double vy=v*sin(yaw);
     // measurement model
    Zsig(0,i)=sqrt(px*px+py*py);
    Zsig(1,i)=atan2(py,px);
    Zsig(2,i)=(px*vx+py*vy)/Zsig(0,i);
  }
 
  //mean predicted measurement
  VectorXd z_pred=VectorXd(n_z);
  z_pred.fill(0);
  for(int i=0;i<2*n_aug_+1;i++){
    z_pred+=weights_(i)*Zsig.col(i);
  }
 //innovation covariance matrix S
 
  MatrixXd S=MatrixXd(n_z,n_z);
  S.fill(0);
  for(int i=0;i<2*n_aug_+1;i++){
    VectorXd z_diff=Zsig.col(i)-z_pred;
    z_diff(1)=tools.NormalizeAngle(z_diff(1));
    S+=weights_(i)*z_diff*z_diff.transpose();
  }
  MatrixXd R=MatrixXd(n_z,n_z);
  R << std_radr_*std_radr_, 0, 0,
      0, std_radphi_*std_radphi_, 0,
      0, 0, std_radrd_*std_radrd_;
  S+=R;
//create vector for incoming radar measurement
VectorXd z=VectorXd(n_z);
z(0)=meas_package.raw_measurements_(0);
z(1)=meas_package.raw_measurements_(1);
z(2)=meas_package.raw_measurements_(2);

//create matrix for cross correlation Tc
MatrixXd tc=MatrixXd(n_x_,n_z);
tc.fill(0);
for(int i=0;i<2*n_aug_+1;i++){
  VectorXd z_diff=Zsig.col(i)-z_pred;
  z_diff(1)=tools.NormalizeAngle(z_diff(1));

  VectorXd x_diff=Xsig_pred_.col(i)-x_;
  x_diff(3)=tools.NormalizeAngle(x_diff(3));

  tc+=weights_(i)*x_diff*z_diff.transpose();
}

MatrixXd K=tc*S.inverse();

VectorXd y=z-z_pred;
y(1)=tools.NormalizeAngle(y(1));

x_+=K*y;
P_-=K*S*K.transpose();
}

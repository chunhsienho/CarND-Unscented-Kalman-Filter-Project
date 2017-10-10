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
  use_radar_ = true;

  // initial state vector
    x_ = VectorXd::Zero(5);

  // initial covariance matrix
    P_ = MatrixXd::Identity(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
    //adjust this to feet the low acceleration
  //std_a_ = 30;
    std_a_=5;
  // Process noise standard deviation yaw acceleration in rad/s^2
  //std_yawdd_ = 30;
    std_yawdd_=2;
    
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
    //other initializations
    
    n_x_=5;
    n_aug_=7;
    lambda=3-n_aug_;
    n_sigmapoint=2*n_aug_+1; // the number of sigma point
    Xsig_pred_=MatrixXd(n_x_,n_sigmapoint);
    
    //lidar
    H_=MatrixXd::Zero(n_sigmapoint);
    H_(0,0)=1.0;
    H_(1,1)=1.0;
    
    //Get the weight for prediction
    //wi=lamda/(lamda+n) for i=0
    weights_=VectorXd(n_sigmapoint);
    weights_(0)=lamda/(lamda+n_aug_);
    for(int i=1;i<n_sigmapoint;i++)
    {
        weights_(i)=0.5/(n_aug_+lambda);
    }
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
   
//if(!user_laser_ && measurement_pack.sensor_type_== MeasurementPackage::Laser)
//   {return;}
  
    if(!is_initialized_)
    {
        x_.fill(0,0);
    }
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
   
 

}

/**
   
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
    MatrixXd Xsig_aug=MatrixXd::Zero(n_aug_,n_sigmapoint);
    GenerateSigmaPoints(Xsig_aug);
    SigmaPointPrediction(Xsig_aug,delta_t);
    PredictMeanAndCovariance();
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
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
    VectorXd z=meas_package.raw_measurement_;
    VectorXd z_pred = H_ * x_;
    VectorXd y = z - z_pred;
    
    MatrixXd Ht = H_.transpose();
    
    MatrixXd S = H_ * P_ * H_.transpose;
    MatrixXd R = MatrixXd::Zero(2,2);
    R(0,0)= std_laspx_ * std_laspx_;
    R(1,1)= std_laspy_ * std_laspy_;
    
    S=S+R;
    MatrixXd K =P_ * H_.transpose * S.inverse();
    x_=x_+K*y
    
    int x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * H_) * P_;
    
    /* calculate the lidar NIS */
    NIS_laser_ = (z -  z_pred).transpose()* S.inverse() * (z - z_pred);
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
    int n_z=3;

    VectorXd z_pred = VectorXd::Zero(n_z);
   
    
    MatrixXd Ht = H_.transpose();
    
    MatrixXd S = MatrixXd::Zero(n_z,n_z);
    MatrixXd Zsig=MatrixXd::Zero(n_z,n_sigmapoint);
    PredictRadarMeasurement(z_pred,S,Zsig);
    VectorXd z=meas_package.raw_measurement_;
    MatrixXd Tc = MatrixXd(n_x, n_z);
    
    //calculate cross correlation matrix
    Tc.fill(0.0);
    for (int i = 0; i < 2 * n_aug + 1; i++) {  //2n+1 simga points
        
        //residual
        VectorXd z_diff = Zsig.col(i) - z_pred;
        //angle normalization
        while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
        while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
        
        // state difference
        VectorXd x_diff = Xsig_pred.col(i) - x;
        //angle normalization
        while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
        while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
        
        Tc = Tc + weights(i) * x_diff * z_diff.transpose();
    }
    //Kalman gain K;
    MatrixXd K = Tc * S.inverse();
    
    //residual
    VectorXd z_diff = z - z_pred;
    
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
    
    //update state mean and covariance matrix
    x = x + K * z_diff;
    P = P - K*S*K.transpose();
    
    /* calculate the radar NIS */
    NIS_radar_ = (z -  z_pred).transpose()* S.inverse() * (z - z_pred);

}
void UKF::GenerateSigmaPoints(MatrixXd& Xsig_aug) {
    
    
    // Initialize augmented state and covariance matrix
    //initial the aug vector
    VectorXd x_aug = VectorXd(n_aug_);
    x_aug.head(n_x_) = x_;
    x_aug(n_x_) = 0.0;
    x_aug(n_x_ + 1) = 0.0;
    //initial the p_aug vector
    MatrixXd P_aug = MatrixXd::Zero(n_aug_, n_aug_);
    P_aug.topLeftCorner(n_x_, n_x_) = P_;
    P_aug(n_x_, n_x_) = std_a_ * std_a_;
    P_aug(n_x_+1, n_x_+1) = std_yawdd_ * std_yawdd_;
    //use the x_aug and caculate the p
    MatrixXd L = P_aug.llt().matrixL();  // square root of P_aug
    L = sqrt(lambda+ n_aug_) * L;
    
    Xsig_aug.col(0) = x_aug;
    for (int i=0; i<n_aug_; i++) {
        Xsig_aug.col(i+1) = x_aug + L.col(i);
        Xsig_aug.col(n_aug_+i+1) = x_aug - L.col(i);
    }
    
    
}

void UKF::SigmaPointPrediction(MatrixXd& Xsig_out,const double delta_t) {
    
    /*******************************************************************************
     * Student part begin
     ******************************************************************************/
    
    //predict sigma points
    for (int i = 0; i< n_sigmapoint; i++)
    {
        //extract values for better readability
        double p_x = Xsig_aug(0,i);
        double p_y = Xsig_aug(1,i);
        double v = Xsig_aug(2,i);
        double yaw = Xsig_aug(3,i);
        double yawd = Xsig_aug(4,i);
        double nu_a = Xsig_aug(5,i);
        double nu_yawdd = Xsig_aug(6,i);
        
        //predicted state values
        double px_p, py_p;
        
        //avoid division by zero
        if (fabs(yawd) > 0.001) {
            px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
            py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
        }
        else {
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
        
        //write predicted sigma point into right column
        Xsig_pred(0,i) = px_p;
        Xsig_pred(1,i) = py_p;
        Xsig_pred(2,i) = v_p;
        Xsig_pred(3,i) = yaw_p;
        Xsig_pred(4,i) = yawd_p;
    }
}
void UKF::PredictMeanAndCovariance() {
    // predicted state mean and covariance using sigmapoints
    x_.fill(0.0);
    for (int i = 0; i < n_sigmapoint; i++) {  //iterate over sigma points
        x_ += weights_(i) * Xsig_pred_.col(i);
    }
    
    // predicted state covariance matrix
    P_.fill(0.0);
    for (int i = 0; i < n_sigmapoint; i++) {  //iterate over sigma points
        // state difference
        VectorXd x_diff = Xsig_pred_.col(i) - x_;
        
        //angle normalization
        x_diff(3) = NormalizeAngle(x_diff(3));
        
        P_ += weights_(i) * x_diff * x_diff.transpose() ;
    }
    
}
void UKF::PredictRadarMeasurement(){
    for (int i = 0; i < n_sigmapoint; i++) {  //2n+1 simga points
        
        // extract values for better readibility
        double p_x = Xsig_pred(0,i);
        double p_y = Xsig_pred(1,i);
        double v  = Xsig_pred(2,i);
        double yaw = Xsig_pred(3,i);
        
        double v1 = cos(yaw)*v;
        double v2 = sin(yaw)*v;
        
        // measurement model
        Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
        Zsig(1,i) = atan2(p_y,p_x);                                 //phi
        Zsig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
    }
    
    //mean predicted measurement
    VectorXd z_pred = VectorXd(n_z);
    z_pred.fill(0.0);
    for (int i=0; i < 2*n_aug+1; i++) {
        z_pred = z_pred + weights(i) * Zsig.col(i);
    }
    
    //measurement covariance matrix S
    MatrixXd S = MatrixXd(n_z,n_z);
    S.fill(0.0);
    for (int i = 0; i < 2 * n_aug + 1; i++) {  //2n+1 simga points
        //residual
        VectorXd z_diff = Zsig.col(i) - z_pred;
        
        //angle normalization
        while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
        while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
        
        S = S + weights(i) * z_diff * z_diff.transpose();
    }
    
    //add measurement noise covariance matrix
    MatrixXd R = MatrixXd(n_z,n_z);
    R <<    std_radr_*std_radr_, 0, 0,
    0, std_radphi_*std_radphi_, 0,
    0, 0,std_radrd_*std_radrd_;
    S = S + R;
}

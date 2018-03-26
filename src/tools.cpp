#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
 VectorXd rmse(4);
  rmse << 0,0,0,0;
	
	//Validate that the estimation vector size is greater than zero and that the ground truth  vector size is equal to the estimation vector size
  if(estimations.size() != ground_truth.size() || estimations.size()==0 ){
    cout << "Invalid data provided"<< endl;
    return rmse;
  }
  //accumulate the squared errors
  for(unsigned int i=0; i<estimations.size();i++){
    VectorXd error= estimations[i]-ground_truth[i];
		//co-efficientwise multiplication
    error= error.array() * error.array();
    rmse+=error;
  }
	//Calculate the mean
  rmse=rmse/estimations.size();
	//Calculate the square root
  rmse=rmse.array().sqrt();
  return rmse;
}
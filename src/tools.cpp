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
	rmse << 1000, 1000, 1000, 1000; // init to "bad" value

									// check the validity of the following inputs:
									//  * the estimation vector size should not be zero
									//  * the estimation vector size should equal ground truth vector size
	if (estimations.size() == 0)
		cout << "CalculateRMSE(): empty estimations vector" << endl;
	else if (estimations.size() != ground_truth.size())
		cout << "CalculateRMSE(): inconsistent vector size" << endl;
	else
	{

		// accumulate squared residuals
		rmse << 0, 0, 0, 0;
		for (int i = 0; i < estimations.size(); ++i) {
			VectorXd diff = estimations[i] - ground_truth[i];
			VectorXd diff_sq = diff.array() * diff.array();
			rmse += diff_sq;
		}

		// calculate the mean
		rmse *= 1.0 / estimations.size();

		// calculate the square root
		rmse = rmse.array().sqrt();
	}

	//return the result
	return rmse;
}
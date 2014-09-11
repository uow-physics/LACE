// Copyright University of Warwick 2014
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// Authors:
// John Back

/*! \file LpcAbsCurve.cc
    \brief  Class that defines an abstract local principal curve
*/

#include "LACE/LpcAbsCurve.hh"

#include "LACE/LpcHitCollection.hh"

#include <iostream>

LpcAbsCurve::LpcAbsCurve(int index,
			 const LpcPoint& startPoint, 
			 const std::vector<LpcPoint>& lpcPoints,
			 const Eigen::MatrixXd& eigenVectors, 
			 const Eigen::VectorXd& cosAngles,
			 const LpcPathLength& lpcPath, 
			 const Eigen::VectorXd& rho, 
			 const Eigen::VectorXd& c0,
			 const std::vector<LpcPoint>& highRhoPoints,
			 const LpcHitCollection* theHits,
			 int flag) :
    type_(LpcAbsCurve::Abstract),
    index_(index),
    startPoint_(startPoint),
    lpcPoints_(lpcPoints),
    eigenVectors_(eigenVectors),
    cosAngles_(cosAngles),
    lpcPath_(lpcPath),
    rho_(rho),
    c0_(c0),
    highRhoPoints_(highRhoPoints),
    theHits_(theHits),
    nDim_(3),
    flag_(flag),
    residuals_(),
    passedPeaks_(),
    peakRanges_()
{
    this->findResiduals();
}

LpcAbsCurve::LpcAbsCurve() :
    type_(LpcAbsCurve::Abstract),
    index_(-1),
    startPoint_(),
    lpcPoints_(),
    eigenVectors_(),
    cosAngles_(),
    lpcPath_(),
    rho_(),
    c0_(),
    highRhoPoints_(),
    theHits_(0),
    nDim_(3),
    residuals_(),
    passedPeaks_(),
    peakRanges_()
{
}

LpcAbsCurve::~LpcAbsCurve()
{
}

void LpcAbsCurve::findResiduals()
{
    
    int nLpc = lpcPoints_.size();
    int nHits(0);
    // Check if the pointer to the hit collection is null
    if (!theHits_) {return;}

    nDim_ = theHits_->getNDimensions();
    nHits = theHits_->getNumberOfHits();

    // Get the unscaled hit co-ordinates
    const Eigen::MatrixXd hitArray = theHits_->getCoords();

    // Get the hit weights
    const Eigen::VectorXd weights = theHits_->getWeights();

    Eigen::MatrixXd lpcArray = Eigen::MatrixXd::Zero(nLpc, nDim_);
    
    for (int i = 0; i < nLpc; i++) {
	LpcPoint thePoint = lpcPoints_[i];
	lpcArray.row(i) = thePoint.getCoords();
    }

    // Arrays for the average residual of the closest data points and
    // the number of closest points (for calculating the average recursively)  
    Eigen::VectorXd lpcResiduals = Eigen::VectorXd::Zero(nLpc);
    Eigen::VectorXd lpcWeightRes = Eigen::VectorXd::Zero(nLpc);
    // Store the number of points as floats, since we need to divide by them
    Eigen::VectorXd lpcNClose = Eigen::VectorXd::Zero(nLpc);
    Eigen::VectorXd lpcWClose = Eigen::VectorXd::Zero(nLpc);

    // Arrays for the nearest lpc residual for each data point
    Eigen::VectorXd hitResiduals = Eigen::VectorXd::Zero(nHits);
    Eigen::VectorXd hitWeightRes = Eigen::VectorXd::Zero(nHits);

    // Array storing the closest lpc point index number for each data point
    Eigen::VectorXi hitNearestLpc = Eigen::VectorXi::Zero(nHits);

    // Here, we calculate two sets of residuals. The first is the average
    // residual of the closest hit points for _each LPC point_.
    // The next is the residual of the closest LPC point for _each hit point_
            
    // For the average residual per LPC point, we loop over all data hit points and find
    // the LPC point that is closest to it. We then recursively add the residual 
    // distance to get the average residual for the matched LPC point.
    // For the residual per hit point, just find the closest LPC point and store
    // the separation distance    

    // Loop over hits
    for (int iH = 0; iH < nHits; iH++) {

	Eigen::VectorXd hitPoint = hitArray.row(iH);
	double weight = weights(iH);
	double minDistSq(0.0);

	int jClosest(0);

	// Loop over lpc points
	for (int jL = 0; jL < nLpc; jL++) {

	    Eigen::VectorXd lpcPoint = lpcArray.row(jL);
	    double distSq = (hitPoint - lpcPoint).squaredNorm();
	    if (jL == 0) {minDistSq = distSq;}

	    if (distSq < minDistSq) {
		minDistSq = distSq;
		jClosest = jL;
	    }

	}

	// For the jClosest lpc point, find the average closest residual
	// distance for each hit point using recursion
	double N = lpcNClose(jClosest);
	double N1 = N + 1.0;
	double minDist = sqrt(minDistSq);

	lpcResiduals(jClosest) = (N*lpcResiduals(jClosest) + minDist)/N1;

	// Keep track of the weights for the jClosest lpc point
	double w = lpcWClose(jClosest);
	double w1 = w + weight;

	if (fabs(w1) > 1e-10) {
	    lpcWeightRes(jClosest) = (w*lpcWeightRes(jClosest) + minDist*weight)/w1;
	}

	lpcNClose(jClosest) = N1;
	lpcWClose(jClosest) = w1;

	hitResiduals(iH) = minDist;
	hitWeightRes(iH) = minDist*weight;
	hitNearestLpc(iH) = jClosest;

    }
	
    // Store the residuals
    residuals_ = LpcResiduals(lpcResiduals, hitResiduals, lpcWeightRes,
			      hitWeightRes, hitNearestLpc);

}

void LpcAbsCurve::print() const
{
    
    std::cout<<"*** Printing the LpcAbsCurve (type = "<<type_
	     <<", index = "<<index_<<") result:"<<std::endl;
    std::cout<<"  Starting point = "<<startPoint_<<std::endl;

    for (size_t i = 0; i < lpcPoints_.size(); i++) {
	
	LpcPoint thePoint = lpcPoints_[i];
	std::cout<<"  LpcPoint is "<<thePoint<<std::endl;

    }

    std::cout<<"  eigenVectors = "<<eigenVectors_<<std::endl;
    std::cout<<"  cosAngles = "<<cosAngles_.transpose()<<std::endl;
    std::cout<<"  rho = "<<rho_.transpose()<<std::endl;
    std::cout<<"  c0 ="<<c0_.transpose()<<std::endl;

    for (size_t i = 0; i < highRhoPoints_.size(); i++) {

	LpcPoint rhoPoint = highRhoPoints_[i];
	std::cout<<"  HighRhoPoint "<<i<<" is "<<rhoPoint<<std::endl;

    }

    std::cout<<"   Pathlengths"<<std::endl;
    std::cout<<"   scaledLambda = "<<lpcPath_.getScaledLambda().transpose()<<std::endl;
    std::cout<<"   scaledLambdaAxes = "<<lpcPath_.getScaledLambdaAxes().transpose()<<std::endl;
    std::cout<<"   Lambda = "<<lpcPath_.getLambda().transpose()<<std::endl;
    std::cout<<"   LambdaAxes = "<<lpcPath_.getLambdaAxes().transpose()<<std::endl;
    std::cout<<"   scaledDeltaLambda = "<<lpcPath_.getScaledDeltaLambda().transpose()<<std::endl;
    std::cout<<"   deltaLambda = "<<lpcPath_.getDeltaLambda().transpose()<<std::endl<<std::endl;

}

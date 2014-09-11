// Copyright University of Warwick 2014
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// Authors:
// John Back

/*! \file LpcBranchAlgorithm.cc
    \brief Class that finds possible branches for the lpc

    This is a C++ implementation of the branch finding algorithm in 
    the R-code package LPCM written by Jochen Einbeck and Ludger Evers
    (http://cran.r-project.org/web/packages/LPCM/index.html)
*/

#include "LACE/LpcBranchAlgorithm.hh"

#include "LACE/LpcBranch.hh"
#include "LACE/LpcBranchCollection.hh"
#include "LACE/LpcCurve.hh"
#include "LACE/LpcHitCollection.hh"
#include "LACE/LpcParameters.hh"
#include "LACE/LpcPoint.hh"

#include <utility>

LpcBranchAlgorithm::LpcBranchAlgorithm(const LpcParameters* theParameters) :
    thePars_(theParameters),
    nDim_(3),
    theHits_(0),
    Xi_(),
    weights_(),
    x0_(),
    startPoint_(),
    theRange_(),
    mu_x_(),
    dx_(),
    E2Vector_(),
    jWeights_(),
    h_(0.05),
    t_(0.05),
    gapt_(0.075),
    iter2_(250),
    iter_(125),
    pen_(2.0),
    rho0_(0.4),
    boundary_(0.005),
    convergence_(1e-6),
    N_(0),
    threshold_(0.0),
    branchLevel_(0),
    functions_(),
    lpcAlgorithm_(theParameters),
    startPoints_(),
    mainCurveIndex_(0)
{
}

LpcBranchAlgorithm::~LpcBranchAlgorithm()
{
}

void LpcBranchAlgorithm::initialise()
{

    // Initialise the lpc parameters
    if (thePars_) {

	h_ = thePars_->getKernelWidth();
	t_ = thePars_->getStepSize();
	iter2_ = thePars_->getNLpcPoints();
	iter_ = iter2_/2;
	pen_ = thePars_->getAnglePenalisation();
	rho0_ = thePars_->getEigenRatio();
	boundary_ = thePars_->getBoundary();
	convergence_ = thePars_->getConvergence();
	gapt_ = thePars_->getBranchGapSize()*t_;
	branchLevel_ = thePars_->getBranchLevel();
    }

}

LpcBranchCollection* LpcBranchAlgorithm::getBranches(const LpcHitCollection* theHits,
						     const std::vector<LpcPoint>& theHighRhoPoints,
						     int mainCurveIndex)
{

    LpcBranchCollection* theBranches  = new LpcBranchCollection();    

    if (!theHits) {return theBranches;}

    nDim_ = theHits->getNDimensions();

    // Check for at least 1 dimension
    if (nDim_ < 1) {return theBranches;}

    this->initialise();

    if (theHighRhoPoints.size() == 0 || branchLevel_ < 1) {return theBranches;}

    theHits_ = theHits;

    // Retrieve the co-ordinate and hit weight information. 
    // These will always be the same for all possible branches
    Xi_ = theHits_->getScaledCoords();
    weights_ = theHits_->getWeights();
    N_ = Xi_.rows();
    theRange_ = theHits_->getRange();  
    // Keep track of the index of the main curve that these branches will be associated with
    mainCurveIndex_ = mainCurveIndex;

    // Select the starting point to initially be from the vector of highRhoPoints
    startPoints_ = theHighRhoPoints;
    // Start the index numbers for the branches at mainCurveIndex+1 (usually = 0+1 = 1)
    int startIndex(mainCurveIndex+1);

    for (int iG = 0; iG < branchLevel_; iG++) {

	std::vector<LpcBranch*> branchVect = this->getGenBranches(startIndex);
	// Add the branches, starting at generation level 1 (hence the +1)
	theBranches->addBranches(iG+1, branchVect);

	// Update the starting index for the next generation
	startIndex += branchVect.size();

	// Also update the new starting points to be the new high rho points
	startPoints_.clear();
	// Loop over the branches we have for this generation
	std::vector<LpcBranch*>::const_iterator bIter;
	for (bIter = branchVect.begin(); bIter != branchVect.end(); ++bIter) {

	    const LpcBranch* theBranch = *bIter;
	    if (theBranch != 0) {
		// Find the high rho points for this branch
		std::vector<LpcPoint> newRhoPoints = theBranch->getHighRhoPoints();

		std::vector<LpcPoint>::const_iterator pIter;
		for (pIter = newRhoPoints.begin(); pIter != newRhoPoints.end(); ++pIter) {
		    // Add the high rho point to the starting points vector
		    startPoints_.push_back(*pIter);

		} // Loop over branch high rho points

	    } // Branch is valid

	} // Loop over branches in current generation

    } // Generation loop

    return theBranches;

}

std::vector<LpcBranch*> LpcBranchAlgorithm::getGenBranches(int startIndex)
{

    std::vector<LpcBranch*> genBranches;

    int nRho = startPoints_.size();
    for (int i = 0; i < nRho; i++) {

	x0_ = startPoints_[i].getScaledCoords();
	int index = startIndex + i;

	LpcBranch* theBranch = this->getBranch(index);

	genBranches.push_back(theBranch);

    }

    return genBranches;

}

LpcBranch* LpcBranchAlgorithm::getBranch(int index)
{
    
    // First, we need to find the new forward position away from the curve
    // We therefore need the 2nd largest eigenvector using the starting
    // position given by the previous highRhoPoint ("x0_").

    // New forward position along the curve. Multiply the stepsize t by
    // the branching gapsize (the total branch step) to move away 
    // from the initial curve
    this->prepareNewStartPoint();

    Eigen::VectorXd new_x1 = mu_x_;
    new_x1 += dx_;

    double kdexVal1 = functions_.kdex(Xi_, new_x1, h_);

    LpcCurve *curve1(0), *curve2(0);

    if (kdexVal1 > threshold_) {

	LpcPoint startPoint1(-1, new_x1);
	startPoint1.unscale(theRange_);

	curve1 = lpcAlgorithm_.getCurve(index, theHits_, startPoint1, 
					jWeights_, E2Vector_, 
					LpcAlgorithm::Forward);

    }

    Eigen::VectorXd new_x2 = mu_x_;
    new_x2 -= dx_;

    double kdexVal2 = functions_.kdex(Xi_, new_x2, h_);

    if (kdexVal2 > threshold_) {

	LpcPoint startPoint2(-1, new_x2);
	startPoint2.unscale(theRange_);

	// The 2nd eigenvector needs to be "reversed" (minus sign)
	curve2 = lpcAlgorithm_.getCurve(index, theHits_, startPoint2,
					jWeights_, -E2Vector_,
					LpcAlgorithm::Backward);

    }

    // From these curves, build the LpcBranch pointer.
    int iter2 = thePars_->getNLpcPoints();
    int iter = iter2/2;
    
    // Set the starting point for the branch result as the initial highRho point
    LpcPoint startPoint(-1, x0_);
    startPoint.unscale(theRange_);

    std::vector<LpcPoint> lpcPoints;
    Eigen::MatrixXd eigenVectors = Eigen::MatrixXd::Zero(iter2, nDim_);
    Eigen::VectorXd cosAngles = Eigen::VectorXd::Ones(iter2);
    Eigen::VectorXd rho = Eigen::VectorXd::Zero(iter2);
    Eigen::VectorXd c0 = Eigen::VectorXd::Ones(iter2);
    std::vector<LpcPoint> highRhoPoints;
    // For the path lengths
    Eigen::VectorXd sLambda = Eigen::VectorXd::Zero(iter2);
    Eigen::MatrixXd sLambdaAxes = Eigen::MatrixXd::Zero(iter2, nDim_);

    if (curve1) {

	std::vector<LpcPoint> lpcP1 = curve1->getLpcPoints();
	Eigen::MatrixXd EVectors1 = curve1->getEigenVectors(); 
	Eigen::VectorXd cAngles1 = curve1->getCosAngles();
	Eigen::VectorXd rho1 = curve1->getRho();
	Eigen::VectorXd c01 = curve1->getc0();

	LpcPathLength lpcPath1 = curve1->getPathLength();
	Eigen::VectorXd sL1 = lpcPath1.getScaledLambda();
	Eigen::MatrixXd sLAxes1 = lpcPath1.getScaledLambdaAxes();

	for (int i = 0; i < iter; i++) {

	    lpcPoints.push_back(lpcP1[i]);
	    eigenVectors.row(i) = EVectors1.row(i);
	    cosAngles(i) = cAngles1(i);
	    rho(i) = rho1(i);
	    c0(i) = c01(i);

	    sLambda(i) = sL1(i);
	    sLambdaAxes.row(i) = sLAxes1.row(i);

	}

	std::vector<LpcPoint> highRP1 = curve1->getHighRhoPoints();
	std::vector<LpcPoint>::const_iterator iter;
	for (iter = highRP1.begin(); iter != highRP1.end(); ++iter) {

	    highRhoPoints.push_back(*iter);

	}

    }

    if (curve2) {

	std::vector<LpcPoint> lpcP2 = curve2->getLpcPoints();
	Eigen::MatrixXd EVectors2 = curve2->getEigenVectors(); 
	Eigen::VectorXd cAngles2 = curve2->getCosAngles();
	Eigen::VectorXd rho2 = curve2->getRho();
	Eigen::VectorXd c02 = curve2->getc0();

	LpcPathLength lpcPath2 = curve2->getPathLength();
	Eigen::VectorXd sL2 = lpcPath2.getScaledLambda();
	Eigen::MatrixXd sLAxes2 = lpcPath2.getScaledLambdaAxes();

	for (int i = iter; i < iter2; i++) {

	    lpcPoints.push_back(lpcP2[i]);
	    eigenVectors.row(i) = EVectors2.row(i);
	    cosAngles(i) = cAngles2(i);
	    rho(i) = rho2(i);
	    c0(i) = c02(i);

	    sLambda(i) = sL2(i);
	    sLambdaAxes.row(i) = sLAxes2.row(i);

	}

	std::vector<LpcPoint> highRP2 = curve2->getHighRhoPoints();
	std::vector<LpcPoint>::const_iterator iter;
	for (iter = highRP2.begin(); iter != highRP2.end(); ++iter) {

	    highRhoPoints.push_back(*iter);

	}

    }


    LpcPathLength lpcPath(sLambda, sLambdaAxes, theRange_);

    LpcBranch* theBranch = new LpcBranch(index, startPoint, lpcPoints, eigenVectors,
					 cosAngles, lpcPath, rho, c0, highRhoPoints, 
					 theHits_, mainCurveIndex_);

    // Delete the two separate curve pointers
    delete curve1;
    delete curve2;

    return theBranch;

}

void LpcBranchAlgorithm::prepareNewStartPoint() 
{

    // Define the array of the kernel weights for the data hits
    Eigen::VectorXd kernelWeights = Eigen::VectorXd::Zero(N_);    

    // Loop over the data points & calculate the kernel weights
    for (int ip = 0; ip < N_; ip++) {

	Eigen::VectorXd hitPoint = Xi_.row(ip);
	double w = weights_(ip);

	kernelWeights(ip) = w*functions_.kernelFunction(hitPoint, x0_, h_);

    }

    // Find and set the weighted local mean
    mu_x_ = functions_.getWeightedMean(Xi_, kernelWeights);

    // Subtract the local mean from the hit positions
    Eigen::MatrixXd mean_sub = functions_.offsetPositions(Xi_, mu_x_);

    // Calculate the symmetric covariance matrix
    Eigen::MatrixXd cov_x = functions_.formCovarianceMatrix(mean_sub, kernelWeights);

    // Get the normalised eigen vectors and values for the covariance matrix
    std::pair<Eigen::VectorXd, Eigen::MatrixXd> eigenPair = 
	functions_.findNormEigenVectors(cov_x);

    Eigen::VectorXd eigenValues = eigenPair.first;
    Eigen::MatrixXd eigenVectors = eigenPair.second;

    // Get the eigenvector with the 2nd (not 1st) largest eigenvalue
    E2Vector_ = eigenVectors.row(1);
    
    // Set the required displacement for the branch starting point
    dx_ = gapt_*E2Vector_;

    // Fill the jWeights array
    jWeights_ = Eigen::VectorXd::Zero(N_);
    for (int ik = 0; ik < N_; ik++) {

	Eigen::VectorXd hitPoint = Xi_.row(ik);
	jWeights_(ik) = 1.0 - functions_.kernelFunction(hitPoint, mu_x_, h_);
	jWeights_(ik) *= weights_(ik);

    }

}

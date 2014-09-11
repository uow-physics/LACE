// Copyright University of Warwick 2014
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// Authors:
// John Back

/*! \file LpcAlgorithm.cc
    \brief Class that follows the scaled data points to find the lpc

    This is a C++ implementation of the main Lpc curve algorithm in 
    the R-code package LPCM written by Jochen Einbeck and Ludger Evers 
    (http://cran.r-project.org/web/packages/LPCM/index.html)
*/

#include "LACE/LpcAlgorithm.hh"

#include "LACE/LpcCurve.hh"
#include "LACE/LpcHitCollection.hh"
#include "LACE/LpcParameters.hh"

#include <utility>

LpcAlgorithm::LpcAlgorithm(const LpcParameters* theParameters) :
    thePars_(theParameters),
    nDim_(3),
    Xi_(),
    weights_(),
    x0_(),
    startPoint_(),
    theRange_(),
    h_(0.05),
    t_(0.05),
    iter2_(250),
    iter_(125),
    pen_(2.0),
    rho0_(0.4),
    boundary_(0.005),
    convergence_(1e-6),
    N_(0),
    functions_()
{
}

LpcAlgorithm::~LpcAlgorithm()
{
}

void LpcAlgorithm::initialise()
{

    // Initialise the lpc parameters. This should be called by getCurve()
    if (thePars_) {
	
	h_ = thePars_->getKernelWidth();
	t_ = thePars_->getStepSize();
	iter2_ = thePars_->getNLpcPoints();
	iter_ = iter2_/2;
	pen_ = thePars_->getAnglePenalisation();
	rho0_ = thePars_->getEigenRatio();
	boundary_ = thePars_->getBoundary();
	convergence_ = thePars_->getConvergence();	

    }

    // Initialise various matrices/vectors
    // The local neighbourhood ("x" or "u") points with eigenvalue ratios > rho0
    highRhoPoints_.clear();

    // The lpc points ("m(u)" or "mu"), equivalent to the local mean
    lpcPoints_ = Eigen::MatrixXd::Zero(iter2_, nDim_);

    // The highest eigenvalue eigenvectors
    eigenVectors_ = Eigen::MatrixXd::Zero(iter2_, nDim_);

    // The cosine between neighbouring eigenvectors
    cosAngles_ = Eigen::VectorXd::Ones(iter2_);

    // The values of the cumulative path length
    lambda_ = Eigen::VectorXd::Zero(iter2_);

    // The matrix of the cumulative path length along each axis (row = point, col = x,y,z..)
    lambdaAxes_ = Eigen::MatrixXd::Zero(iter2_, nDim_);

    // The ratio of the eigenvalues between neighbouring eigenvectors
    rho_ = Eigen::VectorXd::Zero(iter2_);

    // The c0 coefficients used to scale the kernel factor
    c0_ = Eigen::VectorXd::Ones(iter2_);

}

LpcCurve* LpcAlgorithm::getCurve(int curveIndex, const LpcHitCollection* theHits, int direction)
{

    // Get the curve using the starting point and weights from the hit collection directly
    if (!theHits) {return 0;}

    LpcPoint startPoint = theHits->getStartPoint();
    Eigen::VectorXd weights = theHits->getWeights();
    Eigen::VectorXd lastEVect = Eigen::VectorXd::Zero(nDim_);

    LpcCurve* theCurve = this->getCurve(curveIndex, theHits, startPoint, weights,
					lastEVect, direction);

    return theCurve;
    

}

LpcCurve* LpcAlgorithm::getCurve(int curveIndex, const LpcHitCollection* theHits, 
				 const LpcPoint& startPoint, const Eigen::VectorXd& weights, 
				 const Eigen::VectorXd& lastEVect, int direction)
{

    if (!theHits) {return 0;}

    // Initialise parameters and internal arrays
    nDim_ = theHits->getNDimensions();
    // Check for at least 1 dimension
    if (nDim_ < 1) {return 0;}

    this->initialise();

    // Retrieve the co-ordinate and hit weight information
    Xi_ = theHits->getScaledCoords();
    weights_ = weights;
    x0_ = startPoint.getScaledCoords();

    startPoint_ = startPoint;
    theRange_ = theHits->getRange();    

    N_ = Xi_.rows();

    // Find the forward and backward curve sections, as required
    if (direction == LpcAlgorithm::Forward) {

	this->forward(lastEVect);

    } else if (direction == LpcAlgorithm::Backward) {

	this->backward(lastEVect);

    } else {

	this->forward(lastEVect);
	this->backward(lastEVect);

    }

    // Create a vector of the lpc points themselves
    std::vector<LpcPoint> lpcPoints;

    for (int i = 0; i < iter2_; i++) {

	Eigen::VectorXd lp = lpcPoints_.row(i);

	LpcPoint thePoint(i, lp);
	// Also find the unscaled co-ordinates
	thePoint.unscale(theRange_);

	lpcPoints.push_back(thePoint);

    }

    // Store the scaled path lengths. This will also calculate the unscaled values
    LpcPathLength lpcPaths(lambda_, lambdaAxes_, theRange_);

    // Construct the lpc curve pointer, storing all relevant info. This will
    // also automatically calculate and store the lpc-to-hit residuals
    int mainCurveFlag(-1);
    LpcCurve* theCurve = new LpcCurve(curveIndex, startPoint_, lpcPoints, eigenVectors_,
				      cosAngles_, lpcPaths, rho_, c0_, highRhoPoints_,
				      theHits, mainCurveFlag);
    
    // Return the lpc
    return theCurve;

}

void LpcAlgorithm::forward(const Eigen::VectorXd& lastEVect)
{

    int iterVal = iter_ - 1;

    // Set the initial local point
    Eigen::VectorXd x = x0_;

    // Find the forward-going lpc points
    for (int i = iterVal; i >= 0; i--) {

	// Finding lpc point i

	// Calculate the kernel weights for the lpc point
	Eigen::VectorXd kernelWeights = Eigen::VectorXd::Zero(N_);

	// The kernel denominator factor c*h
	double kernel_c = c0_(i)*h_;

	// Loop over the data points & calculate the kernel weights
	for (int ip = 0; ip < N_; ip++) {

	    Eigen::VectorXd hitPoint = Xi_.row(ip);
	    double w = weights_(ip);

	    kernelWeights(ip) = w*functions_.kernelFunction(hitPoint, x, kernel_c);

	}

	// Find the weighted local mean
	Eigen::VectorXd mu_x = functions_.getWeightedMean(Xi_, kernelWeights);

	// Subtract the local mean from the hit positions
	Eigen::MatrixXd mean_sub = functions_.offsetPositions(Xi_, mu_x);

	// Calculate the symmetric covariance matrix
	Eigen::MatrixXd cov_x = functions_.formCovarianceMatrix(mean_sub, kernelWeights);

	// Store the current lpc point
	lpcPoints_.row(i) = mu_x;

	// Update the cumulative path lengths along the lpc
	int j = i + 1;

	if (i == iterVal) {

	    lambda_(i) = 0.0;
	    lambdaAxes_.row(i) = Eigen::VectorXd::Zero(nDim_);

	} else {

	    Eigen::VectorXd prevLpcPoint = lpcPoints_.row(j);
	    Eigen::VectorXd lpcDiff = mu_x - prevLpcPoint;
	    // Get the distance (magnitude only)
	    for (int k = 0; k < lpcDiff.size(); k++) {
		if (lpcDiff(k) < 0.0) {lpcDiff(k) *= -1.0;}
	    }

	    lambda_(i) = lambda_(j) + lpcDiff.norm();
	    Eigen::VectorXd lAxes_j = lambdaAxes_.row(j);
	    lambdaAxes_.row(i) = lAxes_j + lpcDiff;

	}


	// Get the normalised eigen vectors and values for the covariance matrix
	std::pair<Eigen::VectorXd, Eigen::MatrixXd> eigenPair = 
	    functions_.findNormEigenVectors(cov_x);

	Eigen::VectorXd eigenValues = eigenPair.first;
	Eigen::MatrixXd eigenVectors = eigenPair.second;

	// Store the eigenvector with the largest eigenvalue
	Eigen::VectorXd mainEVector = eigenVectors.row(0);
	eigenVectors_.row(i) = mainEVector;

	// Store the ratio of the two largest eigenvalues
	rho_(i) = 0.0;

	double maxEigenValue = eigenValues(0);
	if (fabs(maxEigenValue) > 1e-10 && nDim_ > 1) {
	    rho_(i) = eigenValues(1)/maxEigenValue;
	}

	if (i < iterVal && rho_(i) > rho0_ && rho_(j) < rho0_) {

	    // Create a "high rho" LpcPoint object
	    LpcPoint theHighPoint(i, x);
	    // Obtain its unscaled co-ordinates
	    theHighPoint.unscale(theRange_);
	    // Store in the internal vector
	    highRhoPoints_.push_back(theHighPoint);

	}

	// Calculate the angle between successive eigenvectors. This needs
	// the previous eigenvector
	double lastEVectNorm(0.0);
	if (i == iterVal) {lastEVectNorm = lastEVect.norm();}

	if (i == iterVal && lastEVectNorm > 0.0) {

	    cosAngles_(i) = lastEVect.dot(eigenVectors_.row(i));

	} else if (i < iterVal) {	   

	    cosAngles_(i) = eigenVectors_.row(j).dot(eigenVectors_.row(i));
	}

	// Sign flipping
	if (cosAngles_(i) < 0.0) {eigenVectors_.row(i) *= -1.0;}

	// Angle penalisation factor
	if (pen_ > 0.0) {

	    double a = std::pow(fabs(cosAngles_(i)), pen_);

	    if (i == iterVal && lastEVectNorm > 0.0) {

		Eigen::VectorXd eigenRow_i = eigenVectors_.row(i);
		eigenVectors_.row(i) = a*eigenRow_i + (1.0 - a)*lastEVect;

	    } else if (i < iterVal) {

		Eigen::VectorXd eigenRow_i = eigenVectors_.row(i);
		Eigen::VectorXd eigenRow_j = eigenVectors_.row(j);
		eigenVectors_.row(i) = a*eigenRow_i + (1.0 - a)*eigenRow_j;

	    }

	}

	// Step along the direction of the (penalised) eigenvector
	x = eigenVectors_.row(i)*t_;
	x += mu_x;

	// Check curve termination criteria
	if (i > 0 && i < iterVal) {

	    double conv_denom = fabs(lambda_(i) + lambda_(j));	    
	    double conv_ratio(0.0);
	    if (conv_denom > 1e-10) {
		conv_ratio = fabs(lambda_(i) - lambda_(j))/conv_denom;
	    }

	    if (conv_ratio < convergence_) {break;}

	    // Adjust the boundary conditions for the next lpc point
	    // Since we are iterating downwards, this means the lpc point index
	    // number is decremented by 1
	    int k = i - 1;
	    if (conv_ratio < boundary_) {

		c0_(k) = 0.995*c0_(i);

	    } else {

		c0_(k) = 1.01*c0_(i);
		if (c0_(k) > 1.0) {c0_(k) = 1.0;}

	    }

	}

    } // Loop over i

}

void LpcAlgorithm::backward(const Eigen::VectorXd& lastEVect)
{
    
    int iterVal = iter_;

    // Set the initial local point
    Eigen::VectorXd x = x0_;

    // Find the backward-going lpc points
    for (int i = iterVal; i < iter2_; i++) {

	// Finding lpc point i

	// Calculate the kernel weights for the lpc point
	Eigen::VectorXd kernelWeights = Eigen::VectorXd::Zero(N_);

	// The kernel denominator factor c*h
	double kernel_c = c0_(i)*h_;

	// Loop over the data points & calculate the kernel weights
	for (int ip = 0; ip < N_; ip++) {

	    Eigen::VectorXd hitPoint = Xi_.row(ip);
	    double w = weights_(ip);

	    kernelWeights(ip) = w*functions_.kernelFunction(hitPoint, x, kernel_c);

	}

	// Find the weighted local mean
	Eigen::VectorXd mu_x = functions_.getWeightedMean(Xi_, kernelWeights);

	// Subtract the local mean from the hit positions
	Eigen::MatrixXd mean_sub = functions_.offsetPositions(Xi_, mu_x);

	// Calculate the symmetric covariance matrix
	Eigen::MatrixXd cov_x = functions_.formCovarianceMatrix(mean_sub, kernelWeights);

	// Store the current lpc point
	lpcPoints_.row(i) = mu_x;

	// Update the cumulative path lengths along the lpc
	int j = i - 1;

	Eigen::VectorXd prevLpcPoint = lpcPoints_.row(j);
	// Get the distance (magnitude only)
	Eigen::VectorXd lpcDiff = mu_x - prevLpcPoint;
	for (int k = 0; k < lpcDiff.size(); k++) {
	    if (lpcDiff(k) < 0.0) {lpcDiff(k) *= -1.0;}
	}

	if (i == iterVal) {

	    lambda_(i) = -lpcDiff.norm();
	    lambdaAxes_.row(i) = -lpcDiff;

	} else {

	    lambda_(i) = lambda_(j) - lpcDiff.norm();
	    Eigen::VectorXd lAxes_j = lambdaAxes_.row(j);
	    lambdaAxes_.row(i) = lAxes_j - lpcDiff;

	}

	// Get the normalised eigen vectors and values for the covariance matrix
	std::pair<Eigen::VectorXd, Eigen::MatrixXd> eigenPair = 
	    functions_.findNormEigenVectors(cov_x);

	Eigen::VectorXd eigenValues = eigenPair.first;
	Eigen::MatrixXd eigenVectors = eigenPair.second;

	// Store the eigenvector with the largest eigenvalue
	Eigen::VectorXd mainEVector = eigenVectors.row(0);
	eigenVectors_.row(i) = mainEVector;

	// Store the ratio of the two largest eigenvalues
	rho_(i) = 0.0;

	double maxEigenValue = eigenValues(0);
	if (fabs(maxEigenValue) > 1e-10 && nDim_ > 1) {
	    rho_(i) = eigenValues(1)/maxEigenValue;
	}

	if (i > iterVal && rho_(i) > rho0_ && rho_(j) < rho0_) {

	    // Create a "high rho" LpcPoint object
	    LpcPoint theHighPoint(i, x);
	    // Obtain its unscaled co-ordinates
	    theHighPoint.unscale(theRange_);
	    // Store in the internal vector
	    highRhoPoints_.push_back(theHighPoint);

	}

	// Calculate the angle between successive eigenvectors. This needs
	// the previous eigenvector
	double lastEVectNorm(0.0);
	if (i == iterVal) {lastEVectNorm = lastEVect.norm();}

	if (i == iterVal && lastEVectNorm > 0.0) {

	    cosAngles_(i) = -lastEVect.dot(eigenVectors_.row(i));

	} else if (i > iterVal) {

	    cosAngles_(i) = eigenVectors_.row(j).dot(eigenVectors_.row(i));

	}

	// Sign flipping
	if (cosAngles_(i) < 0.0) {eigenVectors_.row(i) *= -1.0;}

	// Angle penalisation factor
	if (pen_ > 0.0) {

	    double a = std::pow(fabs(cosAngles_(i)), pen_);

	    if (i == iterVal && lastEVectNorm > 0.0) {

		Eigen::VectorXd eigenRow_i = eigenVectors_.row(i);
		eigenVectors_.row(i) = a*eigenRow_i + (1.0 - a)*lastEVect;

	    } else if (i > iterVal) {

		Eigen::VectorXd eigenRow_i = eigenVectors_.row(i);
		Eigen::VectorXd eigenRow_j = eigenVectors_.row(j);
		eigenVectors_.row(i) = a*eigenRow_i + (1.0 - a)*eigenRow_j;
		
	    }

	}

	// Step along the direction of the (penalised) eigenvector
	x = -eigenVectors_.row(i)*t_;
	x += mu_x;

	// Check curve termination criteria
	if (i > iterVal && i < iter2_ - 1) {

	    double conv_denom = fabs(lambda_(i) + lambda_(j));	    
	    double conv_ratio(0.0);
	    if (conv_denom > 1e-10) {
		conv_ratio = fabs(lambda_(i) - lambda_(j))/conv_denom;
	    }

	    if (conv_ratio < convergence_) {break;}

	    // Adjust the boundary conditions for the next lpc point
	    // Since we are iterating upwards, this means the lpc point index
	    // number is incremented by 1
	    int k = i + 1;
	    if (conv_ratio < boundary_) {

		c0_(k) = 0.995*c0_(i);

	    } else {

		c0_(k) = 1.01*c0_(i);
		if (c0_(k) > 1.0) {c0_(k) = 1.0;}

	    }

	}

    } // Loop over i

}


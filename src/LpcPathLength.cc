// Copyright University of Warwick 2014
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// Authors:
// John Back

/*! \file LpcPathLength.cc
    \brief Class used to define a series of pathlengths along the lpc
*/

#include "LACE/LpcPathLength.hh"

LpcPathLength::LpcPathLength() :
    N_(0),
    nDim_(3),
    sL_(),
    sLAxes_(),
    theRange_(Eigen::VectorXd::Ones(3)),
    sdL_(),
    L_(),
    LAxes_(),
    dL_()    
{
}

LpcPathLength::LpcPathLength(const Eigen::VectorXd& sLambda, 
			     const Eigen::MatrixXd& sLambdaAxes,
			     const Eigen::VectorXd& theRange) :
    N_(sLambda.size()),
    nDim_(sLambdaAxes.cols()),
    sL_(sLambda),
    sLAxes_(sLambdaAxes),
    theRange_(theRange),
    sdL_(),
    L_(),
    LAxes_(),
    dL_()    
{

    this->initialise();

}

LpcPathLength::LpcPathLength(const LpcPathLength& other) :
    N_(other.N_),
    nDim_(other.nDim_),
    sL_(other.sL_),
    sLAxes_(other.sLAxes_),
    theRange_(other.theRange_),
    sdL_(other.sdL_),
    L_(other.L_),
    LAxes_(other.LAxes_),
    dL_(other.dL_)    
{
}

LpcPathLength& LpcPathLength::operator = (const LpcPathLength& other)
{
    
    N_ = other.N_;
    nDim_ = other.nDim_;
    sL_ = other.sL_;
    sLAxes_ = other.sLAxes_;
    theRange_ = other.theRange_;
    sdL_ = other.sdL_;
    L_ = other.L_;
    LAxes_ = other.LAxes_;
    dL_ = other.dL_;	

    return *this;
}


void LpcPathLength::initialise()
{
    // Using the co-ordinate range and the scaled pathlengths,
    // find the unscaled pathlengths and also calculate the
    // delta pathlength values such that it starts at 0 and
    // ends at the total pathlength.

    // First unscale the x,y,z components
    LAxes_ = Eigen::MatrixXd::Zero(N_, nDim_);

    int i(0), j(0);
    for (i = 0; i < N_; i++) {

	Eigen::VectorXd sLAxesRow = sLAxes_.row(i);
	LAxes_.row(i) = sLAxesRow.cwiseProduct(theRange_);

    }

    // From these, find the unscaled 3d pathlengths
    L_ = Eigen::VectorXd::Zero(N_);

    int N2 = N_/2;
    int iterVal = N2 - 1;

    // Forward section
    for (i = iterVal; i >= 0; i--) {

	j = i + 1;
	
	if (i == iterVal) {
	    L_(i) = 0.0;
	} else {
	    Eigen::VectorXd diff = this->getDiff(i, j);
	    L_(i) = L_(j) + diff.norm();
	}
    }

    // Backward section
    iterVal = N2;

    for (i = iterVal; i < N_; i++) {

	j = i - 1;

	Eigen::VectorXd diff = this->getDiff(i, j);
	double length = diff.norm();

	if (i == iterVal) {
	    L_(i) = -length;
	} else {
	    L_(i) = L_(j) - length;
	}	
    }

    // Now find the ordered pathlengths, starting at 0
    double max_sL = sL_.maxCoeff();
    sdL_ = Eigen::VectorXd::Constant(N_, max_sL) - sL_;

    double max_L = L_.maxCoeff();
    dL_ = Eigen::VectorXd::Constant(N_, max_L) - L_;

}

Eigen::VectorXd LpcPathLength::getDiff(int i, int j)
{
    Eigen::VectorXd diff;

    if (i >= 0 && i < N_ && j >= 0 && j < N_) {
	diff = LAxes_.row(i) - LAxes_.row(j);
    }
    
    return diff;
}
    
LpcPathLength::~LpcPathLength()
{
}

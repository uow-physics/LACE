// Copyright University of Warwick 2014
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// Authors:
// John Back

/*! \file LpcBranch.cc
    \brief Class that stores the result of a branch of the local principal curve
*/

#include "LACE/LpcBranch.hh"

LpcBranch::LpcBranch(int index, const LpcPoint& startPoint, 
		     const std::vector<LpcPoint>& lpcPoints,
		     const Eigen::MatrixXd& eigenVectors, 
		     const Eigen::VectorXd& cosAngles,
		     const LpcPathLength& lpcPath,
		     const Eigen::VectorXd& rho, 
		     const Eigen::VectorXd& c0,
		     const std::vector<LpcPoint>& highRhoPoints,
		     const LpcHitCollection* theHits,
		     int flag) :
    LpcAbsCurve(index, startPoint, lpcPoints, eigenVectors, cosAngles,
		lpcPath, rho, c0, highRhoPoints, theHits, flag)
{
    type_ = LpcAbsCurve::Branch;
}

LpcBranch::LpcBranch() : LpcAbsCurve()
{
}

LpcBranch::~LpcBranch()
{
}

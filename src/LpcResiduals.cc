// Copyright University of Warwick 2014
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// Authors:
// John Back

/*! \file LpcResiduals.cc
    \brief Class that stores the lpc-to-hit residuals
*/

#include "LACE/LpcResiduals.hh"

LpcResiduals::LpcResiduals(const Eigen::VectorXd& lpcResiduals,
			   const Eigen::VectorXd& hitResiduals,
			   const Eigen::VectorXd& lpcWeightRes,
			   const Eigen::VectorXd& hitWeightRes,
			   const Eigen::VectorXi& hitNearestLpc) :
    lpcResiduals_(lpcResiduals),
    hitResiduals_(hitResiduals),
    lpcWeightRes_(lpcWeightRes),
    hitWeightRes_(hitWeightRes),
    hitNearestLpc_(hitNearestLpc)
{
}
 
LpcResiduals::LpcResiduals() :
    lpcResiduals_(),
    hitResiduals_(),
    lpcWeightRes_(),
    hitWeightRes_(),
    hitNearestLpc_()
{
}

LpcResiduals::LpcResiduals(const LpcResiduals& other) :
    lpcResiduals_(other.lpcResiduals_),
    hitResiduals_(other.hitResiduals_),
    lpcWeightRes_(other.lpcWeightRes_),
    hitWeightRes_(other.hitWeightRes_),
    hitNearestLpc_(other.hitNearestLpc_)
{
}

LpcResiduals& LpcResiduals::operator = (const LpcResiduals& other)
{
    // Assignment operator
    lpcResiduals_ = other.lpcResiduals_;
    hitResiduals_ = other.hitResiduals_;
    lpcWeightRes_ = other.lpcWeightRes_;
    hitWeightRes_ = other.hitWeightRes_;
    hitNearestLpc_ = other.hitNearestLpc_;

    return *this;

}

LpcResiduals::~LpcResiduals()
{
}

std::vector<double> LpcResiduals::getLpcResVector() const
{

    std::vector<double> vect;
    for (int i = 0; i < lpcResiduals_.size(); i++) {
	vect.push_back(lpcResiduals_(i));
    }
    return vect;

}

std::vector<double> LpcResiduals::getHitResVector() const
{

    std::vector<double> vect;
    for (int i = 0; i < hitResiduals_.size(); i++) {
	vect.push_back(hitResiduals_(i));
    }
    return vect;

}

std::vector<double> LpcResiduals::getWeightedLpcResVector() const
{

    std::vector<double> vect;
    for (int i = 0; i < lpcWeightRes_.size(); i++) {
	vect.push_back(lpcWeightRes_(i));
    }
    return vect;

}

std::vector<double> LpcResiduals::getWeightedHitResVector() const
{

    std::vector<double> vect;
    for (int i = 0; i < hitWeightRes_.size(); i++) {
	vect.push_back(hitWeightRes_(i));
    }
    return vect;

}

std::vector<int> LpcResiduals::getHitNearestLpcVector() const
{

    std::vector<int> vect;
    for (int i = 0; i < hitNearestLpc_.size(); i++) {
	vect.push_back(hitNearestLpc_(i));
    }
    return vect;

}

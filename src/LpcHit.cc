// Copyright University of Warwick 2014
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// Authors:
// John Back

/*! \file LpcHit.cc
    \brief Class used to define a hit in the dataset or point cloud
*/

#include "LACE/LpcHit.hh"

#include <iostream>

LpcHit::LpcHit() :
    index_(0),
    coords_(Eigen::VectorXd::Zero(3)),
    weight_(0.0),
    nDim_(3),
    truthIds_()
{
}

LpcHit::LpcHit(int id, const Eigen::VectorXd& coords, double weight) :
    index_(id),
    coords_(coords),
    weight_(weight),
    nDim_(coords.size()),
    truthIds_()
{
}

LpcHit::LpcHit(int id, const Eigen::VectorXd& coords, double weight,
	       const Eigen::VectorXi& truthIds) :
    index_(id),
    coords_(coords),
    weight_(weight),
    nDim_(coords.size()),
    truthIds_(truthIds)
{
}

LpcHit::LpcHit(int id, double x, double y, double z, double weight) :
    index_(id),
    coords_(Eigen::VectorXd::Zero(3)),
    weight_(weight),
    nDim_(3),
    truthIds_()
{
    coords_(0) = x;
    coords_(1) = y;
    coords_(2) = z;   
}

LpcHit::LpcHit(int id, double x, double y, double z, double weight, 
	       const Eigen::VectorXi& truthIds) :
    index_(id),
    coords_(Eigen::VectorXd::Zero(3)),
    weight_(weight),
    nDim_(3),
    truthIds_(truthIds)
{
    coords_(0) = x;
    coords_(1) = y;
    coords_(2) = z;   
}

LpcHit::LpcHit(int id, double x, double y, double weight) :
    index_(id),
    coords_(Eigen::VectorXd::Zero(2)),
    weight_(weight),
    nDim_(2),
    truthIds_()
{
    coords_(0) = x;
    coords_(1) = y;
}

LpcHit::LpcHit(int id, double x, double y, double weight, 
	       const Eigen::VectorXi& truthIds) :
    index_(id),
    coords_(Eigen::VectorXd::Zero(2)),
    weight_(weight),
    nDim_(2),
    truthIds_(truthIds)
{
    coords_(0) = x;
    coords_(1) = y;
}

LpcHit::LpcHit(int id, double x, double weight) :
    index_(id),
    coords_(Eigen::VectorXd::Zero(1)),
    weight_(weight),
    nDim_(1),
    truthIds_()
{
    coords_(0) = x;
}

LpcHit::LpcHit(int id, double x, double weight, 
	       const Eigen::VectorXi& truthIds) :
    index_(id),
    coords_(Eigen::VectorXd::Zero(1)),
    weight_(weight),
    nDim_(1),
    truthIds_(truthIds)
{
    coords_(0) = x;
}


LpcHit::~LpcHit()
{
}

LpcHit::LpcHit(const LpcHit& other) :
    index_(other.index_),
    coords_(other.coords_),
    weight_(other.weight_),
    nDim_(other.nDim_),
    truthIds_(other.truthIds_)
{
}

LpcHit& LpcHit::operator = (const LpcHit& other)
{
    // Assignment operator
    index_ = other.index_;
    coords_ = other.coords_;
    weight_ = other.weight_;
    nDim_ = other.nDim_;
    truthIds_ = other.truthIds_;

    return *this;
}

double LpcHit::getFirstCoord() const
{
    double value(0.0);
    if (nDim_ > 0) {value = coords_(0);}
    return value;

}

double LpcHit::getSecondCoord() const
{
    double value(0.0);
    if (nDim_ > 1) {value = coords_(1);}
    return value;

}

double LpcHit::getThirdCoord() const
{
    double value(0.0);
    if (nDim_ > 2) {value = coords_(2);}
    return value;

}

double LpcHit::getCoord(int index) const
{
    double value(0.0);
    if (index >= 0 && index < nDim_) {value = coords_(index);}
    return value;

}

void LpcHit::print() const
{
    std::cout<<"Hit "<<index_<<" has "<<nDim_<<" co-ordinates: ";
    for (int i = 0; i < nDim_; i++) {
	std::cout<<coords_(i)<<", ";
    }
    std::cout<<"weight = "<<weight_<<std::endl;
}

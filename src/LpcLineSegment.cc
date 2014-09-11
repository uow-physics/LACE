// Copyright University of Warwick 2014
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// Authors:
// John Back

/*! \file LpcLineSegment.cc
    \brief Class that defines the line segments used for vertexing/clustering
*/

#include "LACE/LpcLineSegment.hh"

#include <iostream>

LpcLineSegment::LpcLineSegment() :
    centroid_(),
    direction_(),
    residualMean_(0.0),
    residualRms_(0.0),
    residualCut_(0.0),
    endPoint_()
{
}

LpcLineSegment::LpcLineSegment(const Eigen::VectorXd& centroid, 
			       const Eigen::VectorXd& direction,
			       double residualMean, double residualRms, 
			       double residualCut) :
    centroid_(centroid),
    direction_(direction),
    residualMean_(residualMean),
    residualRms_(residualRms),
    residualCut_(residualCut),
    endPoint_(Eigen::VectorXd::Zero(centroid.size()))
{
}

LpcLineSegment::LpcLineSegment(const LpcLineSegment& other) :
    centroid_(other.centroid_),
    direction_(other.direction_),
    residualMean_(other.residualMean_),
    residualRms_(other.residualRms_),
    residualCut_(other.residualCut_),
    endPoint_(other.endPoint_)
{
}

LpcLineSegment& LpcLineSegment::operator = (const LpcLineSegment& other)
{
    // Assignment operator
    centroid_ = other.centroid_;
    direction_ = other.direction_;
    residualMean_ = other.residualMean_;
    residualRms_ = other.residualRms_;
    residualCut_ = other.residualCut_;
    endPoint_ = other.endPoint_;

    return *this;

}

LpcLineSegment::~LpcLineSegment()
{
}

void LpcLineSegment::setEndPoint(const Eigen::VectorXd& endPoint)
{
    endPoint_ = endPoint;

    // Check the direction vector is pointing to this endpoint
    // from the centroid
    Eigen::VectorXd delta = endPoint - centroid_;
    if (delta.dot(direction_) < 0.0) {
	// Reverse the line direction
	direction_ *= -1.0;
    }

}

void LpcLineSegment::print() const
{
    std::cout<<"LineSegment: C = "<<centroid_.transpose()<<", n = "<<direction_.transpose()
	     <<", V = "<<endPoint_.transpose()<<", rMean = "<<residualMean_
	     <<", rRms = "<<residualRms_<<", rCut = "<<residualCut_<<std::endl;
}

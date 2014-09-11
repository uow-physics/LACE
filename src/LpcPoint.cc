// Copyright University of Warwick 2014
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// Authors:
// John Back

/*! \file LpcPoint.cc
    \brief Class that defines the scaled and unscaled co-ordinates of a point in the LpcAbsCurve
*/

#include "LACE/LpcPoint.hh"

LpcPoint::LpcPoint() :
    index_(-1),
    nDim_(3),
    scaledCoords_(Eigen::VectorXd::Zero(3)),
    coords_(Eigen::VectorXd::Zero(3))    
{
}

LpcPoint::LpcPoint(int index, const Eigen::VectorXd& scaledCoords) :
    index_(index),
    nDim_(scaledCoords.size()),
    scaledCoords_(scaledCoords),
    coords_(Eigen::VectorXd::Zero(nDim_))
{
}

LpcPoint::LpcPoint(int index, const Eigen::VectorXd& scaledCoords,
		   const Eigen::VectorXd& coords) :
    index_(index),
    nDim_(scaledCoords.size()),
    scaledCoords_(scaledCoords),
    coords_(coords)
{
}

LpcPoint::LpcPoint(const LpcPoint& other) :
    index_(other.index_),
    nDim_(other.nDim_),
    scaledCoords_(other.scaledCoords_),
    coords_(other.coords_)
{
}

LpcPoint& LpcPoint::operator = (const LpcPoint& other)
{
    index_ = other.index_;
    nDim_ = other.nDim_;
    scaledCoords_ = other.scaledCoords_;
    coords_ = other.coords_;

    return *this;
}

LpcPoint::~LpcPoint()
{
}

void LpcPoint::unscale(const Eigen::VectorXd& theRange)
{

    // Reverse the scaling factor to get the lpc point in 
    // normal co-ordinates
    if (theRange.size() == nDim_) {

	for (int i = 0; i < nDim_; i++) {
	    coords_(i) = scaledCoords_(i)*theRange(i);
	}

    }

}

std::ostream& operator << (std::ostream& os, const LpcPoint& thePoint)
{

    os << "index " << thePoint.getIndex() << ": scaled = (";
    int i;
    int N = thePoint.getNDimensions() - 1;
    Eigen::VectorXd sCoords = thePoint.getScaledCoords();
    Eigen::VectorXd coords = thePoint.getCoords();

    for (i = 0; i < N; i++) {
	os << sCoords(i) << ", ";
    }

    os << sCoords(N) << "); unscaled = (";

    for (i = 0; i < N; i++) {
	os << coords(i) << ", ";
    }

    os << coords(N) << ")";

    return os;

}

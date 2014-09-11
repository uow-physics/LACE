// Copyright University of Warwick 2014
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// Authors:
// John Back

/*! \file LpcVertex.cc
    \brief Simple class that defines an lpc vertex object
*/

#include "LACE/LpcVertex.hh"

#include <iostream>

LpcVertex::LpcVertex() :
    index_(0),
    coords_(Eigen::VectorXd::Zero(3)),
    curveId_(0),
    branchId_(0)
{
}

LpcVertex::LpcVertex(int index, const Eigen::VectorXd& coords,
		     int curveId, int branchId) :
    index_(index),
    coords_(coords),
    curveId_(curveId),
    branchId_(branchId)
{
}

LpcVertex::LpcVertex(const LpcVertex& other) :
    index_(other.index_),
    coords_(other.coords_),
    curveId_(other.curveId_),
    branchId_(other.branchId_)
{
}

LpcVertex& LpcVertex::operator = (const LpcVertex& other)
{
    index_ = other.index_;
    coords_ = other.coords_;
    curveId_ = other.curveId_;
    branchId_ = other.branchId_;

    return *this;

}

LpcVertex::~LpcVertex()
{
}

void LpcVertex::print() const
{

    std::cout<<"Vertex "<<index_<<" ["<<curveId_
	     <<","<<branchId_<<"] = "<<coords_.transpose()<<std::endl;

}

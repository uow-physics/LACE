// Copyright University of Warwick 2014
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// Authors:
// John Back

/*! \file LpcHitCollection.cc
    \brief Class used to define a collection of hit pointers
*/

#include "LACE/LpcHitCollection.hh"

#include "LACE/LpcFunctions.hh"
#include "LACE/LpcHit.hh"

LpcHitCollection::LpcHitCollection(int nDim) :
    nDim_(nDim),
    theHits_(),
    Xi_(),
    scaledXi_(),
    weights_(),
    theRange_(),
    startPoint_()
{
}

LpcHitCollection::LpcHitCollection(const std::vector<LpcHit*>& theHits, int nDim) :
    nDim_(nDim),
    theHits_(theHits),
    Xi_(),
    scaledXi_(),
    weights_(),
    theRange_(),
    startPoint_()
{
}

LpcHitCollection::~LpcHitCollection()
{
    
    std::vector<LpcHit*>::iterator iter;
    for (iter = theHits_.begin(); iter != theHits_.end(); ++iter) {

	delete (*iter);

    }

}

void LpcHitCollection::process(bool applyScaling)
{
    this->storeCoordWeights();
      
    theRange_ = Eigen::VectorXd::Ones(nDim_);
    if (applyScaling) {this->findCoordRange();}

    this->selectStartPoint();
    this->scaleCoords();
}
 
void LpcHitCollection::storeCoordWeights() 
{

    int nHits = this->getNumberOfHits();
    
    Xi_ = Eigen::MatrixXd::Zero(nHits, nDim_);
    weights_ = Eigen::VectorXd::Zero(nHits);

    Eigen::VectorXd origin = Eigen::VectorXd::Zero(nDim_);

    for (int iH = 0; iH < nHits; iH++) {

	LpcHit* aHit = theHits_[iH];

	// The hit should have a valid pointer
	if (aHit != 0) {

	    // Store the hit position and its weight
	    Xi_.row(iH) = aHit->getCoords();
	    weights_(iH) = aHit->getWeight();

	} else {

	    Xi_.row(iH) = origin;
	    weights_(iH) = 0.0;

	}

    }

}

void LpcHitCollection::findCoordRange()
{

    for (int i = 0; i < nDim_; i++) {

	Eigen::VectorXd theColumn = Xi_.col(i);

	double r = theColumn.maxCoeff() - theColumn.minCoeff();

	// If the range is very small, or zero, then set the effective 
	// scaling factor to 1.0, i.e. no scaling is applied
	if (fabs(r) < 1e-6) {r = 1.0;}

	theRange_(i) = r;

    }

}


void LpcHitCollection::selectStartPoint()
{

    // Select the hit nearest to the weighted centroid as the starting point
    LpcFunctions functions;
    Eigen::VectorXd centroid = functions.getWeightedMean(Xi_, weights_);

    Eigen::VectorXd x0 = Eigen::VectorXd::Zero(nDim_);

    // Find the nearest hit to this centroid location
    double minDistSq = 1e10;
    for (int i = 0; i < Xi_.rows(); i++) {

	Eigen::VectorXd thePoint = Xi_.row(i);
	double distSq = (thePoint - centroid).squaredNorm();
	if (distSq < minDistSq) {
	    minDistSq = distSq;
	    x0 = thePoint;
	}

    }

    Eigen::VectorXd scaledx0 = x0.cwiseQuotient(theRange_);

    startPoint_ = LpcPoint(-1, scaledx0, x0);

}

void LpcHitCollection::scaleCoords()
{

    // Scale the hit co-ordinates with their (non-zero) range
    int nHits = this->getNumberOfHits();

    scaledXi_ = Eigen::MatrixXd::Zero(nHits, nDim_);

    // Scale the co-ordinates
    for (int i = 0; i < nDim_; i++) {
	scaledXi_.col(i) = Xi_.col(i)/theRange_(i);
    }

}

LpcHit* LpcHitCollection::getHit(int hitIndex) const
{

    LpcHit* theHit(0);

    if (hitIndex >= 0 && hitIndex < this->getNumberOfHits()) {
	theHit = theHits_[hitIndex];
    }

    return theHit;

}

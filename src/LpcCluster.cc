// Copyright University of Warwick 2014
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// Authors:
// John Back

/*! \file LpcCluster.cc
    \brief Class that defines an lpc cluster object
*/

#include "LACE/LpcCluster.hh"
#include "LACE/LpcFunctions.hh"
#include "LACE/LpcHit.hh"

#include <algorithm>
#include <iostream>

LpcCluster::LpcCluster() :
    index_(0),
    nDim_(3),
    theHits_(),
    hitResiduals_(),
    curveId_(-1),
    branchId_(0),
    lpcRange_(),
    isShower_(false),
    doUpdate_(true),
    centroid_(),
    PCAxes_(),
    PCLengths_(),
    convexHull_(Eigen::VectorXd::Zero(3))
{
}

LpcCluster::LpcCluster(int index, int nDim, 
		       const std::vector<LpcHit*>& theHits,
		       const std::vector<double>& hitResiduals,
		       int curveId, int branchId, 
		       const LpcBinRange& lpcRange,
		       bool isShower) :
    index_(index),
    nDim_(nDim),
    theHits_(theHits),
    hitResiduals_(hitResiduals),
    curveId_(curveId),
    branchId_(branchId),
    lpcRange_(lpcRange),
    isShower_(isShower),
    doUpdate_(true),
    centroid_(),
    PCAxes_(),
    PCLengths_(),
    convexHull_(Eigen::VectorXd::Zero(nDim))
{
}

LpcCluster::LpcCluster(const LpcCluster& other) :
    index_(other.index_),
    nDim_(other.nDim_),
    theHits_(other.theHits_),
    hitResiduals_(other.hitResiduals_),
    curveId_(other.curveId_),
    branchId_(other.branchId_),
    lpcRange_(other.lpcRange_),
    isShower_(other.isShower_),
    doUpdate_(other.doUpdate_),
    centroid_(other.centroid_),
    PCAxes_(other.PCAxes_),
    PCLengths_(other.PCLengths_),
    convexHull_(other.convexHull_)
{
}

LpcCluster& LpcCluster::operator = (const LpcCluster& other)
{
    // Assignment operator
    index_ = other.index_;
    nDim_ = other.nDim_;
    theHits_ = other.theHits_;
    hitResiduals_ = other.hitResiduals_;
    curveId_ = other.curveId_;
    branchId_ = other.branchId_;
    lpcRange_ = other.lpcRange_;
    isShower_ = other.isShower_;
    doUpdate_ = other.doUpdate_;
    nDim_ = other.nDim_;
    centroid_ = other.centroid_;
    PCAxes_ = other.PCAxes_;
    PCLengths_ = other.PCLengths_;
    convexHull_ = other.convexHull_;

    return *this;

}

LpcCluster::~LpcCluster()
{
    // The cluster does not own the LpcHit pointers! They are owned by the
    // LpcHitCollection stored within the event; the LpcHits are deleted 
    // when the event is deleted.
}

void LpcCluster::addHit(LpcHit* theHit, double hitResidual)
{
    if (theHit) {
	theHits_.push_back(theHit);
	hitResiduals_.push_back(hitResidual);
	// We will need to update the calculation of the centroid and PC axes
	doUpdate_ = true;
    }
}

std::vector<int> LpcCluster::getHitIndices() const
{

    std::vector<int> indices;
    std::vector<LpcHit*>::const_iterator hIter;
    for (hIter = theHits_.begin(); hIter != theHits_.end(); ++hIter) {

	int index(-1);

	LpcHit* theHit = *hIter;
	if (theHit) {index = theHit->getIndex();}

	indices.push_back(index);

    }

    // Ascending order for the indices
    std::sort(indices.begin(), indices.end());
    return indices;

}

Eigen::VectorXd LpcCluster::getCentroid()
{
    this->update();
    return centroid_;

}

Eigen::MatrixXd LpcCluster::getPCAxes()
{
    this->update();
    return PCAxes_;

}

Eigen::VectorXd LpcCluster::getPrincipalAxis()
{
    this->update();

    Eigen::VectorXd principalAxis = Eigen::VectorXd::Zero(nDim_);
    // Get the main principal axis (the first eigenvector)
    if (PCAxes_.rows() > 0) {principalAxis = PCAxes_.row(0);}

    return principalAxis;

}

Eigen::VectorXd LpcCluster::getPCLengths()
{
    this->update();
    return PCLengths_;

}

void LpcCluster::update()
{

    if (!doUpdate_) {return;}

    // First, create an Eigen MatrixXd of the hit co-ordinates.
    // We don't want to store this matrix as a private data member
    // since it "repeats" the information already contained in the
    // vector of LpcHit pointers. The matrix is only required so that
    // we can use calculation methods in LpcFunctions.

    Eigen::MatrixXd hitPositions = this->getHitPositions();

    // Find the cluster centroid using the hit positions
    this->calcCentroid(hitPositions);

    // Find the principal component axes (normalised eigenvectors) and
    // their lengths (eigenvalues)
    this->calcPCAxes(hitPositions);

    doUpdate_ = false;

}

Eigen::MatrixXd LpcCluster::getHitPositions() const
{
    
    // Form an Eigen MatrixXd of the hit co-ordinates 
    // associated to the cluster

    // Number of hits in the cluster
    int nHits = theHits_.size();

    if (nHits == 0) {
	// Return an empty one row matrix
	return Eigen::MatrixXd::Zero(1,nDim_);
    }

    // Create and initialise the required hit position matrix
    Eigen::MatrixXd hitPositions = Eigen::MatrixXd::Zero(nHits, nDim_);

    // Loop over the hits and store the position within each row
    for (int iH = 0; iH < nHits; iH++) {

	LpcHit* theHit = theHits_[iH];
	if (theHit) {

	    hitPositions.row(iH) = theHit->getCoords();

	}

    }

    return hitPositions;

}

std::vector<double> LpcCluster::getHitWeights() const
{

    std::vector<double> weights;

    int nHits = theHits_.size();
    for (int iH = 0; iH < nHits; iH++) {

	LpcHit* theHit = theHits_[iH];
	if (theHit) {
	    weights.push_back(theHit->getWeight());
	}

    }

    return weights;

}

void LpcCluster::calcCentroid(const Eigen::MatrixXd& hitPositions)
{

    LpcFunctions functions;
    centroid_ = functions.getMean(hitPositions);

}

void LpcCluster::calcPCAxes(const Eigen::MatrixXd& hitPositions)
{

    LpcFunctions functions;

    // Mean centred data
    Eigen::MatrixXd meanData = functions.offsetPositions(hitPositions, centroid_);

    // Covariance matrix
    Eigen::MatrixXd covMatrix = functions.formCovarianceMatrix(meanData);

    // Get the eigenvalues and normalised eigenvectors
    std::pair<Eigen::VectorXd, Eigen::MatrixXd> eigenPair = 
	functions.findNormEigenVectors(covMatrix);

    // Set the internal data members
    PCLengths_ = eigenPair.first;
    PCAxes_ = eigenPair.second;

}

void LpcCluster::print()
{
    this->update();
    std::cout<<"Cluster "<<index_<<" ["<<curveId_<<","<<branchId_<<"], C = "
	     <<centroid_.transpose()<<", axes = "<<PCAxes_<<std::endl;

}

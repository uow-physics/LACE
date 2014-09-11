// Copyright University of Warwick 2014
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// Authors:
// John Back

/*! \file LpcProcess.cc
    \brief Class defining the process of finding curves for events
*/

#include "LpcmRec/LpcProcess.hh"

#include "LpcmRec/LpcAlgorithm.hh"
#include "LpcmRec/LpcBranchAlgorithm.hh"
#include "LpcmRec/LpcBranchCollection.hh"
#include "LpcmRec/LpcCluster.hh"
#include "LpcmRec/LpcClusterData.hh"
#include "LpcmRec/LpcCurve.hh"
#include "LpcmRec/LpcEvent.hh"
#include "LpcmRec/LpcFeatures.hh"
#include "LpcmRec/LpcLineClusterAlgorithm.hh"
#include "LpcmRec/LpcParameters.hh"
#include "LpcmRec/LpcPoint.hh"
#include "LpcmRec/LpcShowerAlgorithm.hh"
#include "LpcmRec/LpcSingleClusterAlgorithm.hh"
#include "LpcmRec/LpcVertex.hh"

#include <iostream>

LpcProcess::LpcProcess(const LpcParameters* theParameters) :
    theAlgorithm_(new LpcAlgorithm(theParameters)),
    theBranchAlgorithm_(new LpcBranchAlgorithm(theParameters)),
    theFeatures_(new LpcFeatures(theParameters)),
    theClusterAlgorithm_(0),
    theShowerAlgorithm_(new LpcShowerAlgorithm(theParameters)),
    clustering_(1),
    branchLevel_(0),
    doScaling_(1)
{
    if (theParameters) {

	clustering_ = theParameters->getClustering();
	branchLevel_ = theParameters->getBranchLevel();

	if (clustering_ == 1) {
	    theClusterAlgorithm_ = new LpcLineClusterAlgorithm(theParameters);
	} else if (clustering_ == 2) {
	    theClusterAlgorithm_ = new LpcSingleClusterAlgorithm(theParameters);
	}

	doScaling_ = theParameters->getScalingFlag();

    }
}

LpcProcess::~LpcProcess()
{
    // Clean-up the pointers (in reverse order of their declaration)
    delete theBranchAlgorithm_;
    delete theAlgorithm_;
    delete theFeatures_;
    delete theClusterAlgorithm_;
    delete theShowerAlgorithm_;
}

void LpcProcess::reconstruct(LpcEvent* theEvent)
{

    if (!theEvent) {
	std::cerr<<"Error in LpcProcess::reconstruct. The event is null"<<std::endl;
	return;
    }

    // Set-up the event so that the co-ordinates of the hit collection have been
    // stored in Eigen matrices & scaled; this also finds the initial starting point
    bool applyScaling(true);
    if (doScaling_ == 0) {applyScaling = false;}

    theEvent->setUp(applyScaling);

    // Retrieve the hit collection from the event
    LpcHitCollection* theHits = theEvent->getHitCollection();

    // Do the lpc calculation and obtain the main curve result
    int mainCurveIndex(0);
    LpcCurve* theCurve = theAlgorithm_->getCurve(mainCurveIndex, theHits, LpcAlgorithm::Both);

    if (theCurve != 0) {

	// Find any branches for this main curve
	if (branchLevel_ > 0) {

	    std::vector<LpcPoint> highRhoPoints = theCurve->getHighRhoPoints();

	    LpcBranchCollection* theBranches = 
		theBranchAlgorithm_->getBranches(theHits, highRhoPoints, mainCurveIndex);

	    theCurve->storeBranches(theBranches);

	}
	
 
	// Find the feature points (peaks in 1-|cosAngle|) used for identifying 
	// vertices. This stores the peak ranges within the curve (and its branches)
	theFeatures_->findFeatures(theCurve);

	// Store the curve in the event
	theEvent->storeCurve(theCurve);	

	// Find the vertices and obtain clusters for the main curve and any of its branches
	if (clustering_ != 0) {

	    LpcClusterData clusterData = theClusterAlgorithm_->findClusters(theCurve);

	    // Store the cluster and vertexing information in the event. The event
	    // will then be responsible for the LpcCluster and LpcVertex pointers,
	    // i.e. they will be deleted when the event is deleted
	    std::vector<LpcVertex*> theVertices = clusterData.getVertices();
	    std::vector<LpcCluster*> theClusters = clusterData.getClusters();

	    theEvent->storeVertices(theVertices);
	    theEvent->storeClusters(theClusters);

	    // Determine which clusters could be showers
	    theShowerAlgorithm_->findShowers(theClusters);

	}

    }

}


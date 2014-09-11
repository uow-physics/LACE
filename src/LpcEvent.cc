// Copyright University of Warwick 2014
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// Authors:
// John Back

/*! \file LpcEvent.cc
    \brief Class used to define an event, which contains a collection of hit pointers
*/

#include "LACE/LpcEvent.hh"

#include "LACE/LpcCluster.hh"
#include "LACE/LpcCurve.hh"
#include "LACE/LpcHitCollection.hh"
#include "LACE/LpcVertex.hh"

LpcEvent::LpcEvent(int eventNumber) : 
    eventNumber_(eventNumber),
    theHits_(0),
    theCurves_(),
    theVertices_(),
    theClusters_()
{
}

LpcEvent::LpcEvent(int eventNumber, LpcHitCollection* theHits) :
    eventNumber_(eventNumber),
    theHits_(theHits),
    theCurves_(),
    theVertices_(),
    theClusters_()
{
}

LpcEvent::~LpcEvent()
{
    // Delete the hit collection
    delete theHits_;

    // Delete the LpcCurve pointers
    std::vector<LpcCurve*>::iterator iter;
    for (iter = theCurves_.begin(); iter != theCurves_.end(); ++iter) {
	delete (*iter);
    }

    // Delete the vertex pointers
    std::vector<LpcVertex*>::iterator vIter;
    for (vIter = theVertices_.begin(); vIter != theVertices_.end(); ++vIter) {
	delete (*vIter);
    }

    // Delete the cluster pointers
    std::vector<LpcCluster*>::iterator cIter;
    for (cIter = theClusters_.begin(); cIter != theClusters_.end(); ++cIter) {
	delete (*cIter);
    }

}

void LpcEvent::setUp(bool applyScaling)
{
    // Make sure the hit collection knows about scaled co-ordinates
    // and also find the starting point for the lpc algorithm
    if (theHits_) {theHits_->process(applyScaling);}

}

LpcCurve* LpcEvent::getCurve() const
{
    // Get the first curve in the event
    LpcCurve* theCurve(0);

    if (this->getNLpcCurves() > 0) {
	theCurve = theCurves_[0];
    }

    return theCurve;
}

void LpcEvent::storeClusters(const std::vector<LpcCluster*>& clusterVect)
{
    // Add the LpcCluster pointers to the internal vector
    theClusters_.insert(theClusters_.end(), clusterVect.begin(), clusterVect.end());
}

void LpcEvent::storeVertices(const std::vector<LpcVertex*>& vertexVect)
{
    // Add the LpcVertex pointers to the internal vector
    theVertices_.insert(theVertices_.end(), vertexVect.begin(), vertexVect.end());
}

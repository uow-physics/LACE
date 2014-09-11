// Copyright University of Warwick 2014
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// Authors:
// John Back

/*! \file LpcAbsOutput.cc
    \brief Abstract class to store the results of the lpc algorithm
*/

#include "LACE/LpcAbsOutput.hh"

#include "LACE/LpcCluster.hh"
#include "LACE/LpcCurve.hh"
#include "LACE/LpcHitCollection.hh"
#include "LACE/LpcEvent.hh"
#include "LACE/LpcVertex.hh"

LpcAbsOutput::LpcAbsOutput(const std::string& outputFileName) :
    outputFileName_(outputFileName),
    theEvent_(0),
    theHits_(0),
    nDim_(3)
{
}

LpcAbsOutput::~LpcAbsOutput()
{
}

void LpcAbsOutput::store(LpcEvent* theEvent)
{

    if (!theEvent) {return;}

    theEvent_ = theEvent;
    theHits_ = theEvent->getHitCollection();
    if (theHits_) {nDim_ = theHits_->getNDimensions();}

    this->storeInitialInfo();

    const LpcCurve* theCurve = theEvent->getCurve();
    this->storeCurve(theCurve);

    const std::vector<LpcVertex*> theVertices = theEvent->getVertices();
    this->storeVertices(theVertices);

    const std::vector<LpcCluster*> theClusters = theEvent->getClusters();
    this->storeClusters(theClusters);

    this->storeExtraInfo();

}

// Copyright University of Warwick 2014
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// Authors:
// John Back

/*! \file LpcSingleClusterAlgorithm.cc
    \brief Class that creates simple clusters by just including all hits
*/

#include "LACE/LpcSingleClusterAlgorithm.hh"

#include "LACE/LpcCluster.hh"
#include "LACE/LpcCurve.hh"
#include "LACE/LpcHit.hh"
#include "LACE/LpcHitCollection.hh"
#include "LACE/LpcParameters.hh"
#include "LACE/LpcResiduals.hh"
#include "LACE/LpcVertex.hh"

#include <Eigen/Dense>
#include <vector>

LpcSingleClusterAlgorithm::LpcSingleClusterAlgorithm(const LpcParameters* thePars) : 
    LpcAbsClusterAlgorithm(thePars)
{
}

LpcSingleClusterAlgorithm::~LpcSingleClusterAlgorithm()
{
}

LpcClusterData LpcSingleClusterAlgorithm::findClusters(const LpcCurve* theCurve)
{

    // Check if the curve pointer is valid
    LpcClusterData theResults;

    if (!theCurve) {return theResults;}

    // Consider all hits as one cluster
    const LpcHitCollection* theHits = theCurve->getHitCollection();

    if (!theHits) {return theResults;}

    LpcResiduals lpcResiduals = theCurve->getResiduals();
    Eigen::VectorXd hitResiduals = lpcResiduals.getHitResiduals();

    // Get the vector of LpcHit pointers from the hit collection
    std::vector<LpcHit*> allHits = theHits->getHits();
    // Store the hit-to-nearest-lpc residuals. This will be used for track/shower id
    std::vector<double> allResiduals;

    std::vector<LpcHit*>::const_iterator iter;
    for (iter = allHits.begin(); iter != allHits.end(); ++iter) {

	LpcHit* theHit = *iter;

	if (theHit) {

	    int hitIndex = theHit->getIndex();

	    double hitResidual = hitResiduals(hitIndex);
	    allResiduals.push_back(hitResidual);

	}

    }

    // Create the single LpcCluster pointer
    int index(0), branchId(0);
    int nDim = theHits->getNDimensions();
    int curveId = theCurve->getIndex();
    // Set the cluster lpc point range to include all points
    int nLpcPoints = theCurve->getNLpcPoints();
    LpcBinRange lpcRange(0, nLpcPoints-1);
    bool isShower(false);

    LpcCluster* oneCluster = new LpcCluster(index, nDim, allHits, allResiduals,
					    curveId, branchId, lpcRange, isShower);

    std::vector<LpcVertex*> theVertices;
    std::vector<LpcCluster*> theClusters;
    theClusters.push_back(oneCluster);

    return LpcClusterData(theVertices, theClusters);

}

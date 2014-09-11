// Copyright University of Warwick 2014
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// Authors:
// John Back

/*! \file LpcShowerAlgorithm.cc
    \brief Class that defines the shower algorithm
*/

#include "LACE/LpcShowerAlgorithm.hh"

#include "LACE/LpcCluster.hh"
#include "LACE/LpcParameters.hh"
#include "LACE/LpcConvexHull.hh"

#include <algorithm>
#include <cmath>

LpcShowerAlgorithm::LpcShowerAlgorithm(const LpcParameters* thePars) :
    thePars_(thePars),
    convexHull_(0.0),
    showerRes_(0.0),
    showerResRatio_(0.0),
    showerResFrac_(0.0)
{
    if (thePars_) {
	convexHull_ = thePars_->getConvexHull();
	showerRes_ = thePars->getShowerRes();
	showerResRatio_ = thePars->getShowerResRatio();
	showerResFrac_ = thePars->getShowerResFrac();
    }

}

LpcShowerAlgorithm::~LpcShowerAlgorithm()
{
}

void LpcShowerAlgorithm::findShowers(const std::vector<LpcCluster*>& theClusters)
{

    // Loop over the vector of LpcCluster pointers
    std::vector<LpcCluster*>::const_iterator iter;
    for (iter = theClusters.begin(); iter != theClusters.end(); ++iter) {

	LpcCluster* aCluster = *iter;

	if (aCluster) {

	    // Check if the cluster is a shower (and set its shower boolean)
	    this->identify(aCluster);
	    
	}

    }

}

void LpcShowerAlgorithm::identify(LpcCluster* theCluster)
{

    if (!theCluster) {return;}

    bool result(false);

    LpcConvexHull qHull;
    qHull.findLengths(theCluster);

    Eigen::VectorXd hullLengths = theCluster->getConvexHull();

    // Form the convex hull ratio: 
    // sum of transverse (y,z) lengths/longitudinal length (x)
    double chRatio(0.0);

    int nHull = hullLengths.size();
    double chx(0.0), chy(0.0), chz(0.0);

    if (nHull == 2 || nHull == 3) {

	chx = hullLengths(0);
	chy = hullLengths(1);
	
	if (nHull == 3) {chz = hullLengths(2);}
	if (fabs(chx) > 1e-10) {chRatio = (chy + chz)/chx;}

    }

    // Check if the convex hull ratio is large enough
    if (chRatio > convexHull_) {

	// Retrieve the hit-to-lpc residuals for the cluster (in ascending order)
	std::vector<double> hitResiduals = theCluster->getHitResiduals();
	std::sort(hitResiduals.begin(), hitResiduals.end());

	// Get the maximum hit-to-lpc residual
	int nRes = hitResiduals.size();
	double maxHitRes = hitResiduals[nRes-1];

	// Only consider hit residuals that are larger than
	// showerResRatio*maxHitRes. If the minimum value of this is above
	// showerRes, and at least showerResFrac of hits remain,
	// then the lpc cluster is a shower

	double minResidual = showerResRatio_*maxHitRes;
	double nPassRes(0.0), nTotRes(0.0);
           
	// Find out how many residuals are above minResidual, and how
	// many of these pass showerRes_
 
	for (int i = 0; i < nRes; i++) {

	    double resValue = hitResiduals[i];
	    if (resValue > minResidual) {

		nTotRes += 1.0;

		if (resValue > showerRes_) {nPassRes += 1.0;}

	    }

	} // Loop over the residuals

	// Find the fraction of residuals that meet the requirements
	double fraction(0.0);
	if (nTotRes > 0.0) {fraction = nPassRes/nTotRes;}

	// Check to see if the fraction is above the threshold for a shower
	if (fraction > showerResFrac_) {result = true;}

    }

    // Set the shower boolean decision
    theCluster->setShowerBool(result);

}

// Copyright University of Warwick 2014
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// Authors:
// John Back

/*! \file LpcLineClusterAlgorithm.cc
    \brief Class that defines the (extended) line vertexing algorithm
*/

#include "LACE/LpcLineClusterAlgorithm.hh"

#include "LACE/LpcAbsCurve.hh"
#include "LACE/LpcBranch.hh"
#include "LACE/LpcBranchCollection.hh"
#include "LACE/LpcCluster.hh"
#include "LACE/LpcCurve.hh"
#include "LACE/LpcHit.hh"
#include "LACE/LpcHitCollection.hh"
#include "LACE/LpcParameters.hh"
#include "LACE/LpcVertex.hh"

#include <cmath>
#include <iostream>

LpcLineClusterAlgorithm::LpcLineClusterAlgorithm(const LpcParameters* thePars) : 
    LpcAbsClusterAlgorithm(thePars),
    theHits_(0),
    lineFitter_(),
    nLpcSegment_(2),
    nResSigma_(3.0),
    minVtxResCut_(20.0),
    maxDist_(100.0),
    doClusterMerging_(true),
    doBranchVtx_(1),
    mergeDegAngle_(20.0),
    mergeCosAngle_(0.0),
    functions_(),
    curveId_(0),
    branchId_(0),
    nHitPoints_(0),
    nDim_(3),
    nLpcPoints_(0),
    nLpcPoints1_(0),
    hitPoints_(),
    lpcPoints_(),
    lpcResiduals_(),
    hitResiduals_(),
    hitNearLpc_(),
    usedHits_(),
    nPeaks_(0),
    cosPeakRanges_(),
    lpcRanges_(),
    lineRanges_(),
    lineSegments_(),
    extLineSegments_(),
    mergedOK_(false)
{

    mergeCosAngle_ = fabs(cos(mergeDegAngle_/57.29578));
    if (thePars) {
	minVtxResCut_ = thePars->getMinVtxResCut();
	doBranchVtx_ = thePars->getBranchVtx();
    }

}

LpcLineClusterAlgorithm::~LpcLineClusterAlgorithm()
{
}

LpcClusterData LpcLineClusterAlgorithm::findClusters(const LpcCurve* theCurve)
{

    // Check if the curve pointer is valid
    if (!theCurve) {return LpcClusterData();}

    // The hits are always the same for the main curve or its branches
    theHits_ = theCurve->getHitCollection();

     // The unscaled hit co-ordinates and number of co-ordinate dimensions
    if (theHits_) {
	hitPoints_ = theHits_->getCoords();
	nHitPoints_ = theHits_->getNumberOfHits();
	nDim_ = theHits_->getNDimensions();
    }

    // An array to keep track of which hits have already been used for clusters
    usedHits_ = Eigen::VectorXi::Zero(nHitPoints_);

    LpcClusterData mainResults = this->processCurve(theCurve);

    LpcBranchCollection* branchCollection = theCurve->getBranchCollection();

    if (branchCollection && doBranchVtx_) {

	int nLevels = branchCollection->getNumberLevels();
	    
	// Loop over the number of generations (starts at 1 not 0)
	for (int iL = 1; iL <= nLevels; iL++) {
	    
	    std::vector<LpcBranch*> branches = branchCollection->getBranches(iL);
	    std::vector<LpcBranch*>::const_iterator bIter;
	    for (bIter = branches.begin(); bIter != branches.end(); ++bIter) {
		
		const LpcBranch* theBranch = *bIter;

		if (!theBranch) {continue;}

		LpcClusterData branchResults = this->processCurve(theBranch);
		// Add the branch results to the main results
		mainResults.addData(branchResults);

	    }
	    
	}
	
    } 

    // Check to see if there are any clusters. If not, then create a single 
    // all-encompassing cluster from the list of hits
    int nClusters = mainResults.getNumberClusters();

    if (nClusters < 1) {

	int clusterIndex(0);
	LpcCluster* oneCluster = this->createSingleCluster(clusterIndex);
	// Specify the vertex as the position of the hit with the lowest 
	// value of the first co-ordinate ("x")
	LpcVertex* oneVertex = this->createSingleVertex();
	mainResults.addCluster(oneCluster);
	mainResults.addVertex(oneVertex);

    }

    // Set the indices of the clusters and vertices to represent
    // the order they appear in the vectors
    mainResults.resetIndices();

    return mainResults;

}

LpcClusterData
LpcLineClusterAlgorithm::processCurve(const LpcAbsCurve* theCurve)
{

    LpcClusterData theResults;
    if (!theCurve) {return theResults;}

    int curveIndex = theCurve->getIndex();
    int curveFlag = theCurve->getFlag();
    curveId_ = curveIndex;
    branchId_ = 0;

    if (curveFlag > -1) {
	// We have a branch. Set the curveId to the main curve index number,
	// which is equal to the flag value, and set the branchId to the
	// index number of the branch
	curveId_ = curveFlag;
	branchId_ = curveIndex;
    }
    
    // The vector of lpc points
    lpcPoints_ = theCurve->getLpcPoints();
    nLpcPoints_ = theCurve->getNLpcPoints();
    nLpcPoints1_ = nLpcPoints_ - 1;
	
    // The lpc residuals
    lpcResiduals_ = theCurve->getResiduals();
	
    // Get the Eigen vector containing the lpc residual distances for each hit
    hitResiduals_ = lpcResiduals_.getHitResiduals();

    // Get the Eigen vector containing the closest lpc point integer for each hit
    hitNearLpc_ = lpcResiduals_.getHitNearestLpc();

    // Get the cosine peak ranges for the feature points
    cosPeakRanges_ = theCurve->getCosPeakRanges();
    nPeaks_ = cosPeakRanges_.size();

    // If there are no peaks, stop
    if (nPeaks_ < 1) {return theResults;}

    // Find the lpc ranges that straddle the lpc feature points
    this->findLpcRanges();
    
    // Find the line segments that are near the feature points
    this->findLineRanges();
    this->findLineSegments();
    
    // Find the initial vertices. Matrix has row = vertex, col = co-ordinates (x,y,z)
    Eigen::MatrixXd initVtxPoints = this->findVtxPoints(lineSegments_);
    
    // Extend the line segments, using the initial vertex position
    this->findExtendedLines(initVtxPoints);
    
    // Find the updated vertex positions
    Eigen::MatrixXd extVtxPoints = this->findVtxPoints(extLineSegments_);

    // Create clusters based on the associated hits, together with
    // the nearest hits along the rest of the lpc point range
    std::vector<LpcCluster*> lpcClusters = this->formClusters(extVtxPoints);
    
    // Add any left over hits to the (nearest) cluster(s)
    this->addRemainingHits(lpcClusters);

    // Create the object storing the cluster and vertex pointer information
    LpcClusterData mergedVtxClusters = this->mergeClusters(extVtxPoints, lpcClusters);

    if (mergedOK_) {
	// We have created a new, merged set of LpcCluster pointers.
	// Delete the original LpcCluster pointers
	std::vector<LpcCluster*>::iterator cIter;
	for (cIter = lpcClusters.begin(); cIter != lpcClusters.end(); ++cIter) {
	    delete (*cIter);
	}
    }

    // Return the merged clusters and vertices
    return mergedVtxClusters;

}

void LpcLineClusterAlgorithm::findLpcRanges()
{

    // Find the ranges of the lpc points that straddle the feature points
    // Clear the internal vector of lpc ranges
    lpcRanges_.clear();

    if (nPeaks_ > 0) {

	LpcBinRange firstRange = cosPeakRanges_[0];
	int firstMinBin = firstRange.getMinBin();
	if (firstMinBin > 0) {

	    lpcRanges_.push_back(LpcBinRange(0, firstMinBin));

	}
	
	if (nPeaks_ == 1) {

	    lpcRanges_.push_back(LpcBinRange(firstRange.getMaxBin(), nLpcPoints1_));

	} else {

	    LpcBinRange prevRange = firstRange;
	    LpcBinRange nextRange;

	    for (int iP = 1; iP < nPeaks_; iP++) {
		nextRange = cosPeakRanges_[iP];
		lpcRanges_.push_back(LpcBinRange(prevRange.getMaxBin(), nextRange.getMinBin()));
		prevRange = nextRange;
	    }

	    int nextMaxBin = nextRange.getMaxBin();
	    if (nextMaxBin < nLpcPoints1_) {
		lpcRanges_.push_back(LpcBinRange(nextMaxBin, nLpcPoints1_));
	    }

	} // nPeaks > 1

    } // nPeaks > 0

}

void LpcLineClusterAlgorithm::findLineRanges()
{
    // For each feature point, find the lpc point range numbers that are "nLpcSegment"
    // lpc points away. This allows us to define line segments for vertex finding.
    // The cosine peak ranges should already be ordered.
    lineRanges_.clear();
    
    for (int iP = 0; iP < nPeaks_; iP++) {

	LpcBinRange peakRange = cosPeakRanges_[iP];

	int minBin = peakRange.getMinBin();
	int maxBin = peakRange.getMaxBin();

	int lowP1 = this->checkLpcInt(minBin - nLpcSegment_);
	int highP1 = this->checkLpcInt(minBin);

	int lowP2 = this->checkLpcInt(maxBin);
	int highP2 = this->checkLpcInt(maxBin + nLpcSegment_);

	lineRanges_.push_back(LpcBinRange(lowP1, highP1));
	lineRanges_.push_back(LpcBinRange(lowP2, highP2));

    }

}

void LpcLineClusterAlgorithm::findLineSegments()
{

    // Find the straight line segments that pass through each 
    // of the line ranges that are next to the lpc feature points.
    lineSegments_.clear();

    // For each lpc point range, retrieve the closest set of data hits
    // that have residuals below a maximum threshold.
    // Then, find the centroid position and fit a best straight line
    // segment to these same data hits

    for (size_t iR = 0; iR < lineRanges_.size(); iR++) {

	// Lpc point line range
	LpcBinRange theRange = lineRanges_[iR];
	int startRange = theRange.getMinBin();
	int endRange = theRange.getMaxBin();
	int nRange = endRange - startRange + 1;
	
	// Keep track of the hit positions that are closest to the lpc 
	// curve "line segment". Also store the lpc-to-hit residuals

	// Create vector containers for the hits that are closest to each lpc point
	std::vector<int> matchedHitInt;

	// Loop over the lpc point integers in this line range
	for (int iL = 0; iL < nRange; iL++) {

	    // Lpc point integer
	    int iLpc = startRange + iL;

	    // Loop over the data hits
	    for (int iH = 0; iH < nHitPoints_; iH++) {		

		// Get the integer of the closest lpc point for this hit
		int iHitLpc = hitNearLpc_(iH);

		// See if the lpc point integers match
		if (iHitLpc == iLpc) {

		    // Store the hit number
		    matchedHitInt.push_back(iH);

		}

	    } // Loop over data hits


	} // Loop over the lpc points within the given line range

	// Create an Eigen vector of the matched hit residuals
	int nMatchedHits = matchedHitInt.size();
	Eigen::VectorXd matchedHitRes = Eigen::VectorXd::Zero(nMatchedHits);

	for (int jH = 0; jH < nMatchedHits; jH++) {

	    int hitIndex = matchedHitInt[jH];
	    matchedHitRes(jH) = hitResiduals_(hitIndex);

	}
	
	// Now find the residual average and its root-mean-square (standard deviation)
	std::pair<double, double> resMeanAndRms = functions_.getMeanAndRms(matchedHitRes);
	double residualMean = resMeanAndRms.first;
	double residualRms = resMeanAndRms.second;

	// Define the maximum selection requirement on the hit residuals
	double residualCut = this->getResidualCut(residualMean, residualRms);

	// Loop over the hits and keep those with the required residuals
	std::vector<int> keptHitIndices;
	for (int jH = 0; jH < nMatchedHits; jH++) {

	    int hitInt = matchedHitInt[jH];
	    double hitRes = matchedHitRes(jH);
	    if (hitRes < residualCut) {
		keptHitIndices.push_back(hitInt);
	    }

	}
	
	int nKeptHits = keptHitIndices.size();
	// Define an Eigen MatrixXd that stores the kept hit coords; row = hit, col = x,y,z,..
	Eigen::MatrixXd matchedHitPos = Eigen::MatrixXd::Zero(nKeptHits, nDim_);
	for (int jH = 0; jH < nKeptHits; jH++) {

	    int hitInt = keptHitIndices[jH];
	    matchedHitPos.row(jH) = hitPoints_.row(hitInt);

	}

	// Check to see if there are no kept hits. In this case, just use
	// the lpc points directly
	if (nKeptHits == 0) {

	    // Re-initalise the matchedHitPos matrix with updated dimensions
	    matchedHitPos = Eigen::MatrixXd::Zero(nRange, nDim_);

	    for (int iL = 0; iL < nRange; iL++) {

		// Lpc point integer
		int iLpc = startRange + iL;
		LpcPoint lpcPoint = lpcPoints_[iLpc];
		matchedHitPos.row(iL) = lpcPoint.getCoords();

	    }

	}

	// Find the centroid and unit direction vector for the straight line section
	// that minimises the perpendicular Euclidean distance to each hit
	std::pair<Eigen::VectorXd, Eigen::VectorXd> theLine = lineFitter_.findLine(matchedHitPos);

	// Store the line segment information
	if (residualCut < minVtxResCut_) {residualCut = minVtxResCut_;}
	LpcLineSegment theSegment(theLine.first, theLine.second, residualMean,
				  residualRms, residualCut);
	lineSegments_.push_back(theSegment);

    } // Loop over the line ranges

}

Eigen::MatrixXd 
LpcLineClusterAlgorithm::findVtxPoints(const std::vector<LpcLineSegment>& lineSegments) const
{

    // Find the point of closest approach between neighbouring ling segments,
    // using linear vector algebra.

    // Loop over the line segments
    int nLines = lineSegments.size();
    int nVtx = nLines/2;
    double small(1e-10);

    // Define a matrix of the vertex points: row = vertex, col = co-ordinates (x,y,z)
    Eigen::MatrixXd theVtxPoints = Eigen::MatrixXd::Zero(nVtx, nDim_);

    for (int iL = 0; iL < nVtx; iL++) {

	int iL2 = 2*iL;
	LpcLineSegment line1 = lineSegments[iL2];
	LpcLineSegment line2 = lineSegments[iL2+1];

	Eigen::VectorXd P0 = line1.getCentroid();
	Eigen::VectorXd Q0 = line2.getCentroid();
	Eigen::VectorXd u  = line1.getDirection();
	Eigen::VectorXd v  = line2.getDirection();

	// Various dot products
	Eigen::VectorXd w0 = P0 - Q0;
	double a = u.dot(u);
	double b = u.dot(v);
	double c = v.dot(v);
	double d = u.dot(w0);
	double e = v.dot(w0);

	// Initialise the lengths along the 2 lines which will be "optimised"
	// to find the closest distance of approach
	double sc(0.0), tc(0.0);

	// Denominator factor for sc and tc terms
	double D = a*c - b*b;

	// First check if the two lines are parallel
	if (D < small) {

	    if (b > c && fabs(b) > small) {
		tc = d/b;
	    } else if (fabs(c) > small) {
		tc = e/c;
	    }

	} else {

	    sc = (b*e - c*d)/D;
	    tc = (a*e - b*d)/D;

	}

	// Find the point of closest approach
	Eigen::VectorXd D0 = P0 + Q0;
	Eigen::VectorXd du = sc*u;
	Eigen::VectorXd dv = tc*v;
	Eigen::VectorXd dL = du + dv;

	Eigen::VectorXd vtx = 0.5*(D0 + dL);

	theVtxPoints.row(iL) = vtx;

    }

    return theVtxPoints;

}

void LpcLineClusterAlgorithm::findExtendedLines(const Eigen::MatrixXd& initVtxPoints)
{

    // Extend the line segments by "hoovering up" hits that lie inside
    // cylinders whose axes are defined by the initial line segments 
    // (lineSegments_ containing the centroids and direction unit vectors)

    // Clear the vector of extended line segments
    extLineSegments_.clear();

    // Set the end-points of the given line segments as the initial
    // vertex positions. This will automatically change the sign of the 
    // direction vectors to point to these vertices
    
    int nSegments = lineSegments_.size();
    int iL(0);
    for (iL = 0; iL < nSegments; iL++) {

	Eigen::VectorXd V0 = initVtxPoints.row(iL/2);
	lineSegments_[iL].setEndPoint(V0);

    }

    // Find all hits that are associated to the initial line segments.
    // This is a dictionary where key = line segment number and the
    // value is a vector of integers of the matched hit indices
    bool uniqueHits(false);
    LpcLineIntMap hitMatchDict = this->findAssociatedHits(lineSegments_, uniqueHits);
    
    // Loop over the initial line segments and obtain new line fits
    // that use the extended set of hits
    for (iL = 0; iL < nSegments; iL++) {

	LpcLineSegment theLine = lineSegments_[iL];

	// All matched hits for the extended line
	std::vector<int> extMatchedHits = hitMatchDict[iL];
	int nExtHits = extMatchedHits.size();

	if (nExtHits > 0) {

	    // Define an Eigen MatrixXd that stores the required hit coords; row = hit, col = x,y,z,..
	    Eigen::MatrixXd extHitPositions = Eigen::MatrixXd::Zero(nExtHits, nDim_);
	    for (int kH = 0; kH < nExtHits; kH++) {

		int hitInt = extMatchedHits[kH];
		extHitPositions.row(kH) = hitPoints_.row(hitInt);

	    }
 
	    // Find the centroid and unit direction vector for the extended straight line section
	    // that minimises the perpendicular Euclidean distance to each hit
	    std::pair<Eigen::VectorXd, Eigen::VectorXd> extLine = lineFitter_.findLine(extHitPositions);

	    LpcLineSegment extSegment(extLine.first, extLine.second, theLine.getResidualMean(),
				      theLine.getResidualRms(), theLine.getResidualCut());

	    extLineSegments_.push_back(extSegment);


	} else {

	    // Keep the original line
	    extLineSegments_.push_back(theLine);

	}

    } // Loop over initial line segments
  
}

LpcLineClusterAlgorithm::LpcLineIntMap 
LpcLineClusterAlgorithm::findAssociatedHits(const std::vector<LpcLineSegment>& lineSegments, 
					    bool uniqueHits)
{

    LpcLineIntMap hitMatchDict;

    // The number of line segments we have
    int nLines = lineSegments.size();

    // Create Eigen objects to store the centroids, direction vectors
    // and residual cuts, which are needed as we loop over all hits (we want
    // to avoid re-obtaining these as we iterate through the hit loop)
    Eigen::MatrixXd centroids = Eigen::MatrixXd::Zero(nLines, nDim_);
    Eigen::MatrixXd nVectors = Eigen::MatrixXd::Zero(nLines, nDim_);
    Eigen::VectorXd residualCuts = Eigen::VectorXd::Zero(nLines);

    // Store the dot product of (V - C)dot n, where V = end point (vertex),
    // C = centroid, n = direction vector
    Eigen::VectorXd VCnArray = Eigen::VectorXd::Zero(nLines);
    
    // Retrieve the centroid and direction vector for the line segments.
    // Reverse the direction vector if it is not pointing towards its
    // corresponding endpoint/vertex position

    for (int jL = 0; jL < nLines; jL++) {

	LpcLineSegment lSegment = lineSegments[jL];
	Eigen::VectorXd C = lSegment.getCentroid();
	Eigen::VectorXd n = lSegment.getDirection();
	Eigen::VectorXd V = lSegment.getEndPoint();

	Eigen::VectorXd VmC = V - C;
	double VCn = VmC.dot(n);

	centroids.row(jL) = C;
	if (VCn < 0.0) {
	    // Reverse the direction vector
	    nVectors.row(jL) = -n;
	    VCnArray(jL) = -VCn;
	} else {
	    nVectors.row(jL) = n;
	    VCnArray(jL) = VCn;
	}

	residualCuts(jL) = lSegment.getResidualCut();

    }

    // Loop over all hits and find those that can extend the initial line segments
    for (int iH = 0; iH < nHitPoints_; iH++) {

	// Get the hit co-ordinates
	Eigen::VectorXd P = hitPoints_.row(iH);

	// Find out which line segment this hit is closest to
	// (within its residual cut)

	int matchedLine(-1);
	double minDelta(1e10);

	// Loop over the line segments
	for (int iL = 0; iL < nLines; iL++) {

	    // Get the centroid and direction vector of the line segment
	    Eigen::VectorXd C = centroids.row(iL);
	    Eigen::VectorXd n = nVectors.row(iL);
	    double VCn = VCnArray(iL);

	    // Perpendicular distance of the hit from the line segment.
	    // This makes sure that the point P is pointing towards the end point
	    // of the line V (in the direction of n), but does not go beyond V, i.e.
	    // the distance PC is less than VC
	    double delta = functions_.getPerpLineDist(P, C, n, VCn);

	    // The cut on the residuals
	    double resCut = residualCuts(iL);

	    // Check if the hit is close to the line, and update the minimum delta so far
	    if (delta > -1.0 && delta < resCut && delta < minDelta) {
		matchedLine = iL;
		minDelta = delta;

	    }

	} // Loop over line segments
   
	// Append the hit integer for the given line segment
	if (matchedLine != -1) {

	    hitMatchDict[matchedLine].push_back(iH);
	    if (uniqueHits) {
		usedHits_(iH) = 1;
	    }

	}

    } // Loop over data hits


    return hitMatchDict;

}

std::vector<LpcCluster*> LpcLineClusterAlgorithm::formClusters(const Eigen::MatrixXd& vtxMatrix)
{

    // Form clusters from the hits associated to the various
    // lpc curve sections separated by feature points/vertices
    std::vector<LpcCluster*> theClusters;

    // Set the end points for the extended lines as the new vertex positions
    int iL(0);
    int nExtLines = extLineSegments_.size();
    for (iL = 0; iL < nExtLines; iL++) {

	Eigen::VectorXd V = vtxMatrix.row(iL/2);
	extLineSegments_[iL].setEndPoint(V);

    }

    // Find the matched hits for the extended line segment
    bool uniqueHits(true);
    LpcLineIntMap hitMatchDict1 = this->findAssociatedHits(extLineSegments_, uniqueHits);

    // Also find hits that lie between the initial and extended line segments.
    // There should always be the same number of initial and extended line segments
    std::vector<LpcLineSegment> interLineSegments;
    for (iL = 0; iL < nExtLines; iL++) {

	LpcLineSegment initLine = lineSegments_[iL];
	LpcLineSegment extLine = extLineSegments_[iL];

	Eigen::VectorXd interC = initLine.getCentroid();
	Eigen::VectorXd interV = extLine.getCentroid();
	Eigen::VectorXd interDir = interV - interC;
	double interNorm = interDir.norm();
	if (interNorm > 0.0) {interDir /= interNorm;}

	double resMean = extLine.getResidualMean();
	double resRms = extLine.getResidualRms();
	double resCut = extLine.getResidualCut();

	LpcLineSegment interLine(interC, interDir, resMean, resRms, resCut);
	interLine.setEndPoint(interV);

	interLineSegments.push_back(interLine);

    }
    
    LpcLineIntMap hitMatchDict2 = this->findAssociatedHits(interLineSegments, uniqueHits);
    
    int nVtx = nExtLines/2;
    int nClusters = nVtx + 1;

    // First create lpcClusters from the initial and extended line segments, then
    // look at the hits closest to the lpc curve range and only add those
    // with residuals within the maximum residual cut

    int iC(0), jC(0);

    for (iC = 0; iC < nClusters; iC++) {

	std::vector<int> matchedHits;

	if (iC == 0 || iC == nVtx) {

	    // First cluster only uses the first line segment
	    // Last cluster only uses the last line segment
	    jC = 0;
	    if (iC == nVtx) {jC = nExtLines - 1;}

	    matchedHits = hitMatchDict1[jC];
	    std::vector<int> hitInts2_jC = hitMatchDict2[jC];

	    matchedHits.insert(matchedHits.end(), hitInts2_jC.begin(), hitInts2_jC.end());
	    
	} else {

	    // Each cluster is formed from two line segments
	    jC = 2*iC;
	    int jC1 = jC - 1;

	    matchedHits = hitMatchDict1[jC1];
	    std::vector<int> hitInts1_jC = hitMatchDict1[jC];
	    std::vector<int> hitInts2_jC1 = hitMatchDict2[jC1];
	    std::vector<int> hitInts2_jC = hitMatchDict2[jC];

	    matchedHits.insert(matchedHits.end(), hitInts1_jC.begin(), hitInts1_jC.end());
	    matchedHits.insert(matchedHits.end(), hitInts2_jC1.begin(), hitInts2_jC1.end());
	    matchedHits.insert(matchedHits.end(), hitInts2_jC.begin(), hitInts2_jC.end());

	}

	double residualCut = extLineSegments_[jC].getResidualCut();

	LpcBinRange lpcRange = lpcRanges_[iC];
	this->formClusterHits(matchedHits, lpcRange, residualCut);

	// From the matchedHits integers, store the full hit information in the cluster
	LpcCluster* theCluster = this->createCluster(iC, matchedHits, lpcRange);
	theClusters.push_back(theCluster);

    } // Loop over clusters

    return theClusters;

}

void LpcLineClusterAlgorithm::formClusterHits(std::vector<int>& matchedHits,
					      const LpcBinRange& lpcRange,
					      double residualCut)
{

    // Create a complete list of cluster hits, starting with those 
    // already found from the extended line hits (matchedLineHits). 
    // Add the nearest hits, within a typical residual cut value, for the 
    // rest of the lpc point range. Skip hits that have already been processed

    // Lpc point line range
    int startRange = lpcRange.getMinBin();
    int endRange = lpcRange.getMaxBin();
    int nRange = endRange - startRange + 1;

    // Loop over the lpc point integers in this line range
    for (int iL = 0; iL < nRange; iL++) {

	// The lpc point integer
	int iLpc = startRange + iL;

	// Loop over the data hits
	for (int iH = 0; iH < nHitPoints_; iH++) {

	    // Check to see if the hit is already matched to another line segment
	    if (usedHits_(iH) == 1) {continue;}

	    // Get the integer of the closest lpc point
	    int iHitLpc = hitNearLpc_(iH);

	    // See if the lpc point integers match, and if the residual is small enough
	    if (iHitLpc == iLpc) {

		double hitLpcRes = hitResiduals_(iH);

		if (hitLpcRes < residualCut) {

		    matchedHits.push_back(iH);
		    usedHits_(iH) = 1;

		} // Is the residual small enough?

	    } // Do the lpc point indices match?

	} // Loop over data hits

    } // Lpc point integer loop

}

LpcCluster* LpcLineClusterAlgorithm::createCluster(int index, 
						   const std::vector<int>& matchedHits,
						   const LpcBinRange& lpcRange) const
{

    // From the vector of matched hit indices, retrieve the pointers of
    // LpcHits that we need for the cluster. Also store the residuals
    // for the hits associated to the cluster
    if (!theHits_) {return 0;}

    std::vector<LpcHit*> clusterHits;
    std::vector<double> clusterResiduals;
    
    std::vector<int>::const_iterator iter;
    for (iter = matchedHits.begin(); iter != matchedHits.end(); ++iter) {

	int hitIndex = *iter;

	LpcHit* theHit = theHits_->getHit(hitIndex);
	clusterHits.push_back(theHit);

	double hitRes = hitResiduals_(hitIndex);
	clusterResiduals.push_back(hitRes);

    }

    // Create and return the LpcCluster pointer
    bool isShower(false);
    LpcCluster* theCluster = new LpcCluster(index, nDim_, clusterHits, clusterResiduals,
					    curveId_, branchId_, lpcRange, isShower);
    return theCluster;

}

void LpcLineClusterAlgorithm::addRemainingHits(const std::vector<LpcCluster*>& lpcClusters)
{

    // Add any remaining hits to the clusters by finding the closest distance 
    // to the cluster principal axis

    LpcFunctions functions;
    int nClusters = lpcClusters.size();

    // First, retrieve the current centroid and direction vectors for the clusters
    // before we add any furher hits to them (which will change the centroid/principal axes)

    // Matrix of the centroids: row = centroid, col = co-ordinates x,y,z...
    Eigen::MatrixXd CMatrix = Eigen::MatrixXd::Zero(nClusters, nDim_);

    // Matrix of the direction vectors: row = centroid, col = co-ordinates x,y,z...
    Eigen::MatrixXd nMatrix = Eigen::MatrixXd::Zero(nClusters, nDim_);

    for (int j = 0; j < nClusters; j++) {

	LpcCluster* theCluster = lpcClusters[j];
	if (theCluster) {

	    Eigen::VectorXd C = theCluster->getCentroid();
	    Eigen::VectorXd n = theCluster->getPrincipalAxis();

	    CMatrix.row(j) = C;
	    nMatrix.row(j) = n;

	}

    }

    // First, loop over all of the data hits
    for (int iH = 0; iH < nHitPoints_; iH++) {

	// Process hits that have not yet been used
	if (usedHits_(iH) == 0) {

	    // Get the hit co-ordinates
	    Eigen::VectorXd P =  hitPoints_.row(iH);

	    int clusterInt(-1);
	    double minDist(1e10);

	    // Loop over the clusters to find the closest one
	    for (int iC = 0; iC < nClusters; iC++) {

		Eigen::VectorXd C = CMatrix.row(iC);
		Eigen::VectorXd n = nMatrix.row(iC);
		    
		double dist = functions.getPerpLineDist(P, C, n);

		if (dist < minDist && dist < maxDist_) {
		    minDist = dist;
		    clusterInt = iC;
		}

	    } // Loop over clusters

	    if (clusterInt != -1) {

		// We have found the closest cluster for the hit point.
		// Add the hit and its hit-to-lpc residual to the cluster
		LpcCluster* closestCluster = lpcClusters[clusterInt];
		if (closestCluster) {

		    LpcHit* theHit = theHits_->getHit(iH);
		    double hitResidual = hitResiduals_(iH);
		    closestCluster->addHit(theHit, hitResidual);

		    // Set the used flag for this hit
		    usedHits_(iH) = 1;

		}

	    }

	} // Has the hit been already used?

    } // Loop over hits
 
}

std::vector<LpcVertex*> 
LpcLineClusterAlgorithm::createVertices(const Eigen::MatrixXd& vtxPointMatrix) const
{

    std::vector<LpcVertex*> theVertices;
    
    int nVtx = vtxPointMatrix.rows();
    for (int iV = 0; iV < nVtx; iV++) {

	LpcVertex* vertex = new LpcVertex(iV, vtxPointMatrix.row(iV), curveId_, branchId_);
	theVertices.push_back(vertex);

    }

    return theVertices;
}

LpcClusterData
LpcLineClusterAlgorithm::mergeClusters(const Eigen::MatrixXd& vtxPointMatrix,
				       const std::vector<LpcCluster*>& lpcClusters)
{

    mergedOK_ = false;
    int nVtx = vtxPointMatrix.rows();
    int nClusters = lpcClusters.size();

    std::vector<LpcVertex*> mergedVertices;
    std::vector<LpcCluster*> mergedClusters;

    // Note that we must have nClusters = nVtx + 1 for the rest of the function logic.
    // If we have the wrong number of vertices/clusters, or we don't want merging, just
    // return the unchanged clusters and form vertices from vtxPointMatrix
    if (!doClusterMerging_ || nClusters != nVtx+1) {

	mergedVertices = this->createVertices(vtxPointMatrix);
	// Copy the pointers to the original LpcClusters. This does not
	// create new LpcClusters, but just copies the pointer info
	mergedClusters = std::vector<LpcCluster*>(lpcClusters);
	return LpcClusterData(mergedVertices, mergedClusters);

    }

    // A map of the kept cluster index number and the vector of cluster indices
    // that need to be merged
    LpcLineIntMap keptClusters;
    int nKept(0);

    // Loop over the vertices
    double precision(1.0e-10);

    for (int iV = 0; iV < nVtx; iV++) {

	// Get the clusters on either side of the vertex
	int iV1 = iV + 1;
	LpcCluster* cluster1 = lpcClusters[iV];
	LpcCluster* cluster2 = lpcClusters[iV+1];

	// Get the normalised principal major axis of the two clusters
	Eigen::VectorXd pca1 = Eigen::VectorXd::Zero(nDim_);
	if (cluster1) {pca1 = cluster1->getPrincipalAxis();}

	Eigen::VectorXd pca2 = Eigen::VectorXd::Zero(nDim_);
	if (cluster2) {pca2 = cluster2->getPrincipalAxis();}

	// Calculate the cosine magnitude of the angle between the axes
	double pcaDot = fabs(pca1.dot(pca2));

	bool validPCA(true);
	// Check if any of the PCA vectors is zero
	if (pca1.isZero(precision) || pca2.isZero(precision)) {validPCA = false;}

	// Always keep the first cluster
	if (iV == 0) {keptClusters[nKept].push_back(iV);}

	// Check to see if the angle between the principal axes is large enough
	if (validPCA == true && pcaDot < mergeCosAngle_) {

	    // Keep the vertex and store the 2nd cluster number with a new nKept value
	    LpcVertex* mergedVtx = new LpcVertex(nKept, vtxPointMatrix.row(iV),
						 curveId_, branchId_);
	    mergedVertices.push_back(mergedVtx);
	    nKept += 1;
	    keptClusters[nKept].push_back(iV1);

	} else {

	    // Discard the vertex. Add the 2nd cluster to the pre-existing one
	    // for this vertex (nKept is not incremented)
	    keptClusters[nKept].push_back(iV1);

	}

    }

    int nKeptClusters = keptClusters.size();

    // Loop over the kept cluster map and merge clusters if a kept cluster
    // entry has more than one cluster associated with it
    for (int iK = 0; iK < nKeptClusters; iK++) {

	std::vector<int> mergedHitIds;
	int minLpcRange(0), maxLpcRange(0);

	// For each kept cluster map entry, get the cluster index numbers
	std::vector<int> clusterIds = keptClusters[iK];
	int nIndices = clusterIds.size();

	for (int iC = 0; iC < nIndices; iC++) {

	    int cIndex = clusterIds[iC];
	    LpcCluster* theCluster = lpcClusters[cIndex];
	    if (theCluster) {

		std::vector<int> clusterHitIds = theCluster->getHitIndices();

		// Append the hit indices to the merged vector
		mergedHitIds.insert(mergedHitIds.end(), clusterHitIds.begin(), clusterHitIds.end());

		LpcBinRange binRange = theCluster->getLpcRange();
		if (iC == 0) {minLpcRange = binRange.getMinBin();}
		maxLpcRange = binRange.getMaxBin();		

	    } // Valid cluster pointer

	} // Loop over cluster indices in the kept map entry

	// Create the merged cluster
	LpcCluster* mergedCluster = this->createCluster(iK, mergedHitIds,
							LpcBinRange(minLpcRange, maxLpcRange));

	mergedClusters.push_back(mergedCluster);

    }

    mergedOK_ = true;
    return LpcClusterData(mergedVertices, mergedClusters);

}


int LpcLineClusterAlgorithm::checkLpcInt(int lpcInt) const
{
    // Check whether the lpc point integer is valid    
    int theInt(lpcInt);
    if (theInt < 0) {
	theInt = 0;
    } else if (theInt > nLpcPoints1_) {
	theInt = nLpcPoints1_;
    }

    return theInt;

}

LpcCluster* LpcLineClusterAlgorithm::createSingleCluster(int index) const
{

    // Store all hits as one cluster
    if (!theHits_) {return 0;}

    // Get the vector of LpcHit pointers from the hit collection
    std::vector<LpcHit*> allHits = theHits_->getHits();
    // Store the hit-to-nearest-lpc residuals. This will be used for track/shower id
    std::vector<double> allResiduals;

    std::vector<LpcHit*>::const_iterator iter;
    for (iter = allHits.begin(); iter != allHits.end(); ++iter) {

	LpcHit* theHit = *iter;

	if (theHit) {

	    int hitIndex = theHit->getIndex();

	    double hitResidual = hitResiduals_(hitIndex);
	    allResiduals.push_back(hitResidual);

	}

    }

    bool isShower(false);
    // Set the cluster lpc point range to include all points
    LpcBinRange lpcRange(0, nLpcPoints1_);
    int branchId(0);

    // Create and return the LpcCluster pointer
    LpcCluster* theCluster = new LpcCluster(index, nDim_, allHits, allResiduals,
					    curveId_, branchId, lpcRange, isShower);

    return theCluster;
    

}

LpcVertex* LpcLineClusterAlgorithm::createSingleVertex() const
{

    int minXRow(0);
    double minX(0.0);

    for (int i = 0; i < nHitPoints_; i++) {

	// All hit co-ordinates should have at least one dimension, 
	// so retrieving the first value via (0) should be safe
	Eigen::VectorXd hP = hitPoints_.row(i);
	if (i == 0) {minX = hP(0);}

	if (hP(0) < minX) {
	    minXRow = i;
	    minX = hP(0);
	}

    }

    int index(0), branchId(0);
    Eigen::VectorXd vtxPoint = hitPoints_.row(minXRow);

    LpcVertex* theVertex = new LpcVertex(index, vtxPoint, curveId_, branchId);
    return theVertex;

}

double LpcLineClusterAlgorithm::getResidualCut(double residualMean, double residualRms) const
{

    double residualCut = nResSigma_*residualRms + residualMean;
    return residualCut;

}

void LpcLineClusterAlgorithm::printBinRanges(const std::string& preamble, 
					     const std::vector<LpcBinRange>& binRanges) const
{

    std::cout<<preamble;
    for (size_t j = 0; j < binRanges.size(); j++) {
	LpcBinRange b = binRanges[j];
	std::cout<<"("<<b.getMinBin()<<", "<<b.getMaxBin()<<"); ";
    }
    std::cout<<std::endl;

}

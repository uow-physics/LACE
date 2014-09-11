// Copyright University of Warwick 2014
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// Authors:
// John Back

/*! \file LpcTextOutput.cc
    \brief Class to write out the results of the lpc algorithm into a text file
*/

#include "LACE/LpcTextOutput.hh"

#include "LACE/LpcAbsCurve.hh"
#include "LACE/LpcBranchCollection.hh"
#include "LACE/LpcBranch.hh"
#include "LACE/LpcCluster.hh"
#include "LACE/LpcCurve.hh"
#include "LACE/LpcEvent.hh"
#include "LACE/LpcPathLength.hh"
#include "LACE/LpcPoint.hh"
#include "LACE/LpcResiduals.hh"
#include "LACE/LpcVertex.hh"

#include <Eigen/Dense>
#include <iomanip>

LpcTextOutput::LpcTextOutput(const std::string& outputFileName, int precision) :
    LpcAbsOutput(outputFileName),
    writeData_(outputFileName_.c_str()),
    precision_(precision)
{
    this->initialise();
}

LpcTextOutput::~LpcTextOutput()
{
}

void LpcTextOutput::initialise()
{

    // Set the stream output precision    
    writeData_.precision(precision_);

    // Set the width between output columns
    // = "-x.precision" & "e+yy" & extra space
    width_ = precision_ + 10;

}

void LpcTextOutput::storeInitialInfo()
{

    if (theEvent_) {
	
	writeData_ << "Event " << theEvent_->getEventNumber() << std::endl;

    }

}

void LpcTextOutput::storeCurve(const LpcCurve* mainCurve)
{

    // Store the main curve and all of its branches    
    if (!mainCurve) {return;}

    int curveId = mainCurve->getIndex();
    writeData_ << "MainCurve:Id " << curveId << std::endl;
    this->storeCurveDetails(mainCurve);

    // Retrieve the branches
    LpcBranchCollection* theBranches = mainCurve->getBranchCollection();
    
    if (theBranches) {

	// Get the number of generational levels for the branches
	int nLevels = theBranches->getNumberLevels();

	// Loop over each branch generation
	for (int i = 0; i < nLevels; i++) {

	    int iG = i + 1; // Generation number
	    // Get the vector of branches
	    std::vector<LpcBranch*> branchVect = theBranches->getBranches(iG);

	    // Loop over the branches
	    std::vector<LpcBranch*>::const_iterator iter;
	    for (iter = branchVect.begin(); iter != branchVect.end(); ++iter) {

		const LpcBranch* theBranch = *iter;
		// Write out the information stored in the branch
		if (theBranch) {

		    writeData_ << "Branch:Id " << theBranch->getIndex() << std::endl;
		    writeData_ << "Branch:MainCurveId " << curveId << std::endl;
		    writeData_ << "Branch:GenerationLevel " << iG << std::endl;

		    this->storeCurveDetails(theBranch);

		} // Check if branch pointer exists

	    } // Iteration over branches

	} // Branch generation loop
	
    } // Check for any stored branches

}

void LpcTextOutput::storeCurveDetails(const LpcAbsCurve* theCurve)
{

    // Store either the main curve or branch, following the info in
    // LpcAbsCurve::print()

    if (!theCurve) {return;}

    // Starting point
    LpcPoint startPoint = theCurve->getStartPoint();
    this->printPoint("LpcStartPoint", startPoint);

    // LpcPoints
    std::vector<LpcPoint> lpcPoints = theCurve->getLpcPoints();
    int nLpc = lpcPoints.size();
    writeData_ << "NumberOfLpcPoints " << nLpc << std::endl;

    int i(0);

    for (i = 0; i < nLpc; i++) {
	std::string word("LpcPoint");
	this->printPoint(word, i, lpcPoints[i]);
    }

    // Eigenvectors
    Eigen::MatrixXd eigenVectors = theCurve->getEigenVectors();
    this->printMatrix("LpcEigenVectors", eigenVectors);

    // Cosine angles
    Eigen::VectorXd cosAngles = theCurve->getCosAngles();
    this->printVector("LpcCosineAngles", cosAngles);

    // Rho
    Eigen::VectorXd rho = theCurve->getRho();
    this->printVector("EigenRatioRho", rho);

    // c0
    Eigen::VectorXd c0 = theCurve->getc0();
    this->printVector("c0", c0);

    // HighRhoPoints
    std::vector<LpcPoint> highRhoPoints = theCurve->getHighRhoPoints();
    int nHRho = highRhoPoints.size();

    writeData_ << "HighRhoPoints " << nHRho << std::endl;
    for (i = 0; i < nHRho; i++) {
	std::string word("HighRhoPoint");
	this->printPoint(word, i, highRhoPoints[i]);
    }


    // Pathlengths object
    LpcPathLength pathLength = theCurve->getPathLength();

    // Lambda
    Eigen::VectorXd lambda = pathLength.getLambda();
    this->printVector("LpcLambda", lambda);

    // Lambda axes
    Eigen::MatrixXd lambdaAxes = pathLength.getLambdaAxes();
    this->printMatrix("LpcLambdaAxes", lambdaAxes);

    // Delta lambda
    Eigen::VectorXd deltaLambda = pathLength.getDeltaLambda();
    this->printVector("LpcPathLength", deltaLambda);

    // Residuals object
    LpcResiduals residuals = theCurve->getResiduals();

    // Lpc point residuals averaged over all nearest hits
    Eigen::VectorXd lpcRes = residuals.getLpcResiduals();
    this->printVector("MeanLpcResiduals", lpcRes);

    // Hit-to-nearest-lpc point residuals
    Eigen::VectorXd hitRes = residuals.getHitResiduals();
    this->printVector("HitResiduals", hitRes);

    // Weighted lpc point residuals averaged over all nearest hits
    Eigen::VectorXd wLpcRes = residuals.getWeightedLpcResiduals();
    this->printVector("WeightMeanLpcResiduals", wLpcRes);

    // Weighted hit-to-nearest-lpc point residuals
    Eigen::VectorXd wHitRes = residuals.getWeightedHitResiduals();
    this->printVector("WeightHitResiduals", wHitRes);

    // The index number of the nearest lpc point for each hit
    Eigen::VectorXi hitNearLpc = residuals.getHitNearestLpc();
    this->printVector("HitNearestLpc", hitNearLpc);

    // Feature points (cosine peaks)
    std::vector<int> cosPeakIndices = theCurve->getCosPeakIndices();
    this->printVector("CosPeaks", cosPeakIndices);

}

void LpcTextOutput::storeVertices(const std::vector<LpcVertex*>& theVertices)
{

    int nVertices = theVertices.size();
    writeData_ << "NumberOfVertices " << nVertices << std::endl;

    for (int i = 0; i < nVertices; i++) {

	LpcVertex* vertex = theVertices[i];
	if (vertex) {

	    int vtxIndex = vertex->getIndex();
	    int curveId = vertex->getCurveId();
	    int branchId = vertex->getBranchId();

	    Eigen::VectorXd vtxPoint = vertex->getCoords();

	    writeData_ << "Vertex:Index " << vtxIndex << " "
		       << "Vertex:MainCurveIndex " << curveId << " "
		       << "Vertex:BranchIndex " << branchId << std::endl;

	    this->printCoords("Vertex:Position", vtxPoint);
	    
	}

    }

}

void LpcTextOutput::storeClusters(const std::vector<LpcCluster*>& theClusters)
{

    int nClusters = theClusters.size();
    writeData_ << "NumberOfClusters " << nClusters << std::endl;

    for (int i = 0; i < nClusters; i++) {

	LpcCluster* cluster = theClusters[i];
	if (cluster) {
	
	    int cIndex = cluster->getIndex();

	    int curveId = cluster->getCurveId();
	    int branchId = cluster->getBranchId();

	    LpcBinRange lpcRange = cluster->getLpcRange();

	    writeData_ << "Cluster:Index " << cIndex << std::endl;

	    writeData_ << "Cluster:MainCurveIndex " << curveId << " "
		       << "Cluster:BranchIndex " << branchId << " "
		       << "Cluster:LpcPointRange " << lpcRange.getMinBin()
		       << " " << lpcRange.getMaxBin() << std::endl;

	    Eigen::VectorXd centroid = cluster->getCentroid();
	    this->printCoords("Cluster:Centroid", centroid);

	    Eigen::MatrixXd PCAxes = cluster->getPCAxes();
	    this->printMatrix("Cluster:PrincipalAxes", PCAxes);

	    Eigen::VectorXd convexHull = cluster->getConvexHull();
	    this->printCoords("Cluster:ConvexHull", convexHull);

	    bool isShower = cluster->isAShower();
	    writeData_ << "Cluster:IsShower " << std::boolalpha << isShower << std::endl;

	    std::vector<int> hitIndices = cluster->getHitIndices();
	    this->printVector("Cluster:HitIndices", hitIndices);


	}

    }

}

void LpcTextOutput::storeExtraInfo()
{

    // Add an extra empty line to separate out the next event info
    writeData_ << std::endl;

}

void LpcTextOutput::printPoint(const std::string& preamble, const LpcPoint& thePoint)
{

    Eigen::VectorXd coords = thePoint.getCoords();
    int nDim1 = coords.size() - 1;

    writeData_ << std::left << std::setw(width_) << preamble << " ";

    for (int i = 0; i < nDim1; i++) {

	writeData_ << std::left << std::setw(width_) << coords(i) << " ";

    }

    writeData_ << std::left << std::setw(width_) << coords(nDim1) << std::endl;

}

void LpcTextOutput::printPoint(const std::string& preamble, int index, const LpcPoint& thePoint)
{

    Eigen::VectorXd coords = thePoint.getCoords();
    int nDim1 = coords.size() - 1;

    writeData_ << std::left << std::setw(width_) << preamble;
    writeData_ << " " << std::left << std::setw(width_) << index << " ";

    for (int i = 0; i < nDim1; i++) {

	writeData_ << std::left << std::setw(width_) << coords(i) << " ";

    }

    writeData_ << std::left << std::setw(width_) << coords(nDim1) << std::endl;

}

void LpcTextOutput::printCoords(const std::string& preamble, const Eigen::VectorXd& coords)
{

    int nDim1 = coords.size() - 1;
    writeData_ << preamble << " ";

    for (int i = 0; i < nDim1; i++) {

	writeData_ << std::left << std::setw(width_) << coords(i) << " ";

    }

    writeData_ << std::left << std::setw(width_) << coords(nDim1) << std::endl;

}

void LpcTextOutput::printVector(const std::string& preamble, const Eigen::VectorXd& vector)
{

    int N = vector.size();
    writeData_ << preamble << " " << N << std::endl;

    for (int i = 0; i < N; i++) {

	writeData_ << std::left << std::setw(width_) << vector(i);

	if (i%10 == 9) {
	    writeData_ << std::endl;
	} else {
	    writeData_ << " ";
	}

    }

    if (N%10 != 0) {writeData_ << std::endl;}

}

void LpcTextOutput::printVector(const std::string& preamble, const Eigen::VectorXi& vector)
{

    int N = vector.size();
    writeData_ << preamble << " " << N << std::endl;

    for (int i = 0; i < N; i++) {

	writeData_ << std::left << std::setw(width_) << vector(i);

	if (i%10 == 9) {
	    writeData_ << std::endl;
	} else {
	    writeData_ << " ";
	}

    }

    if (N%10 != 0) {writeData_ << std::endl;}

}

void LpcTextOutput::printMatrix(const std::string& preamble, const Eigen::MatrixXd& matrix)
{

    int nRow = matrix.rows();
    int nCol = matrix.cols();    
    int nCol1 = nCol - 1;

    writeData_ << preamble << " " << nRow << " " << nCol << std::endl;

    for (int i = 0; i < nRow; i++) {
	
	for (int j = 0; j < nCol1; j++) {

	    writeData_ << std::left << std::setw(width_) 
		       << matrix(i, j) << " ";

	}

	writeData_ << std::left << std::setw(width_) 
		   << matrix(i, nCol1) << std::endl;

    }    

}

void LpcTextOutput::printVector(const std::string& preamble, const std::vector<int>& vector)
{

    int N = vector.size();
    writeData_ << preamble << " " << N <<std::endl;

    for (int i = 0; i < N; i++) {

	writeData_ << std::left << std::setw(width_) << vector[i];

    	if (i%10 == 9) {
	    writeData_ << std::endl;
	} else {
	    writeData_ << " ";
	}

    }

    if (N%10 != 0) {writeData_ << std::endl;}

}

#ifdef LPC_USE_ROOT

// Copyright University of Warwick 2014
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// Authors:
// John Back

/*! \file LpcRootOutput.cc
    \brief Class to write out the results of the lpc algorithm into a ROOT file
*/

#include "LACE/LpcRootOutput.hh"

#include "LACE/LpcAbsCurve.hh"
#include "LACE/LpcBranchCollection.hh"
#include "LACE/LpcBranch.hh"
#include "LACE/LpcCluster.hh"
#include "LACE/LpcCurve.hh"
#include "LACE/LpcEvent.hh"
#include "LACE/LpcHit.hh"
#include "LACE/LpcHitCollection.hh"
#include "LACE/LpcPathLength.hh"
#include "LACE/LpcPoint.hh"
#include "LACE/LpcResiduals.hh"
#include "LACE/LpcVertex.hh"

#include "TFile.h"
#include "TString.h"
#include "TTree.h"

#include <Eigen/Dense>
#include <iostream>

LpcRootOutput::LpcRootOutput(const std::string& outputFileName) :
    LpcAbsOutput(outputFileName),
    rootFile_(0),
    lpcTree_(0),
    vtxTree_(0),
    clTree_(0),
    definedTrees_(false),
    eventId_(0), curveId_(0), branchId_(0), branchGen_(0),
    lpcStartPoint_(), nLpc_(0), lpcPoints_(), lpcEigenVectors_(),
    lpcCosAngles_(0), lpcRho_(0), lpcC0_(0), 
    lpcLambda_(0), lpcLambdaAxes_(),
    lpcPathLength_(0), lpcResiduals_(0), wLpcResiduals_(0),
    nHighRho_(0), highRhoPoints_(), nCosPeaks_(0), lpcCosPeaks_(),
    nHits_(0), hitCoords_(), hitResiduals_(0), 
    wHitResiduals_(0), hitNearLpc_(0), nVtx_(0),
    vtxCurveId_(0), vtxBranchId_(0), vtxPoint_(), nCl_(0),
    clIndex_(0), clCurveId_(0), clBranchId_(0), clMinLpc_(0),
    clMaxLpc_(0), centroid_(), axes_(), convexHull_(),
    isAShower_(0), nClHits_(0), clHitIds_(0), 
    clHitCoords_(), clHitWeights_()
{
    this->initialise();
}

void LpcRootOutput::initialise()
{

    // Create the output file
    rootFile_ = TFile::Open(outputFileName_.c_str(), "recreate");

    // Lpc info
    lpcTree_ = new TTree("lpcTree", "lpcTree");
    lpcTree_->SetDirectory(rootFile_);

    // Vertexing info
    vtxTree_ = new TTree("vtxTree", "vtxTree");
    vtxTree_->SetDirectory(rootFile_);

    // Cluster info
    clTree_ = new TTree("clsTree", "clsTree");
    clTree_->SetDirectory(rootFile_);

}

LpcRootOutput::~LpcRootOutput()
{
    // The file and the internal vectors/maps should be deleted
    // when finalise is called, so do nothing here
}

void LpcRootOutput::finalise()
{
    // Write out the trees
    rootFile_->cd();
    lpcTree_->Write();
    vtxTree_->Write();
    clTree_->Write();

    // Close the file
    rootFile_->Close(); 

    // Delete all internal vectors and maps
    this->cleanUp();

}

void LpcRootOutput::cleanUp()
{

    this->cleanUpMap(lpcPoints_);
    this->cleanUpMap(lpcEigenVectors_);
    this->cleanUpVector(lpcCosAngles_);
    this->cleanUpVector(lpcRho_);
    this->cleanUpVector(lpcC0_);
    this->cleanUpVector(lpcLambda_);
    this->cleanUpMap(lpcLambdaAxes_);
    this->cleanUpVector(lpcPathLength_);
    this->cleanUpVector(lpcResiduals_);
    this->cleanUpVector(wLpcResiduals_);

    this->cleanUpMap(highRhoPoints_);

    this->cleanUpVector(lpcCosPeaks_);

    this->cleanUpMap(hitCoords_);
    this->cleanUpVector(hitResiduals_);
    this->cleanUpVector(wHitResiduals_);
    this->cleanUpVector(hitNearLpc_);

    this->cleanUpMap(axes_);
    this->cleanUpVector(clHitIds_);
    this->cleanUpMap(clHitCoords_);
    this->cleanUpVector(clHitWeights_);

}

void LpcRootOutput::cleanUpMap(LpcRootMap& theMap)
{
    
    LpcRootMap::const_iterator iter;
    for (iter = theMap.begin(); iter != theMap.end(); ++iter) {

	std::vector<double>* theVect = iter->second;
	this->cleanUpVector(theVect);

    }

    theMap.clear();
}

template <class T> 
void LpcRootOutput::cleanUpVector(std::vector<T>* theVector)
{

    if (theVector) {
	theVector->clear();
	delete theVector;
    }

}

void LpcRootOutput::storeInitialInfo()
{

    // Check that the event pointer and hit collection pointers are OK
    if (!theEvent_ || !theHits_) {return;}

    // Check if the branches have been set-up
    if (definedTrees_ == false) {this->setupTrees();}

    // Set the event number
    eventId_ = theEvent_->getEventNumber();

}

void LpcRootOutput::setupTrees()
{

    // Return if we have already defined the trees
    if (definedTrees_ == true) {return;}

    int i(0);

    if (lpcTree_) {

	// The event number
	lpcTree_->Branch("eventId", &eventId_, "eventId/I");

	// The main curve index number
	lpcTree_->Branch("curveId", &curveId_, "curveId/I");

	// The branch index number
	lpcTree_->Branch("branchId", &branchId_, "branchId/I");
	
	// The branch generation number
	lpcTree_->Branch("branchGen", &branchGen_, "branchGen/I");

	// The number of lpc points
	lpcTree_->Branch("nLpc", &nLpc_, "nLpc/I");

	// The lpc starting point.
	// Need to reserve the required memory for the lpcStartPoint vector
	lpcStartPoint_.reserve(nDim_);
	for (i = 0; i < nDim_; i++) {
	    TString lpcInitName("lpcInitX"); lpcInitName += i+1;
	    TString lpcInitNameD(lpcInitName); lpcInitNameD += "/D";
	    lpcTree_->Branch(lpcInitName.Data(), &lpcStartPoint_[i], lpcInitNameD.Data());
	}

		
	// The lpc points
	for (i = 0; i < nDim_; i++) {

	    TString lpcPointName("lpcX"); lpcPointName += i+1;
	    lpcPoints_[i] = new std::vector<double>();
	    lpcTree_->Branch(lpcPointName.Data(), &lpcPoints_[i]);

	}
	
	// The lpc principal eigenvector components
	for (i = 0; i < nDim_; i++) {

	    TString lpcEigenName("eigenX"); lpcEigenName += i+1;
	    lpcEigenVectors_[i] = new std::vector<double>();
	    lpcTree_->Branch(lpcEigenName.Data(), lpcEigenVectors_[i]);

	}

	// The lpc cosine angles
	lpcCosAngles_ = new std::vector<double>();
	lpcTree_->Branch("cosAngles", lpcCosAngles_);

	// The lpc eigenratios
	lpcRho_ = new std::vector<double>();
	lpcTree_->Branch("rho", lpcRho_);

	// The lpc c0 values
	lpcC0_ = new std::vector<double>();
	lpcTree_->Branch("c0", lpcC0_);

	// The lpc lambda values
	lpcLambda_ = new std::vector<double>();
	lpcTree_->Branch("lambda", lpcLambda_);

	// The lpc lambda axes
	for (i = 0; i < nDim_; i++) {

	    TString lpcLAxesName("lambdaX"); lpcLAxesName += i+1;
	    lpcLambdaAxes_[i] = new std::vector<double>();
	    lpcTree_->Branch(lpcLAxesName.Data(), lpcLambdaAxes_[i]);
	}

	// Lpc pathlength
	lpcPathLength_ = new std::vector<double>();
	lpcTree_->Branch("pathLength", lpcPathLength_);

	// The average lpc residuals
	lpcResiduals_ = new std::vector<double>();
	lpcTree_->Branch("lpcResiduals", lpcResiduals_);

	// The weighted average lpc residuals
	wLpcResiduals_ = new std::vector<double>();
	lpcTree_->Branch("wLpcResiduals", wLpcResiduals_);

	// The number of high rho points
	lpcTree_->Branch("nHighRho", &nHighRho_, "nHighRho/I");

	// The high rho points
	for (i = 0; i < nDim_; i++) {

	    TString lpcHighRhoName("highRhoX"); lpcHighRhoName += i+1;
	    highRhoPoints_[i] = new std::vector<double>();		
	    lpcTree_->Branch(lpcHighRhoName.Data(), highRhoPoints_[i]);

	}

	// The number of cosine angle peaks
	lpcTree_->Branch("nCosPeaks", &nCosPeaks_, "nCosPeaks/I");

	// The lpc points corresponding to the cosine angle peaks
	lpcCosPeaks_ = new std::vector<int>();
	lpcTree_->Branch("cosPeakId", lpcCosPeaks_);

	// The number of hits
	lpcTree_->Branch("nHits", &nHits_, "nHits/I");

	// The hit co-ordinates
	for (i = 0; i < nDim_; i++) {

	    TString hitName("hitX"); hitName += i+1;
	    hitCoords_[i] = new std::vector<double>();
	    lpcTree_->Branch(hitName.Data(), hitCoords_[i]);

	}

	// The nearest lpc residual for each hit
	hitResiduals_ = new std::vector<double>();
	lpcTree_->Branch("hitResiduals", hitResiduals_);

	// The weighted nearest lpc residual for each hit
	wHitResiduals_ = new std::vector<double>();
	lpcTree_->Branch("wHitResiduals", wHitResiduals_);

	// The nearest lpc point for each hit
	hitNearLpc_ = new std::vector<int>();
	lpcTree_->Branch("hitNearestLpc", hitNearLpc_);

    }

    // Vertexing information
    if (vtxTree_) {

	// The event number
	vtxTree_->Branch("eventId", &eventId_, "eventId/I");

	// The number of vertices
	vtxTree_->Branch("nVtx", &nVtx_, "nVtx/I");

	// The vertex index number
	vtxTree_->Branch("vtxIndex", &vtxIndex_, "vtxIndex/I");

	// The index of the main curve
	vtxTree_->Branch("vtxCurveId", &vtxCurveId_, "vtxCurveId/I");

	// The index of the branch
	vtxTree_->Branch("vtxBranchId", &vtxBranchId_, "vtxBranchId/I");

	// The co-ordinates of the vertex.
	// Reserve the memory for vtxPoint_
	vtxPoint_.reserve(nDim_);
	for (i = 0; i < nDim_; i++) {
	    
	    TString vtxName("vtxX"); vtxName += i+1;
	    TString vtxNameD(vtxName); vtxNameD += "/D";
	    vtxTree_->Branch(vtxName.Data(), &vtxPoint_[i], vtxNameD.Data());

	}

    }

    // Cluster information
    if (clTree_) {

	// The event number
	clTree_->Branch("eventId", &eventId_, "eventId/I");

	// The number of clusters
	clTree_->Branch("nCl", &nCl_, "nCl/I");

	// The index number of the cluster
	clTree_->Branch("clIndex", &clIndex_, "clIndex/I");

	// The index of the main curve for the cluster
	clTree_->Branch("clCurveId", &clCurveId_, "clCurveId/I");

	// The index of the branch for the cluster
	clTree_->Branch("clBranchId", &clBranchId_, "clBranchId/I");

	// The minimum lpc point number for the cluster
	clTree_->Branch("clMinLpc", &clMinLpc_, "clMinLpc/I");

	// The maximum lpc point number for the cluster
	clTree_->Branch("clMaxLpc", &clMaxLpc_, "clMaxLpc/I");

	centroid_.reserve(nDim_);
	for (i = 0; i < nDim_; i++) {

	    // The cluster centroid co-ordinates
	    TString centroidName("clX"); centroidName += i+1;
	    TString centroidNameD(centroidName); centroidNameD += "/D";
	    clTree_->Branch(centroidName.Data(), &centroid_[i], centroidNameD.Data());

	}

	// The number of axes = number of co-ordinate dimensions
	clTree_->Branch("nAxes", &nDim_, "nAxes/I");
	   
	// The principal axes (nDim*nDim)
	for (i = 0; i < nDim_; i++) {

	    TString axesName("clAxisX"); axesName += i+1;
	    TString axesNameD(axesName); axesNameD += "[nAxes]/D";
	    axes_[i] = new std::vector<double>();
	    clTree_->Branch(axesName.Data(), axes_[i]);

	}

	// The convex hull lengths
	convexHull_.reserve(nDim_);
	for (i = 0; i < nDim_; i++) {

	    TString cHullName("clHullX"); cHullName += i+1;
	    TString cHullNameD(cHullName); cHullNameD += "/D";
	    clTree_->Branch(cHullName.Data(), &convexHull_[i], cHullNameD.Data());

	}

	// Is the cluster a shower
	clTree_->Branch("clShower", &isAShower_, "clShower/I");

	// The number of hits
	clTree_->Branch("nClHits", &nClHits_, "nClHits/I");

	// The hit indices
	clHitIds_ = new std::vector<int>();
	clTree_->Branch("clHitId", &clHitIds_);
	
	// The co-ordinates of the hits in the cluster
	for (i = 0; i < nDim_; i++) {

	    TString clHitName("clHitX"); clHitName += i+1;
	    clHitCoords_[i] = new std::vector<double>();
	    clTree_->Branch(clHitName.Data(), clHitCoords_[i]);
	}

	// The hit weights
	clHitWeights_ = new std::vector<double>();
	clTree_->Branch("clHitW", &clHitWeights_);

    }

    definedTrees_ = true;

}

void LpcRootOutput::storeCurve(const LpcCurve* mainCurve)
{

    // Store the main curve and all of its branches    
    if (!mainCurve) {return;}

    this->storeCurveDetails(mainCurve, false, 0);

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

		    //int branchId = theBranch->getIndex();
		    this->storeCurveDetails(theBranch, true, iG);

		} // Check if branch pointer exists

	    } // Iteration over branches

	} // Branch generation loop
	
    } // Check for any stored branches

}

void LpcRootOutput::storeCurveDetails(const LpcAbsCurve* theCurve, bool isABranch,
				      int branchGeneration)
{

    // Store either the main curve or one of its branches, following 
    // the methodology in LpcAbsCurve::print()

    if (!theCurve || !lpcTree_) {return;}

    int curveIndex = theCurve->getIndex();
    int curveFlag = theCurve->getFlag();

    if (isABranch == false) {
	// Main curve
	curveId_ = curveIndex; 
	branchId_ = 0; 
	branchGen_ = 0;
    } else {
	// Branch
	curveId_ = curveFlag; 
	branchId_ = curveIndex; 
	branchGen_ = branchGeneration;
    }

    // The number of lpc points
    nLpc_ = theCurve->getNLpcPoints();
        
    // Starting point
    lpcStartPoint_.clear(); lpcStartPoint_.reserve(nDim_);

    LpcPoint startPoint = theCurve->getStartPoint();
    Eigen::VectorXd x0 = startPoint.getCoords();
    for (int i = 0; i < nDim_; i++) {
	lpcStartPoint_.push_back(x0(i));
    }
    
    // LpcPoints    
    std::vector<LpcPoint> lpcPoints = theCurve->getLpcPoints();     
    
    // Eigenvectors
    Eigen::MatrixXd eigenVectors = theCurve->getEigenVectors();

    for (int i = 0; i < nDim_; i++) {

	std::vector<double>* lpcPointVect = lpcPoints_[i];
	lpcPointVect->clear(); lpcPointVect->reserve(nLpc_);
	    
	std::vector<double>* lpcEigenVect = lpcEigenVectors_[i];
	lpcEigenVect->clear(); lpcEigenVect->reserve(nLpc_);
	    
    }
	
    // Cosine angles
    Eigen::VectorXd cosAngles = theCurve->getCosAngles();
    lpcCosAngles_->clear(); lpcCosAngles_->reserve(nLpc_);

    // Rho
    Eigen::VectorXd rho = theCurve->getRho();
    lpcRho_->clear(); lpcRho_->reserve(nLpc_);

    // c0
    Eigen::VectorXd c0 = theCurve->getc0();
    lpcC0_->clear(); lpcC0_->reserve(nLpc_);
    
    // Pathlengths object
    LpcPathLength pathLength = theCurve->getPathLength();    

    // Lambda
    Eigen::VectorXd lambda = pathLength.getLambda();
    lpcLambda_->clear(); lpcLambda_->reserve(nLpc_);

    // Lambda axes
    Eigen::MatrixXd lambdaAxes = pathLength.getLambdaAxes();
    for (int i = 0; i < nDim_; i++) {

	std::vector<double>* lpcLAxesVect = lpcLambdaAxes_[i];
	lpcLAxesVect->clear(); lpcLAxesVect->reserve(nLpc_);

    }

    // Delta lambda
    Eigen::VectorXd deltaLambda = pathLength.getDeltaLambda();
    lpcPathLength_->clear(); lpcPathLength_->reserve(nLpc_);

    // Residuals object
    LpcResiduals residuals = theCurve->getResiduals();

    // Lpc point residuals averaged over all nearest hits
    Eigen::VectorXd lpcRes = residuals.getLpcResiduals();
    lpcResiduals_->clear(); lpcResiduals_->reserve(nLpc_);

    // Weighted lpc point residuals averaged over all nearest hits
    Eigen::VectorXd wLpcRes = residuals.getWeightedLpcResiduals();
    wLpcResiduals_->clear(); wLpcResiduals_->reserve(nLpc_);

    // Fill the internal vectors with the above information

    for (int i = 0; i < nLpc_; i++) {

	lpcCosAngles_->push_back(cosAngles(i));
	lpcRho_->push_back(rho(i));
	lpcC0_->push_back(c0(i));
	lpcLambda_->push_back(lambda(i));
	lpcPathLength_->push_back(deltaLambda(i));
	lpcResiduals_->push_back(lpcRes(i));
	wLpcResiduals_->push_back(wLpcRes(i));

	Eigen::VectorXd lpcCoord = lpcPoints[i].getCoords();
	Eigen::VectorXd mainEVector = eigenVectors.row(i);

	for (int j = 0; j < nDim_; j++) {

	    lpcPoints_[j]->push_back(lpcCoord(j));
	    lpcEigenVectors_[j]->push_back(mainEVector(j));
	    lpcLambdaAxes_[j]->push_back(lambdaAxes(j));

	}

    }

    // HighRhoPoints
    std::vector<LpcPoint> highRhoPoints = theCurve->getHighRhoPoints();
    nHighRho_ = highRhoPoints.size();

    for (int i = 0; i < nDim_; i++) {

	std::vector<double>* highRhoVect = highRhoPoints_[i];
	highRhoVect->clear(); highRhoVect->reserve(nHighRho_);

    }

    for (int i = 0; i < nHighRho_; i++) {
	
	Eigen::VectorXd highRhoP = highRhoPoints[i].getCoords();

	for (int j = 0; j < nDim_; j++) {
	    highRhoPoints_[j]->push_back(highRhoP(j));
	}

    }

    // The number of cosine angle peaks
    std::vector<int> cosPeakIndices = theCurve->getCosPeakIndices();
    nCosPeaks_ = cosPeakIndices.size();
    lpcCosPeaks_->clear(); lpcCosPeaks_->reserve(nCosPeaks_);

    for (int i = 0; i < nCosPeaks_; i++) {
	lpcCosPeaks_->push_back(cosPeakIndices[i]);
    }    

    // Hit information
    if (theHits_) {

	// The number of hits
	nHits_ = theHits_->getNumberOfHits();

	// The hit co-ordinates
	Eigen::MatrixXd hitCoords = theHits_->getCoords();

	for (int i = 0; i < nDim_; i++) {

	    std::vector<double>* hitVect = hitCoords_[i];
	    hitVect->clear(); hitVect->reserve(nHits_);

	}

	// Hit-to-nearest-lpc point residuals
	Eigen::VectorXd hitRes = residuals.getHitResiduals();
	hitResiduals_->clear(); hitResiduals_->reserve(nHits_);

	// Weighted hit-to-nearest-lpc point residuals
	Eigen::VectorXd wHitRes = residuals.getWeightedHitResiduals();
	wHitResiduals_->clear(); wHitResiduals_->reserve(nHits_);

	// The index number of the nearest lpc point for each hit
	Eigen::VectorXi hitNearLpc = residuals.getHitNearestLpc();
	hitNearLpc_->clear(); hitNearLpc_->reserve(nHits_);

	for (int i = 0; i < nHits_; i++) {

	    hitResiduals_->push_back(hitRes(i));
	    wHitResiduals_->push_back(wHitRes(i));
	    hitNearLpc_->push_back(hitNearLpc(i));

	    Eigen::VectorXd hitPos = hitCoords.row(i);

	    for (int j = 0; j < nDim_; j++) {

		hitCoords_[j]->push_back(hitPos(j));

	    }

	}

    }

    // Store the data in the tree
    lpcTree_->Fill();

}

void LpcRootOutput::storeVertices(const std::vector<LpcVertex*>& theVertices)
{

    if (!vtxTree_) {return;}

    nVtx_ = theVertices.size();

    for (int i = 0; i < nVtx_; i++) {

	LpcVertex* vertex = theVertices[i];
	if (vertex) {
	    	    
	    vtxIndex_ = vertex->getIndex();
	    vtxCurveId_ = vertex->getCurveId();
	    vtxBranchId_ = vertex->getBranchId();

	    Eigen::VectorXd vtxCoord = vertex->getCoords();

	    vtxPoint_.clear(); vtxPoint_.reserve(nDim_);
	    for (int j = 0; j < nDim_; j++) {

		vtxPoint_[j] = vtxCoord(j);
		
	    }

	    // Store the info for this vertex
	    vtxTree_->Fill();	
	    
	} // Check if the vertex pointer exists
	
    } // Loop over vertices

}

void LpcRootOutput::storeClusters(const std::vector<LpcCluster*>& theClusters)
{

    if (!clTree_) {return;}
	
    nCl_ = theClusters.size();

    // Loop over clusters
    for (int iC = 0; iC < nCl_; iC++) {
	
	LpcCluster* cluster = theClusters[iC];
	if (cluster) {
	    
	    clIndex_ = cluster->getIndex();
	    clCurveId_ = cluster->getCurveId();
	    clBranchId_ = cluster->getBranchId();
	    
	    LpcBinRange lpcRange = cluster->getLpcRange();
	    
	    clMinLpc_ = lpcRange.getMinBin();
	    clMaxLpc_ = lpcRange.getMaxBin();
	    
	    Eigen::VectorXd centroid = cluster->getCentroid();
	    Eigen::MatrixXd PCAxes = cluster->getPCAxes();
	    Eigen::VectorXd convexHull = cluster->getConvexHull();

	    centroid_.clear(); centroid_.reserve(nDim_);
	    convexHull_.clear(); convexHull_.reserve(nDim_);
	    for (int i = 0; i < nDim_; i++) {
		
		centroid_.push_back(centroid(i));
		convexHull_.push_back(convexHull(i));
		
		std::vector<double>* theAxis = axes_[i];
		theAxis->clear(); theAxis->reserve(nDim_);

		Eigen::VectorXd axisProjection = PCAxes.col(i);

		for (int j = 0; j < nDim_; j++) {
		    theAxis->push_back(axisProjection(j));
		}
		
	    }
	    
	    isAShower_ = 0;
	    if (cluster->isAShower() == true) {isAShower_ = 1;}
	    
	    // Hit indices
	    std::vector<int> hitIndices = cluster->getHitIndices();
	    nClHits_ = hitIndices.size();
	    clHitIds_->clear(); clHitIds_->reserve(nClHits_);

	    // Hit integer indices and positions
	    Eigen::MatrixXd clHitPosMatrix = cluster->getHitPositions();	    
	    for (int i = 0; i < nDim_; i++) {
		
		std::vector<double>* clHitXVect = clHitCoords_[i];
		clHitXVect->clear(); clHitXVect->reserve(nClHits_);
		
	    }
	    
	    // Vector of hit weights
	    std::vector<double> hitWeights = cluster->getHitWeights();
	    clHitWeights_->clear(); clHitWeights_->reserve(nClHits_);

	    // Loop over the cluster hits
	    for (int i = 0; i < nClHits_; i++) {
		
		// Store the hit index integer
		clHitIds_->push_back(hitIndices[i]);
		
		Eigen::VectorXd clusterCoord = clHitPosMatrix.row(i);
		// Store the hit co-ordinates
		for (int j = 0; j < nDim_; j++) {
		    clHitCoords_[j]->push_back(clusterCoord(j));
		}
		
		// Store the hit weight
		clHitWeights_->push_back(hitWeights[i]);

	    }
	    
	    // Store the info for this cluster
	    clTree_->Fill();	    
	    
	} // Check if cluster pointer exists
	
    } // Cluster loop

}

void LpcRootOutput::storeExtraInfo()
{
}

#endif

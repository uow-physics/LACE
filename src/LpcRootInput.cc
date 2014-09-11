#ifdef LPC_USE_ROOT

// Copyright University of Warwick 2014
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// Authors:
// John Back

/*! \file LpcRootInput.cc
    \brief Class to read in a ROOT file containing a dataset of hits (point cloud)
*/

#include "LACE/LpcRootInput.hh"

#include "LACE/LpcEvent.hh"
#include "LACE/LpcHit.hh"
#include "LACE/LpcHitCollection.hh"

#include "TBranch.h"
#include "TFile.h"
#include "TObjArray.h"
#include "TString.h"
#include "TTree.h"

#include <Eigen/Dense>
#include <iostream>

LpcRootInput::LpcRootInput(const std::string& inputFileName,
			   const std::string& inputTreeName) :
    LpcAbsInput(inputFileName),
    inFile_(0),
    treeName_(inputTreeName),
    inTree_(0),
    gotBranches_(false),
    minEvent_(0),
    maxEvent_(0),
    eventId_(-1),
    nHits_(0),
    hitCoords_(),
    weights_(0)
{
    this->initialise();
}

void LpcRootInput::initialise()
{
    inFile_ = TFile::Open(inputFileName_.c_str(), "read");

    if (!inFile_ || inFile_->IsZombie() || !inFile_->IsOpen()) {
	std::cerr << "Error in LpcRootInput::initialise. Cannot open " 
		  << inputFileName_ << std::endl;
	return;
    }

    inTree_ = dynamic_cast<TTree*>(inFile_->Get(treeName_.c_str()));

    if (!inTree_) {
	std::cerr << "Error in LpcRootInput::initilise. Cannot find the tree "
		  << treeName_ << " in the file "<< inputFileName_ << std::endl;
	return;
    }

    // Check that we have the necessary branches
    this->checkBranches();

}


void LpcRootInput::checkBranches()
{

    nDim_ = 0; nEvents_ = 0; minEvent_ = 0; maxEvent_ = 0; gotBranches_ = false;
    if (!inTree_) {return;}

    // Make sure the tree contains the correct branch structure
    TObjArray* theBranches = inTree_->GetListOfBranches();
    if (!theBranches) {return;}

    const int nBranches = theBranches->GetEntries();

    this->cleanUpMap(hitCoords_);
    this->cleanUpVector(weights_);

    bool gotEvent(false), gotNHits(false), gotCoords(false), gotWeights(false);
    // First, disable all branches, then enable the ones we actually need
    inTree_->SetBranchStatus("*", 0);

    // Keep track of the hitX and weight branch names
    hitXNames_.clear();
    weightName_ = "";

    for (int iB = 0; iB < nBranches; iB++) {

	TBranch* theBranch = dynamic_cast<TBranch*>((*theBranches)[iB]);

	if (theBranch) {

	    // Store the lowercase name of the branch with its index
	    TString branchName = theBranch->GetName();

	    // Keep track of how many hitX branches we have, which
	    // will correspond to the number of spatial dimensions
	    if (branchName.Contains("eventId")) {

		inTree_->SetBranchStatus("eventId", 1);
		inTree_->SetBranchAddress("eventId", &eventId_);
		gotEvent = true;

	    } else if (branchName.Contains("nHits")) {

		inTree_->SetBranchStatus("nHits", 1);
		inTree_->SetBranchAddress("nHits", &nHits_);
		gotNHits = true;
		
	    } else if (branchName.Contains("hitX")) {

		// Create the pointer to the vector of the co-ordinate values
		hitCoords_[nDim_] = new std::vector<double>();
		// Set the branch address
		inTree_->SetBranchStatus(branchName, 1);
		inTree_->SetBranchAddress(branchName, &hitCoords_[nDim_]);

		// Keep track of the branch name for the hit co-ordinate
		hitXNames_.push_back(branchName.Data());

		// Increment the number of spatial dimensions
		nDim_++;

	    } else if (branchName.Contains("hitW")) {

		// Create the pointer to the vector of weights
		weights_ = new std::vector<double>();
		// Set the branch address
		inTree_->SetBranchStatus(branchName, 1);
		inTree_->SetBranchAddress(branchName, &weights_);

		weightName_ = branchName.Data();

		gotWeights = true;

	    }

	}

    }

    if (nDim_ > 0) {gotCoords = true;}

    gotBranches_ = gotEvent && gotNHits && gotCoords && gotWeights;
    std::cout << "LpcRootInput: Got branches = " << std::boolalpha << gotBranches_ << std::endl;

    // If we have all of the required branches, set the number of
    // events as the number of entries in the input tree
    if (gotBranches_) {
	nEvents_ = inTree_->GetEntries();
	minEvent_ = inTree_->GetMinimum("eventId");
	maxEvent_ = inTree_->GetMaximum("eventId");
    }

}

LpcRootInput::~LpcRootInput()
{
    // Do nothing here, since the finalise function should 
    // clean-up everything
}

void LpcRootInput::finalise()
{
    this->cleanUpMap(hitCoords_);
    this->cleanUpVector(weights_);
    if (inFile_) {inFile_->Close();}
}

void LpcRootInput::cleanUpMap(LpcRootMap& theMap)
{
    
    LpcRootMap::const_iterator iter;
    for (iter = theMap.begin(); iter != theMap.end(); ++iter) {

	std::vector<double>* theVect = iter->second;
	this->cleanUpVector(theVect);

    }

    theMap.clear();
}

template <class T> 
void LpcRootInput::cleanUpVector(std::vector<T>* theVector)
{

    if (theVector) {
	theVector->clear();
	delete theVector;
    }

}


LpcEvent* LpcRootInput::getEvent(int eventNo)
{

    // Various checks
    if (!inTree_ || nEvents_ < 1 || eventNo < minEvent_ 
	|| eventNo > maxEvent_) {return 0;}

    // Disable all branches except for eventId and nHits
    inTree_->SetBranchStatus("*", 0);
    inTree_->SetBranchStatus("eventId", 1);
    inTree_->SetBranchStatus("nHits", 1);

    nHits_ = 0;
    eventId_ = 0;

    // Try the tree entry given by the difference of the 
    // required event number with the minimum event number
    int entryIndex = eventNo - minEvent_;

    if (entryIndex >= 0 && entryIndex < nEvents_) {
	inTree_->GetEntry(entryIndex);
    }

    // Check that the event id number is the one we want
    bool gotEvent(false);
    if (eventId_ == eventNo) {gotEvent = true;}

    if (gotEvent == false) {

	// We can't seem to find the event. 
	// Try looping through all entries...
	for (int i = 0; i < nEvents_; i++) {

	    inTree_->GetEntry(i);

	    if (eventId_ == eventNo) {
		// We have the event. Stop the loop
		entryIndex = i;
		gotEvent = true;
		break;
	    }

	}

    }

    LpcEvent* theEvent(0);

    // Create the event pointer
    if (gotEvent == true) {

	// We have the required event
        // Enable the rest of the branches and get the info we need
        this->getTreeData(entryIndex);

	// Create a new hit collection that will be stored 
	// and owned by the LpcEvent
	LpcHitCollection* theHits = new LpcHitCollection(nDim_);

	// Store the hit co-ordinate information in an Eigen MatrixXd object.
	// This enables us to only loop over the branch vectors once.
	// Row = hit, column = x,y,z,... co-ordinate components
	Eigen::MatrixXd coordMatrix = Eigen::MatrixXd::Zero(nHits_, nDim_);

	int iD(0), iH(0);
	// Loop over the number of co-ordinate dimensions ("hitX" branches)
	for (iD = 0; iD < nDim_; iD++) {

	    // Get the vector of the given "hitX" co-ordinate axis
	    std::vector<double>* coordVect = hitCoords_[iD];

	    // Check if this pointer is null. If so, there is a (major) problem
	    // with the ROOT file input, i.e. just stop processing it any further
	    if (!coordVect) {break;}

	    // Check that the size of this vector is equal to the number of hits
	    int nVect = coordVect->size();
	    if (nVect != nHits_) {
		std::cerr << "Error in LpcRootInput::getEvent! "
			  << "The branch for hit co-ordinate dimension "
			  << iD << " has " << nVect
			  << " entries which is not equal to the number of hits "
			  << nHits_ << ". Please check the input file!" << std::endl;
		// Stop further processing the ROOT input file
		break;
	    }
	    
	    // Store the co-ordinate components in the matrix for all hits
	    for (iH = 0; iH < nHits_; iH++) {
		coordMatrix(iH, iD) = (*coordVect)[iH];
	    }

	}

	// We need to check if the pointer to the vector of weights is OK
	bool useWeights(false);
	if (weights_) {
	    int nWeights = weights_->size();
	    if (nWeights == nHits_) {useWeights = true;}
	}

	// Create the hit pointers and add them to the collection
	for (iH = 0; iH < nHits_; iH++) {

	    Eigen::VectorXd hitPos = coordMatrix.row(iH);
	    double w(1.0);
	    if (useWeights) {w = (*weights_)[iH];}

	    LpcHit* aHit = new LpcHit(iH, hitPos, w);
	    theHits->addHit(aHit);
	    
	}

	
	theEvent = new LpcEvent(eventNo, theHits);

    }

    // Return the event pointer
    return theEvent;

}

void LpcRootInput::getTreeData(int entryIndex)
{

    if (!inTree_) {return;}

    for (int iD = 0; iD < nDim_; iD++) {
  
        inTree_->SetBranchStatus(hitXNames_[iD].c_str(), 1);

        // Set-up the vector of the given "hitX" co-ordinate axis
        std::vector<double>* coordVect = hitCoords_[iD];
	coordVect->clear(); coordVect->reserve(nHits_);

    }

    // Set-up the weight vector
    inTree_->SetBranchStatus(weightName_.c_str(), 1);
    weights_->clear(); weights_->reserve(nHits_);

    // Fill the above internal vectors
    inTree_->GetEntry(entryIndex);

}

#endif

// Copyright University of Warwick 2014
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// Authors:
// John Back

/*! \file LpcTextInput.cc
    \brief Class to read in an ascii file containing a dataset of hits (point cloud)
*/

#include "LACE/LpcTextInput.hh"

#include "LACE/LpcEvent.hh"
#include "LACE/LpcHit.hh"
#include "LACE/LpcHitCollection.hh"

#include <Eigen/Dense>

#include <cstdlib>

LpcTextInput::LpcTextInput(const std::string& inputFileName) :
    LpcAbsInput(inputFileName),
    getData_(inputFileName.c_str()),
    whiteSpace_(" "),
    maxNChar_(250),
    hitIndex_(-1)
{
    this->initialise();
}

LpcTextInput::~LpcTextInput()
{
}

void LpcTextInput::initialise()
{
    // Read in the input file and store the number of events

    // Start the input stream at the beginning
    getData_.clear(); // reset any flags
    getData_.seekg(0, std::ios::beg);

    std::pair<int, int> indices = this->getNextEventIndices();
    nEvents_ = indices.first;
    nDim_ = indices.second;

}

LpcEvent* LpcTextInput::getEvent(int eventNo)
{   

    // Check if the event number is valid
    if (eventNo < 0 || eventNo >= nEvents_ || nDim_ < 1) {return 0;}

    // Usually, the events in the input file will be processed
    // in sequential order by the LpcRun class

    // Retrieve the event number and number of hits, which should be
    // on the present line if the input file is formatted correctly
    std::pair<int, int> indices = this->getNextEventIndices();
    int evtIndex = indices.first;
    int nHits = indices.second;
    int iEvt(0), iLine(0);

    if (eventNo < evtIndex) {

	// We have gone beyond the event. Reset the input streamer
	// to the beginning of the file and try again.
	this->initialise();

	// Get the first event index and number of events
	indices = this->getNextEventIndices();
	evtIndex = indices.first;
	nHits = indices.second;

    }

    // Process the lines in the input file until we reach the correct event
    LpcEvent* theEvent(0);

    for (iEvt = evtIndex; iEvt < nEvents_; iEvt++) {

	if (eventNo == iEvt) {

	    // We have the required event

	    // Create a new hit collection that will be stored 
	    // and owned by the LpcEvent
	    LpcHitCollection* theHits = new LpcHitCollection(nDim_);

	    // Loop over the hits and add them to the collection
	    // First, initialise the hitIndex_ to -1. This is always incremented
	    // before we create a new hit, and is assigned as the hit index value
	    hitIndex_ = -1;
	    for (iLine = 0; iLine < nHits; iLine++) {

		LpcHit* aHit = this->getNextHit();
		// Add the hit if we have processed the data line OK.
		// This will also increment hitIndex_ by one
		if (aHit) {theHits->addHit(aHit);}

	    }

	    // Create the new LpcEvent
	    //std::cout<<"Creating new LpcEvent, evtNo = "<<eventNo<<", evtIndex = "
	    //     <<evtIndex<<", nHits = "<<nHits<<std::endl;

	    theEvent = new LpcEvent(eventNo, theHits);
	    break;

	} else if (eventNo > iEvt) {

	    //std::cout<<"Skipping event "<<iEvt<<std::endl;

	    // The event number we need has not yet been reached.
	    // Skip the lines for the current indexed event.
	    for (iLine = 0; iLine < nHits; iLine++) {

		// Skip the remaining lines until we reach the next event
		getData_.ignore(maxNChar_, '\n');

		// Stop the loop if we have failed to read the line correctly
		if (getData_.fail()) {break;}

	    }

	    // Get the next event index and number of events
	    indices = this->getNextEventIndices();
	    evtIndex = indices.first;
	    nHits = indices.second;

	} // The current event number is lower that what we need

    } // Loop over the number of possible events in the file

    return theEvent;

}

std::pair<int, int> LpcTextInput::getNextEventIndices()
{
    std::string lineString("");
    std::getline(getData_, lineString);

    int evtIndex(0), nHits(0);

    if (!getData_.fail()) {

	std::vector<std::string> lineVect = this->split(lineString, whiteSpace_);
    
	if (lineVect.size() == 2) {
	    evtIndex = atoi(lineVect[0].c_str());
	    nHits = atoi(lineVect[1].c_str());
	}
    }

    return std::pair<int,int>(evtIndex, nHits);

}

LpcHit* LpcTextInput::getNextHit()
{
    // Read the next line of hit data
    std::string lineString("");
    std::getline(getData_, lineString);

    // Stop if we have failed to read the line
    if (getData_.fail()) {return 0;}

    // Split up the line according to white spaces
    std::vector<std::string> lineVect = this->split(lineString, whiteSpace_);

    // Create a vector storing the hit information (co-ordinates and weight)
    Eigen::VectorXd coords = Eigen::VectorXd::Zero(nDim_);

    // Convert the separated sub-strings into doubles
    int nV = lineVect.size() - 1;
    for (int i = 0; i < nV; i++) {

	double value = std::atof(lineVect[i].c_str());
	coords(i) = value;

    }

    double weight = std::atof(lineVect[nV].c_str());

    // Increment the hit index
    hitIndex_ += 1;
    // Create the LpcHit pointer, then return it
    LpcHit* theHit = new LpcHit(hitIndex_, coords, weight);

    return theHit;

}

std::vector<std::string> LpcTextInput::split(const std::string& theString,
					     const std::string& splitter) const
{
    // Code from STLplus
    std::vector<std::string> result;

    if (!theString.empty() && !splitter.empty()) {
	
	for (std::string::size_type offset = 0;;) {
	    
	    std::string::size_type found = theString.find(splitter, offset);
	    
	    if (found != std::string::npos) {
		std::string tmpString = theString.substr(offset, found-offset);
		if (tmpString.size() > 0) {result.push_back(tmpString);}
		offset = found + splitter.size();
	    } else {
		std::string tmpString = theString.substr(offset, theString.size()-offset);
		if (tmpString.size() > 0) {result.push_back(tmpString);}
		break;
	    }
	}
    }
    
    return result;
}

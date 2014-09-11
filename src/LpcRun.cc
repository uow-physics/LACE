// Copyright University of Warwick 2014
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// Authors:
// John Back

/*! \file LpcRun.cc
  \brief Class to run the overall lpc algorithm with a given input dataset
*/

#include "LACE/LpcRun.hh"

#include "LACE/LpcAbsInput.hh"
#include "LACE/LpcAbsOutput.hh"
#include "LACE/LpcEvent.hh"
#include "LACE/LpcInputFactory.hh"
#include "LACE/LpcOutputFactory.hh"
#include "LACE/LpcParameters.hh"
#include "LACE/LpcProcess.hh"

#include <iostream>

LpcRun::LpcRun(const std::string& parameterFileName) :
    theParameters_(new LpcParameters(parameterFileName)),
    theProcess_(new LpcProcess(theParameters_)),
    theInput_(0),
    theOutput_(0)
{    
    this->initialise();
}

LpcRun::~LpcRun()
{
    delete theOutput_;
    delete theInput_;
    delete theProcess_;
    delete theParameters_;
}

void LpcRun::initialise() {

    // Use the input factory to select the input reader, which returns
    // a new LpcAbsInput pointer
    LpcInputFactory inputFactory;
    theInput_ = inputFactory.getInputReader(theParameters_);

    // Use the output factory to select the output writer, which returns
    // a new LpcAbsOutput pointer
    LpcOutputFactory outputFactory;
    theOutput_ = outputFactory.getOutputWriter(theParameters_);

}

void LpcRun::run() {

    // Run the lpc algorithm over the required set of data.

    // Have we defined valid input & output pointers in initialise?
    if (!theInput_) {
	std::cerr<<"Input pointer is null. Please check that the input file "
		 <<"has a valid filename and format option..."<<std::endl;
	return;
    } 

    if (!theOutput_) {
	std::cerr<<"Output pointer is null. Please check that the output file "
		 <<"has a valid filename and format option..."<<std::endl;
	return;
    }

    // From the input, get the number of events & spatial dimensions
    int nEvents = theInput_->getNEvents();
    int nDim = theInput_->getNDimensions();

    std::cout<<"Number of events in "<<theInput_->getInputFileName()
	     <<" = "<<nEvents<<std::endl;
    std::cout<<"Number of spatial dimensions for the hits = "<<nDim<<std::endl;

    // Loop over the events
    // First, check to see if we have set (valid) indices of the
    // first and last event to process
    int firstIndex(0), lastIndex(nEvents);
    if (theParameters_) {

	firstIndex = theParameters_->getFirstEvent();
	lastIndex = theParameters_->getLastEvent() + 1;

    }

    // Check for negative values
    if (firstIndex < 0) {firstIndex = 0;} // First event
    if (lastIndex <= 0) {lastIndex = nEvents;} // All events

    // Loop over the events
    for (int iEvt = firstIndex; iEvt < lastIndex; iEvt++) {

	// Get the event from the input
	std::cout<<"Event "<<iEvt<<std::endl;
	LpcEvent* theEvent = theInput_->getEvent(iEvt);

	if (theEvent != 0) {

	    // Reconstruct this event to find curves, vertices & clusters
	    theProcess_->reconstruct(theEvent);

	    // Store the output within the event
	    theOutput_->store(theEvent);

	} else {

	    std::cerr<<"LpcRun: The event "<<iEvt<<" is null"<<std::endl;
	}

	// Delete the event after we have stored/accessed its info
	delete theEvent;

    } // Event loop

    // Finalise the input and output steps
    theOutput_->finalise();
    theInput_->finalise();

    std::cout<<"Results stored in "<<theOutput_->getOutputFileName()<<std::endl;

}


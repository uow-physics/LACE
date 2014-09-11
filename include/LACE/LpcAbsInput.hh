// Copyright University of Warwick 2014
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// Authors:
// John Back

/*! \file LpcAbsInput.hh
    \brief File containing the declaration of the LpcAbsInput class
*/

/*! \class LpcAbsInput
    \brief Pure abstract class to read in a dataset of hits (point cloud)    
*/

#ifndef LPC_ABS_INPUT_HH
#define LPC_ABS_INPUT_HH

#include <string>

class LpcEvent;

class LpcAbsInput {

public:

    //! Constructor
    /*!
      \param [in] inputFileName The name of the input file
    */
    LpcAbsInput(const std::string& inputFileName) :
	inputFileName_(inputFileName),
	nEvents_(0),
	nDim_(3)
    {;}
    
    //! Destructor
    virtual ~LpcAbsInput() {;}

    //! Initialisation function. This must be implemented in the derived classes
    virtual void initialise() = 0;

    //! Finalising function. This must be implemented in the derived classes
    virtual void finalise() = 0;
    
    //! Get the required event. This must be implemented in the derived classes
    /*!
      \param [in] eventNo is the event number to retrieve
      \returns a new pointer to the LpcEvent
    */
    virtual LpcEvent* getEvent(int eventNo) = 0;

    // Accessor functions

    //! The name of the input file
    /*!
      \return the name of the input file
    */
    std::string getInputFileName() const {return inputFileName_;}

    //! The number of events in the input file
    /*!
      \returns the number of events in the input file
    */
    int getNEvents() const {return nEvents_;}

    //! The number of dimensions for the hit co-ordinates
    /*!
      \returns the number of dimensions for the hit co-ordinates
    */
    int getNDimensions() const {return nDim_;}

protected:

    //! The name of the input file
    std::string inputFileName_;

    //! The number of events
    int nEvents_;

    //! The number of dimensions for the hit co-ordinates
    int nDim_;

private:

    //! Private default constructor
    LpcAbsInput();

    //! Private copy constructor
    LpcAbsInput(const LpcAbsInput& other);

    //! Private assignment operator
    LpcAbsInput& operator=(const LpcAbsInput& other);

};

#endif


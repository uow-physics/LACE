// Copyright University of Warwick 2014
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// Authors:
// John Back

/*! \file LpcProcess.hh
    \brief File containing the declaration of the LpcProcess class
*/

/*! \class LpcProcess
    \brief Class containing the processing steps to obtain lpc curves for events
*/

#ifndef LPC_PROCESS_HH
#define LPC_PROCESS_HH

class LpcAbsClusterAlgorithm;
class LpcAlgorithm;
class LpcBranchAlgorithm;
class LpcEvent;
class LpcFeatures;
class LpcParameters;
class LpcShowerAlgorithm;

class LpcProcess {

public:

    //! Constructor
    /*!
      \param [in] theParameters The pointer to the constant parameters
    */
    LpcProcess(const LpcParameters* theParameters);

    //! Destructor
    virtual ~LpcProcess();

    //! Reconstruct the curves, vertices and clusters and store them in the given event
    /*!
      \param [in] theEvent a pointer to the event
    */
    void reconstruct(LpcEvent* theEvent);

protected:

    //! Retrive the hit positions and weights from the event
    void retrieveHitInfo();

    //! Find the starting point for the algorithm
    void selectStartPoint();

    //! Scale the hit positions with their (non-zero) range
    void scalePoints();

    //! Pointer to the lpc algorithm
    LpcAlgorithm* theAlgorithm_;

    //! Pointer to the branching part of the lpc algorithm
    LpcBranchAlgorithm* theBranchAlgorithm_;

    //! Pointer to the features code
    LpcFeatures* theFeatures_;

    //! Pointer to the clustering (and vertexing) algorithm
    LpcAbsClusterAlgorithm* theClusterAlgorithm_;

    //! Pointer to the shower algorithm
    LpcShowerAlgorithm* theShowerAlgorithm_;

private:

    //! Private default constructor
    LpcProcess();

    //! Private copy constructor
    LpcProcess(const LpcProcess& other);

    //! Private assignment operator
    LpcProcess& operator = (const LpcProcess& other);

    //! Integer to specify the clustering and vertexing algorithm
    int clustering_;
    
    //! The number of branching levels
    int branchLevel_;

    //! Decide if scaling is to be applied to the data hits
    int doScaling_;

};

#endif

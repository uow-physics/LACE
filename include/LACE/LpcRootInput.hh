#ifdef LPC_USE_ROOT

// Copyright University of Warwick 2014
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// Authors:
// John Back

/*! \file LpcRootInput.hh
    \brief File containing the declaration of the LpcRootInput class
*/

/*! \class LpcRootInput
    \brief Class to read in a ROOT file containing a dataset of hits (point cloud)
*/

#ifndef LPC_ROOT_INPUT_HH
#define LPC_ROOT_INPUT_HH

#include "LACE/LpcAbsInput.hh"

#include <map>
#include <string>
#include <vector>

class LpcEvent;
class TBranch;
class TFile;
class TTree;

class LpcRootInput : public LpcAbsInput {

public:

    //! Constructor
    /*!
      \param [in] inputFileName The name of the input file
      \param [in] inputTreeName The name of the input tree (Default = "Data")
    */
    LpcRootInput(const std::string& inputFileName,
		 const std::string& inputTreeName = "Data");

    //! Destructor
    virtual ~LpcRootInput();

    //! Initialisation function
    virtual void initialise();

    //! Finalising function
    virtual void finalise();
    
    //! Get the required event
    /*!
      \param [in] eventNo is the event number to retrieve
      \returns a new pointer to the LpcEvent
    */
    virtual LpcEvent* getEvent(int eventNo);

protected:

private:

    //! Private default constructor
    LpcRootInput();

    //! Private copy constructor
    LpcRootInput(const LpcRootInput& other);

    //! Private assignment operator
    LpcRootInput& operator=(const LpcRootInput& other);

    //! Typdef for a map storing a pointer to a vector of doubles. 
    //! The key represents the co-ordinate component, 0 = "x", 1 = "y", etc..
    typedef std::map<int, std::vector<double>* > LpcRootMap;
    
    //! Check the branches
    void checkBranches();

    //! Clean up function for the internal maps
    /*!
      \param [in] theMap The map to clean up
    */
    void cleanUpMap(LpcRootMap& theMap);

    //! Clean up function for the internal vectors
    /*!
      \param [in] theVector The vector pointer to clean up
    */
    template <class T> void cleanUpVector(std::vector<T>* theVector);

    //! Get all of the required tree data for the specific entry index
    /*!
      \param [in] entryIndex The tree entry index number
    */
    void getTreeData(int entryIndex);
  
    //! The ROOT input file pointer
    TFile* inFile_;

    //! The ROOT input tree name
    std::string treeName_;

    //! The ROOT input tree pointer
    TTree* inTree_;

    //! Boolean to specify if we have the required branches
    bool gotBranches_;

    //! The minimum event number
    int minEvent_;

    //! The maximum event number
    int maxEvent_;

    //! The event number
    int eventId_;

    //! The number of hits
    int nHits_;

    //! The hit co-ordinates
    LpcRootMap hitCoords_;

    //! The hit weights
    std::vector<double>* weights_;

    //! The names of the hit branches
    std::vector<std::string> hitXNames_;

    //! The name of the weight branch
    std::string weightName_;

};

#endif

#endif

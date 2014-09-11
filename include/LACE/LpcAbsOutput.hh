// Copyright University of Warwick 2014
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// Authors:
// John Back

/*! \file LpcAbsOutput.hh
    \brief File containing the declaration of the LpcAbsOutput class
*/

/*! \class LpcAbsOutput
    \brief Abstract class to store the results of the lpc algorithm
*/

#ifndef LPC_ABS_OUTPUT_HH
#define LPC_ABS_OUTPUT_HH

#include <string>
#include <vector>

class LpcCluster;
class LpcCurve;
class LpcEvent;
class LpcHitCollection;
class LpcVertex;

class LpcAbsOutput {

public:

    //! Constructor
    /*!
      \param [in] outputFileName The name of the output file
    */
    LpcAbsOutput(const std::string& outputFileName);

    //! Destructor
    virtual ~LpcAbsOutput();

    //! Pure virtual initialisation function
    virtual void initialise() = 0;
    
    //! Finalising function. This must be implemented in the derived classes
    virtual void finalise() = 0;
 
    //! Store all of the results from the event. This calls the pure virtual methods below
    /*!
      \param [in] theEvent The pointer to the event
    */
    void store(LpcEvent* theEvent);

    // Accessor functions

    //! The name of the output file
    /*!
      \return the name of the output file
    */
    std::string getOutputFileName() const {return outputFileName_;}

protected:

    //! The name of the output file
    std::string outputFileName_;

    //! The pointer to the (read-only) current event
    const LpcEvent* theEvent_;

    //! The pointer to the (read-only) hit collection
    const LpcHitCollection* theHits_;

    //! The number of co-ordinate dimensions
    int nDim_;

    //! Store any initial information about the event
    virtual void storeInitialInfo() = 0;

    //! Store the main curve (& its branches). Pure virtual method that must be implemented
    /*!
      \param [in] mainCurve The main lpc curve (read-only access)
    */
    virtual void storeCurve(const LpcCurve* mainCurve) = 0;

    //! Store the vertices
    /*!
      \param [in] theVertices The vector of vertex pointers (read-only access)
    */
    virtual void storeVertices(const std::vector<LpcVertex*>& theVertices) = 0;

    //! Store the clusters
    /*!
      \param [in] theClusters The vector of cluster pointers (read-only access)
    */
    virtual void storeClusters(const std::vector<LpcCluster*>& theClusters) = 0;

    //! Store anything else. Pure virtual function
    virtual void storeExtraInfo() = 0;

private:

    //! Private default constructor
    LpcAbsOutput();

    //! Private copy constructor
    LpcAbsOutput(const LpcAbsOutput& other);

    //! Private assignment operator
    LpcAbsOutput& operator=(const LpcAbsOutput& other);

};

#endif


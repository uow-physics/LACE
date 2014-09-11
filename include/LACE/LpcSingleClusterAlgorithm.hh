// Copyright University of Warwick 2014
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// Authors:
// John Back

/*! \file LpcSingleClusterAlgorithm.hh
    \brief File containing the declaration of the LpcSingleClusterAlgorithm class
*/

/*! \class LpcSingleClusterAlgorithm
    \brief Class that creates one cluster by including all hits
*/

#ifndef LPC_SINGLE_CLUSTER_ALGORITHM_HH
#define LPC_SINGLE_CLUSTER_ALGORITHM_HH

#include "LACE/LpcAbsClusterAlgorithm.hh"

#include "LACE/LpcClusterData.hh"

class LpcCurve;
class LpcParameters;

class LpcSingleClusterAlgorithm : public LpcAbsClusterAlgorithm {

public:

    //! Constructor
    /*!
      \param [in] thePars a pointer to the constant lpc parameters
    */
    LpcSingleClusterAlgorithm(const LpcParameters* thePars);

    //! Destructor
    virtual ~LpcSingleClusterAlgorithm();

    //! Process the main lpc curve (and its branches), finding vertices and clusters
    /*!
      \param [in] theCurve The pointer to the constant main lpc curve (read-only access)
      \returns an LpcClusterData object, which stores vectors of vertex and cluster pointers
    */
    virtual LpcClusterData findClusters(const LpcCurve* theCurve);

protected:

private:

    //! Private default constructor
    LpcSingleClusterAlgorithm();
    
    //! Copy constructor
    LpcSingleClusterAlgorithm(const LpcSingleClusterAlgorithm& other);
    
    //! Assignment operator
    LpcSingleClusterAlgorithm& operator=(const LpcSingleClusterAlgorithm& other);
  
};

#endif

// Copyright University of Warwick 2014
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// Authors:
// John Back

/*! \file LpcAbsClusterAlgorithm.hh
    \brief File containing the declaration of the LpcAbsClusterAlgorithm class
*/

/*! \class LpcAbsClusterAlgorithm
    \brief Class that defines the abstract vertex & cluster algorithm interface
*/

#ifndef LPC_ABS_CLUSTER_ALGORITHM_HH
#define LPC_ABS_CLUSTER_ALGORITHM_HH

#include "LACE/LpcClusterData.hh"

class LpcCurve;
class LpcParameters;

class LpcAbsClusterAlgorithm {

public:

    //! Constructor
    /*!
      \param [in] thePars a pointer to the constant lpc parameters
    */
    LpcAbsClusterAlgorithm(const LpcParameters* thePars) :
	thePars_(thePars) {;}

    //! Destructor
    virtual ~LpcAbsClusterAlgorithm() {;}

    //! Process the main lpc curve (and its branches), finding vertices and clusters.
    //! This must be implemented in any derived class
    /*!
      \param [in] theCurve The pointer to the main lpc curve (read-only)
      \returns an LpcClusterData object, which stores vectors of vertex and cluster pointers
    */
    virtual LpcClusterData findClusters(const LpcCurve* theCurve) = 0;

protected:    

    //! Pointer to the constant lpc parameters
    const LpcParameters* thePars_;

private:

    //! Private default constructor
    LpcAbsClusterAlgorithm();

    //! Copy constructor
    LpcAbsClusterAlgorithm(const LpcAbsClusterAlgorithm& other);

    //! Assignment operator
    LpcAbsClusterAlgorithm& operator=(const LpcAbsClusterAlgorithm& other);
    

};

#endif

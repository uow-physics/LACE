// Copyright University of Warwick 2014
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// Authors:
// John Back

/*! \file LpcConvexHull.hh
    \brief File containing the declaration of the LpcConvexHull class
*/

/*! \class LpcConvexHull
    \brief Class that calculates the "convex hull" extent of a cluster of hits
*/

#ifndef LPC_CONVEX_HULL_HH
#define LPC_CONVEX_HULL_HH

#include "LACE/LpcFunctions.hh"

class LpcCluster;

class LpcConvexHull {

public:

    //! Constructor
    LpcConvexHull();

    //! Destructor
    virtual ~LpcConvexHull();

    //! Find the convex hull lengths for the cluster
    /*!
      \param [in] theCluster The cluster to process
    */
    void findLengths(LpcCluster* theCluster);

protected:

private:

    //! The LpcFunctions object
    LpcFunctions functions_;

};

#endif




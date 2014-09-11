// Copyright University of Warwick 2014
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// Authors:
// John Back

/*! \file LpcShowerAlgorithm.hh
    \brief File containing the declaration of the LpcShowerAlgorithm class
*/

/*! \class LpcShowerAlgorithm
    \brief Class that defines the shower algorithm
*/

#ifndef LPC_SHOWER_ALGORITHM_HH
#define LPC_SHOWER_ALGORITHM_HH

#include <vector>

class LpcCluster;
class LpcParameters;

class LpcShowerAlgorithm {

public:

    //! Constructor
    /*!
      \param [in] thePars a pointer to the constant lpc parameters
    */
    LpcShowerAlgorithm(const LpcParameters* thePars);

    //! Destructor
    virtual ~LpcShowerAlgorithm();

    //! Find out which clusters are potential showers (setting their shower boolean)
    /*!
      \param [in] theClusters The vector of LpcCluster pointers
    */
    void findShowers(const std::vector<LpcCluster*>& theClusters);

    //! Identify the cluster as a shower or track and set its shower boolean
    /*!
      \param [in] theCluster The pointer to the LpcCluster
    */
    void identify(LpcCluster* theCluster);

protected:

    //! Pointer to the constant lpc parameters
    const LpcParameters* thePars_;

private:

    //! Private default constructor
    LpcShowerAlgorithm();

    //! Private copy constructor
    LpcShowerAlgorithm(const LpcShowerAlgorithm& other);

    //! Private assignment operator
    LpcShowerAlgorithm& operator=(const LpcShowerAlgorithm& other);

    //! The convex hull selection cut value
    double convexHull_;

    //! The minimum threshold for hit-to-lpc residuals for showers
    double showerRes_;

    //! The minimum hit-to-lpc residual ratio for showers
    double showerResRatio_;

    //! The minimum fraction of residuals that need to be above the shower threshold
    double showerResFrac_;

};

#endif

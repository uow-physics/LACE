// Copyright University of Warwick 2014
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// Authors:
// John Back

/*! \file LpcBranch.hh
    \brief File containing the declaration of the LpcBranch class
*/

/*! \class LpcBranch
    \brief Class that stores the result of a branch of the main principal curve
*/

#ifndef LPC_BRANCH_HH
#define LPC_BRANCH_HH

#include "LACE/LpcAbsCurve.hh"

#include "LACE/LpcPathLength.hh"
#include "LACE/LpcPoint.hh"

#include <Eigen/Dense>
#include <vector>

class LpcBranch : public LpcAbsCurve {

public:

    //! Constructor storing the local principal curve branch result
    /*!
      \param [in] index The index number of the lpc branch
      \param [in] startPoint The scaled starting position for the branch
      \param [in] lpcPoints The calculated scaled branch points
      \param [in] eigenVectors The largest eigenvectors (row = lpc point, col = x,y,z,...)
      \param [in] cosAngles The cosine of the angles between adjacent eigenvectors
      \param [in] lpcPath The cumulative path length object along the curve
      \param [in] rho The ratio of the eigenvalues (2nd/1st) for each branch point
      \param [in] c0 The adjusted kernel function denominator exponent scale factors
      \param [in] highRhoPoints The scaled local neighbourhood points when rho > rho0 (~0.4)
      \param [in] theHits The pointer to the constant hit collection used to find the curve
      \param [in] flag The status flag set to be the index of the main curve
     */
    LpcBranch(int index, const LpcPoint& startPoint,
	      const std::vector<LpcPoint>& lpcPoints,
	      const Eigen::MatrixXd& eigenVectors, 
	      const Eigen::VectorXd& cosAngles,
	      const LpcPathLength& lpcPath,
	      const Eigen::VectorXd& rho, 
	      const Eigen::VectorXd& c0,
	      const std::vector<LpcPoint>& highRhoPoints,
	      const LpcHitCollection* theHits,
	      int flag);

    //! Empty constructor
    LpcBranch();

    //! Destructor
    virtual ~LpcBranch();

protected:
  
private:

    //! Copy constructor
    LpcBranch(const LpcBranch& other);

    //! Assignment operator
    LpcBranch& operator=(const LpcBranch& other);
};

#endif


// Copyright University of Warwick 2014
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// Authors:
// John Back

/*! \file LpcCurve.hh
    \brief File containing the declaration of the LpcCurve class
*/

/*! \class LpcCurve
    \brief Class that stores the main principal curve, keeping track of any branches.
*/

#ifndef LPC_CURVE_HH
#define LPC_CURVE_HH

#include "LACE/LpcAbsCurve.hh"

#include "LACE/LpcPathLength.hh"
#include "LACE/LpcPoint.hh"

#include <Eigen/Dense>
#include <vector>

class LpcBranchCollection;

class LpcCurve : public LpcAbsCurve {

public:

    //! Constructor storing the local principal curve result
    /*!
      \param [in] index The index number of the lpc curve
      \param [in] startPoint The scaled starting position for the lpc
      \param [in] lpcPoints The calculated scaled lpc points
      \param [in] eigenVectors The largest eigenvectors (row = lpc point, col = x,y,z,...)
      \param [in] cosAngles The cosine of the angles between adjacent eigenvectors
      \param [in] lpcPath The cumulative path length object along the curve
      \param [in] rho The ratio of the eigenvalues (2nd/1st) for each lpc point
      \param [in] c0 The adjusted kernel function denominator exponent scale factors
      \param [in] highRhoPoints The scaled local neighbourhood points when rho > rho0 (~0.4)
      \param [in] theHits The pointer to the constant hit collection used to find the curve
      \param [in] flag An integer status flag (default value of -1 to mean the "main curve")
     */
    LpcCurve(int index, const LpcPoint& startPoint, 
	     const std::vector<LpcPoint>& lpcPoints,
	     const Eigen::MatrixXd& eigenVectors, 
	     const Eigen::VectorXd& cosAngles,
	     const LpcPathLength& lpcPath, 
	     const Eigen::VectorXd& rho, 
	     const Eigen::VectorXd& c0,
	     const std::vector<LpcPoint>& highRhoPoints,
	     const LpcHitCollection* theHits,
	     int flag = -1);

    //! Empty constructor
    LpcCurve();

    //! Destructor
    virtual ~LpcCurve();
  
    // Modifiers

    //! Store the collection of branches for this curve
    /*!
      \param [in] theBranches The branch collection to store
    */
    void storeBranches(LpcBranchCollection* theBranches) {theBranches_ = theBranches;}

    // Accessor methods
    
    //! Get the collection of branches
    /*!
      \returns the pointer to the LpcBranchCollection
    */
    LpcBranchCollection* getBranchCollection() const {return theBranches_;}

    //! Print the main curve as well as branches
    virtual void print() const;

protected:
  
private:

    //! Copy constructor
    LpcCurve(const LpcCurve& other);

    //! Assignment operator
    LpcCurve& operator=(const LpcCurve& other);

    //! The pointer to the collection of possible branches, which this class owns
    LpcBranchCollection* theBranches_;

};

#endif


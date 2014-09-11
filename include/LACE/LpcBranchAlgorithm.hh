// Copyright University of Warwick 2014
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// Authors:
// John Back

/*! \file LpcBranchAlgorithm.hh
    \brief File containing the declaration of the LpcBranchAlgorithm class
*/

/*! \class LpcBranchAlgorithm
    \brief Class that follows the scaled data points to find lpc branches

    This is a C++ implementation of the branch finding algorithm in 
    the R-code package LPCM written by Jochen Einbeck and Ludger Evers
    (http://cran.r-project.org/web/packages/LPCM/index.html)
*/

#ifndef LPC_BRANCH_ALGORITHM_HH
#define LPC_BRANCH_ALGORITHM_HH

#include "LACE/LpcAlgorithm.hh"
#include "LACE/LpcFunctions.hh"
#include "LACE/LpcPoint.hh"

#include <Eigen/Dense>
#include <vector>

class LpcBranch;
class LpcBranchCollection;
class LpcCurve;
class LpcHitCollection;
class LpcParameters;

class LpcBranchAlgorithm {

public:

    //! Constructor
    /*!
      \param [in] thePars a pointer to the constant lpc parameters
    */
    LpcBranchAlgorithm(const LpcParameters* thePars);

    //! Destructor
    virtual ~LpcBranchAlgorithm();
    
 
    //! Get the collection of branches for the lpc
    /*!
      \param [in] theHits the pointer to the constant collection of hits
      \param [in] theHighRhoPoints the vector local neighbourhood points with high eigenvalue ratios
      \param [in] mainCurveIndex The index number of the main curve that is associated to these branches
    */
    LpcBranchCollection* getBranches(const LpcHitCollection* theHits,
				     const std::vector<LpcPoint>& theHighRhoPoints,
				     int mainCurveIndex);
    
protected:
    
    //! Initalise the parameters
    void initialise();

    //! Get the branches for the current generation
    /*!
      \param [in] startIndex The starting index number for the branches (starts at 1)
      \returns a vector of new LpcBranch pointers
    */
    std::vector<LpcBranch*> getGenBranches(int startIndex);

    //! Find the next branch
    /*!
      \param [in] index The index number of the branch (starts at index = 1)
      \returns a new LpcBranch pointer
    */
    LpcBranch* getBranch(int index);

    //! Initialise various private data members for calculating the new forward and backward
    //! starting points away from the main curve, also updating the extra weighting factors
    void prepareNewStartPoint();
    

private:

    //! Private default & copy constructors and assignment operator
    LpcBranchAlgorithm();

    //! Copy constructor
    LpcBranchAlgorithm(const LpcBranchAlgorithm& other);

    //! Assignment operator
    LpcBranchAlgorithm& operator=(const LpcBranchAlgorithm& other);

    //! Pointer to the constant lpc parameters
    const LpcParameters* thePars_;

    //! The number of co-ordinate dimensions
    int nDim_;

    //! Pointer to the constant hit collection
    const LpcHitCollection* theHits_;

    //! The matrix of the scaled positions of all of the hits
    //! Each row = a hit, each column = x, y, z co-ordinate
    Eigen::MatrixXd Xi_;

    //! The 1-row vector of the hit weights
    Eigen::VectorXd weights_;

    //! The scaled starting point for the algorithm
    Eigen::VectorXd x0_;

    //! The starting point as an LpcPoint
    LpcPoint startPoint_;

    //! The co-ordinate range along x, y, and z
    Eigen::VectorXd theRange_;

    //! The local weighted mean
    Eigen::VectorXd mu_x_;

    //! The branching displacement
    Eigen::VectorXd dx_;
    
    //! The local 2nd eigenvector
    Eigen::VectorXd E2Vector_;

    //! The updated weights for the branching sections
    Eigen::VectorXd jWeights_;

    //! The kernel width
    double h_;

    //! The step size
    double t_;

    //! The branch gap step size * the step size
    double gapt_;

    //! The number of lpc points
    int iter2_;

    //! The number of lpc points/2
    int iter_;

    //! The angle penalisation factor
    double pen_;

    //! The ratio of eigenvalues for identifying possible branching points
    double rho0_;

    //! The boundary condition
    double boundary_;

    //! The convergence requirement
    double convergence_;

    //! The number of data points
    int N_;

    //! The threshold level for including branching sections
    double threshold_;

    //! The number of branching levels (generations)
    int branchLevel_;

    //! Object to access the functions within LpcFunctions
    LpcFunctions functions_;

    //! Object to access the main lpc algorithm for use with finding the branch sections
    LpcAlgorithm lpcAlgorithm_;

    //! The starting points for the branching
    std::vector<LpcPoint> startPoints_;

    //! The main curve index that is associated to the branches
    int mainCurveIndex_;

};

#endif

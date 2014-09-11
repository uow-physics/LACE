// Copyright University of Warwick 2014
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// Authors:
// John Back

/*! \file LpcAlgorithm.hh
    \brief File containing the declaration of the LpcAlgorithm class
*/

/*! \class LpcAlgorithm
    \brief Class that follows the scaled data points to find the lpc

    This is a C++ implementation of the main Lpc curve algorithm in 
    the R-code package LPCM written by Jochen Einbeck and Ludger Evers
    (http://cran.r-project.org/web/packages/LPCM/index.html)
*/

#ifndef LPC_ALGORITHM_HH
#define LPC_ALGORITHM_HH

#include "LACE/LpcFunctions.hh"
#include "LACE/LpcPoint.hh"

#include <Eigen/Dense>
#include <vector>

class LpcCurve;
class LpcHitCollection;
class LpcParameters;

class LpcAlgorithm {

public:

    //! Constructor
    /*!
      \param [in] thePars a pointer to the constant lpc parameters
    */
    LpcAlgorithm(const LpcParameters* thePars);

    //! Destructor
    virtual ~LpcAlgorithm();
    
    //! Enumerator specifying the direction of the curve
    enum LpcDirection {Both = 0, Forward = 1, Backward = 2};

    //! Find the local principal curve for the given dataset and return the result
    /*!
      \param [in] curveIndex the index integer number for the curve
      \param [in] theHits a pointer to the constant LpcHitCollection
      \param [in] direction an integer specifying if we want the forward and/or backward parts
      \returns a new pointer to the lpc curve result
    */
    LpcCurve* getCurve(int curveIndex, const LpcHitCollection* theHits, 
		       int direction = LpcAlgorithm::Both);

    //! Find the local principal curve for the given dataset and return the result
    /*!
      \param [in] curveIndex the index integer number for the curve
      \param [in] theHits a pointer to the constant LpcHitCollection
      \param [in] startPoint the starting point which may be different to the centroid of theHits
      \param [in] weights the Eigen VectorXd object of the weights (of each hit) to be used
      \param [in] lastEVect the last, previous eigen vector for branching sections
      \param [in] direction an integer specifying if we want the forward and/or backward parts
      \returns a new pointer to the lpc curve result
    */
    LpcCurve* getCurve(int curveIndex, const LpcHitCollection* theHits, 
		       const LpcPoint& startPoint, const Eigen::VectorXd& weights, 
		       const Eigen::VectorXd& lastEVect, int direction = LpcAlgorithm::Both);

protected:
    
    //! Initialise
    void initialise();

    // ! Find the forward-going part of the lpc
    /*!
      \param [in] lastEVect The last, previous eigenvector for branching sections
    */
    void forward(const Eigen::VectorXd& lastEVect);

    //! Find the backward-going part of the lpc
    /*!
      \param [in] lastEVect The last, previous eigenvector for branching sections
    */
    void backward(const Eigen::VectorXd& lastEVect);

private:

    //! Private default & copy constructors and assignment operator
    LpcAlgorithm();

    //! Copy constructor
    LpcAlgorithm(const LpcAlgorithm& other);

    //! Assignment operator
    LpcAlgorithm& operator=(const LpcAlgorithm& other);

    //! Pointer to the constant lpc parameters
    const LpcParameters* thePars_;

    //! The number of dimensions for each hit
    int nDim_;

    //! The matrix of the scaled positions of all of the hits
    //! Each row = a hit, each column = x, y, z,.. co-ordinate
    Eigen::MatrixXd Xi_;

    //! The 1-row vector of the hit weights
    Eigen::VectorXd weights_;

    //! The scaled starting point for the algorithm
    Eigen::VectorXd x0_;

    //! The starting point as an LpcPoint
    LpcPoint startPoint_;

    //! The co-ordinate ranges (x,y,z,...)
    Eigen::VectorXd theRange_;

    //! The kernel width
    double h_;

    //! The step size
    double t_;

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

    //! The local neighbourhood points ("x" or "u") with eigenvalue ratios > rho0
    std::vector<LpcPoint> highRhoPoints_;

    //! The lpc points ("m(u)" or "mu"), equivalent to the local weighted mean
    Eigen::MatrixXd lpcPoints_;

    //! The principal eigen vectors
    Eigen::MatrixXd eigenVectors_;

    //! The cosine between neighbouring principal eigenvectors
    Eigen::VectorXd cosAngles_;

    //! The values of the cumulative path length
    Eigen::VectorXd lambda_;
 
    //! The values of the ordered cumulative path length, starting at dlambda = 0.
    Eigen::VectorXd dlambda_;

    //! The matrix of the cumulative path length along each axis (row = point, col = x,y,z..)
    Eigen::MatrixXd lambdaAxes_;

    //! The ratio of the eigenvalues between neighbouring eigenvectors
    Eigen::VectorXd rho_;

    //! The c0 coefficients used to scale the kernel factor
    Eigen::VectorXd c0_;

    //! Object to access the functions within LpcFunctions
    LpcFunctions functions_;

};

#endif

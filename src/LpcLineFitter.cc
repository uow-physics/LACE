// Copyright University of Warwick 2014
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// Authors:
// John Back

/*! \file LpcLineFitter.cc
    \brief Class that calculates the regression best fit line through a set of data points
*/

#include "LACE/LpcLineFitter.hh"

#include "LACE/LpcFunctions.hh"

LpcLineFitter::LpcLineFitter()
{
}

LpcLineFitter::~LpcLineFitter()
{
}

std::pair<Eigen::VectorXd, Eigen::VectorXd> 
LpcLineFitter::findLine(const Eigen::MatrixXd& data) const
{

    LpcFunctions functions;

    // Get the non-weighted mean of the data points
    Eigen::VectorXd mean = functions.getMean(data);

    // Form a matrix of the centred data points
    Eigen::MatrixXd meanData = functions.offsetPositions(data, mean);

    // Now find the Singular Value Decomposition (SVD) of the mean data
    // which uses the covariance matrix of this adjusted data. The
    // Rayleigh quotient is minimized and maximized by the eigen vectors of 
    // the covariance matrix that correspond to its smallest and largest 
    // eigenvalues, i.e. min eigen vector = normal vector of best-fit plane,
    // max eigen vector = direction of best fit line. Note that the best line
    // lies on the best plane; they are equivalent.
    
    // Get the (symmetric nDim*nDim) covariance matrix of the mean data
    Eigen::MatrixXd covMatrix = functions.formCovarianceMatrix(meanData);
    
    // Create the Jacobi SVD matrix, with covMatrix = A = U S V*
    Eigen::JacobiSVD<Eigen::MatrixXd> SVD(covMatrix, Eigen::ComputeThinU | Eigen::ComputeThinV);

    // The direction of the best fit line will be the eigenvector corresponding
    // to the largest eigenvalue of the matrix U
    Eigen::MatrixXd U = SVD.matrixU();

    // Get the first column of U, corresponding to the largest eigenvalue, since
    // the Jacobi eigenvalues are automatically in descending order.

    int nDim = data.cols();
    Eigen::VectorXd dirVect = Eigen::VectorXd::Zero(nDim);
    // Ensure U is not empty
    if (U.cols() > 0) {dirVect = U.col(0);}

    // Return the mean and unit direction vector along the regression line
    return std::pair<Eigen::VectorXd, Eigen::VectorXd>(mean, dirVect);

}

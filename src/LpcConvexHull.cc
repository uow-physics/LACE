// Copyright University of Warwick 2014
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// Authors:
// John Back

/*! \file LpcConvexHull.cc
    \brief Class that calculates the "convex hull" extent of a cluster of hits
*/

#include "LACE/LpcConvexHull.hh"

#include "LACE/LpcCluster.hh"

#include <Eigen/Dense>

LpcConvexHull::LpcConvexHull()
{
}

LpcConvexHull::~LpcConvexHull()
{
}

void LpcConvexHull::findLengths(LpcCluster* theCluster)
{

    // Find the convex hull lengths along the spatial dimensions of 
    // the cluster and store them in the cluster pointer. Originally,
    // this method would require the calculation of the convex hull
    // points with QHull. But, the extent of the hull is simply given
    // by the size of the rectangular box that will enclose the
    // cluster "ellipsoid", i.e. the lengths (max-min range) along 
    // each of the principal axes

    if (!theCluster) {return;}

    //Retrieve the matrix of cluster hit positions
    Eigen::MatrixXd hitCoords = theCluster->getHitPositions();
    
    // We need to transform the point co-ordinates to lie along the principal axes.
    // Find the eigenvectors of the covariance matrix of the hit co-ordinates
    Eigen::MatrixXd covMatrix = functions_.getCovarianceMatrix(hitCoords);

    std::pair<Eigen::VectorXd, Eigen::MatrixXd> pcaEigen = 
	functions_.findNormEigenVectors(covMatrix);

    Eigen::MatrixXd eVectors = pcaEigen.second;

    Eigen::MatrixXd eV = eVectors.transpose();

    // Transform the hit co-ordinates. Row = point, col = x,y,z,...
    Eigen::MatrixXd transCoords = hitCoords*eV;

    // Get the lengths along each axis
    int nDim = hitCoords.cols();
    Eigen::VectorXd lengths = Eigen::VectorXd::Zero(nDim);

    for (int i = 0; i < nDim; i++) {

	Eigen::VectorXd coordCol = transCoords.col(i);
	double maxVal = coordCol.maxCoeff();
	double minVal = coordCol.minCoeff();
	double range = fabs(maxVal - minVal);

	lengths(i) = range;

    }

    theCluster->storeConvexHull(lengths);

}


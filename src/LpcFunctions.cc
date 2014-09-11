// Copyright University of Warwick 2014
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// Authors:
// John Back

/*! \file LpcFunctions.cc
    \brief Class containing various functions used by other classes
*/

#include "LACE/LpcFunctions.hh"

LpcFunctions::LpcFunctions()
{
}

LpcFunctions::~LpcFunctions()
{
}

Eigen::VectorXd LpcFunctions::getWeightedMean(const Eigen::MatrixXd& Xi, 
					      const Eigen::VectorXd& weights) const
{
    double totalWeight(0.0);

    int nHits = Xi.rows();
    int nDim = Xi.cols();
    Eigen::VectorXd centroid = Eigen::VectorXd::Zero(nDim);

    for (int i = 0; i < nHits; i++) {

	double w = weights(i);
	totalWeight += w;

	for (int j = 0; j < nDim; j++) {
	    centroid(j) += Xi(i,j)*w;
	}

    }
    
    if (fabs(totalWeight) > 1e-10) {
	centroid /= totalWeight;
    }

    return centroid;

}

Eigen::VectorXd LpcFunctions::getMean(const Eigen::MatrixXd& Xi) const
{

    int nHits = Xi.rows();
    int nDim = Xi.cols();
    Eigen::VectorXd centroid = Eigen::VectorXd::Zero(nDim);

    for (int i = 0; i < nHits; i++) {

	for (int j = 0; j < nDim; j++) {
	    centroid(j) += Xi(i,j);
	}

    }
    
    double N = nHits*1.0;
    if (fabs(N) > 1e-10) {centroid /= N;}

    return centroid;

}


Eigen::MatrixXd LpcFunctions::offsetPositions(const Eigen::MatrixXd& Xi,
					      const Eigen::VectorXd& offset) const
{

    Eigen::MatrixXd offsetXi(Xi);

    for (int i = 0; i < Xi.cols(); i++) {
	// Create a vector with the same number of rows as the Xi matrix
	// but with the constant value of the offset co-ordinate to subtract
	Eigen::VectorXd offsetVect = Eigen::VectorXd::Constant(Xi.rows(), offset(i));

 	offsetXi.col(i) -= offsetVect;
    }

    return offsetXi;

}

double LpcFunctions::kernelFunction(const Eigen::VectorXd& hitPoint,
				    const Eigen::VectorXd& localPoint, 
				    double factor) const
{

    double k(1.0);
    for (int i = 0; i < localPoint.size(); i++) {

	k *= this->kernelPart(hitPoint(i), localPoint(i), factor);

    }

    return k;

}

double LpcFunctions::kernelPart(double x, double u, double factor) const
{

    double expTerm(0.0);
    if (fabs(factor) > 1e-10) {expTerm = (x - u)/factor;}

    double k = exp(-0.5*expTerm*expTerm);
    return k;

}

double LpcFunctions::kdex(const Eigen::MatrixXd& data, 
			  const Eigen::VectorXd& u, double factor) const
{

    double sumKern(0.0);
    int N = data.rows();

    for (int i = 0; i < N; i++) {

	Eigen::VectorXd thePoint = data.row(i);
	sumKern += this->kernelFunction(thePoint, u, factor);

    }

    double kdex(0.0);
    if (N > 0 && fabs(sumKern) > 1e-10) {
	kdex = 1.0/(N*sumKern);
    }

    return kdex;

}

Eigen::MatrixXd LpcFunctions::formCovarianceMatrix(const Eigen::MatrixXd& meanData,
						   const Eigen::VectorXd& weights) const
{

    // Calculate the covariance matrix for the set of points
    // within the (mean subtracted) data matrix (row = hit, col = x,y,z,...)
    // with the vector of weights (each entry = charge or energy of hit)

    // Declare and initialise the covariance matrix.
    // Get the number of co-ordinate variables
    const int nCols = meanData.cols();
    Eigen::MatrixXd covMatrix = Eigen::MatrixXd::Zero(nCols, nCols);

    // Loop over the set of data points (rows) and calculate the
    // top part of the symmetric covariance matrix.
    // The number of co-ordinate dimensions is generalised to nCols,
    // but will usually be equal to 3 (x,y,z)

    for (int ip = 0; ip < meanData.rows(); ip++) {

	double w = weights(ip);
	
	// First, only fill the top right of the symmetric matrix
	for (int ix = 0; ix < nCols; ix++) {

	    double x_i = meanData(ip, ix);

	    for (int iy = ix; iy < nCols; iy++) {

		double y_i = meanData(ip, iy);

		covMatrix(ix, iy) += w*x_i*y_i;

	    }

	}

    }


    // Fill the other part of the symmetric matrix
    for (int ix = 0; ix < nCols; ix++) {

	for (int iy = ix; iy < nCols; iy++) {

	    covMatrix(iy, ix) = covMatrix(ix, iy);

	}

    }

    // Now divide by the sum of the weights
    double sum_weights = weights.sum();
    covMatrix /= sum_weights;

    // Return the covariance matrix
    return covMatrix;

}

Eigen::MatrixXd LpcFunctions::formCovarianceMatrix(const Eigen::MatrixXd& meanData) const
{

    // Number of data points
    int N = meanData.rows();
    // Define unit weights for each point
    Eigen::VectorXd unitWeights = Eigen::VectorXd::Ones(N);

    Eigen::MatrixXd covMatrix = this->formCovarianceMatrix(meanData, unitWeights);

    return covMatrix;    

}

Eigen::MatrixXd LpcFunctions::getCovarianceMatrix(const Eigen::MatrixXd& data,
						  const Eigen::VectorXd& weights) const
{

    // First, create the set of mean centred weighted data
    Eigen::VectorXd mean = this->getWeightedMean(data, weights);
    Eigen::MatrixXd meanData = this->offsetPositions(data, mean);
    Eigen::MatrixXd covMatrix = this->formCovarianceMatrix(meanData, weights);

    return covMatrix;

}

Eigen::MatrixXd LpcFunctions::getCovarianceMatrix(const Eigen::MatrixXd& data) const
{

    // First, create the set of mean centred data
    Eigen::VectorXd mean = this->getMean(data);
    Eigen::MatrixXd meanData = this->offsetPositions(data, mean);
    Eigen::MatrixXd covMatrix = this->formCovarianceMatrix(meanData);

    return covMatrix;

}

std::pair<Eigen::VectorXd, Eigen::MatrixXd> 
LpcFunctions::findNormEigenVectors(const Eigen::MatrixXd& covMatrix) const
{
 
    // Covariance matrix will be symmetric
    const int N = covMatrix.rows();

    Eigen::MatrixXd eigenVectors = Eigen::MatrixXd::Zero(N, N);
    Eigen::VectorXd eigenValues = Eigen::VectorXd::Zero(N);

    // Eigenvectors and eigenvalues
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver(covMatrix);

    // Check if the solver worked OK
    if (eigenSolver.info() == Eigen::Success) {

	Eigen::MatrixXd origEigenVectors = eigenSolver.eigenvectors().real();
	Eigen::VectorXd origEigenValues = eigenSolver.eigenvalues().real();

	// The eigenvectors are in ascending order. Reorder them to be descending.
	// Also arrange the eigenvectors to be along each row 
	// of the returned eigenvector matrix

	for (int i = 0; i < N; i++) {

	    int k = N - i - 1;
	    // Eigenvalues should be in descending order
	    eigenValues(i) = origEigenValues(k);

	    for (int j = 0; j < N; j++) {
		// Each eigenvector = matrix row
		eigenVectors(i, j) = origEigenVectors(j, k);
	    }

	}

    }

    // The eigenvectors should already be normalised.
    // Return the result as a STL pair

    return std::pair<Eigen::VectorXd, Eigen::MatrixXd>(eigenValues, eigenVectors);

}

std::pair<double, double> LpcFunctions::getMeanAndRms(const Eigen::VectorXd& data) const
{

    // Calculate the mean and rms of the set of data recursively
    double mean(0.0), rms(0.0);

    double N(0.0), N1(0.0);

    for (int i = 0; i < data.size(); i++) {

	// The current data value
	double x = data(i);

	// Increment the number of data points
	N += 1.0;
	N1 = N - 1.0;

	double mean1 = mean; // Previous mean

	// Update the mean
	if (fabs(N) > 1e-10) {
	    mean = (N1*mean1 + x)/N;
	}

	// Update the (square) of the rms
	double diff = mean - x;
	if (fabs(N1) > 1e-10) {
	    rms += diff*diff*N/N1;
	}

    }

    // Take the square root and divide by N
    if (fabs(N) > 1e-10) {rms = sqrt(rms/N);}

    return std::pair<double, double>(mean, rms);

}

double LpcFunctions::getPerpLineDist(const Eigen::VectorXd& point,
				     const Eigen::VectorXd& centroid,
				     const Eigen::VectorXd& direction) const
{

    Eigen::VectorXd PC = point - centroid;
    double PCn = PC.dot(direction);
    Eigen::VectorXd PCnn = PCn*direction;

    Eigen::VectorXd d = PC - PCnn;

    double dist = d.norm();
    return dist;

}

double LpcFunctions::getPerpLineDist(const Eigen::VectorXd& point,
				     const Eigen::VectorXd& centroid,
				     const Eigen::VectorXd& direction,
				     double dL) const
{

    double dist(-1.0);

    Eigen::VectorXd PC = point - centroid;
    double PCn = PC.dot(direction);

    // Only continue if the point P is towards the line end-point V 
    // (in the direction of n), but does not go beyond V, i.e. the 
    // distance PC is less than VC)
    if (PCn > 0.0 && PCn < dL) {

	Eigen::VectorXd PCnn = PCn*direction;
	Eigen::VectorXd d = PC - PCnn;
	dist = d.norm();
    }

    return dist;

}

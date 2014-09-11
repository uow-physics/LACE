// Copyright University of Warwick 2014
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// Authors:
// John Back

/*! \file LpcFunctions.hh
    \brief File containing the declaration of the LpcFunctions class
*/

/*! \class LpcFunctions
    \brief Class containing various functions used by other classes
*/

#ifndef LPC_FUNCTIONS_HH
#define LPC_FUNCTIONS_HH

#include <Eigen/Dense>

class LpcFunctions {

public:

    //! Default constructor
    LpcFunctions();

    //! Destructor
    virtual ~LpcFunctions();
    
    //! Weighted mean (centroid) of the hits using Eigen matrices
    /*!
      \param [in] Xi The matrix of hit positions (row = hit, col = x,y,z)
      \param [in] weights The vector of the energy or charge weights of the hits
      \returns The weighted mean (centroid) of the hits as an Eigen vector.
    */
    Eigen::VectorXd getWeightedMean(const Eigen::MatrixXd& Xi, 
				    const Eigen::VectorXd& weights) const;

    //! Non-weighted centroid of the hits using Eigen matrices
    /*!
      \param [in] Xi The matrix of hit positions (row = hit, col = x,y,z)
      \returns The non-weighted centroid of the hits as an Eigen vector.
    */
    Eigen::VectorXd getMean(const Eigen::MatrixXd& Xi) const;


    //! Subtract the Eigen matrix of hit positions with a constant vector offset
    /*!
      \param [in] Xi The matrix of hit positions (row = hit, col = x,y,z)
      \param [in] offset The constant vector position offset
      \returns The updated offset hit positions (row = hit, col = x,y,z)
    */
    Eigen::MatrixXd offsetPositions(const Eigen::MatrixXd& Xi,
				    const Eigen::VectorXd& offset) const;

    //! The kernel function for each hit point
    /*!
      \param [in] hitPoint the data hit point
      \param [in] localPoint the local neighbourhood point u
      \param [in] factor the kernel denominator factor
      \returns the kernel function value
    */
    double kernelFunction(const Eigen::VectorXd& hitPoint,
			  const Eigen::VectorXd& localPoint, 
			  double factor) const;

    //! The kernel part used by the kernel function
    /*!
      \param [in] x The 1D hit co-ordinate variable
      \param [in] u The 1D co-ordinate of the local neighbourhood point
      \param [in] factor the kernel denominator factor
    */
    double kernelPart(double x, double u, double factor) const;

    //! A function used by the LpcBranchAlgorithm
    /*!
      \param [in] data the matrix of zero mean data (row = hit, col = x,y,z)
      \param [in] u the local neighbourhood point
      \param [in] factor the kernel denominator factor
    */
    double kdex(const Eigen::MatrixXd& data, 
		const Eigen::VectorXd& u, double factor) const;

    //! Covariance matrix from zero mean matrix data and vector of weights
    /*!
      \param [in] meanData the matrix of zero mean data (row = hit, col = x,y,z)
      \param [in] weights the vector of hit weights
      \returns the symmetric convariance matrix of the data
    */
    Eigen::MatrixXd formCovarianceMatrix(const Eigen::MatrixXd& meanData,
					 const Eigen::VectorXd& weights) const;

    //! Unweighted ovariance matrix from zero mean matrix data
    /*!
      \param [in] meanData the matrix of zero mean data (row = hit, col = x,y,z)
      \returns the symmetric convariance matrix of the unweighted data
    */
    Eigen::MatrixXd formCovarianceMatrix(const Eigen::MatrixXd& meanData) const;

    //! Covariance matrix from a matrix of data points and vector of weights
    /*!
      \param [in] data the matrix of the data (row = hit, col = x,y,z)
      \param [in] weights the vector of hit weights
      \returns the symmetric convariance matrix of the data
    */
    Eigen::MatrixXd getCovarianceMatrix(const Eigen::MatrixXd& data,
					const Eigen::VectorXd& weights) const;

    //! Covariance matrix from a matrix of data points
    /*!
      \param [in] data the matrix of the data (row = hit, col = x,y,z)
      \returns the symmetric convariance matrix of the data
    */
    Eigen::MatrixXd getCovarianceMatrix(const Eigen::MatrixXd& data) const;


    //! Find the eigenvalues and normalised eigenvectors of the covariance matrix
    /*!
      \param [in] covMatrix The covariance matrix
      \returns a pair of the descending eigenvalues with the normalised eigenvectors (rows)
    */
    std::pair<Eigen::VectorXd, Eigen::MatrixXd>
    findNormEigenVectors(const Eigen::MatrixXd& covMatrix) const;

    //! Calculate the mean and rms of the set of data recursively
    /*!
      \param [in] data The Eigen VectorXd of data points
      \returns the (mean, rms) pair for the set of data
    */
    std::pair<double, double> getMeanAndRms(const Eigen::VectorXd& data) const;

    //! Calculate the perpendicular distance of a point from a line
    /*!
      \param [in] point The point that we want to find the distance for
      \param [in] centroid A point on the line, defined as the centroid
      \param [in] direction The unit direction vector of the line
    */
    double getPerpLineDist(const Eigen::VectorXd& point,
			   const Eigen::VectorXd& centroid,
			   const Eigen::VectorXd& direction) const;

    //! Calculate the perpendicular distance of a point from a line
    //! for points that lie towards the end-point of the line V
    /*!
      \param [in] point The point that we want to find the distance for
      \param [in] centroid A point on the line, defined as the centroid
      \param [in] direction The unit direction vector of the line
      \param [in] dL The distance along the line between the centroid and end-point
    */
    double getPerpLineDist(const Eigen::VectorXd& point,
			   const Eigen::VectorXd& centroid,
			   const Eigen::VectorXd& direction,
			   double dL) const;


protected:
    
private:

};

#endif

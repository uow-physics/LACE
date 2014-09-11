// Copyright University of Warwick 2014
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// Authors:
// John Back

/*! \file LpcPathLength.hh
    \brief File containing the declaration of the LpcPathLength class
*/

/*! \class LpcPathLength
    \brief Class used to define a series of pathlengths along the lpc
*/

#ifndef LPC_PATH_LENGTH_HH
#define LPC_PATH_LENGTH_HH

#include <Eigen/Dense>

class LpcPathLength {

public:

    //! Default, empty constructor
    LpcPathLength();

    //! Constructor using the results of the scaled pathlengths found by LpcFollow
    /*!
      \param [in] sLambda the scaled cumulative 3d path length
      \param [in] sLambdaAxes the scaled cumulative path length along each co-ordinate axis
      \param [in] theRange the scaling range for x,y,z
    */
    LpcPathLength(const Eigen::VectorXd& sLambda, 
		  const Eigen::MatrixXd& sLambdaAxes,
		  const Eigen::VectorXd& theRange);

    //! Copy constructor
    /*!
      \param [in] other The LpcPathLength to copy
    */
    LpcPathLength(const LpcPathLength& other);

    //! Assignment operator using another LpcPathLength
    /*!
      \param [in] other The LpcPathLength used in the assignment
    */
    LpcPathLength& operator=(const LpcPathLength& other);

    //! Destructor
    virtual ~LpcPathLength();


    // Accessors.

    //! Get the number of pathlength segments
    /*!
      \return the number of pathlength segments
    */
    int getN() const {return N_;}

    //! Get the scaled cumulative 3d path length as an Eigen VectorXd
    /*!
      \returns the scaled cumulative 3d path length
    */
    Eigen::VectorXd getScaledLambda() const {return sL_;}

    //! Get the scaled cumulative path length along x as an Eigen VectorXd
    /*!
      \returns the scaled cumulative path length along x
    */
    Eigen::MatrixXd getScaledLambdaAxes() const {return sLAxes_;}

    //! Get the scaling co-ordinate range
    /*!
      \returns the scaling co-ordinate range as an Eigen 3-vector
    */
    Eigen::VectorXd getRange() const {return theRange_;}

    //! Get the ordered (start = 0) scaled cumulative 3d path length as an Eigen VectorXd
    /*!
      \returns the ordered scaled cumulative 3d path length
    */
    Eigen::VectorXd getScaledDeltaLambda() const {return sdL_;}


    //! Get the unscaled cumulative 3d path length as an Eigen VectorXd
    /*!
      \returns the cumulative 3d path length
    */
    Eigen::VectorXd getLambda() const {return L_;}

    //! Get the unscaled cumulative path length along each axis (row = point, col = x,y,z,...)
    /*!
      \returns the unscaled cumulative path length along each axis as an Eigen MatrixXd
    */
    Eigen::MatrixXd getLambdaAxes() const {return LAxes_;}

    //! Get the unscaled ordered (start = 0) cumulative path length as an Eigen VectorXd
    /*!
      \returns the ordered unscaled cumulative path length
    */
    Eigen::VectorXd getDeltaLambda() const {return dL_;}

 
protected:


private:

    //! The number of pathlength segments
    int N_;

    //! The number of co-ordinate dimensions
    int nDim_;

    //! The values of the scaled cumulative path length
    Eigen::VectorXd sL_;
 
    //! The values of the scaled cumulative path length along each co-ordinate axis
    Eigen::MatrixXd sLAxes_;

    //! The values of the co-ordinate scaling range
    Eigen::VectorXd theRange_;

    //! The values of the scaled ordered cumulative path length, starting at 0
    Eigen::VectorXd sdL_;

    //! The values of the unscaled cumulative path length
    Eigen::VectorXd L_;
 
    //! The values of the unscaled cumulative path length along each co-ordinate axis
    Eigen::MatrixXd LAxes_;

    //! The values of the unscaled ordered cumulative path length, starting at 0
    Eigen::VectorXd dL_;

    //! Initialise function
    void initialise();

    //! Helper function to find the path length difference between lpc points
    /*!
      \param [in] i the index of the first lpc point
      \param [in] j the index of the second (neighbouring) lpc point
      \returns an Eigen vector of the x,y,z,... path length between the points
    */
    Eigen::VectorXd getDiff(int i, int j);

};

#endif


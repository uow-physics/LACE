// Copyright University of Warwick 2014
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// Authors:
// John Back

/*! \file LpcPoint.hh
    \brief File containing the declaration of the LpcPoint class
*/

/*! \class LpcPoint
    \brief Class that defines the scaled and unscaled co-ordinates of a point in the LpcAbsCurve
*/

#ifndef LPC_POINT_HH
#define LPC_POINT_HH

#include <Eigen/Dense>
#include <iostream>

class LpcPoint  {

public:

    //! Default, empty constructor
    LpcPoint();

    //! Constructor using an Eigen VectorXd containing the scaled co-ordinates
    /*!
      \param [in] index The lpc point index number
      \param [in] scaledCoords The Eigen vector containing the scaled co-ordinates
    */
    LpcPoint(int index, const Eigen::VectorXd& scaledCoords);

    //! Constructor using an Eigen VectorXd containing the scaled & unscaled co-ordinates
    /*!
      \param [in] index The lpc point index number
      \param [in] scaledCoords The Eigen vector containing the scaled co-ordinates
      \param [in] coords The Eigen vector containing the unscaled co-ordinates
    */
    LpcPoint(int index, const Eigen::VectorXd& scaledCoords,
	     const Eigen::VectorXd& coords);
    
    //! Copy constructor
    /*!
      \param [in] other The LpcPoint to copy
    */
    LpcPoint(const LpcPoint& other);
    
    //! Assignment operator using another LpcPoint
    /*!
      \param [in] other The LpcPoint used for assignment
    */
    LpcPoint& operator=(const LpcPoint& other);

    //! Destructor
    virtual ~LpcPoint();

    //! Find the unscaled co-ordinates using an Eigen vector
    /*!
      \param [in] theRange The scaling range as an Eigen vector
    */
    void unscale(const Eigen::VectorXd& theRange);

    // Accessor functions

    //! Get the index number of the lpc point
    /*!
      \return the index number of the lpc point
    */
    int getIndex() const {return index_;}

    //! Return the number of co-ordinate dimensions
    /*!
      \return the number of co-ordinate dimensions
    */
    int getNDimensions() const {return nDim_;}

    //! Return the scaled co-ordinates as an Eigen vector object
    /*!
      \return the scaled co-ordinates as an Eigen vector
    */
    Eigen::VectorXd getScaledCoords() const {return scaledCoords_;}

    //! Return the unscaled co-ordinates as an Eigen vector object
    /*!
      \return the unscaled co-ordinates as an Eigen vector
    */
    Eigen::VectorXd getCoords() const {return coords_;}

    //! Overload the ostream insertion operator to print out the LpcPoint info
    /*!
      \param [in] os The string stream
      \param [in] thePoint The LpcPoint object to print out
    */
    friend std::ostream& operator << (std::ostream& os, const LpcPoint& thePoint);

protected:

    //! The integer index number
    int index_;

    //! The number of co-ordinate dimensions
    int nDim_;

    //! The scaled co-ordinates
    Eigen::VectorXd scaledCoords_;

    //! The unscaled co-ordinates
    Eigen::VectorXd coords_;

private:

};

#endif

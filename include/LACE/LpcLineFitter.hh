// Copyright University of Warwick 2014
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// Authors:
// John Back

/*! \file LpcLineFitter.hh
    \brief File containing the declaration of the LpcLineFitter class
*/

/*! \class LpcLineFitter
    \brief Class that calculates the regression best fit line through a set of data points
*/

#ifndef LPC_LINE_FITTER_HH
#define LPC_LINE_FITTER_HH

#include <Eigen/Dense>

class LpcLineFitter {

public:

    //! Constructor
    LpcLineFitter();

    //! Destructor
    virtual ~LpcLineFitter();

    //! Find the centroid and unit direction vector for the line
    /*!
      \param [in] data The Eigen MatrixXd of the data (row = hit, col = x,y,z,.. coords)
      \returns an STL pair of the centroid and unit direction vector
    */
    std::pair<Eigen::VectorXd, Eigen::VectorXd> findLine(const Eigen::MatrixXd& data) const;

protected:

private:

};

#endif

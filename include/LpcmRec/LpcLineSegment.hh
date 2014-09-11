// Copyright University of Warwick 2014
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// Authors:
// John Back

/*! \file LpcLineSegment.hh
    \brief File containing the declaration of the LpcLineSegment class
*/

/*! \class LpcLineSegment
    \brief Class that defines the line segments used for vertexing/clustering
*/

#ifndef LPC_LINE_SEGMENT_HH
#define LPC_LINE_SEGMENT_HH

#include <Eigen/Dense>

class LpcLineSegment {

public:

    //! Default, empty constructor
    LpcLineSegment();

    //! Constructor defining all information for a line segment
    /*!
      \param [in] centroid The point on the line corresponding to the centroid of the hits
      \param [in] direction The unit direction along the line
      \param [in] residualMean The mean value of the hit-to-lpc residuals for the segment
      \param [in] residualRms The rms value of the hit-to-lpc residuals for the segment
      \param [in] residualCut The cut on the hit-to-lpc residuals to be used for further analysis
    */
    LpcLineSegment(const Eigen::VectorXd& centroid, const Eigen::VectorXd& direction,
		   double residualMean, double residualRms, double residualCut);

    //! Copy constructor
    /*!
      \param [in] other The LpcLineSegment to copy
    */
    LpcLineSegment(const LpcLineSegment& other);

    //! Assignment operator
    /*!
      \param [in] other The LpcLineSegment used for assignment
    */
    LpcLineSegment& operator=(const LpcLineSegment& other);

    //! Destructor
    virtual ~LpcLineSegment();

    // Modifiers

    //! Set the end point of the line. The direction vector can change sign in
    //! order to point to this from the centroid
    void setEndPoint(const Eigen::VectorXd& endPoint);

    // Accessors

    //! Get the centroid position
    /*!
      \returns the centroid as an Eigen VectorXd
    */
    Eigen::VectorXd getCentroid() const {return centroid_;}

    //! Get the unit direction vector
    /*!
      \returns the unit direction vector as an Eigen VectorXd
    */
    Eigen::VectorXd getDirection() const {return direction_;}

    //! Get the mean of the hit-to-lpc residuals for the segment
    /*!
      \returns the mean of the hit-to-lpc residuals
    */
    double getResidualMean() const {return residualMean_;}

    //! Get the rms of the hit-to-lpc residuals for the segment
    /*!
      \returns the rms of the hit-to-lpc residuals
    */
    double getResidualRms() const {return residualRms_;}

    //! Get the selection cut for the hit-to-lpc residuals for the segment
    /*!
      \returns the selection cut for the hit-to-lpc residuals
    */
    double getResidualCut() const {return residualCut_;}

    //! Get the end point of the line
    /*!
      \returns the end point
    */
    Eigen::VectorXd getEndPoint() const {return endPoint_;}

    //! Print out the line segment information
    void print() const;

protected:


private:

    //! The centroid position
    Eigen::VectorXd centroid_;

    //! The direction vector
    Eigen::VectorXd direction_;

    //! The residual mean
    double residualMean_;

    //! The residual rms
    double residualRms_;

    //! The residual cut
    double residualCut_;

    //! The end point of the line
    Eigen::VectorXd endPoint_;

};

#endif


// Copyright University of Warwick 2014
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// Authors:
// John Back

/*! \file LpcResiduals.hh
    \brief File containing the declaration of the LpcResiduals class
*/

/*! \class LpcResiduals
    \brief Class that stores the lpc-to-hit residuals
*/

#ifndef LPC_RESIDUALS_HH
#define LPC_RESIDUALS_HH

#include <Eigen/Dense>
#include <vector>

class LpcResiduals {

public:

    //! Constructor using residual arrays stored as Eigen vector objects
    /*!
      \param [in] lpcResiduals the average residual over all hits for each lpc point
      \param [in] hitResiduals the residual of the nearest lpc point for each hit
      \param [in] lpcWeightRes the average weighted residual over all hits for each lpc point
      \param [in] hitWeightRes the weighted residual of the nearest lpc point for each hit
      \param [in] hitNearestLpc the nearest lpc point index number for each hit
    */
    LpcResiduals(const Eigen::VectorXd& lpcResiduals,
		 const Eigen::VectorXd& hitResiduals,
		 const Eigen::VectorXd& lpcWeightRes,
		 const Eigen::VectorXd& hitWeightRes,
		 const Eigen::VectorXi& hitNearestLpc);

    //! Default, empty constructor
    LpcResiduals();

    //! Copy constructor
    /*!
      \param [in] other The LpcResiduals object to copy
    */
    LpcResiduals(const LpcResiduals& other);

    //! Assignment operator
    /*!
      \param [in] other The LpcResiduals used for the assignment
    */
    LpcResiduals& operator=(const LpcResiduals& other);

    //! Destructor
    virtual ~LpcResiduals();

    // Accessors

    //! Get the average residual over all hits for each lpc point
    /*!
      \returns an Eigen VectorXd of the average residual with the hits per lpc point
    */
    Eigen::VectorXd getLpcResiduals() const {return lpcResiduals_;}

    //! Get the residual of the nearest lpc point for each hit
    /*!
      \returns an Eigen VectorXd of the residual of the nearest lpc point for each hit
    */
    Eigen::VectorXd getHitResiduals() const {return hitResiduals_;}

    //! Get the average weighted residual over all hits for each lpc point
    /*!
      \returns an Eigen VectorXd of the average weighted residual with the hits per lpc point
    */
    Eigen::VectorXd getWeightedLpcResiduals() const {return lpcWeightRes_;}

    //! Get the weighted residual of the nearest lpc point for each hit
    /*!
      \returns an Eigen VectorXd of the weighted residual of the nearest lpc point for each hit
    */
    Eigen::VectorXd getWeightedHitResiduals() const {return hitWeightRes_;}

    //! Get the nearest lpc point index number for each hit
    /*!
      \returns an Eigen VectorXi of the index number of the nearest lpc point for each hit
    */
    Eigen::VectorXi getHitNearestLpc() const {return hitNearestLpc_;}


    // STL vector versions of the above access functions

    //! Get the average residual over all hits for each lpc point
    /*!
      \returns a STL vector of the average residual with the hits per lpc point
    */
    std::vector<double> getLpcResVector() const;

    //! Get the residual of the nearest lpc point for each hit
    /*!
      \returns a STL vector of the residual of the nearest lpc point for each hit
    */
    std::vector<double> getHitResVector() const;

    //! Get the average weighted residual over all hits for each lpc point
    /*!
      \returns a STL vector of the average weighted residual with the hits per lpc point
    */
    std::vector<double> getWeightedLpcResVector() const;

    //! Get the weighted residual of the nearest lpc point for each hit
    /*!
      \returns a STL vector of the weighted residual of the nearest lpc point for each hit
    */
    std::vector<double> getWeightedHitResVector() const;

    //! Get the nearest lpc point index number for each hit
    /*!
      \returns a STL vector of the index number of the nearest lpc point for each hit
    */
    std::vector<int> getHitNearestLpcVector() const;


protected:

private:

    //! The average residual over all hits for each lpc point
    Eigen::VectorXd lpcResiduals_;

    //! The residual of the nearest lpc point for each hit
    Eigen::VectorXd hitResiduals_;

    //! The average weighted residual over all hits for each lpc point
    Eigen::VectorXd lpcWeightRes_;

    //! The weighted residual of the nearest lpc point for each hit
    Eigen::VectorXd hitWeightRes_;

    //! The nearest lpc point index number for each hit
    Eigen::VectorXi hitNearestLpc_;

};

#endif

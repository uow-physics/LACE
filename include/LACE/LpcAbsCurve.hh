// Copyright University of Warwick 2014
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// Authors:
// John Back

/*! \file LpcAbsCurve.hh
    \brief File containing the declaration of the LpcAbsCurve class
*/

/*! \class LpcAbsCurve
    \brief Class that defines an abstract local principal curve
*/

#ifndef LPC_ABS_CURVE_HH
#define LPC_ABS_CURVE_HH

#include "LACE/LpcBinRange.hh"
#include "LACE/LpcPathLength.hh"
#include "LACE/LpcPoint.hh"
#include "LACE/LpcResiduals.hh"

#include <Eigen/Dense>
#include <vector>

class LpcHitCollection;

class LpcAbsCurve {

public:

    //! Constructor storing the local principal curve result
    /*!
      \param [in] index The index number of the curve
      \param [in] startPoint The scaled starting position for the lpc
      \param [in] lpcPoints The calculated scaled lpc points
      \param [in] eigenVectors The largest eigenvectors (row = lpc point, col = x,y,z,...)
      \param [in] cosAngles The cosine of the angles between adjacent eigenvectors
      \param [in] lpcPath The cumulative path length object along the curve
      \param [in] rho The ratio of the eigenvalues (2nd/1st) for each lpc point
      \param [in] c0 The adjusted kernel function denominator exponent scale factors
      \param [in] highRhoPoints The scaled local neighbourhood points when rho > rho0 (~0.4)
      \param [in] theHits The pointer to the constant hit collection used to find the curve
      \param [in] flag An integer status flag (default value of -1 to mean the "main curve")
     */
    LpcAbsCurve(int index,
		const LpcPoint& startPoint, 
		const std::vector<LpcPoint>& lpcPoints,
		const Eigen::MatrixXd& eigenVectors, 
		const Eigen::VectorXd& cosAngles,
		const LpcPathLength& lpcPath, 
		const Eigen::VectorXd& rho, 
		const Eigen::VectorXd& c0,
		const std::vector<LpcPoint>& highRhoPoints,
		const LpcHitCollection* theHits,
		int flag = -1);

    //! Empty constructor
    LpcAbsCurve();

    //! Destructor
    virtual ~LpcAbsCurve();
  
    //! Enumeration to specify the curve type
    enum LpcCurveType {Abstract = -1, Main = 0, Branch = 1};

    // Accessor methods for retrieving the results in either the original internal
    // storage format or using helper classes, e.g. LpcPoint of LpcPathlength

    //! Get the curve type (LpcCurveType enum)
    /*!
      \return an integer to specify the curve type (LpcCurveType enum)
    */
    int getType() const {return type_;}

    //! Get the curve index number 
    /*!
      \return the curve index number
    */
    int getIndex() const {return index_;}

    //! Get the status flag integer (-1 = main curve; equal to main curve index for branches)
    /*!
      \return the status flag integer
    */
    int getFlag() const {return flag_;}

    //! Get the number of dimensions
    /*!
      \return the number of dimensions
    */
    int getNDimensions() const {return nDim_;}

    //! Get the starting point as an LpcPoint object
    /*!
      \returns the starting point as an LpcPoint object with index = -1
    */
    LpcPoint getStartPoint() const {return startPoint_;}

    //! Get the number of lpc points
    /*!
      \returns the number of lpc points
    */
    int getNLpcPoints() const {return lpcPoints_.size();}

    //! Get the lpc points as a vector of LpcPoint objects
    /*!
      \returns a vector of LpcPoint objects
    */
    std::vector<LpcPoint> getLpcPoints() const {return lpcPoints_;}

    //! Get the principal eigenvectors with the largest eigenvalues of the local covariance matrix
    /*!
      \returns the eigenvectors as an Eigen MatrixXd (row = eigenvector, column = x,y,z,... components)
    */
    Eigen::MatrixXd getEigenVectors() const {return eigenVectors_;}

    //! Get the cosine of the angles between neighbouring principal (largest) eigenvectors
    /*!
      \returns the cosine of the angles as an Eigen VectorXd
    */
    Eigen::VectorXd getCosAngles() const {return cosAngles_;}

    //! Get the cumulative path length info along the curve (scaled and unscaled)
    /*!
      \returns the LpcPathlengths object containing all path length segments
    */
    LpcPathLength getPathLength() const {return lpcPath_;}

    //! Get the ratio of the eigenvalues (2nd/1st) for each lpc point
    /*!
      \returns the ratio of the eigenvalues as an Eigen VectorXd
    */
    Eigen::VectorXd getRho() const {return rho_;}

    //! Get the c0 coefficients used to scale the kernel factor
    /*!
      \returns the c0 coefficients as an Eigen VectorXd
    */
    Eigen::VectorXd getc0() const {return c0_;}

    //! Get the local neighbourhood points ("x" or "u") with eigenvalue ratios > rho0
    /*!
      \returns a stl vector of LpcPoint objects, with both scaled and unscaled co-ordinates
    */
    std::vector<LpcPoint> getHighRhoPoints() const {return highRhoPoints_;}

    //! Get the pointer to the constant hit collection used for this curve
    /*!
      \returns a pointer to the constant hit collection used for the curve
    */
    const LpcHitCollection* getHitCollection() const {return theHits_;}

    //! Get the lpc-to-hit residuals object (use it to access individual residuals)
    /*!
      \returns the lpc residuals object
    */
    LpcResiduals getResiduals() const {return residuals_;}

    //! Print method, which is read-only (no contents can be changed)
    virtual void print() const;

    //! Reset the curve index number
    /*!
      \param [in] index the new index number
    */
    void setIndex(int index) {index_ = index;}

    //! Reset the status flag integer
    /*!
      \param [in] flag The status flag integer
    */
    void setFlag(int flag) {flag_ = flag;}

    //! Store the lpc point peak and range for the 1-|cosAngle| peaks (LpcFeatures class)
    /*!
      \param [in] passedPeaks Vector of the lpc point indices that correspond to the peaks
      \param [in] peakRanges Vector of LpcBinRanges to specify the neighbourhood of each peak
    */
    void storeCosPeaks(const std::vector<int>& passedPeaks, 
		       const std::vector<LpcBinRange>& peakRanges) {
	passedPeaks_ = passedPeaks; peakRanges_ = peakRanges;}

    //! Retrieve the lpc point integers that correspond to the cosine angle peaks
    /*!
      \return the vector of lpc point integers corrsponding to the peaks
    */
    std::vector<int> getCosPeakIndices() const {return passedPeaks_;}

    //! Retrieve the lpc point ranges corresponding to the neighbourhood of the cosine peaks
    /*!
      \return the vector of LpcBinRange objects for each cosine peak range
    */
    std::vector<LpcBinRange> getCosPeakRanges() const {return peakRanges_;}


protected:

    //! Method to calculate the lpc-to-hit residuals and store them in residuals_
    void findResiduals();

    //! Integer specifying the curve type defined by the LpcCurveType enumeration
    int type_;

    //! The index number of the curve
    int index_;

    //! The starting local point of the lpc
    LpcPoint startPoint_;

    //! The lpc points ("m(u)" or "mu"), equivalent to the local weighted mean
    std::vector<LpcPoint> lpcPoints_;

    //! The largest eigenvectors of the covariance matrix (row = lpc point, col = x,y,z)
    Eigen::MatrixXd eigenVectors_;

    //! The cosine between neighbouring principal eigenvectors
    Eigen::VectorXd cosAngles_;

    //! The path lengths along the lpc (scaled and unscaled)
    LpcPathLength lpcPath_;

    //! The ratio of the eigenvalues between neighbouring eigenvectors
    Eigen::VectorXd rho_;

    //! The c0 coefficients used to scale the kernel factor
    Eigen::VectorXd c0_;

    //! The local neighbourhood points ("x" or "u") with eigenvalue ratios > rho0
    std::vector<LpcPoint> highRhoPoints_;

    //! The pointer to the constant used hit collection which is owned by the event
    const LpcHitCollection* theHits_;

    //! The number of dimensions
    int nDim_;

    //! The integer status flag
    int flag_;

    //! The object storing the lpc to hit residuals
    LpcResiduals residuals_;

    //! The vector of lpc point indices for the cosine peaks
    std::vector<int> passedPeaks_;

    //! The vector of LpcBinRanges for the lpc point range around each cosine peak
    std::vector<LpcBinRange> peakRanges_;

private:

    //! Copy constructor
    LpcAbsCurve(const LpcAbsCurve& other);

    //! Assignment operator
    LpcAbsCurve& operator=(const LpcAbsCurve& other);

};

#endif


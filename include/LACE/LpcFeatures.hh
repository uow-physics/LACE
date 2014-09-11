// Copyright University of Warwick 2014
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// Authors:
// John Back

/*! \file LpcFeatures.hh
    \brief File containing the declaration of the LpcFeatures class
*/

/*! \class LpcFeatures
    \brief Class containing the processing steps to find features for a given lpc
*/

#ifndef LPC_FEATURES_HH
#define LPC_FEATURES_HH

#include "LACE/LpcBinRange.hh"

#include <Eigen/Dense>
#include <vector>

class LpcAbsCurve;
class LpcCurve;
class LpcParameters;

class LpcFeatures {

public:

    //! Constructor
    /*!
      \param [in] theParameters a pointer to the constant lpc parameters
    */
    LpcFeatures(const LpcParameters* theParameters);

    //! Destructor
    virtual ~LpcFeatures();

    //! Process the main curve to find all relevant cosine peak ranges
    /*!
      \param [in] theCurve The main curve
    */
    void findFeatures(LpcCurve* theCurve);

    //! Find the cosine peak ranges and store them in the lpc. This is called by process.
    /*!
      \param [in] theCurve The main curve or any one of its branches
    */
    void findCosPeakRanges(LpcAbsCurve* theCurve);

protected:

private:

    //! Private default & copy constructors and assignment operator
    LpcFeatures();

    //! Copy constructor
    LpcFeatures(const LpcFeatures& other);

    //! Assignment operator
    LpcFeatures& operator=(const LpcFeatures& other);

    //! Initialisation function
    void initialise();

    //! Print out a vector of LpcBinRange objects
    /*!
      \param [in] binRanges The given vector of LpcBinRange objects
    */
    void printBinRanges(const std::vector<LpcBinRange>& binRanges) const;

    //! Setup a struct to store (bin, height) information for the 1-|cosAngle| distribution
    /*!
      \struct LpcPeakInfo
      \brief Struct to store (bin, height) information for  the 1-|cosAngle| distribution
    */
    struct LpcPeakInfo {

	//! The bin
	int bin_;
	//! The height
	double height_;

    };

    //! Find the lpc point numbers corresponding to the cosine peaks
    /*!
      \param [in] absCosArray The vector of 1-|cosAngle| values
      \returns the lpc point peak positions of the distribution
    */
    void getCosPeaks(const Eigen::VectorXd& absCosArray);

    //! Find the lpc point numbers corresponding to the cosine peaks using a simple method
    /*!
      \param [in] absCosArray The vector of 1-|cosAngle| values
      \returns the lpc point peak positions of the distribution
    */
    void getSimplePeaks(const Eigen::VectorXd& absCosArray);

    //! Find the lpc point numbers corresponding to the cosine peaks using TSpectrum
    /*!
      \param [in] absCosArray The vector of 1-|cosAngle| values
      \returns the lpc point peak positions of the distribution
    */
    void getTSpectrumPeaks(const Eigen::VectorXd& absCosArray);


    //! Pointer to the constant lpc parameters
    const LpcParameters* thePars_;

    //! The integer specifying the peak finder method
    int peakFinder_;

    //! The minimum threshold cut for 1-|cosAngle|
    double cosAngleCut_;

    //! The minimum fraction for the next peak w.r.t the previous peak
    double minPeakFrac_;

    //! Set the "distance" (squared) limit for checking if lpc peaks "overlap"
    double peakDiffSq_;

    //! Peak range bin differences
    int dBin_, dBin1_;

    //! Vector of LpcPeakInfo (bin, height) pairs
    std::vector<LpcPeakInfo> allPeaks_;

    // ! Predicate for sorting LpcPeakPairs object with peak height in descending order
    struct greater {
	bool operator()(const LpcPeakInfo& a, const LpcPeakInfo& b) const {
	    return a.height_ > b.height_;
	}

    };

};

#endif

// Copyright University of Warwick 2014
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// Authors:
// John Back

/*! \file LpcFeatures.cc
    \brief Class containing the processing steps to find features for a given lpc
*/

#include "LACE/LpcFeatures.hh"

#include "LACE/LpcAbsCurve.hh"
#include "LACE/LpcBranch.hh"
#include "LACE/LpcBranchCollection.hh"
#include "LACE/LpcCurve.hh"
#include "LACE/LpcParameters.hh"
#include "LACE/LpcPoint.hh"

#ifdef LPC_USE_ROOT
#include "TSpectrum.h"
#endif

#include <algorithm>
#include <iostream>

LpcFeatures::LpcFeatures(const LpcParameters* theParameters) :
    thePars_(theParameters),
    peakFinder_(0),
    cosAngleCut_(0.01),
    minPeakFrac_(0.01),
    peakDiffSq_(3.0),
    dBin_(2),
    dBin1_(3),
    allPeaks_()
{
    this->initialise();
}

LpcFeatures::~LpcFeatures()
{
}

void LpcFeatures::initialise()
{

    if (thePars_) {

	// The abs(cosAngle) distribution has spikes pointing downwards.
        // We use 1.0 - abs(cosAngle) to get these peaks pointing upwards.
	// Specify the minimum threshold cut for identifying peaks
	cosAngleCut_ = thePars_->getCosAngleCut();

        // Specify what must be the next peaks minimum height w.r.t the 
        // previous peak, e.g. frac = 0.1 => 2nd peak must have 
        // height >= 10% of 1st peak, etc..
	minPeakFrac_ = thePars_->getMinPeakFrac();

	// Set the distance (squared) limit for checking whether an lpc point 
        // corresponding to a peak is close to other lpc peak positions.
        // In this case, the lpc point ranges are merged, avoiding duplicate
        // vertices that will be essentially the same point
	peakDiffSq_ = thePars_->getPeakDiffSq();

	// Select the peak finder method
	peakFinder_ = thePars_->getPeakFinder();

    }

}

void LpcFeatures::findFeatures(LpcCurve* theCurve)
{

    // Process the main curve and any of its branches
    if (!theCurve) {return;}

    this->findCosPeakRanges(theCurve);

    LpcBranchCollection* branchCollection = theCurve->getBranchCollection();

    if (branchCollection) {

	int nLevels = branchCollection->getNumberLevels();
	    
	// Loop over the number of generations (starts at 1 not 0)
	for (int iL = 1; iL <= nLevels; iL++) {
	    
	    std::vector<LpcBranch*> branches = branchCollection->getBranches(iL);
	    std::vector<LpcBranch*>::iterator bIter;
	    for (bIter = branches.begin(); bIter != branches.end(); ++bIter) {
		
		LpcBranch* theBranch = *bIter;

		this->findCosPeakRanges(theBranch);
		
	    }
	    
	}
	
    }

}

void LpcFeatures::findCosPeakRanges(LpcAbsCurve* theCurve)
{

    // First, get the peak positions for the 1-|cosAngle| distribution
    allPeaks_.clear();
    if (!theCurve) {return;}

    // Get the vector of LpcPoints for the curve
    std::vector<LpcPoint> lpcPoints = theCurve->getLpcPoints();

    // Create a vector that will store the lpc number for each passed peak
    std::vector<int> passedPeaks;

    // Also store the LpcPoints that correspond to a unique peak
    std::vector<LpcPoint> peakCoord;

    // Need to also store a list of the peak bin ranges
    std::vector<LpcBinRange> theRanges;

    // The number of lpc points, and other integers to check for bin boundaries
    int nLpc = theCurve->getNLpcPoints();
    int nLpc1 = nLpc - 1;
    int nLpc2 = nLpc - 2;

    Eigen::VectorXd cosAngles = theCurve->getCosAngles();
    Eigen::VectorXd absCosArray = Eigen::VectorXd::Zero(nLpc);
    for (int i = 0; i < nLpc; i++) {

	absCosArray(i) = 1.0 - fabs(cosAngles(i));

    }

    // Get the vector of (bin, height) pairs for the cosine peaks
    this->getCosPeaks(absCosArray);

    // If no peaks are found, stop 
    if (allPeaks_.size() < 1) {return;}

    // Sort the (bin, height) pairs, where height is in descending order
    std::sort(allPeaks_.begin(), allPeaks_.end(), greater());
  
    int nPassedPeaks(0);
    double prevPeakHeight(0.0);
  
    // Loop over the sorted peaks to find a first pass of the cosine ranges 
    std::vector<LpcPeakInfo>::const_iterator iter;
    for (iter = allPeaks_.begin(); iter != allPeaks_.end(); ++iter) {

	LpcFeatures::LpcPeakInfo peak = *iter;
	int peakBin0 = peak.bin_;
	double peakHeight0 = peak.height_;
	
	// Sometimes the peak is +- 1 bin away. 
	// Check to see which has the largest height
	int peakBin1 = peakBin0 - 1;
	double peakHeight1(0.0);
	if (peakBin1 >= 0) {
	    peakHeight1 = absCosArray(peakBin1);
	}

	int peakBin2 = peakBin0 + 1;
	double peakHeight2(0.0);
	if (peakBin2 < nLpc) {
	    peakHeight2 = absCosArray(peakBin2);
	}

	int peakBin(peakBin0);
	double peakHeight(peakHeight0);
	if (peakHeight1 > peakHeight) {
	    peakBin = peakBin1;
	    peakHeight = peakHeight1;
	}
	if (peakHeight2 > peakHeight) {
	    peakBin = peakBin2;
	    peakHeight = peakHeight2;
	}
  
	// Check that this peak satisfies the minimum selection requirements.
	// Exclude peaks that are right at the boundaries of the lpc curve, since
	// we cannot look at lpc points on either side to find vertices etc..
	if (peakBin > 0 && peakBin < nLpc1 && peakHeight > cosAngleCut_ &&
	    peakHeight > minPeakFrac_*prevPeakHeight) {
                                            
	    int minBin = peakBin - dBin_;
	    int maxBin = peakBin + dBin_;
	    if (minBin < 1) {minBin = 1;}
	    if (maxBin > nLpc2) {maxBin = nLpc2;}
 
	    LpcPoint lpcPoint = lpcPoints[peakBin];
	    Eigen::VectorXd lpcPos = lpcPoint.getCoords();

	    // Check to see if the lpc point co-ordinate of this peak is
	    // very close to another one. If so, merge their ranges by
	    // taking the minimum and maximum limits of the lpc point integers

	    bool gotSamePeak(false);
	    int samePeakInt(-1);

	    int nPeakPoints = peakCoord.size();
	    for (int iP = 0; iP < nPeakPoints; iP++) {

		LpcPoint lpcPeak = peakCoord[iP];
		Eigen::VectorXd peakPos = lpcPeak.getCoords();
		Eigen::VectorXd diffPos = peakPos - lpcPos;
		double distSq = diffPos.squaredNorm();

		if (distSq < peakDiffSq_) {

		    // Peak points match
		    gotSamePeak = true;
		    samePeakInt = iP;

		    LpcBinRange bin0 = theRanges[iP];
		    int minBin0 = bin0.getMinBin();
		    int maxBin0 = bin0.getMaxBin();
		    // Merge the minimum and maximum limits of the lpc point range
		    // by updating the minBin and maxBin values
		    if (minBin0 < minBin) {minBin = minBin0;}
		    if (maxBin0 > maxBin) {maxBin = maxBin0;}
		    // This peak matches a previous one. Stop this loop
		    break;

		}

	    } // Loop over the current peak bins

	    if (gotSamePeak) {

		// Set the range with the updated min/max lpc point limits
		theRanges[samePeakInt] = LpcBinRange(minBin, maxBin);

	    } else {
		// We have a unique peak
		nPassedPeaks += 1;
		peakCoord.push_back(lpcPoint);

		LpcBinRange peakBinRange(minBin, maxBin);
		theRanges.push_back(peakBinRange);

		passedPeaks.push_back(peakBin);

	    }
	
	    // Reset the previous peak height to this peak height
	    prevPeakHeight = peakHeight;

 
	} // Checks on peak bin/height boundaries

    } // Loop over the sorted peak bins

    // Loop over the bin ranges and check to see if any ranges overlap. 
    // If they do, merge them. Also merge if any range < dBin_.
    if (nPassedPeaks > 0 && theRanges.size() > 0) {

	// Sort the bin ranges in ascending (minimum bin) order
	std::sort(theRanges.begin(), theRanges.end());

	// Define the vector of merged bin ranges
	std::vector<LpcBinRange> mergedRanges;
	// Add the first entry in the vector of sorted bin ranges
	mergedRanges.push_back(theRanges[0]);
	
	for (int iR = 1; iR < nPassedPeaks; iR++) {

	    LpcBinRange nextRange = theRanges[iR];
	    int lowNext = nextRange.getMinBin();
	    int highNext = nextRange.getMaxBin();

	    int nMerge1 = mergedRanges.size() - 1;
	    LpcBinRange mergedBinR = mergedRanges[nMerge1];
	    int lowMerge = mergedBinR.getMinBin();
	    int highMerge = mergedBinR.getMaxBin();

	    if (lowNext > lowMerge && highNext < highMerge) {
		// The nextRange is already within the current merged range
		continue;

	    } else if (lowNext > lowMerge && lowNext < highMerge && highNext > highMerge) {
		// Adjust the high end of the current merged range
		mergedRanges[nMerge1].setMaxBin(highNext);

	    } else if ((lowNext - highMerge) < dBin1_ && highNext > highMerge) {
		// The next range is only dBin - 1 units away from the previous one.
		// Adjust the high end of the current merged range
		mergedRanges[nMerge1].setMaxBin(highNext);

	    } else if (lowNext < lowMerge && highNext < highMerge && highNext > lowMerge) {
		// Adjust the low end of the current merged range
		mergedRanges[nMerge1].setMinBin(lowNext);

	    } else if (lowNext < lowMerge && highNext > highMerge) {
		// Adjust the low and high ends of the current merged range
		mergedRanges[nMerge1] = nextRange;

	    } else {
		// Add new merged range interval
		mergedRanges.push_back(nextRange);

	    }
	    
	} // Loop over passed peaks

	// Store the passed peaks and the final peak ranges to the lpc curve
	theCurve->storeCosPeaks(passedPeaks, mergedRanges);


    } // We have at least one passed peak       

}


void LpcFeatures::printBinRanges(const std::vector<LpcBinRange>& binRanges) const
{

    for (size_t j = 0; j < binRanges.size(); j++) {
	LpcBinRange b = binRanges[j];
	std::cout<<"("<<b.getMinBin()<<", "<<b.getMaxBin()<<"); ";
    }
    std::cout<<std::endl;

}

void LpcFeatures::getCosPeaks(const Eigen::VectorXd& absCosArray)
{

    if (peakFinder_ == 1) {

	// Use the simple peak finder method
	this->getSimplePeaks(absCosArray);

    } else if (peakFinder_ == 2) {

	// Use the ROOT TSpectrum method
	this->getTSpectrumPeaks(absCosArray);

    }	

}

void LpcFeatures::getSimplePeaks(const Eigen::VectorXd& absCosArray) 
{

    // Find the peaks using a simple smoothing method by averaging
    // neighbouring "bins" in the 1-|cosAngle| distribution. Peaks
    // are identified if the 2 bins on either side of the smoothed
    // distribution are lower in height

    allPeaks_.clear();

    int N = absCosArray.size();
    if (N < 1) {return;}

    // Create a smoothed series of absCosValues
    Eigen::VectorXd smoothArray = Eigen::VectorXd::Zero(N);
    double maxVal(0.0);
    int N1 = N - 1;
    for (int iL = 1; iL < N1; iL++) {

	double v0 = absCosArray(iL-1);
	double v1 = absCosArray(iL);
	double v2 = absCosArray(iL+1);

	double smoothedVal = 0.5*(v0 + v2) + v1;
	if (smoothedVal > maxVal) {maxVal = smoothedVal;}
	smoothArray(iL) = smoothedVal;

    }

    // Find the smoothed peaks
    double threshold = 0.05*maxVal;
    int nS = smoothArray.size();

    std::vector<int> peakBins;

    for (int iS = 2; iS < nS-2; iS++) {
	    
	double y0  = smoothArray(iS);
	if (y0 < threshold) {continue;}

	double ym2 = smoothArray(iS-2);
	double ym1 = smoothArray(iS-1);
	double yp1 = smoothArray(iS+1);
	double yp2 = smoothArray(iS+2);

	if (ym1 > ym2 && y0 > ym1 && yp1 < y0 && yp2 < yp1) {
	    // We have a peak
	    peakBins.push_back(iS+1);

	} else if (fabs(ym2) < threshold && fabs(yp2) < threshold) {

	    // If the 2nd bins on either side are very small, we 
	    // may have a narrow peak
	    if (y0 > ym1 && yp1 < y0) {
		peakBins.push_back(iS+1);
	    }

	}

    }

    // Sometimes, the smoothing will move the peak by +-1 bin. 
    // Try to correct for this
    int nPeaks = peakBins.size();

    std::vector<int> corrPeakBins;
    for (int iP = 0; iP < nPeaks; iP++) {

	int pBin = peakBins[iP];

	int corrBin(pBin);
	double yValue = absCosArray(pBin);
	double thePeak(yValue);

	if (pBin > 1 && pBin < N1) {

	    int pBin1 = pBin - 1;
	    int pBin2 = pBin + 1;
	    double yp0 = absCosArray(pBin1);
	    double yp2 = absCosArray(pBin2);

	    if (yp0 > yValue) {
		corrBin = pBin1;
		thePeak = yp0;
	    } else if (yp2 > yValue) {
		corrBin = pBin2;
		thePeak = yp2;
	    }
	    
	}

	// Store the lpc peak
	LpcPeakInfo peakInfo;
	peakInfo.bin_ = corrBin;
	peakInfo.height_ = thePeak;
	allPeaks_.push_back(peakInfo);

    }
    

}

void LpcFeatures::getTSpectrumPeaks(const Eigen::VectorXd& absCosArray)
{

    // Use the ROOT TSpectrum class for finding the peaks
    // First, clear the vector of peak positions
    allPeaks_.clear();

    if (absCosArray.size() < 1) {return;}

    // Check that the LPC_ROOT variable is set; otherwise the 
    // rest of the function is ignored

#ifdef LPC_USE_ROOT

    TSpectrum spectrum(20, 10);
    
    float sigma(1.0);
    double threshold(0.50);
    bool bkgndRemove(false), markov(true);
    int deconIter(5), avWindow(1);

    int nLpc = absCosArray.size();

    // Unfortunately, TSpectrum requires float "arrays"...
    float* source = new float[nLpc];
    float* destVector = new float[nLpc];
    int i(0);
    for (i = 0; i < nLpc; i++) {
	source[i] = absCosArray(i);
	destVector[i] = 0.0;
    }

    int nPeaks = spectrum.SearchHighRes(source, destVector, nLpc, sigma, threshold, 
					bkgndRemove, deconIter, markov, avWindow);

    float* xPeaks = spectrum.GetPositionX();

    // Convert xPeaks to the integer values of the lpc point indices
    for (i = 0; i < nPeaks; i++) {

	LpcPeakInfo peakInfo;
	peakInfo.bin_ = int(xPeaks[i] + 0.5);
	peakInfo.height_ = absCosArray(peakInfo.bin_);
	allPeaks_.push_back(peakInfo);

    }

    // Delete the float arrays
    delete [] source;
    delete [] destVector;

#endif

}

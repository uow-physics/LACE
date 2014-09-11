// Copyright University of Warwick 2014
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// Authors:
// John Back

/*! \file LpcLineClusterAlgorithm.hh
    \brief File containing the declaration of the LpcLineClusterAlgorithm class
*/

/*! \class LpcLineClusterAlgorithm
    \brief Class that defines the (extended) line vertexing algorithm
*/

#ifndef LPC_LINE_CLUSTER_ALGORITHM_HH
#define LPC_LINE_CLUSTER_ALGORITHM_HH

#include "LACE/LpcAbsClusterAlgorithm.hh"

#include "LACE/LpcBinRange.hh"
#include "LACE/LpcClusterData.hh"
#include "LACE/LpcFunctions.hh"
#include "LACE/LpcLineFitter.hh"
#include "LACE/LpcLineSegment.hh"
#include "LACE/LpcPoint.hh"
#include "LACE/LpcResiduals.hh"

#include <Eigen/Dense>
#include <map>
#include <string>
#include <vector>

class LpcAbsCurve;
class LpcCurve;
class LpcCluster;
class LpcHitCollection;
class LpcParameters;
class LpcVertex;

class LpcLineClusterAlgorithm : public LpcAbsClusterAlgorithm {

public:

    //! Constructor
    /*!
      \param [in] thePars a pointer to the constant lpc parameters
    */
    LpcLineClusterAlgorithm(const LpcParameters* thePars);

    //! Destructor
    virtual ~LpcLineClusterAlgorithm();

    //! Process the main lpc curve (and its branches), finding vertices and clusters
    /*!
      \param [in] theCurve The pointer to the constant main lpc curve (read-only access)
      \returns an LpcClusterData object, which stores vectors of vertex and cluster pointers
    */
    virtual LpcClusterData findClusters(const LpcCurve* theCurve);

protected:

    //! Various functions and methods
    
    //! Typedef for a map of matched hits for each line segment
    typedef std::map<int, std::vector<int> > LpcLineIntMap;

    //! Process either the main curve or one of its branches
    /*!
      \param [in] theCurve The pointer to the (read-only) main lpc or one of its branches
      \returns the paired vector<LpcVertex*>, vector<LpcCluster*> for the curve
    */
    LpcClusterData processCurve(const LpcAbsCurve* theCurve);

    //! Find the lpc ranges that straddle the lpc feature points.
    //! This function fills the lpcRanges_ vector
    void findLpcRanges();

    //! Find the ranges of the line segments that are near the feature points.
    //! This function fills the lineRanges_ vector
    void findLineRanges();

    //! Find the line segments that are near the feature points.
    //! This function fills the lineSegments_ vector
    void findLineSegments();

    //! Find the vertex points for the given line segments
    /*!
      \param [in] lineSegments The vector of line segments
      \returns an Eigen MatrixXd of the vertices (row = vertex, col = co-ordinates)
    */
    Eigen::MatrixXd findVtxPoints(const std::vector<LpcLineSegment>& lineSegments) const;

    //! Find the extended line segments that are near the feature points.
    //! This function fills the extLineSegments_ vector
    /*!
      \param [in] initVtxPoints A MatrixXd of initial vertices (row = vtx, col = co-ords)
    */
    void findExtendedLines(const Eigen::MatrixXd& initVtxPoints);

    //! Find the hits associated to the initial line segments
    /*!
      \param [in] lineSegments Vector of LpcLineSegments that we want to find the associated hits for
      \param [in] uniqueHits Boolean to specify if any previously used hit is to be ignored
      \returns the LpcLineIntMap dictionary, key = line segment number, value = vector of hit ids
    */
    LpcLineIntMap findAssociatedHits(const std::vector<LpcLineSegment>& lineSegments, 
				     bool uniqueHits = false);

    //! Create clusters using the associated hits, line segments and vertices
    /*!
      \param [in] vtxMatrix An Eigen MatrixXd of the vertices (row = vtx, col = coords)
      \returns a vector of LpcCluster pointers
    */
    std::vector<LpcCluster*> formClusters(const Eigen::MatrixXd& vtxMatrix);

    //! Create a complete list of cluster hits (indices) using extended lines, lpcRange and residual cut
    /*!
      \param [in] matchedHits The vector of hit indices to be considered for the cluster
      \param [in] lpcRange The range of the lpc points to be considered for the cluster
      \param [in] residualCut The maximum allowed residual for a hit to be part of the cluster
    */
    void formClusterHits(std::vector<int>& matchedHits,
			 const LpcBinRange& lpcRange,
			 double residualCut);

    //! Create a new LpcCluster pointer from the required hits
    /*!
      \param [in] index The index integer number of the cluster
      \param [in] matchedHits The vector of hit indices to be used for the cluster
      \param [in] lpcRange The range of lpc points used for the cluster
      \returns a new LpcCluster pointer
    */
    LpcCluster* createCluster(int index, const std::vector<int>& matchedHits, 
			      const LpcBinRange& lpcRange) const;

    //! Create a new single LpcCluster pointer from all available hits
    /*!
      \param [in] index The index integer number of the cluster
      \returns a new LpcCluster pointer
    */
    LpcCluster* createSingleCluster(int index) const;

    //! Create a new single LpcVertex pointer as the position of the hit with the lowest 1st co-ord ("x")
    /*!
      \returns a new LpcVertex pointer
    */
    LpcVertex* createSingleVertex() const;

    //! Add any remaining hits to the clusters by finding the closest 
    //! distance to the cluster principal axis
    /*!
      \param [in] lpcClusters The vector of LpcCluster pointers
    */
    void addRemainingHits(const std::vector<LpcCluster*>& lpcClusters);

    //! Create a vector of LpcVertex pointers containing the vertex information
    /*!
      \param [in] vtxPointMatrix The matrix of the vertex co-ordinates (row = vtx, col = x,y,z..)
      \return a vector of new LpcVertex pointers
    */
    std::vector<LpcVertex*> createVertices(const Eigen::MatrixXd& vtxPointMatrix) const;

    //! Merge clusters
    /*!
      \param [in] vtxPointMatrix The matrix of the vertex co-ordinates (row = vtx, col = x,y,z..)
      \param [in] lpcClusters The vector of original LpcCluster pointers
      \returns the pair of merged vertex and (a copy of the) cluster vectors
    */
    LpcClusterData mergeClusters(const Eigen::MatrixXd& vtxPointMatrix, 
				 const std::vector<LpcCluster*>& lpcClusters);


    //! Check and potentially correct the allowed lpc index number
    /*!
      \param [in] lpcInt The initial lpc index point number
      \returns the index point number within the range [0, nLpcPoints1_]
    */
    int checkLpcInt(int lpcInt) const;

    //! Get the residual cut given the residual mean and rms
    /*!
      \param [in] residualMean The mean of the hit-to-lpc residuals
      \param [in] residualRms The rms of the hit-to-lpc residuals
      \returns the residual cut
    */
    double getResidualCut(double residualMean, double residualRms) const;

    //! Print out a vector of LpcBinRange objects
    /*!
      \param [in] preamble A starting string for the printout line
      \param [in] binRanges The given vector of LpcBinRange objects
    */
    void printBinRanges(const std::string& preamble,
			const std::vector<LpcBinRange>& binRanges) const;


private:

    //! Private default constructor
    LpcLineClusterAlgorithm();

    //! Copy constructor
    LpcLineClusterAlgorithm(const LpcLineClusterAlgorithm& other);

    //! Assignment operator
    LpcLineClusterAlgorithm& operator=(const LpcLineClusterAlgorithm& other);
  
    //! The pointer to the constant hit collection
    const LpcHitCollection* theHits_;

    //! The line fitter used to find the regression line for the set of hits
    LpcLineFitter lineFitter_;

    //! Number of lpc points to consider for each line segment
    int nLpcSegment_;

    //! Number of sigma for residual cut
    double nResSigma_;

    //! The minimum selection cut for lpc-to-hit residuals for the vertexing
    double minVtxResCut_;

    //! Maximum distance to associate any remaining hits with a cluster
    double maxDist_;

    //! Specify if we want to merge neighbouring clusters
    bool doClusterMerging_;

    //! Specify if we want to do vertexing and clustering for branches
    int doBranchVtx_;

    //! Specify the minimum angle for merging clusters
    double mergeDegAngle_;

    //! The cosine of the merging cluster angle
    double mergeCosAngle_;

    //! The LpcFunctions object
    LpcFunctions functions_;

    //! The curve index number
    int curveId_;

    //! The branch index number
    int branchId_;

    //! Hit and lpc point information

    //! The number of hit points 
    int nHitPoints_;

    //! The number of co-ordinate dimensions
    int nDim_;

    //! The number of lpc points
    int nLpcPoints_;

    //! The number of lpc points - 1
    int nLpcPoints1_;

    //! The Eigen matrix of unscaled hit positions
    Eigen::MatrixXd hitPoints_;

    //! The vector of LpcPoints
    std::vector<LpcPoint> lpcPoints_;

    //! The lpc residual object
    LpcResiduals lpcResiduals_;

    //! The hit-to-lpc residuals
    Eigen::VectorXd hitResiduals_;

    //! The closest lpc point integer for each hit
    Eigen::VectorXi hitNearLpc_;

    //! An array to keep track of which hits have already been used for clusters
    Eigen::VectorXi usedHits_;

    //! Internal vectors of various lpc point bin ranges

    //! The number of cosine peaks
    int nPeaks_;

    //! The vector of lpc (min,max) point ranges for each cosine feature peak
    std::vector<LpcBinRange> cosPeakRanges_;

    //! The vector of lpc (min,max) point ranges that straddle the feature points
    std::vector<LpcBinRange> lpcRanges_;

    //! The vector of lpc (min,max) point ranges that are near the feature points
    std::vector<LpcBinRange> lineRanges_;

    //! The vector of line segments
    std::vector<LpcLineSegment> lineSegments_;

    //! The vector of extended line segments
    std::vector<LpcLineSegment> extLineSegments_;

    //! An internal boolean to specify if cluster merging was successful
    bool mergedOK_;

};

#endif

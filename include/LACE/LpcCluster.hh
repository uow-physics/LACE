// Copyright University of Warwick 2014
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// Authors:
// John Back

/*! \file LpcCluster.hh
    \brief File containing the declaration of the LpcCluster class
*/

/*! \class LpcCluster
    \brief Class that defines an lpc cluster object
*/

#ifndef LPC_CLUSTER_HH
#define LPC_CLUSTER_HH

#include "LACE/LpcBinRange.hh"

#include <Eigen/Dense>
#include <vector>

class LpcHit;

class LpcCluster {

public:

    //! Empty constructor
    LpcCluster();

    //! Constructor
    /*!
      \param [in] index The index number of the cluster
      \param [in] nDim The number of co-ordinate dimensions
      \param [in] theHits The vector of the LpcHit pointers that form the cluster
      \param [in] hitResiduals The vector of hit-to-lpc residuals for each cluster hit
      \param [in] curveId The main curve index number for the cluster
      \param [in] branchId The branch curve index number for the cluster (0 = main curve)
      \param [in] lpcRange The lpc point (bin) range used for the cluster
      \param [in] isShower Specify if the cluster is a shower or track (false, by default)
    */
    LpcCluster(int index, int nDim, 
	       const std::vector<LpcHit*>& theHits,
	       const std::vector<double>& hitResiduals,
	       int curveId, int branchId, 
	       const LpcBinRange& lpcRange,
	       bool isShower = false);
    
    //! Copy constructor
    /*!
      \param [in] other The LpcCluster to copy
    */
    LpcCluster(const LpcCluster& other);

    //! Assignment operator
    /*!
      \param [in] other The LpcCluster used for assignment
    */
    LpcCluster& operator=(const LpcCluster& other);

    //! Destructor
    virtual ~LpcCluster();

    // Modifiers

    //! Add another LpcHit. This sets the doUpdate_ flag to true
    /*!
      \param [in] theHit The hit to add
      \param [in] hitResidual The hit-to-lpc residual for the hit
    */
    void addHit(LpcHit* theHit, double hitResidual);

    //! Set the shower boolean
    /*!
      \param [in] result The true/false value for the shower boolean
    */
    void setShowerBool(bool result) {isShower_ = result;}

    //! Change the index number
    /*!
      \param [in] index The index number of the cluster
    */
    void setIndex(int index) {index_ = index;}

    //! Set the convex hull lengths
    /*!
      \param [in] convexHull The vector of convex hull lengths
    */
    void storeConvexHull(const Eigen::VectorXd& convexHull) {convexHull_ = convexHull;}

    // Accessors

    //! Get the index number of the cluster
    /*!
      \return the index integer
    */
    int getIndex() const {return index_;}

    //! Get the number of co-ordinate dimensions
    /*!
      \return the number of dimensions
    */
    int getNDimensions() const {return nDim_;}

    //! Get the list of LpcHits used to form the cluster
    /*!
      \return the vector of LpcHit pointers used for the cluster
    */
    std::vector<LpcHit*> getHits() const {return theHits_;}

    //! Get the number of hits in the cluster
    /*!
      \return the number of hits in the cluster
    */
    int getNumberOfHits() const {return theHits_.size();}

    //! Get the list of LpcHit indices used to form the cluster
    /*!
      \return the vector of LpcHit index integers
    */
    std::vector<int> getHitIndices() const;

    //! Get the weights of the LpcHits used to form the cluster
    /*!
      \return the vector of weights of the hits in the cluster
    */
    std::vector<double> getHitWeights() const;

    //! Get the hit-to-lpc residuals for hits used in the cluster
    /*!
      \return the vector of residuals of hits used in the cluster
    */
    std::vector<double> getHitResiduals() const {return hitResiduals_;}

    //! Get the lpc point range
    /*!
      \return the lpc point range as an LpcBinRange object
    */
    LpcBinRange getLpcRange() const {return lpcRange_;}

    //! Get the main curve index number used for the cluster
    /*!
      \return the main curve index integer
    */
    int getCurveId() const {return curveId_;}

    //! Get the branch curve index number used for the cluster
    /*!
      \return the branch curve index integer
    */
    int getBranchId() const {return branchId_;}

    //! Is the cluster a shower or track?
    /*!
      \return the shower identification boolean
    */
    bool isAShower() const {return isShower_;}

    //! Get the cluster centroid
    /*!
      \return the cluster centroid as an Eigen VectorXd
    */
    Eigen::VectorXd getCentroid();

    //! Get the principle component axes (normalised eigenvectors)
    /*!
      \return the principal component axes as an Eigen MatrixXd (row = axis, col = co-ords)
    */
    Eigen::MatrixXd getPCAxes();

    //! Get the main principle axis (normalised eigenvector)
    /*!
      \return the main principal axes as an Eigen VectorXd
    */
    Eigen::VectorXd getPrincipalAxis();

    //! Get the lengths along the principal component axes (eigenvalues)
    /*!
      \return the lengths along the principal component axes
    */
    Eigen::VectorXd getPCLengths();

    //! Get the convex hull lengths
    /*!
      \return the vector of convex hull lengths
    */
    Eigen::VectorXd getConvexHull() const {return convexHull_;}


    //! Print out the cluster info
    void print();

    //! Form an Eigen MatrixXd of the hit co-ordinates contained in the cluster
    /*!
      \return an Eigen MatrixXd of the hit co-ordinates (row = hit, col = x,y,z,..)
    */
    Eigen::MatrixXd getHitPositions() const;

protected:

    //! Update the centroid and principal component axes if doUpdate_ is true
    void update();

    //! Calculate the cluster centroid
    /*!
      \param [in] hitPositions The Eigen MatrixXd of the hit positions (row = hit, col = x,y,z..)
    */
    void calcCentroid(const Eigen::MatrixXd& hitPositions);

    //! Calculate the cluster principal component axes (directions and lengths)
    /*!
      \param [in] hitPositions The Eigen MatrixXd of the hit positions (row = hit, col = x,y,z..)
    */
    void calcPCAxes(const Eigen::MatrixXd& hitPositions);

private:

    //! The cluster index integer
    int index_;

    //! The number of co-ordinate dimensions
    int nDim_;

    //! The vector of LpcHit pointers, which are owned by 
    //! the LpcHitCollection in LpcEvent
    std::vector<LpcHit*> theHits_;

    //! The vector of hit-to-lpc residuals; same ordering as theHits_
    std::vector<double> hitResiduals_;

    //! The main curve index integer
    int curveId_;

    //! The branch curve index integer
    int branchId_;

    //! The lpc point range
    LpcBinRange lpcRange_;

    //! The shower boolean
    bool isShower_;

    //! Boolean to specify if we need to call update again if
    //! retrieving the cluster centroid and principal components
    bool doUpdate_;

    //! The centroid of the cluster
    Eigen::VectorXd centroid_;

    //! The principal component axes for the cluster (row = axis, col = coords)
    Eigen::MatrixXd PCAxes_;

    //! The principal component lengths for the cluster
    Eigen::VectorXd PCLengths_;

    //! The vector of convex hull lengths
    Eigen::VectorXd convexHull_;

};

#endif

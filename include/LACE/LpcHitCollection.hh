// Copyright University of Warwick 2014
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// Authors:
// John Back

/*! \file LpcHitCollection.hh
    \brief File containing the declaration of the LpcHitCollection class
*/

/*! \class LpcHitCollection
    \brief Class used to define a collection of hit pointers
*/

#ifndef LPC_HIT_COLLECTION_HH
#define LPC_HIT_COLLECTION_HH

#include "LACE/LpcPoint.hh"

#include <Eigen/Dense>
#include <vector>

class LpcHit;

class LpcHitCollection {

public:

    //! Empty constructor containing no hits
    /*!
      \param [in] nDim The number of co-ordinate dimensions (default = 3)
    */
    LpcHitCollection(int nDim = 3);

    //! Constructor using a vector of hit pointers (which it owns)
    /*!
      \param [in] theHits The vector of hit pointers
      \param [in] nDim The number of co-ordinate dimensions (default = 3)
    */
    LpcHitCollection(const std::vector<LpcHit*>& theHits, int nDim = 3);

    //! Destructor
    virtual ~LpcHitCollection();

    //! Add another hit to the end of the collection
    /*!
      \param [in] aHit additional hit to be included at the end of the collection
    */
    void addHit(LpcHit* aHit) {theHits_.push_back(aHit);}

    //! Process the hit collection to create Xi_, scaledXi_, weights_ and x0_ (start point)
    /*!
      \param [in] applyScaling Bool to decide if scaled co-ordinates are used for the Lpc algorithm
    */
    void process(bool applyScaling = true);

    //! Reset the starting point
    /*!
      \param [in] newStartPoint the new starting point represented as an LpcPoint
    */
    void setStartPoint(const LpcPoint& newStartPoint) {startPoint_ = newStartPoint;}

    // Accessor functions

    //! The number of co-ordinate dimensions stored in each hit
    /*!
      \return the number of co-ordinate dimensions
    */
    int getNDimensions() const {return nDim_;}

    //! The vector of hit pointers
    /*!
      \return a stl vector of hit pointers
    */
    std::vector<LpcHit*> getHits() const {return theHits_;}

    //! Retrieve a specific hit in the collection, given the index number
    /*!
      \param [in] hitIndex The integer index number of the hit we need
      \returns the pointer to the LpcHit
    */
    LpcHit* getHit(int hitIndex) const;

    //! The number of hits
    /*!
      \return the number of hits
    */
    int getNumberOfHits() const {return theHits_.size();}

    //! Return the unscaled hit co-ordinates as an Eigen matrix (row = hit, col = x,y,z)
    /*!
      \returns the constant reference to the Eigen matrix of the unscaled co-ordinates
    */
    const Eigen::MatrixXd& getCoords() const {return Xi_;}

    //! Return the scaled hit co-ordinates as an Eigen matrix (row = hit, col = x,y,z)
    /*!
      \returns the constant reference to the Eigen matrix of the scaled co-ordinates
    */
    const Eigen::MatrixXd& getScaledCoords() const {return scaledXi_;}

    //! Return the hit weights as an Eigen vector
    /*!
      \returns the constant reference to the Eigen vector of the hit weights
    */
    const Eigen::VectorXd& getWeights() const {return weights_;}

    //! Retrieve the starting point for the lpc algorithm
    /*!
      \returns an LpcPoint of the starting point
    */
    LpcPoint getStartPoint() const {return startPoint_;}

    //! Retrieve the co-ordinate range that was used to scale the hit co-ordinates
    /*!
      \returns an Eigen vector of the scaling range
    */
    Eigen::VectorXd getRange() const {return theRange_;}

    //! Retrieve the scaled starting point as an Eigen vector
    /*!
      \returns the scaled starting point as an Eigen vector
    */
    Eigen::VectorXd getx0() const {return startPoint_.getScaledCoords();}


protected:

private:

    //! Private copy constructor
    LpcHitCollection(const LpcHitCollection& other);

    //! Private assignment operator
    LpcHitCollection& operator=(const LpcHitCollection& other);
    
    //! Store the hit co-ordinates and weights in Xi_ and weights_
    void storeCoordWeights();

    //! Find the co-ordinate range for scaling
    void findCoordRange();

    //! Select the starting point for the lpc algorithm
    void selectStartPoint();

    //! Store the scaled co-ordinates in scaledXi_
    void scaleCoords();

    //! The number of co-ordinate dimensions
    int nDim_;

    //! The vector containing all of the hits; this class owns the hit pointers.
    std::vector<LpcHit*> theHits_;

    //! An Eigen Matrix storing the original, unscaled hit co-ordinates
    Eigen::MatrixXd Xi_;

    //! An Eigen Matrix storing the scaled hit co-ordinates
    Eigen::MatrixXd scaledXi_;

    //! An Eigen Vector storing the hit weights
    Eigen::VectorXd weights_;

    //! The vector of co-ordinate ranges
    Eigen::VectorXd theRange_;

    //! The starting point for the lpc algorithm (scaled & unscaled co-ordinates)
    LpcPoint startPoint_;

};

#endif

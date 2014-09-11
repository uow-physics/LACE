// Copyright University of Warwick 2014
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// Authors:
// John Back

/*! \file LpcEvent.hh
    \brief File containing the declaration of the LpcEvent class
*/

/*! \class LpcEvent
    \brief Class used to define an event containing hits, curves, clusters and vertices
*/

#ifndef LPC_EVENT_HH
#define LPC_EVENT_HH

#include <vector>

class LpcCluster;
class LpcCurve;
class LpcHitCollection;
class LpcVertex;

class LpcEvent {

public:

    //! Default empty constructor containing no hits
    /*!
      \param [in] eventNumber the index number labelling the event
    */
    LpcEvent(int eventNumber);

    //! Constructor using a vector of hit pointers
    /*!
      \param [in] eventNumber the index number labelling the event
      \param [in] theHits The pointer to the collection of hits
    */
    LpcEvent(int eventNumber, LpcHitCollection* theHits);
    
    //! Destructor
    virtual ~LpcEvent();

    //! Set-up the event so that the hit collection contains scaled co-ords etc.
    /*!
      \param [in] applyScaling Bool to decide if scaled co-ordinates are used for the Lpc algorithm
    */
    void setUp(bool applyScaling = true);


    //! Store (and own) a local principal curve pointer
    /*!
      \param [in] aCurve a local principal curve pointer
    */
    void storeCurve(LpcCurve* aCurve) {theCurves_.push_back(aCurve);}

    //! Store and own the LpcCluster pointers stored in a vector
    /*!
      \param [in] theClusters The vector of LpcCluster pointers
    */
    void storeClusters(const std::vector<LpcCluster*>& theClusters);

    //! Store and own the LpcVertex pointers stored in a vector
    /*!
      \param [in] theVertices The vector of LpcVertex pointers
    */
    void storeVertices(const std::vector<LpcVertex*>& theVertices);

    //! Store (and own) an LpcVertex pointer
    /*!
      \param [in] theVertex the pointer to the LpcVertex
    */
    void storeVertex(LpcVertex* theVertex) {theVertices_.push_back(theVertex);}

    //! Store (and own) an LpcCluster pointer
    /*!
      \param [in] theCluster the pointer to the LpcCluster
    */
    void storeCluster(LpcCluster* theCluster) {theClusters_.push_back(theCluster);}

    // Accessor functions

    //! The hit collection, containing hit co-ordinates and weights
    /*!
      \return the hit collection
    */
    LpcHitCollection* getHitCollection() const {return theHits_;}

    //! The index number for the event
    /*!
      \return the index number labelling the event
    */
    int getEventNumber() const {return eventNumber_;}

    //! Get the list of local principal curves
    /*!
      \returns a vector of local principal curves
    */
    std::vector<LpcCurve*> getCurves() const {return theCurves_;}

    //! Get the number of local principal curves
    /*!
      \returns the number of local principal curves
    */
    int getNLpcCurves() const {return theCurves_.size();}

    //! Return the first local principal curve
    /*!
      \returns the first local principal curve
    */
    LpcCurve* getCurve() const;

    //! Get the list of vertices
    /*!
      \returns a vector of vertices
    */
    std::vector<LpcVertex*> getVertices() const {return theVertices_;}

    //! Get the number of vertices
    /*!
      \returns the number of vertices
    */
    int getNVertices() const {return theVertices_.size();}

    //! Get the list of clusters
    /*!
      \returns a vector of clusters
    */
    std::vector<LpcCluster*> getClusters() const {return theClusters_;}

    //! Get the number of clusters
    /*!
      \returns the number of clusters
    */
    int getNClusters() const {return theClusters_.size();}

protected:

private:

    // Private copy constructor and assignment operator
    
    //! Copy constructor
    LpcEvent(const LpcEvent& other);

    //! Assignment operator
    LpcEvent& operator=(const LpcEvent& other);

    //! The event number index label
    int eventNumber_;

    //! The hit collection
    LpcHitCollection* theHits_;

    //! Vector storing the local principal curve pointers
    std::vector<LpcCurve*> theCurves_;

    //! Vector storing the vertices
    std::vector<LpcVertex*> theVertices_;

    //! Vector storing the clusters
    std::vector<LpcCluster*> theClusters_;

};

#endif

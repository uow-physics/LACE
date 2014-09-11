// Copyright University of Warwick 2014
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// Authors:
// John Back

/*! \file LpcClusterData.hh
    \brief File containing the declaration of the LpcClusterData class
*/

/*! \class LpcClusterData
    \brief Class that is a container of the LpcCluster and LpcVertex pointers from clustering
*/

#ifndef LPC_CLUSTER_DATA_HH
#define LPC_CLUSTER_DATA_HH

#include <vector>

class LpcCluster;
class LpcVertex;

class LpcClusterData {

public:

    //! Empty constructor
    LpcClusterData();

    //! Constructor with vectors of LpcVertex and LpcCluster pointers
    /*!
      \param [in] vertices The vector of LpcVertex pointers
      \param [in] clusters The vector of LpcCluster pointers
    */
    LpcClusterData(const std::vector<LpcVertex*>& vertices,
		   const std::vector<LpcCluster*>& clusters);

    //! Copy constructor
    /*!
      \param [in] other The LpcClusterData object to copy
    */
    LpcClusterData(const LpcClusterData& other);

    //! Assignment operator
    /*!
      \param [in] other The LpcClusterData to assign to
    */
    LpcClusterData& operator = (const LpcClusterData& other);

    //! Destructor
    virtual ~LpcClusterData();

    //! Add another set of LpcClusterData to the present object
    /*!
      \param [in] data More LpcClusterData to add
    */
    void addData(const LpcClusterData& data);

    //! Add additional vertices
    /*!
      \param [in] vertices The vector of LpcVertex pointers to add
    */
    void addVertices(const std::vector<LpcVertex*>& vertices);

    //! Add additional clusters
    /*!
      \param [in] clusters The vector of LpcCluster pointers to add
    */
    void addClusters(const std::vector<LpcCluster*>& clusters);

    //! Add an additional vertex
    /*!
      \param [in] theVertex The pointer to the LpcVertex to add
    */
    void addVertex(LpcVertex* theVertex) {theVertices_.push_back(theVertex);}

    //! Add an additional cluster
    /*!
      \param [in] theCluster The pointer to the LpcCluster to add
    */
    void addCluster(LpcCluster* theCluster) {theClusters_.push_back(theCluster);}

    //! Reset the indices of the vertices and clusters to match their vector ordering
    void resetIndices();

    // Accessors
    
    //! Retrieve the vector of LpcVertex pointers
    /*!
      \return the vector of LpcVertex pointers
    */
    std::vector<LpcVertex*> getVertices() const {return theVertices_;}

    //! Retrieve the vector of LpcCluster pointers
    /*!
      \return the vector of LpcCluster pointers
    */
    std::vector<LpcCluster*> getClusters() const {return theClusters_;}

    //! Get the number of stored vertices
    /*!
      \returns the number of vertices
    */
    int getNumberVertices() const {return theVertices_.size();}

    //! Get the number of stored clusters
    /*!
      \returns the number of clusters
    */
    int getNumberClusters() const {return theClusters_.size();}

protected:

private:

    //! The vector of LpcVertex pointers
    std::vector<LpcVertex*> theVertices_;

    //! The vector of LpcCluster pointers
    std::vector<LpcCluster*> theClusters_;

};

#endif


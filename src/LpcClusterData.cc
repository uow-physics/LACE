// Copyright University of Warwick 2014
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// Authors:
// John Back

/*! \file LpcClusterData.cc
    \brief Class that is a container of the LpcCluster and LpcVertex pointers from clustering
*/

#include "LACE/LpcClusterData.hh"

#include "LACE/LpcCluster.hh"
#include "LACE/LpcVertex.hh"

LpcClusterData::LpcClusterData() :
    theVertices_(),
    theClusters_()
{
}

LpcClusterData::LpcClusterData(const std::vector<LpcVertex*>& vertices,
			       const std::vector<LpcCluster*>& clusters) :
    theVertices_(vertices),
    theClusters_(clusters)
{
}

LpcClusterData::LpcClusterData(const LpcClusterData& other) :
    theVertices_(other.theVertices_),
    theClusters_(other.theClusters_)
{
}

LpcClusterData& LpcClusterData::operator=(const LpcClusterData& other)
{
    theVertices_ = std::vector<LpcVertex*>(other.theVertices_);
    theClusters_ = std::vector<LpcCluster*>(other.theClusters_);
    return *this;
}

LpcClusterData::~LpcClusterData()
{
    // This class does not own the pointers; it's just a convenient
    // way to store the pointer information for vertices and clusters.
}

void LpcClusterData::addData(const LpcClusterData& data)
{

    std::vector<LpcVertex*> moreVertices = data.getVertices();
    std::vector<LpcCluster*> moreClusters = data.getClusters();

    // Add them to the already existing vectors
    this->addVertices(moreVertices);
    this->addClusters(moreClusters);

}

void LpcClusterData::addVertices(const std::vector<LpcVertex*>& vertices)
{    
    theVertices_.insert(theVertices_.end(), vertices.begin(), vertices.end());
}

void LpcClusterData::addClusters(const std::vector<LpcCluster*>& clusters)
{    
    theClusters_.insert(theClusters_.end(), clusters.begin(), clusters.end());
}

void LpcClusterData::resetIndices()
{

    int vIndex(0);
    std::vector<LpcVertex*>::iterator vIter;    
    for (vIter = theVertices_.begin(); vIter != theVertices_.end(); ++vIter) {

	LpcVertex* theVertex = *vIter;
	if (theVertex) {
	    theVertex->setIndex(vIndex);
	    vIndex += 1;
	}

    }

    int cIndex(0);
    std::vector<LpcCluster*>::iterator cIter;    
    for (cIter = theClusters_.begin(); cIter != theClusters_.end(); ++cIter) {

	LpcCluster* theCluster = *cIter;
	if (theCluster) {
	    theCluster->setIndex(cIndex);
	    cIndex += 1;
	}

    }

}

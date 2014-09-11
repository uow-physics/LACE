// Copyright University of Warwick 2014
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// Authors:
// John Back

/*! \file LpcBranchCollection.cc
    \brief Class that stores a series of branches for the main principal curve
*/

#include "LACE/LpcBranchCollection.hh"

#include "LACE/LpcBranch.hh"

#include <iostream>

LpcBranchCollection::LpcBranchCollection() :
    theMap_()
{
}

LpcBranchCollection::~LpcBranchCollection()
{
    // Delete the LpcBranch pointers stored in the map
    LpcBranchMap::iterator iter;
    for (iter = theMap_.begin(); iter != theMap_.end(); ++iter) {

	std::vector<LpcBranch*> theBranches = iter->second;

	// Iterate over the LpcBranch pointers in this vector and delete them
	std::vector<LpcBranch*>::iterator bIter;

	for (bIter = theBranches.begin(); bIter != theBranches.end(); ++bIter) {

	    delete (*bIter);

	}

    }

}

void LpcBranchCollection::addBranch(int generation, LpcBranch* theBranch)
{

    if (theMap_.find(generation) == theMap_.end()) {

	// The map does not contain the generation number. Create a STL
	// vector and store the branch in the map with the generation index number
	std::vector<LpcBranch*> theBranches;
	theBranches.push_back(theBranch);

	theMap_[generation] = theBranches;

    } else {

	// Append the branch pointer to the vector for the given map entry
	theMap_[generation].push_back(theBranch);

	std::vector<LpcBranch*> theBranches = theMap_[generation];

    }

}

void LpcBranchCollection::addBranches(int generation, std::vector<LpcBranch*>& theBranches)
{

    std::vector<LpcBranch*>::iterator iter;
    for (iter = theBranches.begin(); iter != theBranches.end(); ++iter) {

	this->addBranch(generation, *iter);
    }

}

int LpcBranchCollection::getNumberLevels() const
{

    int nLevels(0);
    // Check that the map is not empty
    if (theMap_.size() > 0) {
	nLevels = theMap_.rbegin()->first;
    }

    return nLevels;

}

std::vector<LpcBranch*> LpcBranchCollection::getBranches(int generation)
{

    std::vector<LpcBranch*> theBranches;

    if (theMap_.find(generation) != theMap_.end()) {

	theBranches = theMap_[generation];

    }

    return theBranches;
}

void LpcBranchCollection::print()
{
    int nLevels = this->getNumberLevels();

    std::cout<<"Printing the branch collection; nLevels = "<<nLevels<<std::endl;

    for (int i = 0; i < nLevels; i++) {

	int iL = i + 1;

	// Get the branches for this level
	std::vector<LpcBranch*> theBranches = this->getBranches(iL);
	int nB = theBranches.size();
	
	std::cout<<"Level "<<iL<<" has "<<nB<<" branches"<<std::endl;
	
	for (int j = 0; j < nB; j++) {
	    LpcBranch* aBranch = theBranches[j];
	    if (aBranch) {
		std::cout<<"Printing branch "<<aBranch->getIndex()<<std::endl;
		aBranch->print();
	    }
	}
	
    }
    
}

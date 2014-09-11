// Copyright University of Warwick 2014
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// Authors:
// John Back

/*! \file LpcBranchCollection.hh
    \brief File containing the declaration of the LpcBranchCollection class
*/

/*! \class LpcBranchCollection
    \brief Class that stores a series of branches for the main principal curve
*/

#ifndef LPC_BRANCH_COLLECTION_HH
#define LPC_BRANCH_COLLECTION_HH

#include <map>
#include <vector>

class LpcBranch;

class LpcBranchCollection {

public:

    //! Default, empty constructor
    LpcBranchCollection();

    //! Destructor
    virtual ~LpcBranchCollection();    

    //! Typedef for the map storing a vector of branches for the given genealogy
    typedef std::map<int, std::vector<LpcBranch*> > LpcBranchMap;
  
    //! Add a branch with the given genealogy
    /*!
      \param [in] level The level or generation index number of the branch (starts at 1)
      \param [in] theBranch The LpcBranch pointer
    */
    void addBranch(int level, LpcBranch* theBranch);
    
    //! Add a vector of branches with the given genealogy
    /*!
      \param [in] level The level or generation index number of the branches (starts at 1)
      \param [in] theBranches The vector of LpcBranch pointers to add
    */
    void addBranches(int level, std::vector<LpcBranch*>& theBranches);

    // Accessors

    //! Get the number of levels = generations
    /*!
      \return the number of levels = generations
    */
    int getNumberLevels() const;

    //! Get the vector of LpcBranch pointers for the given generation
    /*!
      \param [in] generation The generation index number of the set of branches (starts at 1)
      \returns a STL vector of the LpcBranch pointers for the given generation
    */    
    std::vector<LpcBranch*> getBranches(int generation = 1);

    //! Print out the stored branch information
    void print();

protected:
  
private:

    //! Copy constructor
    LpcBranchCollection(const LpcBranchCollection& other);

    //! Assignment operator
    LpcBranchCollection& operator=(const LpcBranchCollection& other);

    //! The map storing the genealogy of the branches
    LpcBranchMap theMap_;

};

#endif


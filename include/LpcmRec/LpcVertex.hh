// Copyright University of Warwick 2014
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// Authors:
// John Back

/*! \file LpcVertex.hh
    \brief File containing the declaration of the LpcVertex class
*/

/*! \class LpcVertex
    \brief Simple class that defines an lpc vertex object
*/

#ifndef LPC_VERTEX_HH
#define LPC_VERTEX_HH

#include <Eigen/Dense>

class LpcVertex {

public:

    //! Default, empty constructor
    LpcVertex();

    //! Main constructor
    /*!
      \param [in] index The index number of the vertex
      \param [in] coords The co-ordinates of the vertex location
      \param [in] curveId The main curve index number used to find the vertex
      \param [in] branchId The branch curve index number used to find the vertex
    */
    LpcVertex(int index, const Eigen::VectorXd& coords, 
	      int curveId = 0, int branchId = 0);

    //! Copy constructor
    /*!
      \param [in] other The vertex to copy
    */
    LpcVertex(const LpcVertex& other);

    //! Assigment operator
    /*!
      \param [in] other The vertex used for the assignment
    */
    LpcVertex& operator=(const LpcVertex& other);

    //! Destructor
    virtual ~LpcVertex();

    // Modifiers

    //! Reset the vertex index number
    /*!
      \param [in] index The vertex index number
    */
    void setIndex(int index) {index_ = index;}

    // Accessors

    //! Get the vertex index number
    /*!
      \returns the vertex index integer
    */
    int getIndex() const {return index_;}

    //! Get the vertex location
    /*!
      \returns an Eigen VectorXd of the vertex co-ordinates
    */
    Eigen::VectorXd getCoords() const {return coords_;}

    //! Get the vertex co-ordinate dimension
    /*!
      \returns the number of co-ordinate dimensions
    */
    int getNDimensions() const {return coords_.size();}

    //! Get the main curve id associated to the vertex
    /*!
      \return the main curve index integer
    */
    int getCurveId() const {return curveId_;}

    //! Get the branch curve id associated to the vertex
    /*!
      \return the branch curve index integer
    */
    int getBranchId() const {return branchId_;}

    //! Print out the vertex info
    void print() const;

protected:


private:

    //! The vertex index number
    int index_;

    //! The vertex co-ordinates
    Eigen::VectorXd coords_;

    //! The curve index number used to find the vertex
    int curveId_;

    //! The branch index number used to find the vertex
    int branchId_;

};

#endif

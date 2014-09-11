// Copyright University of Warwick 2014
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// Authors:
// John Back

/*! \file LpcHit.hh
    \brief File containing the declaration of the LpcHit class
*/

/*! \class LpcHit
    \brief Class used to define a hit in the dataset or point cloud

    Class used to define a hit in the dataset or point cloud. Contains
    co-ordinate and weight (charge/energy) information. Also contains the 
    list of truth particles that was used to create the hit "voxel"
*/

#ifndef LPC_HIT_HH
#define LPC_HIT_HH

#include <Eigen/Dense>

class LpcHit {

public:

    //! Default, empty constructor. Set the number of dimensions to 3
    /*! (x,y,z) = (0,0,0), no charge nor energy deposition, no truth ids
    */
    LpcHit();

    //! Constructor using a vector of co-ordinates; use this for nDim > 3.
    /*!
      \param [in] index the integer index of the hit
      \param [in] coords Vector containing the co-ordinates (x,y,z,..)
      \param [in] weight The weight assigned to the hit, e.g. charge or deposited energy
    */
    LpcHit(int index, const Eigen::VectorXd& coords, double weight);

    //! Constructor using a vector of co-ordinates with truth information
    /*!
      \param [in] index the integer index of the hit
      \param [in] coords Vector containing the co-ordinates (x,y,z,..)
      \param [in] weight The weight assigned to the hit, e.g. charge or deposited energy
      \param [in] truthIds The list of truth ids used to create the (MC) hit
    */
    LpcHit(int index, const Eigen::VectorXd& coords, double weight,
	   const Eigen::VectorXi& truthIds);

    //! Constructor for 3D hits, no truth information
    /*!
      \param [in] index the integer index of the hit
      \param [in] x The x co-ordinate
      \param [in] y The y co-ordinate
      \param [in] z The z co-ordinate
      \param [in] weight The weight assigned to the hit, e.g. charge or deposited energy
    */
    LpcHit(int index, double x, double y, double z, double weight);

    //! Constructor for 3D hits with truth information
    /*!
      \param [in] index the integer index of the hit
      \param [in] x The x co-ordinate
      \param [in] y The y co-ordinate
      \param [in] z The z co-ordinate
      \param [in] weight The weight assigned to the hit, e.g. charge or deposited energy
      \param [in] truthIds The list of truth ids used to create the (MC) hit
    */
    LpcHit(int index, double x, double y, double z, double weight,
	   const Eigen::VectorXi& truthIds);

    //! Constructor for 2D hits (x & y co-ordinates only), no truth information
    /*!
      \param [in] index the integer index of the hit
      \param [in] x The x co-ordinate
      \param [in] y The y co-ordinate
      \param [in] weight The weight assigned to the hit, e.g. charge or deposited energy
    */
    LpcHit(int index, double x, double y, double weight);

    //! Constructor for 2D hits (x & y co-ordinates only) with truth information
    /*!
      \param [in] index the integer index of the hit
      \param [in] x The x co-ordinate
      \param [in] y The y co-ordinate
      \param [in] weight The weight assigned to the hit, e.g. charge or deposited energy
      \param [in] truthIds The list of truth ids used to create the (MC) hit      
    */
    LpcHit(int index, double x, double y, double weight, 
	   const Eigen::VectorXi& truthIds);

    //! Constructor for 1D hits (x co-ordinate only), no truth information
    /*!
      \param [in] index the integer index of the hit
      \param [in] x The x co-ordinate
      \param [in] weight The weight assigned to the hit, e.g. charge or deposited energy
    */
    LpcHit(int index, double x, double weight);

    //! Constructor for 1D hits (x co-ordinate only) with  truth information
    /*!
      \param [in] index the integer index of the hit
      \param [in] x The x co-ordinate
      \param [in] weight The weight assigned to the hit, e.g. charge or deposited energy
      \param [in] truthIds The list of truth ids used to create the (MC) hit      
    */
    LpcHit(int index, double x, double weight, const Eigen::VectorXi& truthIds);

    //! Destructor
    virtual ~LpcHit();

    //! Copy constructor
    /*!
      \param [in] other The LpcHit to copy
    */
    LpcHit(const LpcHit& other);

    //! Assignment operator
    /*!
      \param [in] other The LpcHit used in the assignment
    */
    LpcHit& operator=(const LpcHit& other);
    
    // Accessor functions in order to retrieve id number, co-ordinates, Q, truthIds

    //! The index index number
    /*!
      \return the index index number
    */
    int getIndex() const {return index_;}

    //! The first co-ordinate, usually x
    /*!
      \return the first co-ordinate, usually x
    */
    double getFirstCoord() const;

    //! The second co-ordinate, usually y
    /*!
      \return the second co-ordinate, usually y
    */
    double getSecondCoord() const;

    //! The third co-ordinate, usually z
    /*!
      \return the third co-ordinate, usually z
    */
    double getThirdCoord() const;

    //! The vector of the generalised, n-dimensional position
    /*!
      \return the vector of the co-ordinates of the hit
    */
    Eigen::VectorXd getCoords() const {return coords_;}

    //! The co-ordinate value given by index (0,1,2,...,nDim-1)
    /*!
      \return the co-ordinate value for the given index
    */
    double getCoord(int index) const;

    //! The weight of the hit, usually the charge or deposited energy
    /*!
      \return the weight of the hit
    */
    double getWeight() const {return weight_;}

    //! The vector of truth integer ids
    /*!
      \return the vector of truth integer ids
    */
    Eigen::VectorXi getTruthIds() const {return truthIds_;}

    //! The number of co-ordinate dimensions
    /*!
      \return the number of co-ordinate dimensions
    */
    int nDimensions() const {return nDim_;}

    //! Print function
    void print() const;

protected:

private:

    //! The index index number
    int index_;

    //! Internal representation of the co-ordinates of the hit
    //! Note that this may include more than 3 dimensions (x,y,z,...)
    Eigen::VectorXd coords_;

    //! The hit weight, usually the charge or deposited energy
    double weight_;

    //! The number of co-ordinate dimensions (set by the constructor)
    int nDim_;

    //! The vector of truth id integers
    Eigen::VectorXi truthIds_;
    
};

#endif

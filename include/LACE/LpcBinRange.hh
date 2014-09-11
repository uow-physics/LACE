// Copyright University of Warwick 2014
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// Authors:
// John Back

/*! \file LpcBinRange.hh
    \brief File containing the declaration of the LpcBinRange class
*/

/*! \class LpcBinRange
    \brief Simple class to define a (min,max) range of lpc point indices

*/

#ifndef LPC_BIN_RANGE_HH
#define LPC_BIN_RANGE_HH

class LpcBinRange {

public:

    //! Empty constructor
    LpcBinRange() : minBin_(0), maxBin_(0) {}

    //! Main constructor
    LpcBinRange(int minBin, int maxBin) : minBin_(minBin),
					  maxBin_(maxBin) {}

    //! Destructor
    virtual ~LpcBinRange() {}

    //! Get the minimum bin entry
    /*!
      \return the minimum lpc point index
    */
    int getMinBin() const {return minBin_;}

    //! Get the maximum bin entry
    /*!
      \return the maximum lpc point index
    */
    int getMaxBin() const {return maxBin_;}

    //! Set the minimum bin entry
    /*!
      \param [in] minBin The minimum lpc point index
    */
    void setMinBin(int minBin) {minBin_ = minBin;}

    //! Set the maximum bin entry
    /*!
      \param [in] maxBin The maximum lpc point index
    */
    void setMaxBin(int maxBin) {maxBin_ = maxBin;}

    //! Less than operator. Only check and compare the minBin values
    bool operator < (const LpcBinRange& other) const 
                    {return minBin_ < other.minBin_;}

protected:

private:

    //! Minimum bin
    int minBin_;

    //! Maximum bin
    int maxBin_;

};

#endif


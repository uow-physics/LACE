// Copyright University of Warwick 2014
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// Authors:
// John Back

/*! \file LpcOutputFactory.hh
    \brief File containing the declaration of the LpcOutputFactory class
*/

/*! \class LpcOutputFactory
    \brief Class to select the required LpcAbsOutput implementation
*/

#ifndef LPC_OUTPUT_FACTORY_HH
#define LPC_OUTPUT_FACTORY_HH

class LpcAbsOutput;
class LpcParameters;

class LpcOutputFactory {

public:

    //! Constructor
    LpcOutputFactory();

    //! Destructor
    virtual ~LpcOutputFactory();

    //! Create the specific LpcAbsOutput implementation we need
    /*!
      \param [in] thePars The pointer to the constant LpcParameters
      \return a new LpcAbsOutput pointer to the specific output implementation
    */
    LpcAbsOutput* getOutputWriter(const LpcParameters* thePars) const;

protected:

private:

};

#endif


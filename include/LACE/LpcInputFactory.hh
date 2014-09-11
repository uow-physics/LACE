// Copyright University of Warwick 2014
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// Authors:
// John Back

/*! \file LpcInputFactory.hh
    \brief File containing the declaration of the LpcInputFactory class
*/

/*! \class LpcInputFactory
    \brief Class to select the required LpcAbsInput implementation
*/

#ifndef LPC_INPUT_FACTORY_HH
#define LPC_INPUT_FACTORY_HH

class LpcAbsInput;
class LpcParameters;

class LpcInputFactory {

public:

    //! Constructor
    LpcInputFactory();

    //! Destructor
    virtual ~LpcInputFactory();

    //! Create the specific LpcAbsInput implementation we need
    /*!
      \param [in] thePars The pointer to the constant LpcParameters
      \return a new LpcAbsInput pointer to the specific input implementation
    */
    LpcAbsInput* getInputReader(const LpcParameters* thePars) const;

protected:

private:

};

#endif


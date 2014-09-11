// Copyright University of Warwick 2014
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// Authors:
// John Back

/*! \file LpcInputFactory.cc
    \brief Class to select the required LpcAbsInput implementation
*/

#include "LACE/LpcInputFactory.hh"

#include "LACE/LpcAbsInput.hh"
#include "LACE/LpcTextInput.hh"
#include "LACE/LpcRootInput.hh"
#include "LACE/LpcParameters.hh"

#include <locale>
#include <string>

LpcInputFactory::LpcInputFactory()
{
}

LpcInputFactory::~LpcInputFactory()
{
}

LpcAbsInput* LpcInputFactory::getInputReader(const LpcParameters* thePars) const
{

    // Try to select the required input reader
    LpcAbsInput* theInput(0);

    // Get the name and format of the input file
    std::string inFileName("lpcInput.txt");
    std::string inFormat("text");

    if (thePars) {
	inFileName = thePars->getInputFileName();
	inFormat = thePars->getInputFormat();
    }

    for (size_t i = 0; i < inFormat.size(); i++) {
	inFormat[i] = std::tolower(inFormat[i], std::locale());
    }

    if (inFormat == "root") {

#ifdef LPC_USE_ROOT
	std::string treeName("Data");
	theInput = new LpcRootInput(inFileName, treeName);
#endif

    } else {

	theInput = new LpcTextInput(inFileName);

    }

    return theInput;

}

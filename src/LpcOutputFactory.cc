// Copyright University of Warwick 2014
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// Authors:
// John Back

/*! \file LpcOutputFactory.cc
    \brief Class to select the required LpcAbsOutput implementation
*/

#include "LACE/LpcOutputFactory.hh"

#include "LACE/LpcAbsOutput.hh"
#include "LACE/LpcTextOutput.hh"
#include "LACE/LpcRootOutput.hh"
#include "LACE/LpcParameters.hh"

#include <locale>
#include <string>

LpcOutputFactory::LpcOutputFactory()
{
}

LpcOutputFactory::~LpcOutputFactory()
{
}

LpcAbsOutput* LpcOutputFactory::getOutputWriter(const LpcParameters* thePars) const
{

    // Try to select the required output writer
    LpcAbsOutput* theOutput(0);

    // Get the name and format of the output file
    std::string outFileName("lpcOutput.txt");
    std::string outFormat("text");

    if (thePars) {
	outFileName = thePars->getOutputFileName();
	outFormat = thePars->getOutputFormat();
    }

    for (size_t i = 0; i < outFormat.size(); i++) {
	outFormat[i] = std::tolower(outFormat[i], std::locale());
    }

    if (outFormat == "root") {

#ifdef LPC_USE_ROOT
	theOutput = new LpcRootOutput(outFileName);
#endif

    } else {

	theOutput = new LpcTextOutput(outFileName);

    }

    return theOutput;

}

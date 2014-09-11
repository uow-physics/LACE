// Copyright University of Warwick 2014
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// Authors:
// John Back

/*! \file LpcTextInput.hh
    \brief File containing the declaration of the LpcTextInput class
*/

/*! \class LpcTextInput
    \brief Class to read in an ascii file containing a dataset of hits (point cloud)
*/

#ifndef LPC_TEXT_INPUT_HH
#define LPC_TEXT_INPUT_HH

#include "LACE/LpcAbsInput.hh"

#include <fstream>
#include <string>
#include <utility>
#include <vector>

class LpcEvent;
class LpcHit;

class LpcTextInput : public LpcAbsInput {

public:

    //! Constructor
    /*!
      \param [in] inputFileName The name of the input file
    */
    LpcTextInput(const std::string& inputFileName);

    //! Destructor
    virtual ~LpcTextInput();

    //! Initialisation function
    virtual void initialise();

    //! Finalising function
    virtual void finalise() {getData_.close();}
    
    //! Get the required event
    /*!
      \param [in] eventNo is the event number to retrieve
      \returns a new pointer to the LpcEvent
    */
    virtual LpcEvent* getEvent(int eventNo);

protected:

private:

    //! Private default constructor
    LpcTextInput();

    //! Private copy constructor
    LpcTextInput(const LpcTextInput& other);

    //! Private assignment operator
    LpcTextInput& operator=(const LpcTextInput& other);

    //! The input string stream
    std::ifstream getData_;

    //! The white space delimiter for splitting up each data hit line
    std::string whiteSpace_;

    //! The maximum readable character size of a given line of text
    const int maxNChar_;

    //! The index number of the current hit to be created
    int hitIndex_;

    //! Reads the line containing the event number and number of hits
    std::pair<int, int> getNextEventIndices();

    //! Form a hit from the next line of input data
    /*!
      \returns a new LpcHit pointer
    */
    LpcHit* getNextHit();

    //! Split a string using a delimiter, such as " " or ","
    /*!
      \param [in] theString the string to be split
      \param [in] splitter the delimiter for splitting the string
      \returns a vector of strings that correspond to the split string segments
    */
    std::vector<std::string> split(const std::string& theString,
				   const std::string& splitter) const;

};

#endif


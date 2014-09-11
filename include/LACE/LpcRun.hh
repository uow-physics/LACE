// Copyright University of Warwick 2014
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// Authors:
// John Back

/*! \file LpcRun.hh
    \brief File containing the declaration of the LpcRun class
*/

/*! \class LpcRun
    \brief Class to run the overall lpc algorithm with a given input dataset
*/

#ifndef LPC_RUN_HH
#define LPC_RUN_HH

#include <string>

class LpcAbsInput;
class LpcAbsOutput;
class LpcParameters;
class LpcProcess;

class LpcRun {

public:

    //! Enumeration to specify the file format
    //enum FileFormat {Text = 0};

    //! Constructor
    /*!
      \param [in] parameterFileName The name of the parameter file
    */
    LpcRun(const std::string& parameterFileName);

    //! Destructor
    virtual ~LpcRun();

    //! Run the lpc algorithm on the provided dataset
    void run();

protected:

private:

    //! Private default constructor
    LpcRun();

    //! Private copy constructor
    LpcRun(const LpcRun& other);

    //! Private assignment operator
    LpcRun& operator=(const LpcRun& other);

    //! Intialisation function called by the constructor
    void initialise();

    //! The parameter filename
    std::string parameterFileName_;

    //! The pointer to the lpc parameters
    LpcParameters* theParameters_;

    //! The name of the input data file
    std::string inFileName_;

    //! The name of the output data file
    std::string outFileName_;

    //! A pointer to the lpc processing steps
    LpcProcess* theProcess_;

    //! A pointer to the input class
    LpcAbsInput* theInput_;

    //! A pointer to the output class
    LpcAbsOutput* theOutput_;

};

#endif

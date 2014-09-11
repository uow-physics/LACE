// Copyright University of Warwick 2014
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// Authors:
// John Back

/*! \file LpcTextOutput.hh
    \brief File containing the declaration of the LpcTextOutput class
*/

/*! \class LpcTextOutput
    \brief Class to write out the results of the lpc algorithm into a text file
*/

#ifndef LPC_TEXT_OUTPUT_HH
#define LPC_TEXT_OUTPUT_HH

#include "LACE/LpcAbsOutput.hh"

#include <Eigen/Dense>
#include <fstream>
#include <string>

class LpcAbsCurve;
class LpcCluster;
class LpcCurve;
class LpcPoint;
class LpcVertex;

class LpcTextOutput : public LpcAbsOutput {

public:

    //! Constructor
    /*!
      \param [in] outputFileName The name of the output file
      \param [in] precision The number of decimal places precision
    */
    LpcTextOutput(const std::string& outputFileName, int precision = 6);

    //! Destructor
    virtual ~LpcTextOutput();

    //! Initialisation function
    virtual void initialise();
    
    //! Finalising function
    virtual void finalise() {writeData_.close();}
    
protected:

    //! Store any initial information about the event
    virtual void storeInitialInfo();

    //! Store the main curve and its branches. Pure virtual method that must be implemented
    /*!
      \param [in] mainCurve The main lpc curve (read-only access)
    */
    virtual void storeCurve(const LpcCurve* mainCurve);

    //! Store the details of either the main curve or a branch
    /*!
      \param [in] theCurve The main curve or a branch
    */
    void storeCurveDetails(const LpcAbsCurve* theCurve);

    //! Store the vertices
    /*!
      \param [in] theVertices The vector of vertex pointers (read-only access)
    */
    virtual void storeVertices(const std::vector<LpcVertex*>& theVertices);

    //! Store the clusters
    /*!
      \param [in] theClusters The vector of cluster pointers (read-only access)
    */
    virtual void storeClusters(const std::vector<LpcCluster*>& theClusters);

    //! Store anything else
    virtual void storeExtraInfo();

private:

    //! Private default constructor
    LpcTextOutput();

    //! Private copy constructor
    LpcTextOutput(const LpcTextOutput& other);

    //! Private assignment operator
    LpcTextOutput& operator=(const LpcTextOutput& other);

    //! The output stream
    std::ofstream writeData_;

    //! The stream precision
    int precision_;

    //! The separation width for the output columns
    int width_;

    //! Print out a LpcPoint
    /*!
      \param [in] preamble The initial string/words for the print line
      \param [in] thePoint The LpcPoint to print (unscaled co-ordinates)
    */
    void printPoint(const std::string& preamble, const LpcPoint& thePoint);

    //! Print out a LpcPoint
    /*!
      \param [in] preamble The initial string for the print line
      \param [in] index An integer index
      \param [in] thePoint The LpcPoint to print (unscaled co-ordinates)
    */
    void printPoint(const std::string& preamble, int index, const LpcPoint& thePoint);

    //! Print out the co-oridnates using an Eigen VectorXd
    /*!
      \param [in] preamble The initial string for the print line
      \param [in] coords The co-ordinates to print
    */
    void printCoords(const std::string& preamble, const Eigen::VectorXd& coords);

    //! Print out Eigen VectorXd object
    /*!
      \param [in] preamble The initial string for the print line
      \param [in] vector The Eigen VectorXd object
    */
    void printVector(const std::string& preamble, const Eigen::VectorXd& vector);

    //! Print out Eigen VectorXi object
    /*!
      \param [in] preamble The initial string for the print line
      \param [in] vector The Eigen VectorXi object
    */
    void printVector(const std::string& preamble, const Eigen::VectorXi& vector);

    //! Print out Eigen MatrixXd object
    /*!
      \param [in] preamble The initial string for the print line
      \param [in] matrix The Eigen MatrixXd object
    */
    void printMatrix(const std::string& preamble, const Eigen::MatrixXd& matrix);
    
    //! Print out a STL vector of integers
    /*!
      \param [in] preamble The initial string for the print line
      \param [in] vector The vector of integers
    */
    void printVector(const std::string& preamble, const std::vector<int>& vector);


};

#endif


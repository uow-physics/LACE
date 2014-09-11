// Copyright University of Warwick 2014
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// Authors:
// John Back

/*! \file LpcParameters.hh
    \brief File containing the declaration of the LpcParameters class
*/

/*! \class LpcParameters
    \brief Class used to define the parameters for the lpc algorithm
*/

#ifndef LPC_PARAMETERS_HH
#define LPC_PARAMETERS_HH

#include <string>

class LpcParameters {

public:

    //! Default constructor
    /*! infile = "input.txt", outfile = "output.txt", informat = text, outformat = text,
        apply scaled co-ordinates = true,
        scaled kernel width = 0.05, scaled step size = 0.05, number of lpc points = 250,
        angle penalisation exponent factor = 2.0, ratio of eigenvalues for branching = 0.4,
	boundary condition = 0.005, convergence requirement = 1e-6, branch = 0, gapSize = 1.5,
	peakFinderMethod = 1 (TSpectrum if ROOT exists), otherwise = 0 (Simple method),
	cosAngleCut = 0.01, minPeakFrac = 0.01, peakDiffSq = 3, firstEvt = -1, lastEvt = -1,
	clustering = 1, doBranchVertexing = 1, minVtxResCut = 20.0, convexHull = 0.12,
	showerRes = 20.0, showerResRatio = 0.3, showerResFrac = 0.9
    */
    LpcParameters();

    //! Constructor using a text file containing the list of parameters
    /*!
      \param [in] parameterFileName name of the text file containing the parameters
    */
    LpcParameters(const std::string& parameterFileName);

    //! Destructor
    virtual ~LpcParameters();

    //! Copy constructor
    /*! 
       \param[in] rhs are the parameters to be copied
    */
    LpcParameters(const LpcParameters& rhs);

    // Setting parameter values using functions

    //! Set the default values
    void setDefaultValues();

    //! Reset the parameters given a new parameter file
    /*!
      \param [in] newParamFileName the new parameter file name
    */
    void resetParameters(const std::string& newParamFileName);

    //! Set the name of the input data file
    /*!
      \param [in] inFileName The name of the input data file
    */
    void setInputFileName(const std::string& inFileName) {inFileName_ = inFileName;}

    //! Set the name of the output data file
    /*!
      \param [in] outFileName The name of the output data file
    */
    void setOutputFileName(const std::string& outFileName) {outFileName_ = outFileName;}

    //! Set the input file format
    /*!
      \param [in] inputFormat The input file format type
    */
    void setInputFormat(const std::string& inputFormat) {inputFormat_ = inputFormat;}

    //! Set the output file format
    /*!
      \param [in] outputFormat The output file format integer
    */
    void setOutputFormat(const std::string& outputFormat) {outputFormat_ = outputFormat;}

    //! Set the flag for scaled co-ordinates: 0 = off, 1 = on (default)
    /*!
      \param [in] flag The flag for turning off (0) or on (1 = default) co-ordinate scaling
    */
    void setScalingFlag(int flag) {scalingFlag_ = flag;}

    //! Set the scaled kernel width
    /*!
        \param [in] kernelWidth is the scaled width of the kernel function
    */
    void setKernelWidth(double kernelWidth) {kernelWidth_ = kernelWidth;}

    //! Set the scaled step size
    /*!
        \param [in] stepSize is the scaled size of the steps for the algorithm
    */
    void setStepSize(double stepSize) {stepSize_ = stepSize;}

    //! Set the number of lpc points
    /*!
        \param [in] nPoints is the number of lpc points
    */
    void setNLpcPoints(int nPoints) {nPoints_ = nPoints;}

    //! Set the angle penalisation exponent factor
    /*!
        \param [in] penalisation is the angle penalisation exponent factor
    */
    void setAnglePenalisation(double penalisation) {penalisation_ = penalisation;}

    //! Set the eigenvalue ratio cut for identifying possible branching points
    /*!
        \param [in] eigenRatio is the minimum ratio of eigenvalues for identifying
	possible branching points
    */
    void setEigenRatio(double eigenRatio) {eigenRatio_ = eigenRatio;}

    //! Set the boundary condition
    /*!
        \param [in] boundary is the boundary condition value
    */
    void setBoundary(double boundary) {boundary_ = boundary;}

    //! Set the convergence requirement
    /*!
        \param [in] convergence is the convergence requirement
    */
    void setConvergence(double convergence) {convergence_ = convergence;}

    //! Set the number of possible branching level (generations)
    /*!
      \param [in] branchLevel The number of possible branching levels (generations)
    */
    void setBranchLevel(int branchLevel) {branchLevel_ = branchLevel;}

    //! Set the branch gap size
    /*!
      \param [in] gapSize The scaled branch gap size
    */
    void setBranchGapSize(double gapSize) {gapSize_ = gapSize;}

    //! Set the peak finding method (0 = none, 1 = simple, 2 = TSpectrum)
    /*!
      \param [in] method The integer specifying the peak finding algorithm
    */
    void setPeakFinder(int method) {peakFinder_ = method;}

    //! Set the minimum threshold cut for 1-|cosAngle|
    /*!
      \param [in] cosAngleCut The minimum threshold value of 1-|cosAngle|
    */
    void setCosAngleCut(double cosAngleCut) {cosAngleCut_ = cosAngleCut;}

    //! Set minimum fraction for the next peak w.r.t the previous peak
    /*!
      \param [in] minPeakFrac The minimum peak fraction
    */
    void setMinPeakFrac(double minPeakFrac) {minPeakFrac_ = minPeakFrac;}

    //! Set the "distance" (squared) limit for checking if lpc peaks "overlap"
    /*!
      \param [in] peakDiffSq The squared-distance limit for cosine peak overlaps
    */
    void setPeakDiffSq(double peakDiffSq) {peakDiffSq_ = peakDiffSq;}

    //! Set the first event to process
    /*!
      \param [in] evtNumber The index of the first event
    */
    void setFirstEvent(int evtNumber) {firstEvent_ = evtNumber;}

    //! Set the last event to process
    /*!
      \param [in] evtNumber The index of the last event
    */
    void setLastEvent(int evtNumber) {lastEvent_ = evtNumber;}

    //! Set the integer to of the clustering algorithm (0 = none, 1 = line, 2 = one)
    /*!
      \param [in] index The integer number of the cluster algorithm
    */
    void setClustering(int index) {clustering_ = index;}

    //! Set the integer to specify if we want vertexing of branches to be done (off = 0, on != 0)
    /*!
      \param [in] enable Specify if we want to turn on/off branch vertexing
    */
    void setBranchVtx(int enable) {doBranchVtx_ = enable;}

    //! Set the minimum selection cut value allowed for hit-to-lpc residuals for vertexing
    /*!
      \param [in] minVtxResCut The minimum selection cut value for residuals in vertexing
    */    
    void setMinVtxResCut(double minVtxResCut) {minVtxResCut_ = minVtxResCut;}

    //! Set the minimum value of the convex hull ratio L_transverse/L_longitudinal for showers
    /*!
      \param [in] convexHull The minimum value of the convex hull ratio for showers
    */
    void setConvexHull(double convexHull) {convexHull_ = convexHull;}

    //! Set the minimum threshold for hit-to-lpc residuals for showers
    /*!
      \param [in] showerRes The minimum threshold for hit-to-lpc residuals for showers
    */
    void setShowerRes(double showerRes) {showerRes_ = showerRes;}

    //! Set the minimum hit-to-lpc residual ratio for showers
    /*!
      param [in] showerResRatio The residual ratio for showers
    */
    void setShowerResRatio(double showerResRatio) {showerResRatio_ = showerResRatio;}

    //! Set the minimum fraction of residuals that need to be above the shower threshold
    /*!
      \param [in] showerResFrac The minimum fraction of residuals above the shower threshold
    */
    void setShowerResFrac(double showerResFrac) {showerResFrac_ = showerResFrac;}


    // Accessor functions in order to retrieve parameters

    //! The name of the input data file
    /*!
      \return the name of the input data file
    */
    std::string getInputFileName() const {return inFileName_;}

    //! The name of the output data file
    /*!
      \return the name of the output data file
    */
    std::string getOutputFileName() const {return outFileName_;}

    //! The format type of the input
    /*!
      \return the format type of the input
    */
    std::string getInputFormat() const {return inputFormat_;}

    //! The format type of the output
    /*!
      \return the format type of the output
    */
    std::string getOutputFormat() const {return outputFormat_;}

    //! Get the flag for scaled co-ordinates: 0 = off, 1 = on (default)
    /*!
      \return the flag for turning off (0) or on (1 = default) co-ordinate scaling
    */
    int getScalingFlag() const {return scalingFlag_;}
 

    //! The scaled kernel width
    /*!
        \return the scaled kernel width
    */
    double getKernelWidth() const {return kernelWidth_;}

    //! The scaled step size
    /*!
        \return the scaled step size
    */
    double getStepSize() const {return stepSize_;}

    //! The number of lpc points
    /*!
        \return the number of lpc points
    */
    int getNLpcPoints() const {return nPoints_;}

    //! The angle penalisation factor
    /*!
        \return the angle penalisation factor
    */
    double getAnglePenalisation() const {return penalisation_;}

    //! The eigenvalue ratio for identifying possible branching points
    /*!
        \return the eigenvalue ratio for identifying possible branching points
    */
    double getEigenRatio() const {return eigenRatio_;}

    //! The boundary condition
    /*!
        \return the boundary condition
    */
    double getBoundary() const {return boundary_;}

    //! The convergence requirement
    /*!
        \return the convergence requirement
    */
    double getConvergence() const {return convergence_;}

    //! The number of branching level generations
    /*!
      \return the number of possible branching levels (generations)
    */
    int getBranchLevel() const {return branchLevel_;}

    //! The step size gap for branches
    /*!
      \return the step size gap for branches
    */
    double getBranchGapSize() const {return gapSize_;}

    //! The peak finder method selection (0 = none, 1 = simple, 2 = TSpectrum)
    /*!
      \return the integer specifying the peak finder selection method
    */
    int getPeakFinder() const {return peakFinder_;}

    //! The minimum threshold cut for 1-|cosAngle|
    /*!
      \return the minimum threshold cut for 1-|cosAngle|
    */
    double getCosAngleCut() const {return cosAngleCut_;}

    //! The minimum fraction for the next peak w.r.t the previous peak
    /*!
      \return the minimum peak fraction
    */
    double getMinPeakFrac() const {return minPeakFrac_;}

    //! The "distance" (squared) limit for checking if lpc peaks "overlap"
    /*!
      \return the squared-distance limit for cosine peak overlaps
    */
    double getPeakDiffSq() const {return peakDiffSq_;}

    //! The first event index
    /*!
      \return the first event number
    */
    int getFirstEvent() const {return firstEvent_;}

    //! The last event index
    /*!
      \return the last event number
    */
    int getLastEvent() const {return lastEvent_;}

    //! The integer index of the cluster algorithm to be used
    /*!
      \return the integer to of the clustering algorithm (0 = none, 1 = line, 2 = one)
    */
    int getClustering() const {return clustering_;}

    //! Check to see if vertexing for branches is turned on or off
    /*!
      \return the integer to specify if branch vertexing is off (0) or on (!=0)
    */
    int getBranchVtx() const {return doBranchVtx_;}

    //! The minimum selection cut value allowed for hit-to-lpc residuals for vertexing
    /*!
      \return the minimum hit-to-lpc residual selection cut for vertexing
    */    
    double getMinVtxResCut() const {return minVtxResCut_;}

    //! The minimum value of the convex hull ratio (Ly+Lz)/Lx for showers
    /*!
      \return the minimum convex hull ratio for showers
    */
    double getConvexHull() const {return convexHull_;}

    //! The minimum threshold for hit-to-lpc residuals for showers
    /*!
      \return the minimum threshold for hit-to-lpc residuals for showers
    */
    double getShowerRes() const {return showerRes_;}
    
    //! The minimum hit-to-lpc residual ratio for showers
    /*!
      \return the minimum hit-to-lpc residual ratio for showers
    */
    double getShowerResRatio() const {return showerResRatio_;}

    //! The minimum fraction of residuals that need to be above the shower threshold
    /*!
      \return the minimum fraction of residuals that need to be above the shower threshold
    */
    double getShowerResFrac() const {return showerResFrac_;}
    

    //! Print out the parameters
    void print() const;

protected:

private:

    //! Internal function to set the parameters with an input file
    /*!
      \param [in] parFile the parameter file
    */
    void readParameterFile(const std::string& parFile);

    //! The input file name
    std::string inFileName_;

    //! The output file name
    std::string outFileName_;

    //! The input file format type
    std::string inputFormat_;

    //! The output file format type
    std::string outputFormat_;

    //! The scaling control flag
    int scalingFlag_;

    //! The kernel width
    double kernelWidth_;

    //! The step size
    double stepSize_;

    //! The number of lpc points
    int nPoints_;

    //! The angle penalisation factor
    double penalisation_;

    //! The ratio of eigenvalues for identifying possible branching points
    double eigenRatio_;

    //! The boundary condition
    double boundary_;

    //! The convergence requirement
    double convergence_;

    //! The number of possible branching levels
    int branchLevel_;

    //! The branch gap step size
    double gapSize_;

    //! The peak finder method integer
    int peakFinder_;

    //! The minimum threshold cut for 1-|cosAngle|
    double cosAngleCut_;

    //! The minimum fraction for the next peak w.r.t the previous peak
    double minPeakFrac_;

    //! Set the "distance" (squared) limit for checking if lpc peaks "overlap"
    double peakDiffSq_;

    //! First event index
    int firstEvent_;

    //! Last event index
    int lastEvent_;

    //! The integer index of the cluster algorithm to be used
    int clustering_;

    //! Turn on/off branch vertexing
    int doBranchVtx_;

    //! Specify the minimum line segment lpc-to-hit residual cut for vertexing
    double minVtxResCut_;

    //! The convex hull ratio for showers
    double convexHull_;

    //! The minimum threshold for hit-to-lpc residuals for showers
    double showerRes_;

    //! The minimum hit-to-lpc residual ratio for showers
    double showerResRatio_;

    //! The minimum fraction of residuals that need to be above the shower threshold
    double showerResFrac_;

};

#endif

// Copyright University of Warwick 2014
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// Authors:
// John Back

/*! \file LpcParameters.cc
    \brief Class used to define parameters for the lpc algorithm
*/

#include "LACE/LpcParameters.hh"

#include <cstdlib>
#include <fstream>
#include <locale>
#include <iostream>

LpcParameters::LpcParameters() :
    inFileName_("lpcInput.txt"),
    outFileName_("lpcOutput.txt"),
    inputFormat_("text"),
    outputFormat_("text"),
    scalingFlag_(1),
    kernelWidth_(0.0),
    stepSize_(0.0),
    nPoints_(0),
    penalisation_(0.0),
    eigenRatio_(0.0),
    boundary_(0.0),
    convergence_(0.0),
    branchLevel_(0),
    gapSize_(0.0),
    peakFinder_(0),
    cosAngleCut_(0.0),
    minPeakFrac_(0.0),
    peakDiffSq_(0.0),
    firstEvent_(0),
    lastEvent_(-1),
    clustering_(0),
    doBranchVtx_(0),
    minVtxResCut_(0.0),
    convexHull_(0.0),
    showerRes_(0.0),
    showerResRatio_(0.0),
    showerResFrac_(0.0)
{

    this->setDefaultValues();

}

LpcParameters::LpcParameters(const std::string& parameterFileName) 
{
    // Retrieve the parameters using a text file containing
    // "name value" lines

    // First, use the default values
    this->setDefaultValues();

    // Now update any required parameters
    this->readParameterFile(parameterFileName);

}

void LpcParameters::resetParameters(const std::string& parameterFileName)
{
    // Retrieve the parameters using a text file containing
    // "name value" lines
    this->readParameterFile(parameterFileName);
}

void LpcParameters::readParameterFile(const std::string& parameterFileName) 
{

    // If the parameterFile is empty, use the default parameter values
    if (parameterFileName.size() == 0) {
	std::cout<<"Using the default parameter values."<<std::endl;
	this->setDefaultValues();
	return;
    } else {
	std::cout<<"Reading the lpc parameters from the file "
		 <<parameterFileName<<std::endl;
    }

    std::ifstream getPar(parameterFileName.c_str());

    if (getPar.fail()) {
	std::cerr<<"Could not open the parameter file "
		 <<parameterFileName<<"; using default values."<<std::endl;
	this->setDefaultValues();
	this->print();
	return;
    }

    std::string whiteSpace("\t");
    
    while (getPar.good()) {

	if (getPar.peek() == '\n') {

	    // Finish reading the line
	    char c;
	    getPar.get(c);

	    // Stop while loop if we have reached the end of the file
	    if (getPar.eof()) {break;}

	} else if (getPar.peek() == '#') {

	    // Skip the comment line
	    getPar.ignore(1000, '\n');
	    getPar.putback('\n');

	    // Stop while loop if we have reached the end of the file
	    if (getPar.eof()) {break;}

	} else {

	    // Read the parameter line
	    std::string parName(""), parValue("");
	    getPar >> parName >> parValue;

	    // Use lowercase for the parameter name
	    for (size_t i = 0; i < parName.size(); i++) {
		parName[i] = std::tolower(parName[i], std::locale());
	    }

	    // Find which parameter we have and set its value
	    if (parName == "infile") {

		this->setInputFileName(parValue);

	    } else if (parName == "outfile") {

		this->setOutputFileName(parValue);

	    } else if (parName == "informat") {

		this->setInputFormat(parValue);

	    } else if (parName == "outformat") {

		this->setOutputFormat(parValue);

	    } else if (parName == "usescaling") {

		int flag = atoi(parValue.c_str());
		this->setScalingFlag(flag);

	    } else if (parName == "kernelwidth") {

		double kernelWidth = atof(parValue.c_str());
		this->setKernelWidth(kernelWidth);

	    } else if (parName == "stepsize") {

		double stepSize = atof(parValue.c_str());
		this->setStepSize(stepSize);

	    } else if (parName == "npoints") {

		int nPoints = atoi(parValue.c_str());
		this->setNLpcPoints(nPoints);

	    } else if (parName == "penalisation") {

		double penalisation = atof(parValue.c_str());
		this->setAnglePenalisation(penalisation);

	    } else if (parName == "eigenratio") {

		double eigenRatio = atof(parValue.c_str());
		this->setEigenRatio(eigenRatio);

	    } else if (parName == "boundary") {

		double boundary = atof(parValue.c_str());
		this->setBoundary(boundary);

	    } else if (parName == "convergence") {

		double convergence = atof(parValue.c_str());
		this->setConvergence(convergence);

	    } else if (parName == "branchlevel") {

		int branchLevel = atoi(parValue.c_str());
		this->setBranchLevel(branchLevel);

	    } else if (parName == "gapsize") {

		double gapSize = atof(parValue.c_str());
		this->setBranchGapSize(gapSize);

	    } else if (parName == "cosanglecut") {

		double cosAngleCut = atof(parValue.c_str());
		this->setCosAngleCut(cosAngleCut);

	    } else if (parName == "minpeakfrac") {

		double minPeakFrac = atof(parValue.c_str());
		this->setMinPeakFrac(minPeakFrac);

	    } else if (parName == "peakdiffsq") {

		double peakDiffSq = atof(parValue.c_str());
		this->setPeakDiffSq(peakDiffSq);

	    } else if (parName == "firstevent") {

		int firstEvent = atoi(parValue.c_str());
		this->setFirstEvent(firstEvent);
		
	    } else if (parName == "lastevent") {

		int lastEvent = atoi(parValue.c_str());
		this->setLastEvent(lastEvent);
		
	    } else if (parName == "clustering") {

		int clustering = atoi(parValue.c_str());
		this->setClustering(clustering);

	    }  else if (parName == "dobranchvtx") {

		int doBranchVtx = atoi(parValue.c_str());
		this->setBranchVtx(doBranchVtx);

	    } else if (parName == "minvtxrescut") {

		double minVtxResCut = atof(parValue.c_str());
		this->setMinVtxResCut(minVtxResCut);

	    } else if (parName == "convexhull") {

		double convexHull = atof(parValue.c_str());
		this->setConvexHull(convexHull);

	    } else if (parName == "showerres") {

		double showerRes = atof(parValue.c_str());
		this->setShowerRes(showerRes);

	    } else if (parName == "showerresratio") {

		double showerResRatio = atof(parValue.c_str());
		this->setShowerResRatio(showerResRatio);

	    } else if (parName == "showerresfrac") {

		double showerResFrac = atof(parValue.c_str());
		this->setShowerResFrac(showerResFrac);

	    } else if (parName == "peakfinder") {

		int peakFinder = atoi(parValue.c_str());
		this->setPeakFinder(peakFinder);

	    }

	    // Stop while loop if we have reached the end of the file
	    if (getPar.eof()) {break;}

	}

    }

    this->print();

}

LpcParameters::~LpcParameters() 
{
}

LpcParameters::LpcParameters(const LpcParameters& rhs) 
{

    this->setInputFileName(rhs.getInputFileName());
    this->setOutputFileName(rhs.getOutputFileName());
    this->setInputFormat(rhs.getInputFormat());
    this->setOutputFormat(rhs.getOutputFormat());
    this->setScalingFlag(rhs.getScalingFlag());
    this->setKernelWidth(rhs.getKernelWidth());
    this->setStepSize(rhs.getStepSize());
    this->setNLpcPoints(rhs.getNLpcPoints());
    this->setAnglePenalisation(rhs.getAnglePenalisation());
    this->setEigenRatio(rhs.getEigenRatio());
    this->setBoundary(rhs.getBoundary());
    this->setConvergence(rhs.getConvergence());
    this->setBranchLevel(rhs.getBranchLevel());
    this->setBranchGapSize(rhs.getBranchGapSize());
    this->setPeakFinder(rhs.getPeakFinder());
    this->setCosAngleCut(rhs.getCosAngleCut());
    this->setMinPeakFrac(rhs.getMinPeakFrac());
    this->setPeakDiffSq(rhs.getPeakDiffSq());
    this->setFirstEvent(rhs.getFirstEvent());
    this->setLastEvent(rhs.getLastEvent());
    this->setClustering(rhs.getClustering());
    this->setBranchVtx(rhs.getBranchVtx());
    this->setMinVtxResCut(rhs.getMinVtxResCut());
    this->setConvexHull(rhs.getConvexHull());
    this->setShowerRes(rhs.getShowerRes());
    this->setShowerResRatio(rhs.getShowerResRatio());
    this->setShowerResFrac(rhs.getShowerResFrac());
}


void LpcParameters::setDefaultValues()
{
    // Set the default values for the parameters
    inFileName_ = "lpcInput.txt";
    outFileName_ = "lpcOutput.txt";
    inputFormat_ = "text";
    outputFormat_ = "text";

    scalingFlag_ = 1;
    kernelWidth_ = 0.05;
    stepSize_ = 0.05;
    nPoints_ = 250;
    penalisation_ = 2.0;
    eigenRatio_ = 0.4;
    boundary_ = 0.005;
    convergence_ = 1e-6;

    branchLevel_ = 0;
    gapSize_ = 1.5;
    cosAngleCut_ = 0.01;
    minPeakFrac_ = 0.01;
    peakDiffSq_ = 3.0;

    firstEvent_ = 0;
    lastEvent_ = -1;
    clustering_ = 1;
    doBranchVtx_ = 0;
    minVtxResCut_ = 20.0;

    convexHull_ = 0.12;
    showerRes_ = 20.0;
    showerResRatio_ = 0.3;
    showerResFrac_ = 0.9;

    // Simple peak finder method
    peakFinder_ = 1;

#ifdef LPC_USE_ROOT
    // Use the TSpectrum method
    peakFinder_ = 2;
#endif

}

void LpcParameters::print() const
{
    // Print out the lpc parameters
    std::cout<<"The lpc parameters are:"<<std::endl;
    std::cout<<"  input file name = "<<inFileName_<<", format type = "<<inputFormat_<<std::endl;
    std::cout<<"  output file name = "<<outFileName_<<", format type = "<<outputFormat_<<std::endl;
    if (scalingFlag_ == 0) {
	std::cout<<"  Co-ordinate scaling is disabled (scale factor = 1)"<<std::endl;
    } else {
	std::cout<<"  Co-ordinate scaling is enabled (scale factor = range)"<<std::endl;
    }

    std::cout<<"  scaled kernel width = "<<kernelWidth_<<std::endl;
    std::cout<<"  scaled step size = "<<stepSize_<<std::endl;
    std::cout<<"  number of Lpc points = "<<nPoints_<<std::endl;
    std::cout<<"  angle penalisation factor = "<<penalisation_<<std::endl;
    std::cout<<"  ratio of eigenvalues for branching = "<<eigenRatio_<<std::endl;
    std::cout<<"  boundary condition = "<<boundary_<<std::endl;
    std::cout<<"  convergence criterion = "<<convergence_<<std::endl;
    std::cout<<"  branching generation level = "<<branchLevel_<<std::endl;
    std::cout<<"  branching scaled gap step size = "<<gapSize_<<std::endl;
    std::cout<<"  peak finding method = "<<peakFinder_<<std::endl;
    std::cout<<"  cosine angle feature threshold cut = "<<cosAngleCut_<<std::endl;
    std::cout<<"  minimum fraction for next cosine angle peak = "<<minPeakFrac_<<std::endl;
    std::cout<<"  distance (sq) limit for overlapping cosine peaks = "<<peakDiffSq_<<std::endl;
    std::cout<<"  first event index = "<<firstEvent_<<std::endl;
    std::cout<<"  last event index = "<<lastEvent_<<std::endl;
    std::cout<<"  selected clustering & vertexing algorithm = "<<clustering_<<std::endl;
    std::cout<<"  enable branch vertexing = "<<doBranchVtx_<<std::endl;
    std::cout<<"  minimum selection cut for hit-to-lpc residuals for vertexing = "
	<<minVtxResCut_<<std::endl;
    std::cout<<"  shower convex hull ratio = "<<convexHull_<<std::endl;
    std::cout<<"  minimum threshold for shower residuals = "<<showerRes_<<std::endl;
    std::cout<<"  minimum residual ratio for showers = "<<showerResRatio_<<std::endl;
    std::cout<<"  minimum fraction of shower residuals above threshold = "<<showerResFrac_<<std::endl;


}


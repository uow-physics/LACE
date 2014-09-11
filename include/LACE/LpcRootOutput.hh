#ifdef LPC_USE_ROOT

// Copyright University of Warwick 2014
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// Authors:
// John Back

/*! \file LpcRootOutput.hh
    \brief File containing the declaration of the LpcRootOutput class
*/

/*! \class LpcRootOutput
    \brief Class to write out the results of the lpc algorithm into a ROOT file
*/

#ifndef LPC_ROOT_OUTPUT_HH
#define LPC_ROOT_OUTPUT_HH

#include "LACE/LpcAbsOutput.hh"

#include <Eigen/Dense>
#include <map>
#include <string>
#include <vector>

class LpcAbsCurve;
class LpcCluster;
class LpcCurve;
class LpcVertex;

class TFile;
class TTree;

class LpcRootOutput : public LpcAbsOutput {

public:

    //! Constructor
    /*!
      \param [in] outputFileName The name of the output file
    */
    LpcRootOutput(const std::string& outputFileName);

    //! Destructor
    virtual ~LpcRootOutput();

    //! Initialisation function
    virtual void initialise();
    
    //! Finalising function
    virtual void finalise();
    
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
      \param [in] isABranch Boolean to specify if the curve is a branch or not
      \param [in] branchGeneration The generation number for the branch
    */
    void storeCurveDetails(const LpcAbsCurve* theCurve, bool isABranch = false, 
			   int branchGeneration = 0);

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
    LpcRootOutput();

    //! Private copy constructor
    LpcRootOutput(const LpcRootOutput& other);

    //! Private assignment operator
    LpcRootOutput& operator=(const LpcRootOutput& other);

    //! Typdef for a map storing a pointer to a vector of doubles. 
    //! The key represents the co-ordinate component
    typedef std::map<int, std::vector<double>* > LpcRootMap;

    //! Set-up the output Tree
    void setupTrees();
    
    //! Clean up the ROOT file and the Trees
    void cleanUp();

    //! Clean up function for the internal maps
    /*!
      \param [in] theMap The map to clean up
    */
    void cleanUpMap(LpcRootMap& theMap);

    //! Clean up function for the internal vectors
    /*!
      \param [in] theVector The vector pointer to clean up
    */
    template <class T> void cleanUpVector(std::vector<T>* theVector);

    //! The output ROOT file
    TFile* rootFile_;

    //! The lpc info ROOT tree
    TTree* lpcTree_;

    //! The vertex ROOT tree
    TTree* vtxTree_;

    //! The cluster ROOT tree
    TTree* clTree_;

    //! Boolean to indicate if the branches of the ROOT trees have been defined
    bool definedTrees_;

    //! Boolean to indicate if the lpc branches have been initialise
    bool initLpcBranches_;

    //! Boolean to indicate if the vertex branches have been initialise
    bool initVtxBranches_;

    //! Boolean to indicate if the cluster branches have been initialise
    bool initClBranches_;

    //! Event number
    int eventId_;

    //! Curve index number
    int curveId_;

    //! Branch index number (= 0 for the main curve)
    int branchId_;

    //! Branch generation number (= 0 for the main curve)
    int branchGen_;

    //! Lpc start point: vector of co-ordinate components (x,y,z,...)
    std::vector<double> lpcStartPoint_;

    //! The number of lpc points
    int nLpc_;

    //! Lpc points
    LpcRootMap lpcPoints_;

    //! Lpc principal eigenvectors (i.e. with the largest eigenvalue)
    LpcRootMap lpcEigenVectors_;

    //! Lpc cosine angles
    std::vector<double>* lpcCosAngles_;

    //! Lpc eigenratios
    std::vector<double>* lpcRho_;

    //! Lpc c0 values
    std::vector<double>* lpcC0_;

    //! Lpc lambda values
    std::vector<double>* lpcLambda_;

    //! Lpc lambda axes
    LpcRootMap lpcLambdaAxes_;

    //! Lpc pathlength
    std::vector<double>* lpcPathLength_;

    //! The average lpc residuals over all hits
    std::vector<double>* lpcResiduals_;

    //! The weighted average lpc residuals over all hits
    std::vector<double>* wLpcResiduals_;

    //! Number of high rho points
    int nHighRho_;

    //! The high rho points
    LpcRootMap highRhoPoints_;

    //! The number of cosine angle peaks
    int nCosPeaks_;

    //! The lpc point indices corresponding to the cosine angle peaks
    std::vector<int>* lpcCosPeaks_;

    //! The number of hits in the event
    int nHits_;

    //! The hit co-ordinates
    LpcRootMap hitCoords_;

    //! The nearest lpc residual for each hit
    std::vector<double>* hitResiduals_;

    //! The weighted nearest lpc residual for each hit
    std::vector<double>* wHitResiduals_;

    //! The nearest lpc point for each hit
    std::vector<int>* hitNearLpc_;

    //! Vertexing info
    //! The number of vertices
    int nVtx_;

    //! The index number of the vertex
    int vtxIndex_;

    //! The index of the main curve for the vertex
    int vtxCurveId_;

    //! The index of the branch for the vertex
    int vtxBranchId_;

    //! The co-ordinates of the vertex
    std::vector<double> vtxPoint_;

    //! Clustering info
    //! The number of clusters
    int nCl_;

    //! The index number of the cluster
    int clIndex_;

    //! The index of the main curve for the cluster
    int clCurveId_;

    //! The index of the branch for the cluster
    int clBranchId_;

    //! The minimum lpc point number for the cluster
    int clMinLpc_;

    //! The maximum lpc point number for the cluster
    int clMaxLpc_;

    //! The centroid of the cluster
    std::vector<double> centroid_;

    //! The principal axes of the cluster
    LpcRootMap axes_;

    //! Convex hull lengths
    std::vector<double> convexHull_;

    //! Is a shower integer (0 = no, 1 = yes)
    int isAShower_;

    //! The number of hits in the cluster
    int nClHits_;

    //! The indices of the hits in the cluster
    std::vector<int>* clHitIds_;

    //! The co-ordinates of the hits in the cluster
    LpcRootMap clHitCoords_;

    //! The weights of the hits in the cluster
    std::vector<double>* clHitWeights_;

};

#endif

#endif

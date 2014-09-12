// Macro to create the convex hull ratio plots for monoenergetic
// particle samples. Only use the single clustering algorithm with no branches
// i.e. each event has only one lpc and one cluster

#include <algorithm>
#include <string>
#include <vector>

using std::cout;
using std::endl;
using std::string;
using std::vector;

TCanvas* theCanvas = new TCanvas("theCanvas", "", 900, 700);
gROOT->SetStyle("Plain");
theCanvas->UseCurrentStyle();

TH1D* muHist = 0;
TH1D* eHist = 0;

void plotConvexHull() {

    createPlots("rootFiles/muon0_5GeVLpcOutput.root", 
		"rootFiles/elec0_5GeVLpcOutput.root",
		"gifFiles/convHull0_5GeV.gif");

    createPlots("rootFiles/muon1_5GeVLpcOutput.root", 
		"rootFiles/elec1_5GeVLpcOutput.root",
		"gifFiles/convHull1_5GeV.gif");


}

void createPlots(string muonFile, string elecFile, string plotFile) {

    double minResCut = 20.0;
    double minResRatio = 0.30;
    double resCutFrac = 0.90;
    double convHullCut = 0.12;

    if (muHist) {delete muHist;}
    if (eHist) {delete eHist;}

    muHist = getConvHullHist(muonFile, "muHist", minResRatio, minResCut, resCutFrac);
    muHist->SetLineColor(kBlue);
    muHist->SetLineWidth(2);

    eHist = getConvHullHist(elecFile, "eHist", minResRatio, minResCut, resCutFrac);
    eHist->SetLineColor(kRed);
    eHist->SetLineWidth(2);      
    
    TLegend theLegend(0.6, 0.6, 0.8, 0.7);
    theLegend.SetTextSize(0.045);
    theLegend.SetFillColor(0);
    theLegend.AddEntry(muHist, "Muon", "lp");
    theLegend.AddEntry(eHist, "Electron", "lp");

    theCanvas->Clear();
    gStyle->SetOptStat(0);
    theCanvas->UseCurrentStyle();
  
    muHist->Draw();
    eHist->Draw("same");

    theLegend.Draw("same");

    // Muon histogram usually has highest peak
    double yMax = muHist->GetMaximum()*1.05;
    cout<<"yMax = "<<yMax<<endl;
    TLine theLine(convHullCut, 0.0, convHullCut, yMax);
    theLine.SetLineStyle(3);
    theLine.SetLineWidth(2);
    theLine.Draw();

    theCanvas->Print(plotFile.c_str());

}

TH1D* getConvHullHist(const string& inFileName, const string& histName, 
		      double minResRatio, double minResCut, double resCutFrac) {

    cout<<"Creating histogram for "<<inFileName<<endl;

    // Get the convex hull ratio histogram after applying a selection cut
    // on the lpc residuals: residual > minResRatio*maxResidual && > minResCut
    // and fraction of residuals passing cut > resCutFrac.

    int nR = 150;
    double rMin = 0.0;
    double rMax = 0.75;

    TH1D* theHist = new TH1D(histName.c_str(), "", nR, rMin, rMax);
    theHist->SetDirectory(0);

    theHist->SetXTitle("Convex hull ratio d_{T}/d_{L}");
    theHist->GetYaxis()->SetTicks("+");
    theHist->GetYaxis()->SetLabelOffset(-0.03);
    theHist->SetTitleOffset(1.25, "Y");
    theHist->GetXaxis()->CenterTitle(kTRUE);
    theHist->SetLabelSize(0.045, "X");
    theHist->SetLabelSize(0.05, "Y");
    theHist->SetTitleSize(0.045, "X");
    
    TFile* theFile = TFile::Open(inFileName.c_str());

    TTree* lpcTree = dynamic_cast<TTree*>(theFile->Get("lpcTree"));    
    int lpcEventId;
    vector<double>* hitResiduals = 0;
    lpcTree->SetBranchAddress("eventId", &lpcEventId);
    lpcTree->SetBranchAddress("hitResiduals", &hitResiduals);

    int nLpcEntries = lpcTree->GetEntries();

    TTree* clsTree = dynamic_cast<TTree*>(theFile->Get("clsTree"));
    int clsEventId;
    double clHullX1, clHullX2, clHullX3;
    clsTree->SetBranchAddress("eventId", &clsEventId);
    clsTree->SetBranchAddress("clHullX1", &clHullX1);
    clsTree->SetBranchAddress("clHullX2", &clHullX2);
    clsTree->SetBranchAddress("clHullX3", &clHullX3);

    int nClsEntries = clsTree->GetEntries();

    // Both trees should have the same number of entries
    if (nLpcEntries != nClsEntries) {
	cout<<"Error. nPid = "<<nPid<<" is not equal to nCls = "<<nCls<<endl;
	return 0;
    }

    for (int i = 0; i < nLpcEntries; i++) {

	lpcTree->GetEntry(i);
	clsTree->GetEntry(i);

	if (lpcEventId != clsEventId) {
	    cout<<"Error. LpcEventId "<<lpcEventId<<" != clsEventId "<<clsEventId<<endl;
	    continue;
	}

	// Sort the vector of hit residuals for this single lpc "cluster"
	std::sort(hitResiduals->begin(), hitResiduals->end());

	// Get the maximum residual
	int nResiduals = hitResiduals->size();
	double maxResidual = (*hitResiduals)[nResiduals-1];

	// Only consider hit residuals that are larger than minResidual
	double minResidual = minResRatio*maxResidual;

	double nPassRes(0.0), nTotRes(0.0);
	for (int j = 0; j < nResiduals; j++) {

	    double resValue = (*hitResiduals)[j];
	    if (resValue > minResidual) {	  
		nTotRes += 1.0;
		if (resValue > minResCut) {nPassRes += 1.0;}
	    }
	    
	} // Residual loop

	// Get the fraction of the remaining residuals that were larger than the cut
	double fraction(0.0);
	if (nTotRes > 0.0) {fraction = nPassRes/nTotRes;}

	// Check if this fraction passes the fraction cut
	bool passed(false);
	if (fraction > resCutFrac) {passed = true;}

	if (passed == true) {

	    // We have passed the selection for the lpc residuals
	    // Get the convex hull variables for this cluster
	    double ratio(0.0);
	    // x is along the main principal axis
	    if (fabs(clHullX1) > 1e-10) {ratio = (clHullX2 + clHullX3)/clHullX1;}

	    theHist->Fill(ratio);

	}

    }

    // Normalise the histogram
    double scale(1.0);
    double integral = theHist->Integral();
    if (integral > 0.0) {scale = 1.0/integral;}
    
    theHist->Scale(scale);
    
    theFile->Close();

    return theHist;
  
}

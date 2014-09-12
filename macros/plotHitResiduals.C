// Macro to create the plots of hit-to-lpc residual (>30% of max res).
// For each plot, show the residual cut at 20 mm = 2 cm

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

void plotHitResiduals() {

    createPlot("rootFiles/muon0_5GeVLpcOutput.root", 
	       "rootFiles/elec0_5GeVLpcOutput.root", 
	       "gifFiles/resHist0_5GeV.gif");

    createPlot("rootFiles/muon1_5GeVLpcOutput.root", 
	       "rootFiles/elec1_5GeVLpcOutput.root",
	       "gifFiles/resHist1_5GeV.gif");

}

void createPlot(string muonFile, string elecFile, string plotFile) {

    double minResCut = 20.0;
    double minResRatio = 0.30;
    double resCut = 0.12;
    
    if (muHist) {delete muHist;}
    if (eHist) {delete eHist;}

    muHist = getResHist(muonFile, "muHist", minResRatio);
    muHist->SetLineColor(kBlue);
    muHist->SetLineWidth(2);

    eHist = getResHist(elecFile, "eHist", minResRatio);
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
    
    theCanvas->Print(plotFile.c_str());

    // Find the fraction of clusters with residuals < cutValues
    double muIntegral = muHist->Integral();

    double resCuts[2] = {50.0, 120.0};
    for (int k = 0; k < 2; k++) {

	double resCut = resCuts[k];
	int muBin = getBin(muHist, resCut);
	double muFrac = muHist->Integral(1, muBin)/muIntegral;

	double eIntegral = eHist->Integral();
	int eBin = getBin(eHist, resCut);
	double eFrac = eHist->Integral(1, eBin)/eIntegral;

	cout<<"Residual cut = "<<resCut<<endl;
	cout<<"Fraction for muon passing cut = "<<1.0-muFrac<<endl;
	cout<<"Fraction for electron passing cut = "<<1.0-eFrac<<endl;

    }

}

int getBin(TH1* theHist, double x) {

    // Find the bin number for the given x value
    if (!theHist) {return 0;}

    TAxis* xAxis = theHist->GetXaxis();
    double xMin = xAxis->GetXmin();
    double xMax = xAxis->GetXmax();
    int nX = theHist->GetNbinsX();
    double dx = (xMax - xMin)/(nX*1.0);
    
    int theBin = int(((x - xMin)/dx) + 0.5);
    return theBin;
    
}

TH1D* getResHist(const string& inFileName, const string& histName, 
		 double minResRatio) {

    cout<<"Creating histogram for "<<inFileName<<endl;

    // Get the hit-to-lpc residual for residuals > 30% of maximum residual
    int nR = 120;
    double rMin = 0.0;
    double rMax = 600.0;

    TH1D* theHist = new TH1D(histName.c_str(), "", nR, rMin, rMax);
    theHist->SetDirectory(0);

    theHist->SetXTitle("Hit-to-lpc residuals #delta r' larger than 30% of #delta r_{max} (mm)");

    theHist->GetYaxis()->SetTicks("+");
    theHist->GetYaxis()->SetLabelOffset(-0.03);
    theHist->SetTitleOffset(1.25, "Y");
    theHist->GetXaxis()->CenterTitle(kTRUE);
    theHist->SetTitleSize(0.045, "X");
    theHist->SetLabelSize(0.045, "X");
    theHist->SetLabelSize(0.05, "Y");

    TFile* theFile = TFile::Open(inFileName.c_str(), "read");

    TTree* lpcTree = dynamic_cast<TTree*>(theFile->Get("lpcTree"));

    vector<double>* hitResiduals = 0;
    lpcTree->SetBranchAddress("hitResiduals", &hitResiduals);

    int nEntries = lpcTree->GetEntries();

    for (int i = 0; i < nEntries; i++) {

	// Get the residual info
	lpcTree->GetEntry(i);

	// Sort the vector of hit residuals for this lpc cluster
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

		// We have a residual distance that is large enough
		theHist->Fill(resValue);

	    }

	} // Loop over residuals

    }

    // Normalise the histogram
    double scale(1.0);
    double integral = theHist->Integral();
    if (integral > 0.0) {scale = 1.0/integral;}
    
    theHist->Scale(scale);
    
    theFile->Close();

    return theHist;

}

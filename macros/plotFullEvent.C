TCanvas* theCanvas = new TCanvas("theCanvas", "", 900, 700);
gROOT->SetStyle("Plain");
gStyle->SetOptStat(0);
gStyle->SetOptFit(0);
theCanvas->UseCurrentStyle();

TFile* theFile = 0;
TH3D* theEvent = 0;

#include <vector>
using std::vector;

void plotFullEvent(string inFileName, int eventNo, 
		   string plotName = "fullEvent.gif",
		   double xMin = 0.0, double xMax = 1000.0,
		   double yMin = -500.0, double yMax = 1500.0, 
		   double zMin = -500.0, double zMax = 1500.0, 
		   double xText1 = 0.4, double yText1 = 0.7, 
		   double xText2 = 0.5, double yText2 = 0.2)
{
    
    if (theFile) {theFile->Close();}

    theFile = TFile::Open(inFileName.c_str(), "read");
    if (!theFile) {return;}

    // Create the 3D histogram for the event display
    int nX = 100;
    int nY = 100;
    int nZ = 100;
    
    if (theEvent) {delete theEvent;}
    TH3D* theEvent = new TH3D("theEvent", "", nX, xMin, xMax, nY, yMin, yMax, nZ, zMin, zMax);
    theEvent->SetDirectory(0);

    theEvent->SetXTitle("x (mm)");
    theEvent->SetTitleOffset(1.25, "X");
    theEvent->SetYTitle("y (mm)");
    theEvent->SetZTitle("z (mm)");
    theEvent->SetTitleOffset(1.45, "X");
    theEvent->SetTitleOffset(1.75, "Y");
    theEvent->SetTitleOffset(1.35, "Z");
    theEvent->GetXaxis()->CenterTitle(kTRUE);
    theEvent->GetYaxis()->CenterTitle(kTRUE);
    theEvent->GetZaxis()->CenterTitle(kTRUE);
    theEvent->GetYaxis()->SetNdivisions(506);
    
    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    theCanvas->UseCurrentStyle();
    
    theEvent->Draw();
    
    TTree* lpcTree = dynamic_cast<TTree*>(theFile->Get("lpcTree"));
  
    // Plot the data hits first
    lpcTree->SetMarkerStyle(kFullDotSmall);
    lpcTree->SetMarkerSize(1.0);
    lpcTree->SetMarkerColor(kBlack);
    TString cutString("eventId=="); cutString += eventNo;
    
    lpcTree->Draw("hitX3:hitX2:hitX1", cutString.Data(), "same");
  
    // Plot the lpc points next
    lpcTree->SetMarkerStyle(kOpenCircle);
    lpcTree->SetMarkerSize(1.0);
    //lpcTree->SetMarkerSize(1.25);
    lpcTree->SetMarkerColor(kBlack);
    lpcTree->Draw("lpcX3:lpcX2:lpcX1", cutString.Data(), "same");

    // Plot the clusters
    TTree* clsTree = dynamic_cast<TTree*>(theFile->Get("clsTree"));
    int clEventId, nCl, clIndex;
    clsTree->SetBranchAddress("eventId", &clEventId);
    clsTree->SetBranchAddress("nCl", &nCl);
    clsTree->SetBranchAddress("clIndex", &clIndex);

    const int nClEntries = clsTree->GetEntries();

    int clColour[6] = {kRed, kBlue, kMagenta+1, kAzure+8, kGreen+2, kOrange-9};

    int i, j;

    for (i = 0; i < nClEntries; i++) {

	clsTree->GetEntry(i);
	if (clEventId == eventNo) {

	    //cout<<"nClusters = "<<nCl<<endl;

	    clsTree->SetMarkerStyle(5);
	    clsTree->SetMarkerSize(0.4);

	    clsTree->SetMarkerColor(clColour[clIndex%6]);
	    
	    TString clsCut("eventId=="); clsCut += clEventId;
	    clsCut += "&&clIndex=="; clsCut += clIndex;
	    clsTree->Draw("clHitX3:clHitX2:clHitX1", clsCut.Data(), "same");
	    
	} else if (clEventId > eventNo) {
	    // Stop processing the TTree when the event number 
	    // is above the one we want
	    break;
	}
	
    }
        
    // Replot the lpc points
    lpcTree->SetMarkerStyle(kOpenCircle);
    lpcTree->SetMarkerSize(1.0);
    //lpcTree->SetMarkerSize(1.25);
    lpcTree->SetMarkerColor(kBlack);
    lpcTree->Draw("lpcX3:lpcX2:lpcX1", cutString.Data(), "same");
    
    TTree* vtxTree = dynamic_cast<TTree*>(theFile->Get("vtxTree"));
  
    // Plot the vertices
    vtxTree->SetMarkerStyle(kFullSquare);
    vtxTree->SetMarkerSize(1.5);
    vtxTree->SetMarkerColor(kGreen+2); 
    vtxTree->Draw("vtxX3:vtxX2:vtxX1", cutString.Data(), "same");
  
    // Highlight lpc cosine angle feature points
    const int nLpcEntries = lpcTree->GetEntries();
    int lpcEventId, nCosPeaks;
    vector<int>* cosPeakId = 0;
    lpcTree->SetBranchAddress("eventId", &lpcEventId);
    lpcTree->SetBranchAddress("nCosPeaks", &nCosPeaks);
    lpcTree->SetBranchAddress("cosPeakId", &cosPeakId);

    for (i = 0; i < nLpcEntries; i++) {
      
	lpcTree->GetEntry(i);

	if (lpcEventId == eventNo) {

	    // Highlight lpc feature points
	    for (j = 0; j < nCosPeaks; j++) {
		int feature = (*cosPeakId)[j];
		lpcTree->SetMarkerColor(kOrange+7);
		lpcTree->SetMarkerStyle(kFullCircle);
		lpcTree->SetMarkerSize(1.25);
		TString varString("lpcX3["); varString += feature; varString += "]:";
		varString += "lpcX2["; varString += feature; varString += "]:";
		varString += "lpcX1["; varString += feature; varString += "]";
		lpcTree->Draw(varString.Data(), cutString.Data(), "same");
	    }

	} else if (lpcEventId > eventNo) {
	    break;
	}
      
    }

    theCanvas->Print(plotName.c_str());

}

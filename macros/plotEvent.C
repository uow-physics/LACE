// ROOT macro script to plot an "event". Data points are represented
// by black dots, lpc points by open red circles.
// To run the code within ROOT:
// .L plotEvent.C
// plotEvent(eventNumber, inFileName)
// where eventNumber is the event number integer you want to see
// and inFileName is the string of the ROOT file containing output
// from the lpc algorithm.

// Create a canvas to plot the information
TCanvas* theCanvas = new TCanvas("theCanvas", "", 900, 700);

#include <string>
using std::string;

// Number of bins
int nX = 100;
int nY = 100;
int nZ = 100;

TH3D* theEvent = 0;

void plotEvent(string inFileName = "rootFiles/mup770MeVLpcOutput.root", int eventNo = 6,
	       string plotFile = "gifFiles/mupEvent.gif",
	       double xMin = 0.0, double xMax = 1000.0, double yMin = -500.0,
	       double yMax = 1500.0, double zMin = -500.0, double zMax = 1500.0) {

    // Open the file
    TFile* theFile = TFile::Open(inFileName.c_str(), "read");

    // Create the 3D histogram for the event display
    if (theEvent) {delete theEvent;}

    theEvent = new TH3D("theEvent", "", nX, xMin, xMax, nY, yMin, yMax, nZ, zMin, zMax);
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
    
    // Set the canvas style
    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    theCanvas->UseCurrentStyle();
    
    // Draw the axes
    theEvent->Draw();

    // Get the lpcTree containing the hit and lpc information
    TTree* lpcTree = dynamic_cast<TTree*>(theFile->Get("lpcTree"));
  
    // Plot the data hits first
    lpcTree->SetMarkerStyle(kFullDotSmall);
    lpcTree->SetMarkerSize(1.0);
    lpcTree->SetMarkerColor(kBlack);
    TString cutString("eventId=="); cutString += eventNo;
    
    lpcTree->Draw("hitX3:hitX2:hitX1", cutString.Data(), "same");
    
    // Plot the lpc points
    lpcTree->SetMarkerStyle(kOpenCircle);
    lpcTree->SetMarkerSize(1.0);
    lpcTree->SetMarkerColor(kRed);
    lpcTree->Draw("lpcX3:lpcX2:lpcX1", cutString.Data(), "same");
 
    // Create a png file of the data & lpc points
    theCanvas->Print(plotFile.c_str());
    
    // Close the input file
    theFile->Close();

}

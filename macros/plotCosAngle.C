// Macro to plot the cosine angle distribution between neighbouring
// main eigenvectors along the lpc
#include <vector>
#include <string>

using std::string;
using std::vector;

TCanvas* theCanvas = new TCanvas("theCanvas", "", 900, 700);
TLine* cutLine = 0;
gROOT->SetStyle("Plain");
gStyle->SetOptStat(0);
gStyle->SetOptFit(0);
theCanvas->UseCurrentStyle();

void plotCosAngle() {

    createPlot(6, "rootFiles/mup770MeVLpcOutput.root", 
	       "gifFiles/mupCosAngle1.gif", 0.35, true);

    createPlot(408, "rootFiles/mup770MeVLpcOutput.root", 
	       "gifFiles/mupCosAngle2.gif", 0.07, true);

}

void createPlot(int eventNo, string inFileName, string plotFileName, double maxValue, 
		bool drawCut = true) {

    cout<<"Running plotCosAngle for "<<inFileName<<", event = "<<eventNo<<endl;

    TFile* theFile = TFile::Open(inFileName.c_str(), "read");
    TTree* lpcTree = dynamic_cast<TTree*>(theFile->Get("lpcTree"));
    
    int eventId, nLpc, branchId;
    vector<double>* cosAngles = 0;

    lpcTree->SetBranchAddress("eventId", &eventId);
    lpcTree->SetBranchAddress("nLpc", &nLpc);
    lpcTree->SetBranchAddress("cosAngles", &cosAngles);
    lpcTree->SetBranchAddress("branchId", &branchId);
    lpcTree->SetBranchStatus("*", 0);
    lpcTree->SetBranchStatus("eventId", 1);
    lpcTree->SetBranchStatus("nLpc", 1);
    lpcTree->SetBranchStatus("cosAngles", 1);
    lpcTree->SetBranchStatus("branchId", 1);

    int nEntries = lpcTree->GetEntries();

    TGraph* theGraph = new TGraph();
    theGraph->SetMarkerStyle(kFullCircle);
    theGraph->SetMarkerSize(0.8);
    theGraph->SetMarkerColor(kBlack);
    
    double yMax = 0.0;
    double x_yMax = 0.0;
    
    for (int i = 0; i < nEntries; i++) {
	
	lpcTree->GetEntry(i);
	
	if (eventId == eventNo && branchId == 0) {

	    // We have the event
	    // Loop over the lpc points and store 1 - |cosphi|
	    cout<<"nLpcPoints = "<<nLpc<<endl;

	    for (int j = 0; j < nLpc; j++) {

		double cPhi = (*cosAngles)[j];
		double yValue = 1.0 - fabs(cPhi);
		double xValue = j*1.0;
		theGraph->SetPoint(theGraph->GetN(), xValue, yValue);
	
		if (yValue > yMax) {
		    yMax = yValue;
		    x_yMax = xValue;
		}
		
	    }
	    
	}

    }

    double yHistMax = yMax;
    //yHistMax = ceil(yHistMax*10.0)/10.0;
    yHistMax = maxValue;
    cout<<"For graph: yMax = "<<yMax<<", x_yMax = "<<x_yMax<<", yHistMax = "<<yHistMax<<endl;  
    
    double xMax = nLpc*1.0 + 1.0;
    //double xMax = 60.0;
    
    theCanvas->Clear();

    TH2D* cosPlot = new TH2D("cosPlot", "", 2, 0.0, xMax, 2, 0.0, yHistMax);
    cosPlot->SetDirectory(0);
    cosPlot->SetXTitle("Lpc point number");
    cosPlot->SetYTitle("1 - |cos#phi|");
    cosPlot->GetXaxis()->CenterTitle(kTRUE);
    cosPlot->GetYaxis()->CenterTitle(kTRUE);
    cosPlot->SetLabelSize(0.045, "X");
    cosPlot->SetLabelSize(0.045, "Y");
    cosPlot->SetTitleOffset(1.05, "Y");
    cosPlot->SetTitleSize(0.05, "X");
    cosPlot->SetTitleSize(0.05, "Y");
    
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    
    cosPlot->Draw();
    theGraph->Draw("lp");
    theCanvas->Update();
    
    // Draw horizontal line for selection cut (0.01)
    if (drawCut == true) {
	if (cutLine) {delete cutLine;}
	TLine* cutLine = new TLine(0.0, 0.01, nLpc*1.0, 0.01);
	cutLine->SetLineStyle(kDotted);
	cutLine->SetLineWidth(2);
	cutLine->Draw();
    }
    
    theCanvas->Print(plotFileName.c_str());
    
    theFile->Close();
    
}

// Macro to plot the number of clusters
#include <sstream>
#include <string>
#include <iomanip>
#include <iostream>
#include <vector>

using std::vector;
using std::string;

TCanvas* theCanvas = new TCanvas("theCanvas", "", 1000, 400);
gROOT->SetStyle("Plain");
theCanvas->UseCurrentStyle();

void plotNClusters() {

    createPlots("rootFiles/mup770MeV.root", 
		"rootFiles/mup770MeVLpcOutput.root",
		"gifFiles/mupNCl.gif", 0.9);

    createPlots("rootFiles/ep770MeV.root", 
		"rootFiles/ep770MeVLpcOutput.root",
		"gifFiles/epNCl.gif", 0.6);

}

void createPlots(string inputFile, string outputFile,
		 string plotName, double plotMax = 0.9) {

    TFile* inFile = TFile::Open(inputFile.c_str(), "read");
    TTree* inTree = dynamic_cast<TTree*>(inFile->Get("Data"));
  
    int eventId, nHits;
    vector< vector<double> >* truthId = 0;
    inTree->SetBranchAddress("eventId", &eventId);
    inTree->SetBranchAddress("nHits", &nHits);
    inTree->SetBranchAddress("truthId", &truthId);
    
    vector<int> passedEvt;
    // Find which events have at least 25 proton hits (via MC truth)
    int i, j;
    
    int nEntries = inTree->GetEntries();
    int nProtonCut = 25;
    int nPassedEvents = 0;

    int nEvents = inTree->GetMaximum("eventId") + 1;

    for (i = 0; i < nEntries; i++) {

	inTree->GetEntry(i);
    
	int nProtonHits = 0;

	// Loop over all hits in this event
	for (j = 0; j < nHits; j++) {

	    // Get the vector of truth ids for this hit
	    vector<double> hitIdVect = (*truthId)[j];
	    // Get the first entry and see if we have a proton
	    if (hitIdVect.size() > 0) {
		if (int(hitIdVect[0]) == 2) {nProtonHits += 1;}
	    }
	    
	}

	if (nProtonHits >= nProtonCut) {
	    //cout<<"Event "<<i<<" passed proton cut: "<<nProtonHits<<endl;
	    passedEvt.push_back(1);
	    nPassedEvents += 1;
	} else {
	    //cout<<"Event "<<i<<" failed proton cut: "<<nProtonHits<<endl;
	    passedEvt.push_back(0);
	}
	
    }

    inFile->Close();

    // Now plot the number of clusters for the events 
    // with/without the minimum proton track length
    TFile* outFile = TFile::Open(outputFile.c_str(), "read");
    TTree* clsTree = dynamic_cast<TTree*>(outFile->Get("clsTree"));

    int nCl(0);
    clsTree->SetBranchAddress("eventId", &eventId);
    clsTree->SetBranchAddress("nCl", &nCl);

    nEntries = clsTree->GetEntries();

    double minVal = 0.0;
    double maxVal = clsTree->GetMaximum("nCl") + 1.0;
    cout<<"maxVal = "<<maxVal<<endl;
    int nBins = int(maxVal + 0.01);

    TH1D* nCHist = new TH1D("nCHist", "", nBins, minVal, maxVal);
    TH1D* nCAllHist = new TH1D("nCAllHist", "", nBins, minVal, maxVal);
    nCHist->SetDirectory(0);
    nCAllHist->SetDirectory(0);

    double nPassCut(0.0), nPassAll(0.0);

    int index(-1), oldEvent(-1);
    bool sameEvent(false);

    while (index < nEntries) {

	index++;
	clsTree->GetEntry(index);

	if (index == 0) {
	    // Initialise the old event index
	    oldEvent = eventId;
	} else if (eventId == oldEvent) {
	    // We have the same event. Skip this entry
	    continue;
	} else if (eventId != oldEvent) {
	    // We have the next new event; update the oldEvent index
	    oldEvent = eventId;
	}	

	nPassAll += 1.0;

	// Fill the number of clusters for all events
	nCAllHist->Fill(nCl*1.0);

	// Check the proton track is "long enough"
	if (passedEvt[eventId] == 1) {
	    nPassCut += 1.0;

	    // Fill the number of clusters for events that pass
	    // the proton selection
	    nCHist->Fill(nCl*1.0);

	}

    }

    // Print entries, mean, rms, under/overflow, integral
    gStyle->SetOptFit(0);
    gStyle->SetOptStat("emr");
    theCanvas->Clear();
    
    float small = 0.01;
    gStyle->SetPadBorderMode(0);
    gStyle->SetFrameBorderMode(0);
    
    TLatex* text = new TLatex();
    text->SetNDC(kTRUE);
    text->SetTextSize(0.05);
    
    theCanvas->cd(1);
    
    TPad* newPad = gPad;
    newPad->Divide(2,1);
    
    newPad->cd(1);
    nCHist->Scale(1.0/nCHist->Integral());
    nCHist->SetMaximum(plotMax);
    
    gPad->SetRightMargin(small);
    cout<<"nCHist integral = "<<nCHist->Integral()<<endl;
    nCHist->Draw();
    
    gStyle->SetStatX(0.9);
    gStyle->SetStatY(0.9);
    gStyle->SetStatW(0.25);
    gStyle->SetStatH(0.25);
    gStyle->SetTextSize(2);
    
    theCanvas->UseCurrentStyle();
    
    nCHist->GetXaxis()->CenterTitle(kTRUE);
    nCHist->GetYaxis()->CenterTitle(kTRUE);
    
    int NDiv = nBins;
    
    nCHist->GetXaxis()->SetNdivisions(NDiv);
    nCHist->SetTitleSize(0.0505, "X");
    nCHist->SetLabelSize(0.0525, "X");
    nCHist->SetTitleSize(0.0525, "Y");
    nCHist->SetLabelSize(0.0525, "Y");
    nCHist->SetXTitle("Number of clusters");
    nCHist->SetYTitle("Fraction of events");
    nCHist->Draw();
    
    newPad->cd(2);
    
    cout<<"nCAllHist integral = "<<nCAllHist->Integral()<<endl;
    nCAllHist->Scale(1.0/nCAllHist->Integral());
    nCAllHist->SetMaximum(plotMax);
    
    nCAllHist->GetXaxis()->SetNdivisions(NDiv);
    nCAllHist->GetXaxis()->CenterTitle(kTRUE);
    nCAllHist->GetYaxis()->CenterTitle(kTRUE);
    nCAllHist->SetTitleSize(0.0505, "X");
    nCAllHist->SetLabelSize(0.0525, "X");
    nCAllHist->SetTitleSize(0.0525, "Y");
    nCAllHist->SetLabelSize(0.0525, "Y");
    nCAllHist->SetXTitle("Number of clusters");
    nCAllHist->SetYTitle("Fraction of events");
    
    nCAllHist->Draw();
    
    theCanvas->Print(plotName.c_str());
    
    outFile->Close();

}

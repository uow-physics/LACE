// Macro to plot the vertex distributions
#include <sstream>
#include <string>
#include <iomanip>
#include <iostream>
#include <vector>

using std::vector;
using std::string;

TCanvas* theCanvas = new TCanvas("theCanvas", "", 1200, 800);
gROOT->SetStyle("Plain");
theCanvas->UseCurrentStyle();

// Define vertex offsets
double xOffset = 0.0;
double yOffset = 0.0;
double zOffset = 0.0;

TFile* inFile = 0;

void plotVertex() {

    createPlots("rootFiles/mup770MeV.root", 
		"rootFiles/mup770MeVLpcOutput.root",
		"gifFiles/mupVtx.gif", 20.0, 40, false);

    createPlots("rootFiles/ep770MeV.root", 
		"rootFiles/ep770MeVLpcOutput.root",
		"gifFiles/epVtx.gif", 20.0, 40, false);

}

void createPlots(string inputFile, string outputFile, string plotName, 
		 double axisLimit = 20.0, int nBins = 40, bool only1Vtx = false) {

    if (inFile) {inFile->Close();}

    inFile = TFile::Open(inputFile.c_str(), "read");
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

    // Now plot the vertex distributions for the events with the minimum
    // proton track length
    TFile* outFile = TFile::Open(outputFile.c_str(), "read");
    TTree* vtxTree = dynamic_cast<TTree*>(outFile->Get("vtxTree"));

    int nVtx(0);
    double vtxX1, vtxX2, vtxX3;

    vtxTree->SetBranchAddress("eventId", &eventId);
    vtxTree->SetBranchAddress("nVtx", &nVtx);
    vtxTree->SetBranchAddress("vtxX1", &vtxX1);
    vtxTree->SetBranchAddress("vtxX2", &vtxX2);
    vtxTree->SetBranchAddress("vtxX3", &vtxX3);

    nEntries = vtxTree->GetEntries();

    double minVal = -axisLimit;
    double maxVal = axisLimit;

    TH1D* xHist = new TH1D("xHist", "", nBins, minVal, maxVal);
    xHist->SetDirectory(0);
    xHist->GetXaxis()->CenterTitle(kTRUE);
    xHist->SetTitleSize(0.05, "X");
    xHist->SetLabelSize(0.045, "X");
    xHist->SetLabelSize(0.045, "Y");
    xHist->SetXTitle("x vertex (mm)");
    
    TH1D* yHist = new TH1D("yHist", "", nBins, minVal, maxVal);
    yHist->SetDirectory(0);
    yHist->GetXaxis()->CenterTitle(kTRUE);
    yHist->SetXTitle("y vertex (mm)");
    yHist->SetTitleSize(0.05, "X");
    yHist->SetLabelSize(0.045, "X");
    yHist->SetLabelSize(0.045, "Y");
    
    TH1D* zHist = new TH1D("zHist", "", nBins, minVal, maxVal);
    zHist->SetDirectory(0);
    zHist->GetXaxis()->CenterTitle(kTRUE);
    zHist->SetXTitle("z vertex (mm)");
    zHist->SetTitleSize(0.05, "X");
    zHist->SetLabelSize(0.045, "X");
    zHist->SetLabelSize(0.045, "Y");

    // Also plot the vertex distance for all events (log plot)
    int nDBins = 120;
    TH1D* dHist = new TH1D("dHist", "", nDBins, 0.0, 1200.0);
    dHist->SetDirectory(0);
    dHist->SetLineColor(kBlack);
    dHist->GetXaxis()->CenterTitle(kTRUE);
    dHist->GetYaxis()->CenterTitle(kTRUE);
    dHist->SetXTitle("Vertex distance (mm)");
    
    TString dYTitle("Fraction of events/(");
    char tmpString[100];
    std::sprintf(tmpString, "%.1f", dHist->GetXaxis()->GetBinWidth(1));
    dYTitle += tmpString; dYTitle += " mm)";
    dHist->SetYTitle(dYTitle.Data());

    dHist->SetTitleSize(0.05, "X");
    dHist->SetTitleSize(0.05, "Y");
    dHist->SetLabelSize(0.045, "X");
    dHist->SetLabelSize(0.045, "Y");
    
    //ofstream writeData("problemVtx.txt");

    double d2cmFrac(0.0);
    int index(-1);

    // Increment the tree index until we 
    // reach the number of entries
    while (index < nEntries) {

	index++;
	vtxTree->GetEntry(index);

	// Skip if the event has no vertex
	if (nVtx < 1) {continue;}

	// Check the proton track is "long enough"
	if (passedEvt[eventId] == 0) {
	    //cout<<"Skipping event "<<eventId<<endl;
	    continue;
	}

	if (only1Vtx == true && nVtx != 1) {
	    nPassedEvents -= 1;
	    continue;
	}

	// Loop through any remaining vertices for the given
	// event number and keep the one with the minimum x value
	double xBest(vtxX1), yBest(vtxX2), zBest(vtxX3);

	for (j = 1; j < nVtx; j++) {

	    index++;
	    vtxTree->GetEntry(index);

	    if (vtxX1 < xBest) {
		xBest = vtxX1;
		yBest = vtxX2;
		zBest = vtxX3;
	    }

	}

	double dx = xBest - xOffset;
	double dy = yBest - yOffset;
	double dz = zBest - zOffset;

	xHist->Fill(dx);
	yHist->Fill(dy);
	zHist->Fill(dz);
	
	double vtxDist = sqrt(dx*dx + dy*dy + dz*dz);
	dHist->Fill(vtxDist);

	if (vtxDist > 20.0) {d2cmFrac += 1.0;}

	//if (fabs(xBest) > 500.0) {
	//  writeData<<eventId<<" "<<nVtx<<endl;
	//}

    }

    gStyle->SetOptFit(111);
    // Print entries, mean, rms, under/overflow, integral
    //gStyle->SetOptStat("eMRuoi");
    //gStyle->SetOptStat("MR");
    gStyle->SetOptStat(0);
    theCanvas->Clear();
    theCanvas->Divide(2,2);
    
    double xHeight = xHist->GetMaximum();
    double yHeight = yHist->GetMaximum();
    double zHeight = zHist->GetMaximum();
    
    TF1* xFun = new TF1("xFun", vtxFun, minVal, maxVal, 5);
    xFun->SetLineColor(kBlack);
    double xMean = xHist->GetMean(1);
    double xSigma = xHist->GetRMS(1);
    double r = 0.02;
    double nSigma = 3.0;
    xFun->SetParameters(xHeight, xMean, xSigma, r, nSigma*xSigma);
    SetParNames(xFun);
    
    TF1* yFun = new TF1("yFun", vtxFun, minVal, maxVal, 5);
    yFun->SetLineColor(kBlack);
    double yMean = yHist->GetMean(1);
    double ySigma = yHist->GetRMS(1);
    yFun->SetParameters(yHeight, yMean, ySigma, r, nSigma*ySigma);
    SetParNames(yFun);
    
    TF1* zFun = new TF1("zFun", vtxFun, minVal, maxVal, 5);
    zFun->SetLineColor(kBlack);
    double zMean = zHist->GetMean(1);
    double zSigma = zHist->GetRMS(1);
    zFun->SetParameters(zHeight, zMean, zSigma, r, nSigma*zSigma);
    SetParNames(zFun);
    
    TLatex* text = new TLatex();
    text->SetNDC(kTRUE);
    text->SetTextSize(0.05);
    
    theCanvas->cd(1);
    gPad->SetLogy(0);

    xHist->Fit("xFun", "L");
    double xSigma = getSigma(xFun);
    double xSigmaErr = getSigmaErr(xFun);
    
    xHist->Draw();
    TString xString = getSigmaString(xSigma, xSigmaErr);
    text->DrawLatex(0.1, 0.95, xString.Data());
    
    theCanvas->cd(2);
    yHist->Fit("yFun", "L");
    double ySigma = getSigma(yFun);
    double ySigmaErr = getSigmaErr(yFun);
    
    yHist->Draw();
    TString yString = getSigmaString(ySigma, ySigmaErr);
    text->DrawLatex(0.1, 0.95, yString.Data());
    
    theCanvas->cd(3);
    zHist->Fit("zFun", "L");
    double zSigma = getSigma(zFun);
    double zSigmaErr = getSigmaErr(zFun);
    
    zHist->Draw();
    TString zString = getSigmaString(zSigma, zSigmaErr);
    text->DrawLatex(0.1, 0.95, zString.Data());
    
    cout<<"Effective widths = "<<xSigma<<", "<<ySigma<<", "<<zSigma<<endl;
    double xEff = xHist->Integral()/xHist->GetEntries();
    double yEff = yHist->Integral()/yHist->GetEntries();
    double zEff = zHist->Integral()/zHist->GetEntries();
    cout<<"Efficiencies for "<<inputFile<<" are "<<xEff<<", "<<yEff<<", "<<zEff<<endl;
    
    // Vertex distance plot, using a log scale
    theCanvas->cd(4);
    gPad->SetLogy();
    // Normalise the histograms
    double ND = dHist->Integral();
    dHist->Scale(1.0/ND);
    dHist->Draw();

    theCanvas->Print(plotName.c_str());
    
    cout<<"Fraction of events passing proton cut = "
	<<(nPassedEvents*100.0)/(nEvents*1.0)<<"%"<<endl;

    outFile->Close();

}

TString getSigmaString(double sigma, double sigmaErr) {

  std::ostringstream str;
  str.precision(2);
  str << std::fixed << "#sigma_{eff} = " << sigma 
      << " #pm " << sigmaErr << " mm"<<endl;

  TString word = TString(str.str().c_str());

  return word;

}

TString getFracString(TH1* theHist, int nPassedEvents) {

  if (theHist == 0) {return TString("");}

  double frac = theHist->Integral()/(nPassedEvents*1.0);
  //double frac = theHist->Integral()/theHist->GetEntries();

  std::ostringstream str;
  str.precision(1);
  str << std::fixed << "; Events = " << frac*100 << "%"<<endl;

  TString word = TString(str.str().c_str());

  return word;


}

void SetParNames(TF1* fitFun) {

  fitFun->SetParName(0, "N_{0}");
  fitFun->SetParName(1, "#mu");
  fitFun->SetParName(2, "#sigma_{1}");
  fitFun->SetParName(3, "r");
  fitFun->SetParName(4, "#sigma_{2}");

}

double getSigma(TF1* fitFun) {

  if (fitFun == 0) {return 0.0;}

  double s1 = fitFun->GetParameter(2);
  double r1 = fitFun->GetParameter(3);
  double s2 = fitFun->GetParameter(4);
  double sigma = s1 + r1*s2;
  //cout<<"Sigma = "<<sigma<<" from function "<<fitFun->GetName()<<endl;

  return sigma;

}

double getSigmaErr(TF1* fitFun) {

  if (fitFun == 0) {return 0.0;}

  double r1 = fitFun->GetParameter(3);
  double s2 = fitFun->GetParameter(4);

  double s1Err = fitFun->GetParError(2);
  double r1Err = fitFun->GetParError(3);
  double s2Err = fitFun->GetParError(4);

  double sigmaErrSq = s1Err*s1Err;
  sigmaErrSq += r1*r1*s2Err*s2Err;
  sigmaErrSq += s2*s2*r1Err*r1Err;

  double sigmaErr = sqrt(sigmaErrSq);

  return sigmaErr;

}


double vtxFun(double* x, double* par) {

  double N = par[0];
  double mu = par[1];
  double sigma = par[2];
  double r1 = par[3];
  double sigma2 = par[4];

  if (sigma < 0.0 || sigma2 < 0.0) {return 0.0;}
  if (sigma2 < sigma) {return 0.0;}

  double xVal = x[0];

  double dx = xVal - mu;
  double expTerm1 = 0.5*dx*dx/(sigma*sigma);
  double expTerm2 = 0.5*dx*dx/(sigma2*sigma2);

  double fun = N*(exp(-expTerm1) + r1*exp(-expTerm2));

  return fun;

}

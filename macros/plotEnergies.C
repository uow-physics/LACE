// Macro to plot the quenched energy distributions for mup and ep events

#include <algorithm>
#include <iostream>
#include <map>
#include <string>
#include <vector>

TCanvas* theCanvas = new TCanvas("theCanvas", "", 1200, 800);
gROOT->SetStyle("Plain");
theCanvas->UseCurrentStyle();

TH1D* trueE1Hist = 0;
TH1D* trueE2Hist = 0;
TH1D* recoE1Hist = 0;
TH1D* recoE2Hist = 0;
TH1D* fracE1Hist = 0;
TH1D* fracE2Hist = 0;
TF1* fun1 = 0;

void plotEnergies() {

    createPlots("rootFiles/ep770MeV.root", "rootFiles/ep770MeVLpcOutput.root", 
		"gifFiles/epEnergies.gif", "Electron", 600.0, 500.0, true);

    createPlots("rootFiles/mup770MeV.root", "rootFiles/mup770MeVLpcOutput.root", 
		"gifFiles/mupEnergies.gif", "Muon", 600.0, 500.0, false);

}

void createPlots(std::string inFileName, std::string outFileName,
		 std::string plotFile, std::string particleName,
		 double maxE1, double maxE2, bool doFit1 = true) {

    if (trueE1Hist) {delete trueE1Hist;}
    if (trueE2Hist) {delete trueE2Hist;}
    if (recoE1Hist) {delete recoE1Hist;}
    if (recoE2Hist) {delete recoE2Hist;}
    if (fracE1Hist) {delete fracE1Hist;}
    if (fracE2Hist) {delete trueE2Hist;}

    char label1[1024];
    sprintf(label1, "Entries/(%.1f MeV)", maxE1*0.01);
    char label2[1024];
    sprintf(label2, "Entries/(%.1f MeV)", maxE2*0.01);

    // Histogram of the true energy for particle 1
    trueE1Hist = new TH1D("trueE1Hist", "", 100, 0.0, maxE1);
    trueE1Hist->SetDirectory(0);
    TString xLabel(particleName.c_str()); xLabel += " E (MeV)";
    trueE1Hist->SetXTitle(xLabel);
    trueE1Hist->SetYTitle(label1);
    trueE1Hist->SetLineStyle(kDashed);
    trueE1Hist->SetLineWidth(2);
    trueE1Hist->GetXaxis()->CenterTitle(true);
    trueE1Hist->GetYaxis()->CenterTitle(true);

    // Histogram of the true energy for particle 2 (proton)
    trueE2Hist = new TH1D("trueE2Hist", "", 100, 0.0, maxE2);
    trueE2Hist->SetDirectory(0);
    trueE2Hist->SetXTitle("Proton E (MeV)");
    trueE2Hist->SetYTitle(label2);
    trueE2Hist->SetLineStyle(kDashed);
    trueE2Hist->SetLineWidth(2);
    trueE2Hist->GetXaxis()->CenterTitle(true);
    trueE2Hist->GetYaxis()->CenterTitle(true);

    // Histogram of the reconstructed energy of particle 1
    recoE1Hist = new TH1D("recoE1Hist", "", 100, 0.0, maxE1);
    recoE1Hist->SetXTitle(xLabel);
    recoE1Hist->SetYTitle(label1);
    recoE1Hist->SetDirectory(0);
    recoE1Hist->GetXaxis()->CenterTitle(true);
    recoE1Hist->GetYaxis()->CenterTitle(true);

    // Histogram of the reconstructed energy of particle 2 (proton)
    recoE2Hist = new TH1D("recoE2Hist", "", 100, 0.0, maxE2);
    recoE2Hist->SetDirectory(0);
    recoE2Hist->SetXTitle("Proton E (MeV)");
    recoE2Hist->SetYTitle(label2);
    recoE2Hist->GetXaxis()->CenterTitle(true);
    recoE2Hist->GetYaxis()->CenterTitle(true);

    // Fractional difference between true and reconstructed energy for particle 1
    fracE1Hist = new TH1D("fracE1Hist", "", 100, 0.0, 0.0);
    TString fxLabel(particleName.c_str()); fxLabel += " E difference (%)";
    fracE1Hist->SetXTitle(fxLabel);
    fracE1Hist->SetDirectory(0);
    fracE1Hist->GetXaxis()->CenterTitle(true);
    fracE1Hist->GetYaxis()->CenterTitle(true);

    // Fractional difference between true and reconstructed energy for particle 2
    fracE2Hist = new TH1D("fracE2Hist", "", 50, -5.0, 5.0);
    fracE2Hist->SetXTitle("Proton E difference (%)");
    fracE2Hist->SetDirectory(0);
    fracE2Hist->GetXaxis()->CenterTitle(true);
    fracE2Hist->GetYaxis()->CenterTitle(true);

    // Open the lpc output file to get the cluster energies
    TFile* outFile = TFile::Open(outFileName.c_str(), "read");
    TTree* outTree = dynamic_cast<TTree*>(outFile->Get("clsTree"));
    int nEntries = outTree->GetEntries();

    int eventId, nCl, nClHits;
    std::vector<int>* clHitId = 0;
    std::vector<double>* clHitW = 0;
    outTree->SetBranchAddress("eventId", &eventId);
    outTree->SetBranchAddress("nCl", &nCl);
    outTree->SetBranchAddress("nClHits", &nClHits);
    outTree->SetBranchAddress("clHitId", &clHitId);
    outTree->SetBranchAddress("clHitW", &clHitW);

    // Also open the input file which will contain the primary energy/id info
    // Note that each entry in the input tree = eventId from the lpc output
    TFile* inFile = TFile::Open(inFileName.c_str(), "read");
    TTree* inTree = dynamic_cast<TTree*>(inFile->Get("Data"));
    int passed;
    double E0, pE;
    std::vector<int>* primeIdVect = 0;

    inTree->SetBranchAddress("E0", &E0);
    inTree->SetBranchAddress("pE", &pE);
    inTree->SetBranchAddress("passed", &passed);
    inTree->SetBranchAddress("primaryId", &primeIdVect);

    int i, iH;
    // Loop over all clusters in all events
    for (i = 0; i < nEntries; i++) {

	// Get the reconstructed cluster info
	outTree->GetEntry(i);
	
	// Get the corresponding primary event info
	inTree->GetEntry(eventId);

	// Do not process further if the event does not pass the
	// proton range cut
	if (passed == 0) {continue;}

	// Get the main truth id for this cluster. 
	// Loop over all hits and find the most common truth id
	int nId1(0), nId2(0);

	// Also find the reconstructed total cluster energy
	double recoE(0.0);

	for (iH = 0; iH < nClHits; iH++) {

	    // Get the hit index number
	    int hitId = (*clHitId)[iH];
	    // Find and store the primary id for this hit
	    int primeId = (*primeIdVect)[hitId];
	    
	    if (primeId == 1) {
		nId1++;
	    } else if (primeId == 2) {
		nId2++;
	    }

	    recoE += (*clHitW)[iH];

	}

	int mainId(0);
	if (nId1 > 0 && nId1 > nId2) {
	    mainId = 1;
	} else if (nId2 > 0 && nId2 > nId1) {
	    mainId = 2;
	}

	if (mainId == 1 && recoE > 0.0 && E0 > 0.0) {

	    // Electron or muon
	    double dEFrac = (recoE - E0)/E0;
	    recoE1Hist->Fill(recoE);
	    fracE1Hist->Fill(dEFrac*100.0);	    

	} else if (mainId == 2 && recoE > 0.0 && pE > 0.0) {

	    // Proton
	    double dEFrac = (recoE - pE)/pE;
	    recoE2Hist->Fill(recoE);
	    fracE2Hist->Fill(dEFrac*100.0);

	}

    }

    // Also create the true energy histograms
    int N = inTree->GetEntries();
    for (i = 0; i < N; i++) {

	inTree->GetEntry(i);

	if (passed == 1) {

	    trueE1Hist->Fill(E0);
	    trueE2Hist->Fill(pE);

	}

    }

    // Plot the histograms
    theCanvas->Clear();
    theCanvas->Divide(2,2);
    gStyle->SetOptFit(1111);
    gStyle->SetOptStat(0);
    gStyle->SetStatX(0.9);
    gStyle->SetStatY(0.9);
    theCanvas->UseCurrentStyle();

    theCanvas->cd(1);
    if (trueE1Hist->GetMaximum() > recoE1Hist->GetMaximum()) {
        trueE1Hist->Draw();
        recoE1Hist->Draw("same");
    } else {
        recoE1Hist->Draw();
        trueE1Hist->Draw("same");
    }
    gPad->Update();

    theCanvas->cd(2);
    if (trueE2Hist->GetMaximum() > recoE2Hist->GetMaximum()) {
        trueE2Hist->Draw();
        recoE2Hist->Draw("same");
    } else {
        recoE2Hist->Draw();
        trueE2Hist->Draw("same");
    }
    gPad->Update();

    theCanvas->cd(3);
    if (fun1) {delete fun1;}
    if (doFit1) {
	TF1* fun1 = new TF1("eFun", gaussQuad, -40.0, 40.0, 6);
	fun1->SetParameters(50.0, 0.0, 2.0, 0.0, 0.1, 0.0);
	fun1->SetParNames("N_{0}", "#mu", "#sigma", "a", "b", "c");    
	fracE1Hist->Fit(fun1, "", "", -40.0, 40.0);
    }
    fracE1Hist->Draw();
    gPad->Update();

    theCanvas->cd(4);
    fracE2Hist->Draw();
    gPad->Update();

    theCanvas->Print(plotFile.c_str());

    inFile->Close();
    outFile->Close();

}

double gaussQuad(double* xVar, double* par) {

    double x = xVar[0];

    double N0 = par[0];
    double mu = par[1];
    double sigma = par[2];

    double a = par[3];
    double b = par[4];
    double c = par[5];

    double qGauss = 0.0;

    if (N0 < 1e-10) {
        return 0.0;
    }

    double dx = x - mu;
    double arg = dx*dx/(sigma*sigma);

    if (fabs(arg) < 30.0) {
        qGauss = N0*exp(-0.5*arg);
    }

    qGauss += (a*x + b)*x + c;

    return qGauss;

}

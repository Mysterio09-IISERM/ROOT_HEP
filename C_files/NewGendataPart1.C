
    #include "TFile.h"
    #include "TTree.h"
    #include "TH1F.h"
    #include "TCanvas.h"
    #include "TROOT.h"
    #include "TRandom.h"
    #include "TStyle.h"
    #include "TF1.h"
    #include "TMath.h"

   // #include "RooRealVar.h"
    //#include "RooDataHist.h"
    //#include "RooBreitWigner.h"
    //#include "RooExponential.h"
    //#include "RooAddPdf.h"
    //#include "RooPlot.h"
    //#include "RooArgSet.h"
    //#include "RooGlobalFunc.h"

    ///using namespace RooFit;


void NewGendataPart1(){

    //UNIQUE ID
    TString uid =
        TString::Format(
            "%d",
            gRandom->Integer(10000000)
        );

    //OPEN ROOT FILE

    TFile *f = new TFile("../Root_files/mixdata.root", "READ");
    TTree *t = (TTree*)f->Get("data");

    //CREATE HISTOGRAM

    TH1F *hmass = new TH1F("hmass"+uid, "Mass Distribution; Mass [GeV]; Counts", 100, 60, 120);
    t->Draw("mass >> hmass"+uid, "", "goff");
    TH1F *henergy = new TH1F("henergy"+uid, "Energy Distribution; Energy [GeV]; Counts", 100, 0, 100);
    t->Draw("energy >> henergy"+uid, "", "goff");

    henergy->SetDirectory(0);
    hmass->SetDirectory(0);
    f->Close();
    delete f;

    // Create 2 Canvases to display the histograms

    TCanvas *cmass = new TCanvas("cmass"+uid, "Mass Distribution", 800, 600);
    TCanvas *cenergy = new TCanvas("cenergy"+uid, "Energy Distribution", 800, 600);
    
    cmass->Divide(3,1);
    cenergy->Divide(3,1);

    // Display raw mass distribution

    cmass->cd(1);
    
    TH1F *hrawmass = (TH1F*)hmass->Clone("hrawmass"+uid);
    hrawmass->SetLineColor(kBlue);
    hrawmass->Draw("E1");

    // Fitting mass distribution

    TF1 *ffitmass = new TF1("ffitmass"+uid, "breitwigner(0) + expo(3)", 60, 120);

    ffitmass->SetParameters(
        hmass->GetMaximum(),  // [0] amplitude
        91.2,                 // [1] mass initial guess
        2.5,                  // [2] width initial guess
        5.0,                  // [3] expo constant
        -0.05);               // [4] expo slope

    ffitmass->SetParLimits(0, 0, 100000);  // amplitude > 0
    ffitmass->SetParLimits(1, 85, 95);     // mass peak observed btw 85 and 95 GeV
    ffitmass->SetParLimits(2, 0.1, 10);    // width > 0

    cmass->cd(2);
    hmass->Fit(ffitmass, "RS0");
    hmass->SetTitle("Combined Fit; Mass [GeV]; Counts");
    hmass->SetLineColor(kBlue);
    hmass->Draw("E1");
    ffitmass->SetLineColor(kRed);
    ffitmass->Draw("same");

    TF1 *fbkgmass = new TF1("fbkgmass"+uid, "expo", 60, 120);
    fbkgmass->SetParameters(ffitmass->GetParameter(3), ffitmass->GetParameter(4));
    fbkgmass->SetLineColor(kGreen);
    fbkgmass->Draw("same");

    // Extract signal component
    cmass->cd(3);
    TH1F *hsigmass = (TH1F*)hmass->Clone("hsigmass"+uid);
    hsigmass->SetTitle("Extracted Signal Data; Mass [GeV]; Counts");
    for(int i=1; i<=hsigmass->GetNbinsX(); i++){
        Double_t x = hsigmass->GetBinCenter(i);
        Double_t bkg = fbkgmass->Eval(x); 
        hsigmass->SetBinContent(i, hsigmass->GetBinContent(i) - bkg);
        Double_t err = hsigmass->GetBinError(i);
        hsigmass->SetBinError(i, err); 
    }

    hsigmass->SetLineColor(kBlue);
    hsigmass->Draw("E1");
    TF1 *fsigmass = new TF1("fsigmass"+uid, "breitwigner", 60, 120);
    fsigmass->SetParameters(
        ffitmass->GetParameter(0),  // amplitude
        ffitmass->GetParameter(1),  // mass
        ffitmass->GetParameter(2)); // width
    fsigmass->SetLineColor(kRed);
    fsigmass->Draw("same");

    cmass->Update();


    //Display raw energy distribution

    cenergy->cd(1);
    TH1F *hrawenergy = (TH1F*)henergy->Clone("hrawenergy"+uid);
    hrawenergy->SetLineColor(kBlue);
    hrawenergy->Draw("E1");

    // Fitting energy distribution

    TF1 *ffitenergy = new TF1("ffitenergy"+uid, "gaus(0) + pol2(3)", 0, 100);
    ffitenergy->SetParameters(
        henergy->GetMaximum(),  // [0] amplitude
        50.0,                 // [1] mean initial guess
        10.0,                  // [2] sigma initial guess
        100.0,                // [3] pol2 constant term
        -1.0,                 // [4] pol2 linear term
        0.01);                // [5] pol2 quadratic term
    ffitenergy->SetParLimits(0, 0, 100000);  // amplitude > 0
    ffitenergy->SetParLimits(1, 40, 60);     //mean observed btw 40 and 60 GeV
    ffitenergy->SetParLimits(2, 0.1, 50); // sigma > 0

    cenergy->cd(2);
    henergy->Fit(ffitenergy, "RS0");
    henergy->SetTitle("Combined Fit; Energy [GeV]; Counts");
    henergy->SetLineColor(kBlue);
    henergy->Draw("E1");
    ffitenergy->SetLineColor(kRed);
    ffitenergy->Draw("same");
    TF1 *fbkgenergy = new TF1("fbkgenergy"+uid, "pol2", 0, 100);
    fbkgenergy->SetParameters(ffitenergy->GetParameter(3), ffitenergy->GetParameter(4), ffitenergy->GetParameter(5));
    fbkgenergy->SetLineColor(kGreen);
    fbkgenergy->Draw("same");

    //Extract signal component

    cenergy->cd(3);
    TH1F *hsigenergy = (TH1F*)henergy->Clone("hsigenergy"+uid);
    hsigenergy->SetTitle("Extracted Signal Data; Energy [GeV]; Counts");
    for(int i=1; i<=hsigenergy->GetNbinsX(); i++){
        Double_t x = hsigenergy->GetBinCenter(i);
        Double_t bkg = fbkgenergy->Eval(x); 
        hsigenergy->SetBinContent(i, hsigenergy->GetBinContent(i) - bkg);
        Double_t err = hsigenergy->GetBinError(i);
        hsigenergy->SetBinError(i, err); 
    }

    hsigenergy->SetLineColor(kBlue);
    hsigenergy->Draw("E1");
    TF1 *fsigenergy = new TF1("fsigenergy"+uid, "gaus", 0, 100);
    fsigenergy->SetParameters(
        ffitenergy->GetParameter(0),  // amplitude
        ffitenergy->GetParameter(1),  // mean
        ffitenergy->GetParameter(2)); // sigma
    fsigenergy->SetLineColor(kRed);
    fsigenergy->Draw("same");

    cenergy->Update();

    cmass->SaveAs("../Png_files/NewGendataPart1_mass.png");
    cenergy->SaveAs("../Png_files/NewGendataPart1_energy.png");

}
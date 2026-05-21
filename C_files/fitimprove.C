void fitimprove(){
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TRandom.h"
#include "TStyle.h"
#include "TF1.h"
#include "TMath.h"

#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooBreitWigner.h"
#include "RooExponential.h"
#include "RooAddPdf.h"
#include "RooPlot.h"
#include "RooArgSet.h"
    //clean up old histograms and canvases
    if (gROOT->FindObject("c1")) delete gROOT->FindObject("c1");
    if (gROOT->FindObject("c2")) delete gROOT->FindObject("c2");

    for (int i = 0; i < 4; i++) {
        TString fname = Form("fSig%d", i);
        if (gROOT->FindObject(fname)) delete gROOT->FindObject(fname);
        TString fbname = Form("fBkg%d", i);
        if (gROOT->FindObject(fbname)) delete gROOT->FindObject(fbname);
        TString hsname = Form("hSig%d", i);
        if (gROOT->FindObject(hsname)) delete gROOT->FindObject(hsname);
        TString hbname = Form("hBkg%d", i);
        if (gROOT->FindObject(hbname)) delete gROOT->FindObject(hbname);
    }


    if (gROOT->FindObject("outFile")) delete gROOT->FindObject("outFile");
    if (gROOT->FindObject("outFileB")) delete gROOT->FindObject("outFileB");
    if (gROOT->FindObject("fgauss")) delete gROOT->FindObject("fgauss");

    //open the file with the data and getting the trees

    TFile *F1 = new TFile("../Root_files/tmva_class_example.root", "READ");
    TTree *t1 = (TTree*)F1->Get("TreeS");
    TTree *t2 = (TTree*)F1->Get("TreeB");

    //create variables to hold the data and set the branches

    TString var[4] = {"var1", "var2", "var3", "var4"};

    //create histograms to hold the data

    TH1F *hSig[4], *hBkg[4];
    for ( int i=0; i<4 ; i++ ){

        //Signal

        TString hnameS=Form("hSig%d",i);
        TString htitleS=Form("Signal %s ; %s ; Count", var[i].Data(), var[i].Data());
        hSig[i] = new TH1F(hnameS, htitleS, 50, -10, 10);
        t1->Draw(var[i]+">>"+hnameS);
        hSig[i]->SetDirectory(0);

        //Background
        TString hnameB=Form("hBkg%d",i);
        TString htitleB=Form("Background %s ; %s ; Count", var[i].Data(), var[i].Data());
        hBkg[i] = new TH1F(hnameB, htitleB, 50, -10, 10);
        t2->Draw(var[i]+">>"+hnameB,"weight");
        hBkg[i]->SetDirectory(0);        

    }

    F1->Close();

    TCanvas *c1 = new TCanvas("c1", "Signal Data", 800, 600);
    TCanvas *c2 = new TCanvas("c2", "Background Data", 800, 600);

    c1->Divide(2,2);
    c2->Divide(2,2);

    for ( int i=0; i<4 ; i++ ){
        c1->cd(i+1);
        hSig[i]->Draw("HIST E1");
        c2->cd(i+1);
        hBkg[i]->Draw("HIST E1");
    }

    c1->Update();
    c1->SaveAs("../Png_files/SignalDataNoFit.png");
    c2->Update();
    c2->SaveAs("../Png_files/BackgroundDataNoFit.png");

    //output into .txt file for Signal fitting results

    ofstream outFile("../Txt_files/fitimprove_output_Signal.txt");

    //fitting Signal Data
    outFile << "= Signal Fitting Results =" << endl;

    for( int i=0; i<4;i++){

        c1->cd(i+1);

        TF1 *fgauss = new TF1(Form("fSig%d",i), "gaus", -10, 10);
        fgauss->SetParameters(hSig[i]->GetMaximum(), hSig[i]->GetMean(), hSig[i]->GetRMS());
        outFile << "\n--- Signal " << var[i] << " ---" << endl;
        outFile << "Initial: Amp=" << hSig[i]->GetMaximum()
                << " \nMean=" << hSig[i]->GetMean()
                << " \nSigma=" << hSig[i]->GetRMS() << endl;

        hSig[i]->Fit(fgauss, "RSQ");

        outFile << "\nFitted: Amp=" << fgauss->GetParameter(0)
                << " \nFitted Mean=" << fgauss->GetParameter(1)
                << " \nFitted Sigma=" << fgauss->GetParameter(2)
                << " \nFitted Chi^2/NDF=" << fgauss->GetChisquare() / fgauss->GetNDF()
                << " \nFit Probablity=" << fgauss->GetProb() << endl;
            

        fgauss->SetLineColor(kRed);
        fgauss->Draw("same");
    }
    c1->Update();
    c1->SaveAs("../Png_files/SignalDataWithFit.png");
    outFile.close();

    //output into .txt file for Background fitting results

    ofstream outFileB("../Txt_files/fitimprove_output_Background.txt");

    //fitting Weighted Background Data
    outFileB << "= Background Fitting Results =" << endl;

    for (int i = 0; i < 4; i++) {
    c2->cd(i+1);

    TF1 *fbkg;

    if (i == 0 || i == 3) {
    fbkg = new TF1(Form("fBkg%d", i), "gaus(0) + gaus(3)", -10, 10);
    fbkg->SetParameters(
        300,                   // amp1 -- narrow peak height
        hBkg[i]->GetMean(),   // mean1
        0.3,                  // sigma1 -- VERY narrow
        150,                  // amp2 -- wide base height
        hBkg[i]->GetMean(),   // mean2
        1.5);                 // sigma2 -- MUCH wider

    fbkg->SetParLimits(0, 0, 10000);   // amp1 > 0
    fbkg->SetParLimits(2, 0.05, 1.0);  // sigma1 constrained narrow
    fbkg->SetParLimits(3, 0, 10000);   // amp2 > 0
    fbkg->SetParLimits(5, 1.0, 5.0);   // sigma2 constrained wide
}
     else if (i == 1) {
        // var2 -- Gaussian works fine
        fbkg = new TF1(Form("fBkg%d", i), "gaus", -10, 10);
        fbkg->SetParameters(hBkg[i]->GetMaximum(),
                            hBkg[i]->GetMean(),
                            hBkg[i]->GetRMS());
} else {
    fbkg = new TF1(Form("fBkg%d", i), "gaus", -10, 10);
    fbkg->SetParameters(hBkg[i]->GetMaximum(),
                        hBkg[i]->GetMean(),
                        hBkg[i]->GetRMS());
}



    outFileB << "\n--- Background " << var[i] << " ---" << endl;
    outFileB << "Initial: Amp=" << hBkg[i]->GetMaximum()
             << "\nMean=" << hBkg[i]->GetMean()
             << "\nSigma=" << hBkg[i]->GetRMS() << endl;

// use WL for weighted histograms
hBkg[i]->Fit(fbkg, "RWSQ");

    outFileB << "Fitted Amp=" << fbkg->GetParameter(0)
             << "\nFitted Mean=" << fbkg->GetParameter(1)
             << "\nFitted Sigma=" << fbkg->GetParameter(2)
             << "\nFitted Chi2/NDF=" << fbkg->GetChisquare()/fbkg->GetNDF()
             << "\nFit Prob=" << fbkg->GetProb() << endl;

    fbkg->SetLineColor(kRed);
    fbkg->Draw("same");
}

    c2->Update();
    c2->SaveAs("../Png_files/BackgroundDataWithFit.png");
    outFileB.close();

}
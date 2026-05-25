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
    #include "RooGlobalFunc.h"
    #include "RooChebychev.h"
    #include <fstream>

    using namespace RooFit;

void Mixedfitdata5(){

    //OPEN ROOT FILE

    TFile *f = new TFile("../Root_files/mixed_fit_data.root", "READ");
    TTree *t = (TTree*)f->Get("data");

    //CREATE HISTOGRAMS

    TH1F *h1 = new TH1F("h1", "Generated Distribution; Var; Counts", 100, 10.0, 120.0);
    t->Draw("twopeak_bkg >>h1", "", "goff");
    h1->SetDirectory(0);
    f->Close();
    delete f;

    // ROOFIT VARIABLE

    RooRealVar var("twopeak_bkg", "Variable", 10.0, 120.0);

    // IMPORT DATA

    RooDataHist data("data", "dataset", var, Import(*h1));

    //SIGNAL PARAMETERS

    RooRealVar mean1("mean1", "Mean 1", 35.0, 30.0, 40.0);
    RooRealVar sigma1("sigma1", "Sigma 1", 2.5, 0.1, 10.0);
    RooRealVar mean2("mean2", "Mean 2", 75.0, 70.0, 80.0);
    RooRealVar width2("width2", "Width 2", 5.0, 0.1, 10.0);

    // BACKGROUND PARAMETERS

    RooRealVar tau("tau", "Tau", -0.5, -1.0, 0.0);

    // SIGNAL AND BACKGROUND PDFS

    RooGaussian signal1("signal1", "Signal PDF", var, mean1, sigma1);
    RooBreitWigner signal2("signal2", "Signal PDF 2", var, mean2, width2);
    RooExponential background("background", "Background PDF", var, tau);

    // YEILD VARIABLES

    RooRealVar nsig1("nsig1", "Signal Yield 1", 50000, 10000, 150000);
    RooRealVar nsig2("nsig2", "Signal Yield 2", 50000, 10000, 150000);
    RooRealVar nbkg("nbkg", "Background Yield", 50000, 10000, 150000);

    // COMBINED MODEL

    RooAddPdf model("model", "Signal + Background Model", RooArgList(signal1, signal2, background), RooArgList(nsig1, nsig2, nbkg));

    //FIT MODEL TO DATA

    model.fitTo(data);

    //CANVAS CREATION

    TCanvas *cmain = new TCanvas("cmain", "Fit Result", 1800, 800);

    cmain->Divide(3,1);

    // RAW DISTRIBUTION

    cmain->cd(1);
    RooPlot *frame1 = var.frame();
    data.plotOn(frame1, Name("data"));
    frame1->SetTitle("Raw Distribution");
    frame1->Draw();

    TLegend *leg1 = new TLegend(0.6, 0.85, 0.9, 0.9);
    leg1->AddEntry(frame1->findObject("data"), "Data", "ep");
    leg1->Draw();


    // FITTED DISTRIBUTION

    cmain->cd(2);
    RooPlot *frame2 = var.frame();
    data.plotOn(frame2, Name("data2"));
    model.plotOn(frame2, Name("fullfit"));
    model.plotOn(frame2, Components(background), LineColor(kGreen), Name("bkg"));
    model.plotOn(frame2, Components(signal1), LineColor(kRed), Name("sig1"));
    model.plotOn(frame2, Components(signal2), LineColor(kMagenta), Name("sig2"));
    frame2->SetTitle("Fitted Distribution");
    frame2->Draw();

    TLegend *leg2 = new TLegend(0.6, 0.65, 0.9, 0.9);
    leg2->AddEntry(frame2->findObject("data2"), "Data", "ep");
    leg2->AddEntry(frame2->findObject("fullfit"), "Full Fit","l");
    leg2->AddEntry(frame2->findObject("sig1"), "Signal PDF 1", "l");
    leg2->AddEntry(frame2->findObject("sig2"), "Signal PDF 2", "l");
    leg2->AddEntry(frame2->findObject("bkg"), "Background PDF", "l");
    leg2->Draw();

    // EXTRACTED SIGNAL COMPONENT

    cmain->cd(3);
    TH1F *hsig = (TH1F*)h1->Clone("hsig");
    hsig->SetTitle("Extracted Signal Data; Var; Counts");

    RooArgSet normSet(var);

    for (int i=1; i<=hsig->GetNbinsX(); i++){
        double x = hsig->GetBinCenter(i);
        double bw = hsig->GetBinWidth(i);
        var.setVal(x);
        double bkg_density = background.getVal(normSet);
        double bkg_counts = bkg_density * nbkg.getVal() * bw;
        double y = hsig->GetBinContent(i);
        hsig->SetBinContent(i, y - bkg_counts);
    }

    hsig->SetLineColor(kBlue);
    hsig->SetMarkerStyle(20);
    

    TF1 *fsig1 = new TF1("fsig1", "gaus", 10.0, 50.0);
    double scale1 = nsig1.getVal() * hsig->GetBinWidth(1)/ (sigma1.getVal() * sqrt(2 * TMath::Pi()));
    double scale2 = nsig2.getVal() * hsig->GetBinWidth(1);
    fsig1->SetParameters(scale1, mean1.getVal(), sigma1.getVal());
    TF1 *fsig2 = new TF1("fsig2", "breitwigner", 50.0, 120.0);
    fsig2->SetParameters(scale2, mean2.getVal(), width2.getVal());
    TF1 *ftotal = new TF1("ftotal", "gaus(0) + breitwigner(3)", 10.0, 120.0);
    ftotal->SetParameters(scale1, mean1.getVal(), sigma1.getVal(), scale2, mean2.getVal(), width2.getVal());
    fsig1->SetLineColor(kRed);
    fsig2->SetLineColor(kMagenta);
    ftotal->SetLineColor(kBlack);
    fsig1->SetLineStyle(2);
    fsig2->SetLineStyle(2);
    ftotal->SetLineStyle(0);
    gStyle->SetOptStat(0);
    hsig->Draw("E1");
    ftotal->Draw("same");
    fsig1->Draw("same");
    fsig2->Draw("same");
   

    TLegend *leg3 = new TLegend(0.6, 0.75, 0.9, 0.9);
    leg3->AddEntry(hsig, "Signal Data", "ep");
    leg3->AddEntry(ftotal, "Fitted Total Signal", "l");
    leg3->AddEntry(fsig1, "Fitted Signal 1", "l");
    leg3->AddEntry(fsig2, "Fitted Signal 2", "l");
    leg3->Draw();

    cmain->Modified();
    cmain->Update();
    gPad->Update();

    cmain->SaveAs("../Png_files/mixed_fit5.png");

    ofstream out("../Txt_files/mixed_fit5_results.txt");

    out << "Fit Results:" << endl;
    out << "Signal Yield 1: " << nsig1.getVal() << " ± " << nsig1.getError() << endl;
    out << "Mean1: " << mean1.getVal() << " ± " << mean1.getError() << endl;
    out << "Width1: " << sigma1.getVal() << " ± " << sigma1.getError() << endl;
    out << "Signal Yield 2: " << nsig2.getVal() << " ± " << nsig2.getError() << endl;
    out << "Mean2: " << mean2.getVal() << " ± " << mean2.getError() << endl;
    out << "Width2: " << width2.getVal() << " ± " << width2.getError() << endl;
    out << "Background Yield: " << nbkg.getVal() << " ± " << nbkg.getError() << endl;

    out.close();
}
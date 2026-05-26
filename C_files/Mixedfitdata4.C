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
    #include "RooCBShape.h"
    #include <fstream>

    using namespace RooFit;

void Mixedfitdata4(){

    //OPEN ROOT FILE

    TFile *f = new TFile("../Root_files/mixed_fit_data.root", "READ");
    TTree *t = (TTree*)f->Get("data");

    //CREATE HISTOGRAMS

    TH1F *h1 = new TH1F("h1", "Generated Distribution; Var; Counts", 100, 0.0, 100.0);
    t->Draw("asym_flat >>h1", "", "goff");
    h1->SetDirectory(0);
    f->Close();
    delete f;

    // ROOFIT VARIABLE

    RooRealVar var("asym_flat", "Variable", 0.0, 100.0);

    // IMPORT DATA

    RooDataHist data("data", "dataset", var, Import(*h1));

    //SIGNAL PARAMETERS

    RooRealVar mean("mean", "Mean", 60.0, 55.0, 65.0);
    RooRealVar sigma("sigma", "Sigma", 2.5, 0.1, 10.0);
    RooRealVar alpha("alpha", "Alpha", 2.0, 0.1, 5.0);
    RooRealVar n("n", "N", 5.0, 0.1, 20.0);

    // BACKGROUND PARAMETERS

    // (No additional parameters needed for the Chebychev background in this case)

    // SIGNAL AND BACKGROUND PDFS

    RooCBShape signal("signal", "Signal PDF", var, mean, sigma, alpha, n);
    RooPolynomial background("background", "Background PDF", var, RooArgList());

    // YEILD VARIABLES

    RooRealVar nsig("nsig", "Signal Yield", 75000, 10000, 150000);
    RooRealVar nbkg("nbkg", "Background Yield", 75000, 10000, 150000);

    // COMBINED MODEL

    RooAddPdf model("model", "Signal + Background Model", RooArgList(signal, background), RooArgList(nsig, nbkg));

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
    model.plotOn(frame2, Components(signal), LineColor(kRed), Name("sig"));
    frame2->SetTitle("Fitted Distribution");
    frame2->Draw();

    TLegend *leg2 = new TLegend(0.6, 0.65, 0.9, 0.9);
    leg2->AddEntry(frame2->findObject("data2"), "Data", "ep");
    leg2->AddEntry(frame2->findObject("fullfit"), "Full Fit","l");
    leg2->AddEntry(frame2->findObject("sig"), "Signal PDF", "l");
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
    
/*
    TF1 *fsig = new TF1("fsig", "crystalball", 0.0, 100.0);
    double scale = nsig.getVal() * hsig->GetBinWidth(1);
    fsig->SetParameters(scale, mean.getVal(), sigma.getVal(), alpha.getVal(), n.getVal());
    fsig->SetLineColor(kRed);
    fsig->SetLineStyle(2);
    gStyle->SetOptStat(0);
    hsig->Draw("E1");
    fsig->Draw("same");
   

    TLegend *leg3 = new TLegend(0.6, 0.75, 0.9, 0.9);
    leg3->AddEntry(hsig, "Signal Data", "ep");
    leg3->AddEntry(fsig, "Fitted Signal", "l");
    leg3->Draw();

    cmain->Modified();
    cmain->Update();
    gPad->Update();*/

    RooPlot *frame3 = var.frame();

    RooDataHist sigData("sigData",
                        "Background Subtracted",
                        var,
                        Import(*hsig));

    sigData.plotOn(frame3,
                Name("sigdata"));

    signal.plotOn(frame3,
                LineColor(kRed),
                LineStyle(kDashed),
                Normalization(nsig.getVal(), RooAbsReal::NumEvent),
                Name("sigpdf"));

    frame3->SetTitle("Extracted Signal");
    frame3->Draw();

    cmain->SaveAs("../Png_files/mixed_fit4.png");

    ofstream out("../Txt_files/mixed_fit4_results.txt");

    out << "Fit Results:" << endl;
    out << "Signal Yield: " << nsig.getVal() << " ± " << nsig.getError() << endl;
    out << "Mean: " << mean.getVal() << " ± " << mean.getError() << endl;
    out << "Width: " << sigma.getVal() << " ± " << sigma.getError() << endl;
    out << "Alpha: " << alpha.getVal() << " ± " << alpha.getError() << endl;
    out << "N: " << n.getVal() << " ± " << n.getError() << endl;
    out << "Background Yield: " << nbkg.getVal() << " ± " << nbkg.getError() << endl;

    out.close();
}
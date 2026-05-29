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

void Newdata2(){

    //OPEN ROOT FILE

    TFile *f = new TFile("../Root_files/Newdata.root", "READ");
    TTree *t = (TTree*)f->Get("data");

    //CREATE HISTOGRAMS

    TH1F *h1 = new TH1F("h1", "Generated Distribution; Var; Counts", 800, 20.0, 80.0);
    t->Draw("overlap_signal >>h1", "", "goff");
    h1->SetDirectory(0);
    f->Close();
    delete f;

    // ROOFIT VARIABLE

    RooRealVar var("overlap_signal", "Variable", 20.0, 80.0);

    // IMPORT DATA

    RooDataHist data("data", "dataset", var, Import(*h1));

    //SIGNAL PARAMETERS

    RooRealVar mean1("mean1", "Mean", 50.0, 45.0, 55.0);
    RooRealVar sigma1("sigma1", "Sigma", 0.5, 0.1, 5.0);
    RooRealVar mean2("mean2", "Mean", 55.0, 50.0, 60.0);
    RooRealVar sigma2("sigma2", "Sigma", 2.5, 1.0, 7.0);

    RooRealVar frac("frac", "Fraction of Signal 1", 0.5, 0.0, 1.0);

    // BACKGROUND PARAMETERS

    RooRealVar tau("tau", "Expo Coeff", -0.05, -1.0, -0.001);

    // SIGNAL AND BACKGROUND PDFS

    RooGaussian signal1("signal1", "Signal 1 PDF", var, mean1, sigma1);
    RooGaussian signal2("signal2", "Signal 2 PDF", var, mean2, sigma2);

    RooAddPdf signal("signal", "Signal PDF", RooArgList(signal1, signal2),frac);
    RooExponential background("background", "Background PDF", var, tau);

    // YEILD VARIABLES

    RooRealVar nsig("nsig", "Signal Yield", 100000, 10000, 500000);
    RooRealVar nbkg("nbkg", "Background Yield", 100000, 10000, 500000);

    // COMBINED MODEL

    RooAddPdf model("model", "Signal + Background Model", RooArgList(signal, background), RooArgList(nsig, nbkg));

    //FIT MODEL TO DATA

    model.fitTo(data);

    //CANVAS CREATION

    TCanvas *cmain = new TCanvas("cmain", "Fit Result", 1600, 1000);

    cmain->Divide(3,1);

    // RAW DISTRIBUTION

    cmain->cd(1);
    RooPlot *frame1 = var.frame();
    data.plotOn(frame1, Name("data"));
    frame1->SetTitle("Raw Distribution");
    frame1->Draw();

    TLegend *leg1 = new TLegend(0.1, 0.85, 0.4, 0.9);
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
        double sig_counts = y - bkg_counts;
        hsig->SetBinContent(i, sig_counts >0 ? sig_counts : 0);
    }

    RooDataHist signalData("signalData", "Signal Data", var, Import(*hsig));
    RooPlot *frame3 = var.frame();
    signalData.plotOn(frame3, Name("sigdata"));
    signal.plotOn(frame3, Name("sigpdf"), LineColor(kRed), LineStyle(kDashed), Normalization(hsig->Integral(),RooAbsReal::NumEvent));
    signal1.plotOn(frame3, Name("sig1"), LineColor(kBlue), LineStyle(kDotted), Normalization(hsig->Integral() * frac.getVal(),RooAbsReal::NumEvent));
    signal2.plotOn(frame3, Name("sig2"), LineColor(kMagenta), LineStyle(kDotted), Normalization(hsig->Integral() * (1 - frac.getVal()),RooAbsReal::NumEvent));
    frame3->SetTitle("Extracted Signal Component");
    frame3->Draw();

    TLegend *leg3 = new TLegend(0.6, 0.75, 0.9, 0.9);
    leg3->AddEntry(frame3->findObject("sigdata"), "Signal Data", "ep");
    leg3->AddEntry(frame3->findObject("sigpdf"), "Fitted Signal", "l");
    leg3->AddEntry(frame3->findObject("sig1"), "Signal 1 Component", "l");
    leg3->AddEntry(frame3->findObject("sig2"), "Signal 2 Component", "l");
    leg3->Draw();

    cmain->Modified();
    cmain->Update();
    gPad->Update();

    cmain->SaveAs("../Png_files/Newdata2.png");

    ofstream out("../Txt_files/Newdata2_results.txt");

    out << "Fit Results:" << endl;
    out << "Signal Yield: " << nsig.getVal() << " ± " << nsig.getError() << endl;
    out << "Mean1: " << mean1.getVal() << " ± " << mean1.getError() << endl;
    out << "Mean2: " << mean2.getVal() << " ± " << mean2.getError() << endl;
    out << "Sigma1: " << sigma1.getVal() << " ± " << sigma1.getError() << endl;
    out << "Sigma2: " << sigma2.getVal() << " ± " << sigma2.getError() << endl;
    out << "Background Yield: " << nbkg.getVal() << " ± " << nbkg.getError() << endl;
    out << "Fraction of Signal 1: " << frac.getVal() << " ± " << frac.getError() << endl;

    out.close();
}
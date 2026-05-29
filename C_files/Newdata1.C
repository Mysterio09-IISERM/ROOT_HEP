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

void Newdata1(){

    //OPEN ROOT FILE

    TFile *f = new TFile("../Root_files/Newdata.root", "READ");
    TTree *t = (TTree*)f->Get("data");

    //CREATE HISTOGRAMS

    TH1F *h1 = new TH1F("h1", "Generated Distribution; Var; Counts", 800, 80.0, 100.0);
    t->Draw("narrow_signal >>h1", "", "goff");
    h1->SetDirectory(0);
    f->Close();
    delete f;

    // ROOFIT VARIABLE

    RooRealVar var("narrow_signal", "Variable", 80.0, 100.0);

    // IMPORT DATA

    RooDataHist data("data", "dataset", var, Import(*h1));

    //SIGNAL PARAMETERS

    RooRealVar mean("mean", "Mean", 91.0, 89.0, 93.0);
    RooRealVar sigma("sigma", "Sigma", 0.5, 0.1, 3.0);

    // BACKGROUND PARAMETERS

    //RooRealVar tau("tau", "Expo Coeff", -0.5, -1.0, 0.0); Initial guess of exponential slope
    RooRealVar c0("c0", "Cheby Coeff 0", 0.5, -1.0, 1.0);
    RooRealVar c1("c1", "Cheby Coeff 1", 0.0, -1.0, 1.0);

    // SIGNAL AND BACKGROUND PDFS

    RooGaussian signal("signal", "Signal PDF", var, mean, sigma);
    //RooExponential background("background", "Background PDF", var, tau); Initial exponential background model
    RooChebychev background("background", "Background PDF", var, RooArgList(c0, c1));

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
    frame3->SetTitle("Extracted Signal Component");
    frame3->Draw();

    TLegend *leg3 = new TLegend(0.6, 0.75, 0.9, 0.9);
    leg3->AddEntry(frame3->findObject("sigdata"), "Signal Data", "ep");
    leg3->AddEntry(frame3->findObject("sigpdf"), "Fitted Signal", "l");
    leg3->Draw();

    cmain->Modified();
    cmain->Update();
    gPad->Update();

    cmain->SaveAs("../Png_files/Newdata1.png");

    ofstream out("../Txt_files/Newdata1_results.txt");

    out << "Fit Results:" << endl;
    out << "Signal Yield: " << nsig.getVal() << " ± " << nsig.getError() << endl;
    out << "Mean: " << mean.getVal() << " ± " << mean.getError() << endl;
    out << "Sigma: " << sigma.getVal() << " ± " << sigma.getError() << endl;
    out << "Background Yield: " << nbkg.getVal() << " ± " << nbkg.getError() << endl;

    out.close();

    cout << "Total entries in h1: " << h1->GetEntries() << endl;
    cout << "Total entries in range: " << h1->Integral() << endl;
}
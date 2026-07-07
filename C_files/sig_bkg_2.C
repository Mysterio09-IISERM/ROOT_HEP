#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TRandom.h"
#include "TStyle.h"
#include "TF1.h"
#include "TMath.h"
#include "TSystem.h"
#include "TLine.h"
#include "TLegend.h"

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

void sig_bkg_2(){

    //OPEN ROOT FILE

    gSystem->RedirectOutput("../logs/sig_bkg_2.log", "w");

    TFile *f = new TFile("../Root_files/sig_bkg_2.root", "READ");
    TTree *t = (TTree*)f->Get("tree");

    // ROOFIT VARIABLE

    RooRealVar var("x", "Variable", 0.0, 100.0);

    // IMPORT DATA

    RooDataSet data("data", "dataset",RooArgList(var), Import(*t));
    f->Close();
    delete f;

    cout << "Entries: " << data.numEntries() << endl;

    //SIGNAL PARAMETERS

    RooRealVar mean("mean", "Mean", 30.0, 0.0, 50.0);
    RooRealVar sigma("sigma", "Sigma", 0.7, 0.1, 100.0);
    mean.setVal(44.984);
    sigma.setVal(3.0528);
    mean.setConstant(kTRUE);
    sigma.setConstant(kTRUE);

    // BACKGROUND PARAMETERS

    RooRealVar tau("tau", "Expo Coeff", -0.05, -5.0, 0.0);
    tau.setVal(-0.037282);
    tau.setConstant(kTRUE);

    // SIGNAL AND BACKGROUND PDFS

    RooGaussian signal("signal", "Signal PDF", var, mean, sigma);
    RooExponential background("background", "Background PDF", var, tau);

    // YIELD VARIABLES

    RooRealVar nsig("nsig", "Signal Yield", 10000, 0, 500000);
    RooRealVar nbkg("nbkg", "Background Yield", 10000, 0, 500000);

    // COMBINED MODEL

    RooAddPdf model("model", "Signal + Background Model", RooArgList(signal, background), RooArgList(nsig, nbkg));

    //FIT MODEL TO DATA

    RooFitResult* fitres = model.fitTo(data, Save(), Extended(kTRUE), PrintLevel(1));

    int nbins = (int)sqrt(data.numEntries());
    cout << "Using " << nbins << " bins for display" << endl;

    // ============================================================
    // CANVAS CREATION -- single column: main fit (top) + pull (bottom)
    // ============================================================

    TCanvas *cmain1 = new TCanvas("cmain1", "Fit Result", 900, 700);

    TPad *pad1 = new TPad("pad1", "", 0, 0.28, 1, 1);
    pad1->SetBottomMargin(0.02);
    pad1->Draw();

    TPad *pad2 = new TPad("pad2", "", 0, 0, 1, 0.28);
    pad2->SetTopMargin(0.06);
    pad2->SetBottomMargin(0.35);
    pad2->Draw();

    // ---------------- MAIN FIT PLOT ----------------

    pad1->cd();
    RooPlot *frame2 = var.frame(Bins(nbins));

    data.plotOn(frame2, Name("data2"), Binning(nbins));

    model.plotOn(frame2, Components(background), LineColor(kGreen+2),
                 LineStyle(kDashed), Name("bkg"));
    model.plotOn(frame2, Components(signal), LineColor(kRed),
                 LineWidth(2), Name("sig"));
    model.plotOn(frame2, LineColor(kBlue), LineWidth(2), Name("fullfit"));

    data.plotOn(frame2, Name("data"), Binning(nbins));

    double nparams = fitres->floatParsFinal().getSize();
    double chi2ndf = frame2->chiSquare("fullfit", "data2", nparams);

    frame2->SetTitle("B^{+} Mass Fit");
    frame2->GetYaxis()->SetTitle("Events / ( 0.0025 )");
    frame2->GetXaxis()->SetLabelSize(0);
    frame2->GetXaxis()->SetTitleSize(0);
    frame2->Draw();

    TLegend *leg2 = new TLegend(0.62, 0.65, 0.88, 0.88);
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);
    leg2->AddEntry(frame2->findObject("data"), "Data", "ep");
    leg2->AddEntry(frame2->findObject("fullfit"), "Full Fit", "l");
    leg2->AddEntry(frame2->findObject("sig"), "Signal PDF", "l");
    leg2->AddEntry(frame2->findObject("bkg"), "Background PDF", "l");
    leg2->Draw();

    // ---------------- PULL DISTRIBUTION (bottom) ----------------

    pad2->cd();
    RooPlot *frame3 = var.frame(Bins(nbins));
    frame3->SetMinimum(-5);
    frame3->SetMaximum(5);

    RooHist *pullHist = frame2->pullHist("data2", "fullfit");

    // Convert pull points into a filled bar-style histogram (like a TH1)
    TH1F *hpullbar = new TH1F("hpullbar", "", nbins, var.getMin(), var.getMax());
    for (int i = 0; i < pullHist->GetN(); i++) {
        double x, y;
        pullHist->GetPoint(i, x, y);
        hpullbar->SetBinContent(hpullbar->FindBin(x), y);
    }
    hpullbar->SetFillColor(kAzure+7);
    hpullbar->SetLineColor(kAzure+7);
    hpullbar->SetLineWidth(1);
    hpullbar->SetTitle(";Mass (GeV/c^{2});Pull");
    hpullbar->GetYaxis()->SetTitle("Pull");
    hpullbar->GetYaxis()->SetNdivisions(505);
    hpullbar->GetXaxis()->SetTitle("Mass (GeV/c^{2})");
    hpullbar->SetMinimum(-5);
    hpullbar->SetMaximum(5);
    hpullbar->SetStats(0);
    hpullbar->Draw("HIST");

    TLine *zero = new TLine(var.getMin(), 0, var.getMax(), 0);
    zero->SetLineColor(kBlack);
    zero->Draw("same");

    cmain1->Modified();
    cmain1->Update();
    gPad->Update();

    cmain1->SaveAs("../Png_files/sig_bkg_2_results.png");

    //Finding Upper Limit of signal data

    cout << "Before PLC: " << nsig.getVal() << endl;

    RooArgSet poi(nsig);
    RooStats::ProfileLikelihoodCalculator plc(data, model, poi);
    plc.SetConfidenceLevel(0.90);
    RooStats::LikelihoodInterval *interval = plc.GetInterval();
    double upper_limit = interval->UpperLimit(nsig);
    cout << "Upper limit on nsig at 90% CL: "
         << upper_limit << " events" << endl;

    TCanvas *c = new TCanvas("c", "Upper Limit", 800, 600);
    RooStats::LikelihoodIntervalPlot plot(interval);
    plot.SetRange(0, upper_limit * 2);
    plot.Draw();

    TLine *cl_line = new TLine(upper_limit, 0, upper_limit, 5);
    cl_line->SetLineColor(kRed);
    cl_line->SetLineStyle(2);
    cl_line->Draw("same");

    cout << "\n=== Upper Limit Result ===" << endl;
    cout << "nsig upper limit (90% CL) = " << upper_limit << endl;
    cout << "After PLC: " << nsig.getVal() << endl;

    delete interval;

    ofstream out("../Txt_files/sig_bkg_2_results.txt");

    out << "Fit Results:\n" << endl;

    out << "Status: " << fitres->status() << endl;
    out << "CovQual: " << fitres->covQual() << endl;
    out << "EDM: " << fitres->edm() << endl;
    out << "MinNLL: " << fitres->minNll() << endl;
    out << "Chi2/NDF: " << chi2ndf << endl;

    fitres->Print("v");
    fitres->correlationMatrix().Print();

    out << "\nSignal Yield: " << nsig.getVal() << " ± " << nsig.getError() << endl;
    out << "Mean: " << mean.getVal() << " ± " << mean.getError() << endl;
    out << "Sigma: " << sigma.getVal() << " ± " << sigma.getError() << endl;
    out << "\nBackground Yield: " << nbkg.getVal() << " ± " << nbkg.getError() << endl;
    out << "Tau: " << tau.getVal() << " ± " << tau.getError() << endl;

    out.close();

    gSystem->RedirectOutput(0);
}
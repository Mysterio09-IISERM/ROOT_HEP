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

void Complexdata3(){

    //OPEN ROOT FILE

    gSystem->RedirectOutput("../logs/Complexdata3.log", "w");

    TFile *f = new TFile("../Root_files/ComplexData.root", "READ");
    TTree *t = (TTree*)f->Get("data");

    //CREATE HISTOGRAMS

    TH1F *h1 = new TH1F("h1", "Generated Distribution; Var; Counts", 1000, 3000.0, 4000.0);
    t->Draw("var3 >>h1", "weight3", "goff");
    h1->SetDirectory(0);
    f->Close();
    delete f;

    // ROOFIT VARIABLE

    RooRealVar var("var3", "Variable", 3000.0, 4000.0);

    // IMPORT DATA

    RooDataHist data("data", "dataset", var, Import(*h1));

    //SIGNAL PARAMETERS

    RooRealVar mean1("mean1", "Mean1", 3100.0, 3060.0, 3140.0);
    //RooRealVar width1("width1", "Width1", 40.0, 0.1, 60.0);
    RooRealVar sigma1("sigma1", "Sigma1", 0.5, 0.1, 100.0);
    RooRealVar mean2("mean2", "Mean2", 3690.0, 3650.0, 3730.0);
    //RooRealVar width2("width2", "Width2", 20.0, 0.1, 50.0);
    RooRealVar sigma2("sigma2", "Sigma2", 0.5, 0.1, 100.0);
    RooRealVar mean3("mean3", "Mean3", 3870.0, 3830.0, 3910.0);
    //RooRealVar width3("width3", "Width3", 15.0, 0.1, 40.0);
    RooRealVar sigma3("sigma3", "Sigma3", 0.5, 0.1, 100.0);

    // BACKGROUND PARAMETERS

    //RooRealVar c0("c0", "c0", 0.5, -1.0, 1.0);
    //RooRealVar c1("c1", "c1", -0.5, -1.0, 1.0);
    RooRealVar a1("a1", "a1", 0.5, 0.0, 10.0);
    RooRealVar a2("a2", "a2", 0.0, 0.0, 10.0);

    // SIGNAL AND BACKGROUND PDFS

    RooGaussian signal1("signal1", "Signal 1", var, mean1, sigma1);
    RooGaussian signal2("signal2", "Signal 2", var, mean2, sigma2);
    RooGaussian signal3("signal3", "Signal 3", var, mean3, sigma3);
    //RooChebychev background("background", "Background", var, RooArgList(c0, c1));
    RooBernstein background("background", "Background", var, RooArgList(a1, a2));

    // YEILD VARIABLES

    RooRealVar nsig1("nsig1", "Signal Yield 1", 100000, 0, 500000);
    RooRealVar nsig2("nsig2", "Signal Yield 2", 100000, 0, 500000);
    RooRealVar nsig3("nsig3", "Signal Yield 3", 100000, 0, 500000);
    RooRealVar nbkg("nbkg", "Background Yield", 100000, 0, 500000);

    // COMBINED MODEL

    RooAddPdf model("model", "Signal + Background Model", RooArgList(signal1, signal2, signal3, background), RooArgList(nsig1, nsig2, nsig3, nbkg));

    //FIT MODEL TO DATA

    RooFitResult* fitres = model.fitTo(data, Save(), SumW2Error(kTRUE));

    //CANVAS CREATION

    TCanvas *cmain = new TCanvas("cmain", "Fit Result", 1600, 900);

    cmain->Divide(3,1);

    // RAW DISTRIBUTION

    cmain->cd(1);
    RooPlot *frame1 = var.frame();
    data.plotOn(frame1, Name("data"),DataError(RooAbsData::SumW2));
    frame1->SetTitle("Raw Distribution");
    frame1->Draw();

    TLegend *leg1 = new TLegend(0.1, 0.85, 0.4, 0.9);
    leg1->AddEntry(frame1->findObject("data"), "Data", "ep");
    leg1->Draw();


    // FITTED DISTRIBUTION

    cmain->cd(2);
    RooPlot *frame2 = var.frame();
    data.plotOn(frame2, Name("data2"), DataError(RooAbsData::SumW2));
    model.plotOn(frame2, Name("fullfit"));
    model.plotOn(frame2, Components(background), LineColor(kGreen), Name("bkg"));
    model.plotOn(frame2, Components(RooArgList(signal1, signal2, signal3)), LineColor(kRed), Name("sig"));
    double chi2ndf = frame2->chiSquare("fullfit", "data2", 12);
    frame2->SetTitle("Fitted Distribution");
    frame2->Draw();

    TLegend *leg2 = new TLegend(0.3, 0.65, 0.6, 0.9);
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
    //signal.plotOn(frame3, Name("sigpdf"), LineColor(kRed), LineStyle(kDashed), Normalization(hsig->Integral(),RooAbsReal::NumEvent));
    signal1.plotOn(frame3, Name("sig1"), LineColor(kBlue), LineStyle(kDotted), Normalization(nsig1.getVal(),RooAbsReal::NumEvent));
    signal2.plotOn(frame3, Name("sig2"), LineColor(kGreen), LineStyle(kDashDotted), Normalization(nsig2.getVal(),RooAbsReal::NumEvent));
    signal3.plotOn(frame3, Name("sig3"), LineColor(kMagenta), LineStyle(kSolid), Normalization(nsig3.getVal(),RooAbsReal::NumEvent));
    frame3->SetTitle("Extracted Signal Component");
    frame3->Draw();

    TLegend *leg3 = new TLegend(0.3, 0.75, 0.6, 0.9);
    leg3->AddEntry(frame3->findObject("sigdata"), "Signal Data", "ep");
    //leg3->AddEntry(frame3->findObject("sigpdf"), "Fitted Signal", "l");
    leg3->AddEntry(frame3->findObject("sig1"), "Signal 1", "l");
    leg3->AddEntry(frame3->findObject("sig2"), "Signal 2", "l");
    leg3->AddEntry(frame3->findObject("sig3"), "Signal 3", "l");
    leg3->Draw();

    cmain->Modified();
    cmain->Update();
    gPad->Update();

    cmain->SaveAs("../Png_files/ComplexData3.png");

    ofstream out("../Txt_files/ComplexData3_results.txt");

    out << "Fit Results:\n" << endl;

    out << "Status: " << fitres->status() << endl;
    out << "CovQual: " << fitres->covQual() << endl;
    out << "EDM: " << fitres->edm() << endl;
    out << "MinNLL: " << fitres->minNll() << endl;
    out << "Chi2/NDF: " << chi2ndf << endl;

    fitres->Print("v");

    out << "\nSignal 1 Yield: " << nsig1.getVal() << " ± " << nsig1.getError() << endl;
    out << "Signal 2 Yield: " << nsig2.getVal() << " ± " << nsig2.getError() << endl;
    out << "Signal 3 Yield: " << nsig3.getVal() << " ± " << nsig3.getError() << endl;
    out << "Mean1: " << mean1.getVal() << " ± " << mean1.getError() << endl;
    out << "Mean2: " << mean2.getVal() << " ± " << mean2.getError() << endl;
    out << "Mean3: " << mean3.getVal() << " ± " << mean3.getError() << endl;
    out << "Sigma1: " << sigma1.getVal() << " ± " << sigma1.getError() << endl;
    out << "Sigma2: " << sigma2.getVal() << " ± " << sigma2.getError() << endl;
    out << "Sigma3: " << sigma3.getVal() << " ± " << sigma3.getError() << endl;

    out << "Background Yield: " << nbkg.getVal() << " ± " << nbkg.getError() << endl;

    out.close();

    gSystem->RedirectOutput(0);
}
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
    #include "RooGaussian.h"
    #include "RooGenericPdf.h"
    #include "RooFitResult.h"
    #include "TPad.h"
    #include "TLatex.h"
    #include <fstream>

    using namespace RooFit;

void Complexdata4(){

    //OPEN ROOT FILE

    gSystem->RedirectOutput("../logs/Complexdata4.log", "w");

    TFile *f = new TFile("../Root_files/ComplexData.root", "READ");
    TTree *t = (TTree*)f->Get("data");

    //CREATE HISTOGRAMS

    TH1F *h1 = new TH1F("h1", "Generated Distribution; Var; Counts", 700, -10.0, 10.0);
    t->Draw("var4 >>h1", "weight4", "goff");
    h1->SetDirectory(0);
    f->Close();
    delete f;

    // ROOFIT VARIABLE

    RooRealVar var("var4", "Variable", -10.0, 10.0);

    // IMPORT DATA

    RooDataHist data("data", "dataset", var, Import(*h1));

    //SIGNAL PARAMETERS

    RooRealVar mean1("mean1", "Mean1", 0.0, -2.0, 2.0);
    RooRealVar sigma1("sigma1", "Sigma1", 1.0, 0.1, 5.0);
    RooRealVar sigma2("sigma2", "Sigma2", 2.75, 0.1, 10.0);
    RooRealVar frac("frac", "Fraction of Signal 1", 0.5, 0.0, 1.0);
    sigma2.setVal(2.5);
    sigma2.setConstant(true);
    

    // BACKGROUND PARAMETERS

    RooRealVar tau("tau", "Tau", 1.0, 0.1, 10.0);
    RooRealVar meanbkg("meanbkg", "Meanbkg", 1.0, -2.0, 2.0);
    RooRealVar sigmabkg("sigmabkg", "SigmaBkg", 1.0, 0.1, 5.0);
    RooRealVar fracbkg("fracbkg", "Fraction of Background 1", 0.5, 0.0, 1.0);


    // SIGNAL AND BACKGROUND PDFS

    RooGaussian signal1("signal1", "Signal PDF", var, mean1, sigma1);
    RooGaussian signal2("signal2", "Signal PDF 2", var, mean1, sigma2);
    RooAddPdf signal("signal", "Combined Signal PDF", RooArgList(signal1, signal2),frac);
    RooGenericPdf background1("background1", "exp(-TMath::Abs(var4)/tau)", RooArgList(var, tau));
    RooGaussian background2("background2", "Gaussian Background", var, meanbkg, sigmabkg);
    RooAddPdf background("background", "Combined Background PDF", RooArgList(background1, background2), fracbkg);



    // YEILD VARIABLES

    RooRealVar nsig("nsig", "Signal Yield", 100000, 0, 500000);
    RooRealVar nbkg("nbkg", "Background Yield", 100000, 0, 500000);

    // COMBINED MODEL

    RooAddPdf model("model", "Signal + Background Model", RooArgList(signal, background), RooArgList(nsig, nbkg));

    //FIT MODEL TO DATA

    RooFitResult* fitres = model.fitTo(data,Save(),SumW2Error(kTRUE));

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
    model.plotOn(frame2, Components(signal), LineColor(kRed), Name("sig"));
    model.plotOn(frame2, Components(background1), LineColor(kBlack), LineStyle(kDashed), Name("bkg1"));
    model.plotOn(frame2, Components(background2), LineColor(kMagenta), LineStyle(kDotted), Name("bkg2"));

    double chi2ndf = frame2->chiSquare("fullfit", "data2", 10);
    frame2->SetTitle("Fitted Distribution");
    frame2->Draw();

    TLegend *leg2 = new TLegend(0.1, 0.65, 0.4, 0.9);
    leg2->AddEntry(frame2->findObject("data2"), "Data", "ep");
    leg2->AddEntry(frame2->findObject("fullfit"), "Full Fit","l");
    leg2->AddEntry(frame2->findObject("sig"), "Signal PDF", "l");
    leg2->AddEntry(frame2->findObject("bkg"), "Background PDF", "l");
    leg2->AddEntry(frame2->findObject("bkg1"), "Exp Background", "l");
    leg2->AddEntry(frame2->findObject("bkg2"), "Gaussian Background", "l");
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
    //signal.plotOn(frame3, Name("sig"), LineColor(kBlue), LineStyle(kDotted), Normalization(nsig.getVal(),RooAbsReal::NumEvent));
    frame3->SetTitle("Extracted Signal Component");
    frame3->Draw();

    TLegend *leg3 = new TLegend(0.1, 0.75, 0.4, 0.9);
    leg3->AddEntry(frame3->findObject("sigdata"), "Signal Data", "ep");
    leg3->AddEntry(frame3->findObject("sigpdf"), "Fitted Signal", "l");
    //leg3->AddEntry(frame3->findObject("sig"), "Signal", "l");
    leg3->Draw();

    cmain->Modified();
    cmain->Update();
    gPad->Update();

    cmain->SaveAs("../Png_files/ComplexData4.png");

    ofstream out("../Txt_files/ComplexData4_results.txt");

    out << "Fit Results:\n";

    out << "Status: " << fitres->status() << endl;
    out << "CovQual: " << fitres->covQual() << endl;
    out << "EDM: " << fitres->edm() << endl;
    out << "MinNLL: " << fitres->minNll() << endl;
    out << "Chi2/NDF: " << chi2ndf << endl;

    fitres->Print("v");

    out << "\nSignal Yield: "<< nsig.getVal()<< " ± "<< nsig.getError()<< endl;
    out << "Background Yield: "<< nbkg.getVal()<< " ± "<< nbkg.getError()<< endl;
    out << "\nSignal Mean: "<< mean1.getVal()<< " ± "<< mean1.getError()<< endl;
    out << "Sigma1: "<< sigma1.getVal()<< " ± "<< sigma1.getError()<< endl;
    out << "Sigma2: "<< sigma2.getVal()<< " ± "<< sigma2.getError()<< endl;
    out << "Signal Fraction: "<< frac.getVal()<< " ± "<< frac.getError()<< endl;

    out.close();

    gSystem->RedirectOutput(0);
}
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

void ComplexFitSamples1(){

    //OPEN ROOT FILE

    gSystem->RedirectOutput("../logs/ComplexFitSamples1.log", "w");

    TFile *f = new TFile("../Root_files/ComplexFitSamples.root", "READ");
    TTree *t = (TTree*)f->Get("data");

    //CREATE HISTOGRAMS

    TH1F *h1 = new TH1F("h1", "Generated Distribution; Var; Counts", 700, 110.0, 160.0);
    t->Draw("var1 >>h1", "", "goff");
    h1->SetDirectory(0);
    f->Close();
    delete f;

    // ROOFIT VARIABLE

    RooRealVar var("var1", "Variable", 110.0, 160.0);

    // IMPORT DATA

    RooDataHist data("data", "dataset", var, Import(*h1));

    //SIGNAL PARAMETERS

    RooRealVar mean1("mean1", "Mean of Signal 1", 125.0, 120.0, 130.0);
    RooRealVar sigma1("sigma1", "Sigma of Signal 1", 2.0, 0.1, 10.0);
    RooRealVar alpha1("alpha1", "Alpha of Crystal Ball 1", 1.5, 0.1, 5.0);
    RooRealVar n1("n1", "N of Crystal Ball 1", 2.0, 0.1, 10.0);

    RooRealVar mean2("mean2", "Mean of Signal 2", 140.0, 130.0, 150.0);
    RooRealVar width2("width2", "Width of Signal 2", 7.0, 0.1, 25.0);

    RooRealVar mean3("mean3", "Mean of Signal 3", 125.0, 120.0, 130.0);
    RooRealVar sigma3("sigma3", "Sigma of Signal 3", 5.0, 0.1, 10.0);

    RooRealVar frac1("frac1", "Fraction of Signal 1", 0.5, 0.0, 1.0);
    RooRealVar frac2("frac2", "Fraction of Signal 2", 0.3, 0.0, 1.0);


    // BACKGROUND PARAMETERS

    RooRealVar a1("a1", "a1", 1.2, 0.0, 100.0);
    RooRealVar a2("a2", "a2", 0.01, 0.0, 100.0);

    // SIGNAL AND BACKGROUND PDFS

    RooCBShape signal1("signal1", "Signal Component 1", var, mean1, sigma1, alpha1, n1);
    RooBreitWigner signal2("signal2", "Signal Component 2", var, mean2, width2);
    RooGaussian signal3("signal3", "Signal Component 3", var, mean3, sigma3);
    RooAddPdf signal("signal", "Total Signal", RooArgList(signal1, signal2, signal3), RooArgList(frac1, frac2));
    RooBernstein background("background", "Background Component", var, RooArgList(a1, a2));
   
    // YEILD VARIABLES

    RooRealVar nsig("nsig", "Signal Yield", 100000, 0, 500000);
    RooRealVar nbkg("nbkg", "Background Yield", 100000, 0, 500000);

    // COMBINED MODEL

    RooAddPdf model("model", "Signal + Background Model", RooArgList(signal, background), RooArgList(nsig, nbkg));

    //FIT MODEL TO DATA

    RooFitResult* fitres = model.fitTo(data, Save(), SumW2Error(kTRUE));

    //CANVAS CREATION

    TCanvas *cmain1 = new TCanvas("cmain1", "Fit Result", 1600, 900);

    cmain1->Divide(3,1);

    // RAW DISTRIBUTION

    cmain1->cd(1);
    RooPlot *frame1 = var.frame();
    data.plotOn(frame1, Name("data"),DataError(RooAbsData::SumW2));
    frame1->SetTitle("Raw Distribution");
    frame1->Draw();

    TLegend *leg1 = new TLegend(0.6, 0.85, 0.9, 0.9);
    leg1->AddEntry(frame1->findObject("data"), "Data", "ep");
    leg1->Draw();


    // FITTED DISTRIBUTION

    cmain1->cd(2);
    RooPlot *frame2 = var.frame();
    data.plotOn(frame2, Name("data2"), DataError(RooAbsData::SumW2));
    model.plotOn(frame2, Name("fullfit"));
    model.plotOn(frame2, Components(background), LineColor(kGreen), Name("bkg"));
    model.plotOn(frame2, Components(signal), LineColor(kRed), Name("sig"));
    double chi2ndf = frame2->chiSquare("fullfit", "data2", 10);
    frame2->SetTitle("Fitted Distribution");
    frame2->Draw();

    TLegend *leg2 = new TLegend(0.6, 0.65, 0.9, 0.9);
    leg2->AddEntry(frame2->findObject("data2"), "Data", "ep");
    leg2->AddEntry(frame2->findObject("fullfit"), "Full Fit","l");
    leg2->AddEntry(frame2->findObject("sig"), "Signal PDF", "l");
    leg2->AddEntry(frame2->findObject("bkg"), "Background PDF", "l");
    leg2->Draw();

    // EXTRACTED SIGNAL COMPONENT

    cmain1->cd(3);
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
    signal1.plotOn(frame3, Name("sig1"), LineColor(kBlue), LineStyle(kDotted), Normalization(hsig->Integral()*frac1.getVal(),RooAbsReal::NumEvent));
    signal2.plotOn(frame3, Name("sig2"), LineColor(kGreen), LineStyle(kDashed), Normalization(hsig->Integral()*frac2.getVal(),RooAbsReal::NumEvent));
    signal3.plotOn(frame3, Name("sig3"), LineColor(kMagenta), LineStyle(kDashed), Normalization(hsig->Integral()*(1 - frac1.getVal() - frac2.getVal()),RooAbsReal::NumEvent));
    frame3->SetTitle("Extracted Signal Component");
    frame3->Draw();

    TLegend *leg3 = new TLegend(0.6, 0.75, 0.9, 0.9);
    leg3->AddEntry(frame3->findObject("sigdata"), "Signal Data", "ep");
    leg3->AddEntry(frame3->findObject("sigpdf"), "Fitted Signal", "l");
    leg3->AddEntry(frame3->findObject("sig1"), "Signal Component 1", "l");
    leg3->AddEntry(frame3->findObject("sig2"), "Signal Component 2", "l");
    leg3->AddEntry(frame3->findObject("sig3"), "Signal Component 3", "l");
    leg3->Draw();


    TCanvas *cmain2 = new TCanvas("cmain2", "Pull Distribution and Pull Gauss", 1600, 900);
    cmain2->Divide(1,2);
    cmain2->cd(1);

    RooPlot *frame4 = var.frame();
    frame4->SetMinimum(-5);
    frame4->SetMaximum(5);
    RooHist *pullHist = frame2->pullHist("data2", "fullfit");
    for (int i=0; i<pullHist->GetN(); i++) 
    {
        pullHist->SetPointEYhigh(i,0);
        pullHist->SetPointEYlow(i,0);
    }
    frame4->addPlotable(pullHist, "P");
    frame4->SetTitle("Pull Distribution; Var; Pull");
    frame4->Draw();

    TLine *zero = new TLine(var.getMin(), 0, var.getMax(), 0);
    TLine *plus = new TLine(var.getMin(), 3, var.getMax(), 3);
    TLine *minus = new TLine(var.getMin(),-3, var.getMax(),-3);
    zero->SetLineColor(kRed);
    plus->SetLineColor(kBlue);
    minus->SetLineColor(kBlue);
    plus->SetLineStyle(2);
    minus->SetLineStyle(2);
    zero->Draw("same");
    plus->Draw("same");
    minus->Draw("same");

    
    cmain2->cd(2);
    TH1F *hpull = new TH1F("hpull", "Pull Values; Pull; Bins",20, -5, 5);
    for (int i = 0; i < pullHist->GetN(); i++) {
        double x, y;
        pullHist->GetPoint(i, x, y);
        hpull->Fill(y);
    }
    hpull->Fit("gaus");

    cmain1->Modified();
    cmain1->Update();
    cmain2->Modified();
    cmain2->Update();
    gPad->Update();

    cmain1->SaveAs("../Png_files/ComplexFitSamples1_results.png");
    cmain2->SaveAs("../Png_files/ComplexFitSamples1_pull.png");

    ofstream out("../Txt_files/ComplexFitSamples1_results.txt");

    out << "Fit Results:\n" << endl;

    out << "Status: " << fitres->status() << endl;
    out << "CovQual: " << fitres->covQual() << endl;
    out << "EDM: " << fitres->edm() << endl;
    out << "MinNLL: " << fitres->minNll() << endl;
    out << "Chi2/NDF: " << chi2ndf << endl;

    fitres->Print("v");
    fitres->correlationMatrix().Print();

    out << "\nSignal Yield: "<< nsig.getVal()<< " ± "<< nsig.getError()<< endl;
    out << "Mean 1: "<< mean1.getVal()<< " ± "<< mean1.getError()<< endl;
    out << "Sigma 1: "<< sigma1.getVal()<< " ± "<< sigma1.getError()<< endl;
    out << "Alpha 1: "<< alpha1.getVal()<< " ± "<< alpha1.getError()<< endl;
    out << "N 1: "<< n1.getVal()<< " ± "<< n1.getError()<< endl;
    out << "Mean 2: "<< mean2.getVal()<< " ± "<< mean2.getError()<< endl;
    out << "Width 2: "<< width2.getVal()<< " ± "<< width2.getError()<< endl;
    out << "Fraction of Signal 1: "<< frac1.getVal()<< " ± "<< frac1.getError()<< endl;
    out << "Fraction of Signal 2: "<< frac2.getVal()<< " ± "<< frac2.getError()<< endl;
    out << "Background Yield: " << nbkg.getVal() << " ± " << nbkg.getError() << endl;

    out.close();
    

    gSystem->RedirectOutput(0);
}
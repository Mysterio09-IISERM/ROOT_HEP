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

void Gendata3_2(){

    //OPEN ROOT FILE

    gSystem->RedirectOutput("../logs/Gendata3_2.log", "w");

    TFile *f = new TFile("../Root_files/GenData3.root", "READ");
    TTree *t = (TTree*)f->Get("data");

    //CREATE HISTOGRAMS

    TH1F *h1 = new TH1F("h1", "Generated Distribution; Var; Counts", 800, 20.0, 70.0);
    t->Draw("var2 >>h1", "", "goff");
    h1->SetDirectory(0);
    f->Close();
    delete f;

    // ROOFIT VARIABLE

    RooRealVar var("var2", "Variable", 20.0, 70.0);

    // IMPORT DATA

    RooDataHist data("data", "dataset", var, Import(*h1));

    //SIGNAL PARAMETERS

    RooRealVar mean1("mean1", "Mean1", 40.0, 35.0, 45.0);
    RooRealVar sigma1("sigma1", "Sigma1", 5.0, 0.1, 12.0);
    RooRealVar mean2("mean2", "Mean2", 52.0, 45.0, 60.0);
    RooRealVar sigma2("sigma2", "Sigma2", 7.0, 0.1, 15.0);
    RooRealVar frac("frac", "Fraction of Signal 1", 0.5, 0.0, 1.0);

    // BACKGROUND PARAMETERS

    RooRealVar tau("tau", "Tau", -0.1, -5.0, 0.0);

    // SIGNAL AND BACKGROUND PDFS

    RooGaussian signal1("signal1", "Signal 1 PDF", var, mean1, sigma1);
    RooGaussian signal2("signal2", "Signal 2 PDF", var, mean2, sigma2);
    RooAddPdf signal("signal", "Combined Signal PDF", RooArgList(signal1, signal2),frac);
    RooExponential background("background", "Background PDF", var, tau);

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
    signal1.plotOn(frame3, Name("sig1"), LineColor(kBlue), LineStyle(kDotted), Normalization(hsig->Integral()*frac.getVal(),RooAbsReal::NumEvent));
    signal2.plotOn(frame3, Name("sig2"), LineColor(kGreen), LineStyle(kDashDotted), Normalization(hsig->Integral()*(1-frac.getVal()),RooAbsReal::NumEvent));
    frame3->SetTitle("Extracted Signal Component");
    frame3->Draw();

    TLegend *leg3 = new TLegend(0.6, 0.75, 0.9, 0.9);
    leg3->AddEntry(frame3->findObject("sigdata"), "Signal Data", "ep");
    leg3->AddEntry(frame3->findObject("sigpdf"), "Fitted Signal", "l");
    leg3->AddEntry(frame3->findObject("sig1"), "Signal 1", "l");
    leg3->AddEntry(frame3->findObject("sig2"), "Signal 2", "l");
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

    cmain1->SaveAs("../Png_files/Gendata3_2_results.png");
    cmain2->SaveAs("../Png_files/Gendata3_2_pull.png");

    ofstream out("../Txt_files/Gendata3_2_results.txt");

    out << "Fit Results:\n" << endl;

    out << "Status: " << fitres->status() << endl;
    out << "CovQual: " << fitres->covQual() << endl;
    out << "EDM: " << fitres->edm() << endl;
    out << "MinNLL: " << fitres->minNll() << endl;
    out << "Chi2/NDF: " << chi2ndf << endl;

    fitres->Print("v");

    out << "\nSignal Yield: "<< nsig.getVal()<< " ± "<< nsig.getError()<< endl;
    out << "Mean1: "<< mean1.getVal()<< " ± "<< mean1.getError()<< endl;
    out << "Sigma1: "<< sigma1.getVal()<< " ± "<< sigma1.getError()<< endl;
    out << "Mean2: "<< mean2.getVal()<< " ± "<< mean2.getError()<< endl;
    out << "Sigma2: "<< sigma2.getVal()<< " ± "<< sigma2.getError()<< endl;
    out << "Background Yield: " << nbkg.getVal() << " ± " << nbkg.getError() << endl;
    out << "Fraction of Signal 1: " << frac.getVal() << " ± " << frac.getError() << endl;

    out.close();
    

    gSystem->RedirectOutput(0);
}
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

void sig_bkg_4(){

    //OPEN ROOT FILE

    gSystem->RedirectOutput("../logs/sig_bkg_4.log", "w");

    TFile *f = new TFile("../Root_files/sig_bkg_4.root", "READ");
    TTree *t = (TTree*)f->Get("tree");

    //CREATE HISTOGRAMS

    //TH1F *h1 = new TH1F("h1", "Generated Distribution; Var; Counts", 150, 0.0, 100.0);
    //t->Draw("x >>h1", "", "goff");
    //h1->SetDirectory(0);
    //f->Close();
    //delete f;

    // ROOFIT VARIABLE

    RooRealVar var("x", "Variable", 0.0, 100.0);

    // IMPORT DATA

    //RooDataHist data("data", "dataset", var, Import(*h1));
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

    
    // YEILD VARIABLES

    RooRealVar nsig("nsig", "Signal Yield", 1000, 0, 500000);
    RooRealVar nbkg("nbkg", "Background Yield", 1000, 0, 500000);
    

    // COMBINED MODEL

    RooAddPdf model("model", "Signal + Background Model", RooArgList(signal, background), RooArgList(nsig, nbkg));

    //FIT MODEL TO DATA

    RooFitResult* fitres = model.fitTo(data, Save(), Extended(kTRUE), PrintLevel(1));

    int nbins = (int)sqrt(data.numEntries());
    cout << "Using " << nbins << " bins for display" << endl;

    //CANVAS CREATION

    TCanvas *cmain1 = new TCanvas("cmain1", "Fit Result", 1600, 900);
    TPad *pad1 = new TPad("pad1","",0,0.25,1,1);
    pad1->SetBottomMargin(0.02);
    pad1->Draw();

    TPad *pad2 = new TPad("pad2","",0,0,1,0.25);
    pad2->SetTopMargin(0.05);
    pad2->SetBottomMargin(0.30);
    pad2->Draw();

    pad1->Divide(2,1);
    pad2->Divide(2,1);

    // RAW DISTRIBUTION

    pad1->cd(1);
    RooPlot *frame1 = var.frame(Bins(nbins));
    data.plotOn(frame1, Name("data"),Binning(nbins));
    frame1->SetTitle("Raw Distribution");
    frame1->Draw();

    TLegend *leg1 = new TLegend(0.1, 0.85, 0.4, 0.9);
    leg1->AddEntry(frame1->findObject("data"), "Data", "ep");
    leg1->Draw();


    // FITTED DISTRIBUTION

    pad1->cd(2);
    RooPlot *frame2 = var.frame(Bins(nbins));
    data.plotOn(frame2, Name("data2"), Binning(nbins));
    model.plotOn(frame2, Name("fullfit"));
    model.plotOn(frame2, Components(background), LineColor(kGreen), Name("bkg"));
    model.plotOn(frame2, Components(signal), LineColor(kRed), Name("sig"));
    double nparams = fitres->floatParsFinal().getSize();
    double chi2ndf = frame2->chiSquare("fullfit", "data2", nparams);
    frame2->SetTitle("Fitted Distribution");
    frame2->Draw();

    TLegend *leg2 = new TLegend(0.6, 0.65, 0.9, 0.9);
    leg2->AddEntry(frame2->findObject("data2"), "Data", "ep");
    leg2->AddEntry(frame2->findObject("fullfit"), "Full Fit","l");
    leg2->AddEntry(frame2->findObject("sig"), "Signal PDF", "l");
    leg2->AddEntry(frame2->findObject("bkg"), "Background PDF", "l");
    leg2->Draw();

    // EXTRACTED SIGNAL COMPONENT

    pad2->cd(1);
    RooPlot *frame3 = var.frame(Bins(nbins));
    frame3->SetMinimum(-5);
    frame3->SetMaximum(5);
    RooHist *pullHist = frame2->pullHist("data2", "fullfit");   
    for (int i=0; i<pullHist->GetN(); i++) 
    {
        pullHist->SetPointEYhigh(i,0);
        pullHist->SetPointEYlow(i,0);
    } 

    frame3->addPlotable(pullHist, "P");
    frame3->SetTitle("Pull Distribution; Var; Pull");
    frame3->Draw();

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

    pad2->cd(2);
    TH1F *hpull = new TH1F("hpull", "Pull Values; Pull; Bins",20, -5, 5);
    for (int i = 0; i < pullHist->GetN(); i++) {
        double x, y;
        pullHist->GetPoint(i, x, y);
        hpull->Fill(y);
    }
    hpull->Fit("gaus");

    cmain1->Modified();
    cmain1->Update();
    gPad->Update();

    cmain1->SaveAs("../Png_files/sig_bkg_4_results.png");

    ofstream out("../Txt_files/sig_bkg_4_results.txt");

    out << "Fit Results:\n" << endl;

    out << "Status: " << fitres->status() << endl;
    out << "CovQual: " << fitres->covQual() << endl;
    out << "EDM: " << fitres->edm() << endl;
    out << "MinNLL: " << fitres->minNll() << endl;
    out << "Chi2/NDF: " << chi2ndf << endl;

    fitres->Print("v");
    fitres->correlationMatrix().Print();

    out << "\nSignal Yield: "<< nsig.getVal()<< " ± "<< nsig.getError()<< endl;
    out << "Mean: " << mean.getVal() << " ± " << mean.getError() << endl;
    out << "Sigma: " << sigma.getVal() << " ± " << sigma.getError() << endl;
    out << "\nBackground Yield: "<< nbkg.getVal()<< " ± "<< nbkg.getError()<< endl;
    out << "Tau: " << tau.getVal() << " ± " << tau.getError() << endl;

    out.close();
    

    gSystem->RedirectOutput(0);
}
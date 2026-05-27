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

    using namespace RooFit;

void Mixedfitdata7(){

    //OPEN ROOT FILE

    TFile *f = new TFile("../Root_files/mixed_fit_data.root", "READ");
    TTree *t = (TTree*)f->Get("data");

    //CREATE HISTOGRAMS

    TH1F *h1 = new TH1F("h1", "Generated Distribution; Var; Counts", 100, 40.0, 80.0);
    t->Draw("tiny_signal >>h1", "", "goff");
    h1->SetDirectory(0);
    f->Close();
    delete f;

    // ROOFIT VARIABLE

    RooRealVar var("tiny_signal", "Variable", 40.0, 80.0);

    // IMPORT DATA

    RooDataHist data("data", "dataset", var, Import(*h1));

    //SIGNAL PARAMETERS

    RooRealVar mean("mean", "Mean", 60.0, 55.0, 65.0);
    RooRealVar sigma("sigma", "Sigma", 0.5, 0.1, 10.0);

    // BACKGROUND PARAMETERS

    RooRealVar tau("tau", "Exponential Tau", -0.05, -1.0, 0.0);

    // SIGNAL AND BACKGROUND PDFS

    RooGaussian signal("signal", "Signal PDF", var, mean, sigma);
    RooExponential background("background", "Background PDF", var, tau);

    // YEILD VARIABLES

    RooRealVar nsig("nsig", "Signal Yield", 10000, 0, 100000);
    RooRealVar nbkg("nbkg", "Background Yield", 40000, 0, 100000);

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

    TLegend *leg2 = new TLegend(0.1, 0.65, 0.4, 0.9);
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

    RooDataHist signalData("signalData", "Signal Data", var, Import(*hsig));
    RooPlot *frame3 = var.frame();
    signalData.plotOn(frame3, Name("sigdata"));
    signal.plotOn(frame3, Name("sigpdf"), LineColor(kRed), LineStyle(kDashed), Normalization(nsig.getVal(),RooAbsReal::NumEvent));
    frame3->SetTitle("Extracted Signal Component");
    frame3->Draw();

    /*hsig->SetLineColor(kBlue);
    hsig->SetMarkerStyle(20);
    hsig->Draw("E1");

    TF1 *fsig = new TF1("fsig", "gaus", 0, 200);
    double scale = nsig.getVal() * hsig->GetBinWidth(1)/(sigma.getVal() * sqrt(2*TMath::Pi()));
    fsig->SetParameters(scale, mean.getVal(), sigma.getVal());
    fsig->SetLineColor(kRed);
    fsig->SetLineWidth(2);
    gStyle->SetOptStat(0);
    fsig->Draw("same");*/

    TLegend *leg3 = new TLegend(0.1, 0.75, 0.4, 0.9);
    leg3->AddEntry(frame3->findObject("sigdata"), "Signal Data", "ep");
    leg3->AddEntry(frame3->findObject("sigpdf"), "Fitted Signal", "l");
    leg3->Draw();

    cmain->Modified();
    cmain->Update();
    gPad->Update();

    cmain->SaveAs("../Png_files/mixed_fit7.png");

    ofstream out("../Txt_files/mixed_fit7_results.txt");

    out << "Fit Results:" << endl;
    out << "Signal Yield: " << nsig.getVal() << " ± " << nsig.getError() << endl;
    out << "Mean: " << mean.getVal() << " ± " << mean.getError() << endl;
    out << "Sigma: " << sigma.getVal() << " ± " << sigma.getError() << endl;
    out << "Background Yield: " << nbkg.getVal() << " ± " << nbkg.getError() << endl;

    out.close();
}
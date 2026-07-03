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
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TLine.h"

#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooBreitWigner.h"
#include "RooExponential.h"
#include "RooAddPdf.h"
#include "RooPlot.h"
#include "RooArgSet.h"
#include "RooGlobalFunc.h"
#include "RooChebychev.h"
#include "RooCBShape.h"
#include "RooBernstein.h"
#include "RooFitResult.h"
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/LikelihoodIntervalPlot.h"
#include "RooStats/LikelihoodInterval.h"
#include <fstream>

using namespace RooFit;

void task3_4() {

    gSystem->RedirectOutput("../logs/task1_2.log", "w");

    TFile *f = TFile::Open("../Root_files/Bplus_to_JpsiKplus_signal.root");
    TTree *t = (TTree*)f->Get("Events");

    Double_t  mu1_px, mu1_py, mu1_pz, mu1_E;
    Double_t mu2_px, mu2_py, mu2_pz, mu2_E;
    Double_t K_px, K_py, K_pz, K_E;
    Double_t InvM;
    Double_t jpsi_InvM;

    t->SetBranchAddress("mu1_px", &mu1_px);
    t->SetBranchAddress("mu1_py", &mu1_py);
    t->SetBranchAddress("mu1_pz", &mu1_pz);
    t->SetBranchAddress("mu1_E",  &mu1_E);

    t->SetBranchAddress("mu2_px", &mu2_px);
    t->SetBranchAddress("mu2_py", &mu2_py);
    t->SetBranchAddress("mu2_pz", &mu2_pz);
    t->SetBranchAddress("mu2_E",  &mu2_E);

    t->SetBranchAddress("K_px", &K_px);
    t->SetBranchAddress("K_py", &K_py);
    t->SetBranchAddress("K_pz", &K_pz);
    t->SetBranchAddress("K_E",  &K_E);

    t->SetBranchAddress("jpsi_InvM", &jpsi_InvM);
    t->SetBranchAddress("InvM", &InvM);

    TH1F *h = new TH1F("h", "K+ mass;M_{K+} (GeV);Events", 300, 0.9, 1.5);

    Long64_t nentries = t->GetEntries();

    RooRealVar var("var", "Variable", 4.0, 6.5);
    RooArgSet vars(var);
    RooDataSet data("data", "dataset", vars);

    // -----------------------------------------------------------------
    // Output tree storing the reconstructed J/psi four-momentum for
    // every event that passes the m2 > 0 cut, so it can be reused
    // downstream (e.g. combining with the kaon to reconstruct B+)
    // without having to redo the mu1+mu2 sum every time.
    // -----------------------------------------------------------------
    TFile *fout = new TFile("../Root_files/K+_momenta.root", "RECREATE");
    TTree *ktree = new TTree("KTree", "Reconstructed K+ four-momentum");

    Double_t k_px, k_py, k_pz, k_E, k_mass, mass_diff;

    ktree->Branch("k_px",   &k_px,   "k_px/D");
    ktree->Branch("k_py",   &k_py,   "k_py/D");
    ktree->Branch("k_pz",   &k_pz,   "k_pz/D");
    ktree->Branch("k_E",    &k_E,    "k_E/D");
    ktree->Branch("k_mass", &k_mass, "k_mass/D");
    ktree->Branch("mass_diff", &mass_diff, "mass_diff/D");


    for (Long64_t i = 0; i < nentries; i++) {

        t->GetEntry(i);

        Double_t px = mu1_px + mu2_px + K_px;
        Double_t py = mu1_py + mu2_py + K_py;
        Double_t pz = mu1_pz + mu2_pz + K_pz;
        Double_t E  = mu1_E  + mu2_E + K_E;

        Double_t m2 = E*E - (px*px + py*py + pz*pz);

        if (m2 > 0) {
            Double_t m = sqrt(m2);

            h->Fill(m);
            var = m;
            data.add(vars);

            k_px   = px;
            k_py   = py;
            k_pz   = pz;
            k_E    = E;
            k_mass = m;
            mass_diff = m - InvM;
            ktree->Fill();
        }
    }

    ktree->Write();
    fout->Close();
    delete fout;

    TCanvas *c = new TCanvas();
    h->Draw();
    h->SetDirectory(0);
    f->Close();
    delete f;

    c->SaveAs("../Png_files/task3.png");

    //RooDataHist data("data", "dataset", var, Import(*h));

    RooRealVar mean("mean", "Mean", 3.1, 3.05, 3.15);
    RooRealVar sigma("sigma", "Sigma", 0.01, 0.001, 0.1);
    RooRealVar alpha("alpha", "Alpha", 1.0, 0.0, 10.0);
    RooRealVar nCB("nCB", "NCB", 1.0, 0.0, 10.0);
    RooRealVar c0("c0", "Background Coefficient", 1.0, 0.0, 100.0);
    RooRealVar c1("c1", "Background Coefficient", 0.0, -0.0, 100.0);

    RooCBShape signal("signal", "Signal PDF", var, mean, sigma, alpha, nCB);
    RooBernstein background("background", "Background PDF", var, RooArgList(c0, c1));

    RooRealVar nsig("nsig", "Number of Signal Events", 1000, 0, 10000);
    RooRealVar nbkg("nbkg", "Number of Background Events", 1000, 0, 10000);

    RooAddPdf model("model", "Signal + Background Model", RooArgList(signal, background), RooArgList(nsig, nbkg));

    RooFitResult* fitres = model.fitTo(data, Save(), Extended(kTRUE), PrintLevel(1));

    int nbins = 60;
    double xlo = var.getMin();
    double xhi = var.getMax();
    double binw = (xhi - xlo) / nbins;

    TCanvas *cmain = new TCanvas("cmain", "Fit Result", 1600, 600);

    TPad *pad1 = new TPad("pad1", "Fit Result", 0, 0.25, 1, 1);
    pad1->SetBottomMargin(0.02);
    pad1->Draw();

    TPad *pad2 = new TPad("pad2", "Pull Distribution", 0, 0, 1, 0.25);
    pad2->SetTopMargin(0.05);
    pad2->SetBottomMargin(0.30);
    pad2->Draw();

    pad1->Divide(2, 1);
    pad2->Divide(2, 1);

    pad1->cd(1);

    RooPlot *frame1 = var.frame(Bins(nbins));
    data.plotOn(frame1, Name("data"), Binning(nbins));
    frame1->SetTitle("Raw Distribution");
    frame1->Draw();

    TLegend *leg1 = new TLegend(0.6, 0.85, 0.9, 0.9);
    leg1->AddEntry(frame1->findObject("data"), "Data", "ep");
    leg1->Draw();

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
    leg2->AddEntry(frame2->findObject("fullfit"), "Full Fit", "l");
    leg2->AddEntry(frame2->findObject("sig"), "Signal PDF", "l");
    leg2->AddEntry(frame2->findObject("bkg"), "Background PDF", "l");
    leg2->Draw();

    // -----------------------------------------------------------------
    // Manual bin-integrated pull calculation.
    // IntegrateBins() is broken in this ROOT 6.40-rc1 build, so instead
    // of relying on plotOn()'s pointwise curve evaluation + pullHist(),
    // we directly integrate the fitted PDF over each bin's width and
    // compare that to the observed bin count. This avoids the
    // curvature-induced pull oscillation caused by comparing a
    // bin-integrated data count to a point-sampled model curve.
    // -----------------------------------------------------------------

    TH1F* hdata = (TH1F*) data.createHistogram("hdata", var, Binning(nbins, xlo, xhi));

    TGraphAsymmErrors* pullGraph = new TGraphAsymmErrors(nbins);
    TH1F *hpull = new TH1F("hpull", "Pull Values; Pull; Bins", 20, -5, 5);

    for (int i = 0; i < nbins; i++) {
        double lo = xlo + i * binw;
        double hi = lo + binw;

        var.setRange("binRange", lo, hi);

        RooAbsReal* sigInt = signal.createIntegral(var, NormSet(var), Range("binRange"));
        RooAbsReal* bkgInt = background.createIntegral(var, NormSet(var), Range("binRange"));

        double expected = nsig.getVal() * sigInt->getVal() + nbkg.getVal() * bkgInt->getVal();
        double observed = hdata->GetBinContent(i + 1);

        double err = (observed > 0) ? sqrt(observed) : sqrt(expected);
        double pull = (err > 0) ? (observed - expected) / err : 0.0;

        double center = lo + binw / 2.0;
        pullGraph->SetPoint(i, center, pull);
        pullGraph->SetPointError(i, binw / 2.0, binw / 2.0, 0, 0);

        hpull->Fill(pull);

        delete sigInt;
        delete bkgInt;
    }

    pad2->cd(1);
    RooPlot *frame3 = var.frame(Bins(nbins));
    frame3->SetMinimum(-5);
    frame3->SetMaximum(5);
    frame3->SetTitle("Pull Distribution; Var; Pull");
    frame3->Draw();
    pullGraph->SetMarkerStyle(20);
    pullGraph->SetMarkerSize(0.6);
    pullGraph->Draw("P SAME");

    TLine *zero = new TLine(var.getMin(), 0, var.getMax(), 0);
    TLine *plus = new TLine(var.getMin(), 3, var.getMax(), 3);
    TLine *minus = new TLine(var.getMin(), -3, var.getMax(), -3);
    zero->SetLineColor(kRed);
    plus->SetLineColor(kBlue);
    minus->SetLineColor(kBlue);
    plus->SetLineStyle(2);
    minus->SetLineStyle(2);
    zero->Draw("same");
    plus->Draw("same");
    minus->Draw("same");

    pad2->cd(2);
    hpull->Fit("gaus");

    cmain->Modified();
    cmain->Update();
    gPad->Update();

    cmain->SaveAs("../Png_files/task4_results.png");

    /*cout << "Before PLC: "<< nsig.getVal() << endl;

    float propernsig = nsig.getVal();
    float propernsigError = nsig.getError();

    RooArgSet poi(nsig);
    RooStats::ProfileLikelihoodCalculator plc(data, model, poi);
    plc.SetConfidenceLevel(0.90);
    RooStats::LikelihoodInterval* interval = plc.GetInterval();
    double upperLimit = interval->UpperLimit(nsig);
    double lowerLimit = interval->LowerLimit(nsig);

    cout << "Upper Limit of Signal Events at 90% CL: " << upperLimit << endl;
    cout << "Lower Limit of Signal Events at 90% CL: " << lowerLimit << endl;

    cout << "After PLC: " << nsig.getVal() << endl;

    delete interval;*/

    ofstream out("../Txt_files/task3_4_results.txt");

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
    out << "Alpha: " << alpha.getVal() << " ± " << alpha.getError() << endl;
    out << "N: " << nCB.getVal() << " ± " << nCB.getError() << endl;

    out.close();

    gSystem->RedirectOutput(0);
}
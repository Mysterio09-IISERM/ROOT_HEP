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


void NewGendataPart2(){

    //UNIQUE ID
    TString uid =TString::Format("%d",gRandom->Integer(10000000));

    //OPEN ROOT FILE

    TFile *f = new TFile("../Root_files/mixdata.root", "READ");
    TTree *t = (TTree*)f->Get("data");

    //CREATE HISTOGRAMS

    TH1F *hmass = new TH1F("hmass"+uid, "Mass Distribution; Mass [GeV]; Counts", 100, 60, 120);
    t->Draw(Form("mass >>%s", hmass->GetName()), "", "goff");
    TH1F *henergy = new TH1F("henergy"+uid, "Energy Distribution; Energy [GeV]; Counts", 100, 0, 100);
    t->Draw(Form("energy >>%s", henergy->GetName()), "", "goff");
    hmass->SetDirectory(0);
    henergy->SetDirectory(0);
    f->Close();
    delete f;

    //ROOFIT VARIABLE

    RooRealVar mass("mass", "Mass [GeV]", 60, 120);
    RooRealVar energy("energy", "Energy [GeV]", 0, 100);

    // IMPORT DATA

    RooDataHist data_mass("data_mass"+uid, "dataset", mass, Import(*hmass));
    RooDataHist data_energy("data_energy"+uid, "dataset", energy, Import(*henergy));

    // SIGNAL PARAMETERS

    //Mass
    RooRealVar mass_mean("mass_mean"+uid, "Z Mass", 90.0, 85.0, 95.0);
    RooRealVar mass_width("mass_width"+uid, "Width", 2.0, 0.1, 10.0);
    //Energy
    RooRealVar energy_mean("energy_mean"+uid, "Energy Mean", 50.0, 0.0, 100.0);
    RooRealVar energy_sigma("energy_sigma"+uid, "Energy Sigma", 10.0, 0.1, 20.0);

    // SIGNAL PDF

    RooBreitWigner signal_mass("signal_mass"+uid, "Signal Mass PDF", mass, mass_mean, mass_width);
    RooGaussian signal_energy("signal_energy"+uid, "Signal Energy PDF", energy, energy_mean, energy_sigma);

    // BACKGROUND PARAMETERS

    RooRealVar bkg_slope_mass("bkg_slope_mass"+uid, "Background Slope Mass", -0.05, -2.0, 0.0);
    RooRealVar c0(TString("c0")+uid, "Cheby c0", 0.0, -1.0, 1.0);
    RooRealVar c1(TString("c1")+uid, "Cheby c1", 0.0, -1.0, 1.0);

    // BACKGROUND PDF

    RooExponential bkg_mass("bkg_mass"+uid, "Background Mass PDF", mass, bkg_slope_mass);
    RooChebychev bkg_energy("bkg_energy"+uid, "Background Energy PDF", energy, RooArgList(c0, c1));

    // YEILD VARIABLES

    RooRealVar nsig_mass("nsig_mass"+uid, "Signal Yield", 10000, 0, 100000);
    RooRealVar nsig_energy("nsig_energy"+uid, "Signal Yield", 40000, 0, 100000);
    RooRealVar nbkg_mass("nbkg_mass"+uid, "Background Yield", 40000, 0, 100000);
    RooRealVar nbkg_energy("nbkg_energy"+uid, "Background Yield", 60000, 0, 100000);

    // COMBINED MODELS

    RooAddPdf model_mass("model_mass"+uid, "Signal + Background Mass Model", RooArgList(signal_mass, bkg_mass), RooArgList(nsig_mass, nbkg_mass));
    RooAddPdf model_energy("model_energy"+uid, "Signal + Background Energy Model", RooArgList(signal_energy, bkg_energy), RooArgList(nsig_energy, nbkg_energy));

    //FITTING

    model_mass.fitTo(data_mass);
    model_energy.fitTo(data_energy);

    //CANVAS CREATION

    TCanvas *C1 = new TCanvas( "c1"+uid, "Mass Fit", 1800, 600);
    TCanvas *C2 = new TCanvas( "c2"+uid, "Energy Fit", 1800, 600);

    C1->Divide(3,1);
    C2->Divide(3,1);

    // RAW MASS AND ENERGY DISTRIBUTION

    C1->cd(1);
    RooPlot *frame_mass1 = mass.frame();
    data_mass.plotOn(frame_mass1, Name("data_mass"));
    frame_mass1->SetTitle("Raw Mass Distribution");
    frame_mass1->Draw();

    TLegend *leg_mass1 = new TLegend(0.6, 0.85, 0.9, 0.9);
    leg_mass1->AddEntry(frame_mass1->findObject("data_mass"), "Data", "ep");
    leg_mass1->Draw();

    C2->cd(1);
    RooPlot *frame_energy1 = energy.frame();
    data_energy.plotOn(frame_energy1, Name("data_energy"));
    frame_energy1->SetTitle("Raw Energy Distribution");
    frame_energy1->Draw();

    TLegend *leg_energy1 = new TLegend(0.6, 0.85, 0.9, 0.9);
    leg_energy1->AddEntry(frame_energy1->findObject("data_energy"), "Data", "ep");
    leg_energy1->Draw();

    // FITTED MASS AND ENERGY DISTRIBUTION

    C1->cd(2);
    RooPlot *frame_mass2 = mass.frame();
    data_mass.plotOn(frame_mass2, Name("data_mass2"));
    model_mass.plotOn(frame_mass2, Name("fullfit_mass"));
    model_mass.plotOn(frame_mass2, Components(bkg_mass), LineColor(kGreen), Name("bkg_mass"));
    model_mass.plotOn(frame_mass2, Components(signal_mass), LineColor(kRed), Name("sig_mass"));
    frame_mass2->SetTitle("Mass Fit");
    frame_mass2->Draw();

    TLegend *leg_mass2 = new TLegend(0.6, 0.65, 0.9, 0.9);
    leg_mass2->AddEntry(frame_mass2->findObject("data_mass2"), "Data", "ep");
    leg_mass2->AddEntry(frame_mass2->findObject("fullfit_mass"), "Full Fit", "l");
    leg_mass2->AddEntry(frame_mass2->findObject("sig_mass"), "Signal, Mass PDF", "l");
    leg_mass2->AddEntry(frame_mass2->findObject("bkg_mass"), "Background, Mass PDF", "l");
    leg_mass2->Draw();

    C2->cd(2);
    RooPlot *frame_energy2 = energy.frame();
    data_energy.plotOn(frame_energy2, Name("data_energy2"));
    model_energy.plotOn(frame_energy2, Name("fullfit_energy"));
    model_energy.plotOn(frame_energy2, Components(bkg_energy), LineColor(kGreen), Name("bkg_energy"));
    model_energy.plotOn(frame_energy2, Components(signal_energy), LineColor(kRed), Name("sig_energy"));
    frame_energy2->SetTitle("Energy Fit");
    frame_energy2->Draw();

    TLegend *leg_energy2 = new TLegend(0.6, 0.65, 0.9, 0.9);
    leg_energy2->AddEntry(frame_energy2->findObject("data_energy2"), "Data", "ep");
    leg_energy2->AddEntry(frame_energy2->findObject("fullfit_energy"), "Full Fit", "l");
    leg_energy2->AddEntry(frame_energy2->findObject("sig_energy"), "Signal, Energy PDF", "l");
    leg_energy2->AddEntry(frame_energy2->findObject("bkg_energy"), "Background, Energy PDF", "l");
    leg_energy2->Draw();

    // EXTRACTED SIGNAL COMPONENT FOR MASS AND ENERGY

    C1->cd(3);
    TH1F *hsig_mass = (TH1F*)hmass->Clone("hsig_mass"+uid);
    hsig_mass->SetTitle("Extracted Signal Data; Mass [GeV]; Counts");

    RooArgSet normSet_mass(mass);

    for(int i=1; i<=hsig_mass->GetNbinsX(); i++){
        double x = hsig_mass->GetBinCenter(i);
        double bw = hsig_mass->GetBinWidth(i);
        mass.setVal(x);
        double bkg_mass_density = bkg_mass.getVal(normSet_mass);
        double bkg_mass_counts = bkg_mass_density * nbkg_mass.getVal() * bw;
        double y = hsig_mass->GetBinContent(i);
        hsig_mass->SetBinContent(i, y - bkg_mass_counts);
    }

    hsig_mass->SetLineColor(kBlue);
    hsig_mass->SetMarkerStyle(20);

    gPad->Update();
    TPaveStats *st_mass = (TPaveStats*)hsig_mass->FindObject("stats");
    if (st_mass) {
        st_mass->SetName("stats");
        st_mass->SetLabel("Signal Mass Data");

    }

    TF1 *fsig_mass = new TF1("fsig_mass"+uid, "breitwigner", 60, 120);
    double scale_mass = nsig_mass.getVal() * hsig_mass->GetBinWidth(1);
    fsig_mass->SetParameters(scale_mass, mass_mean.getVal(), mass_width.getVal());
    fsig_mass->SetLineColor(kRed);
    fsig_mass->SetLineWidth(2);
    gStyle->SetOptStat(0);
    hsig_mass->Draw("E1");
    fsig_mass->Draw("same");

    TLegend *leg_mass3 = new TLegend(0.6, 0.75, 0.9, 0.9);
    leg_mass3->AddEntry(hsig_mass, "Mass Signal Data", "ep");
    leg_mass3->AddEntry(fsig_mass, "Mass Signal Fit", "l");
    leg_mass3->Draw();

    C2->cd(3);
    TH1F *hsig_energy = (TH1F*)henergy->Clone("hsig_energy"+uid);
    hsig_energy->SetTitle("Extracted Signal Data; Energy [GeV]; Counts");

    RooArgSet normSet_energy(energy);

    for(int i=1; i<=hsig_energy->GetNbinsX(); i++){
        double x = hsig_energy->GetBinCenter(i);
        double bw = hsig_energy->GetBinWidth(i);
        energy.setVal(x);
        double bkg_energy_density = bkg_energy.getVal(normSet_energy);
        double bkg_energy_counts = bkg_energy_density * nbkg_energy.getVal() * bw;
        double y = hsig_energy->GetBinContent(i);
        hsig_energy->SetBinContent(i, y - bkg_energy_counts);
    }

    hsig_energy->SetLineColor(kBlue);
    hsig_energy->SetMarkerStyle(20);

    gPad->Update();
    TPaveStats *st_energy = (TPaveStats*)hsig_energy->FindObject("stats");
    if (st_energy) {
        st_energy->SetName("stats");
        st_energy->SetLabel("Signal Energy Data");

    }

    TF1 *fsig_energy = new TF1("fsig_energy"+uid, "gaus", 0, 100);
    double sigma = energy_sigma.getVal();
    double binwidth = hsig_energy->GetBinWidth(1);
    double yield = nsig_energy.getVal();

    // Convert total yield -> Gaussian peak height
    double amplitude = yield * binwidth / (sigma * sqrt(2*TMath::Pi()));
    fsig_energy->SetParameters(amplitude,energy_mean.getVal(),sigma);
    fsig_energy->SetLineColor(kRed);
    fsig_energy->SetLineWidth(2);
    gStyle->SetOptStat(0);
    hsig_energy->Draw("E1");
    fsig_energy->Draw("same");

    TLegend *leg_energy3 = new TLegend(0.6, 0.75, 0.9, 0.9);
    leg_energy3->AddEntry(hsig_energy, "Energy Signal Data", "ep");
    leg_energy3->AddEntry(fsig_energy, "Energy Signal Fit", "l");
    leg_energy3->Draw();


    // UPDATE CANVASES

    C1->Modified();
    C1->Update();
    C2->Modified();
    C2->Update();
    gPad->Update();

    ofstream out("../Txt_files/fit_results.txt");

    out << "Mass Fit Results:" << endl;
    out << "Signal Yield: " << nsig_mass.getVal() << " ± " << nsig_mass.getError() << endl;
    out << "Mass Mean: " << mass_mean.getVal() << " ± " << mass_mean.getError() << " GeV" << endl;
    out << "Mass Width: " << mass_width.getVal() << " ± " << mass_width.getError() << " GeV" << endl;
    out << endl;
    out << endl;
    out << "Energy Fit Results:" << endl;
    out << "Signal Yield: " << nsig_energy.getVal() << " ± " << nsig_energy.getError() << endl;
    out << "Energy Mean: " << energy_mean.getVal() << " ± " << energy_mean.getError() << " GeV" << endl;
    out << "Energy Sigma: " << energy_sigma.getVal() << " ± " << energy_sigma.getError() << " GeV" << endl;

    out.close();

    C1->SaveAs("../Png_files/mass_fit.png");
    C2->SaveAs("../Png_files/energy_fit.png");


}
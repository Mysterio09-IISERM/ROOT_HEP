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

using namespace RooFit;

void Newtest(){


    // =====================================================
    // UNIQUE ID
    // =====================================================

    TString uid =
        TString::Format(
            "%d",
            gRandom->Integer(1000000)
        );

    // =====================================================
    // OPEN ROOT FILE
    // =====================================================

    TFile *file =
        new TFile(
            "../Root_files/generatedata.root",
            "READ"
        );

    TTree *t =
        (TTree*)file->Get("data");

    // =====================================================
    // CREATE HISTOGRAM
    // =====================================================

    TH1F *h =
        new TH1F(
            "h"+uid,
            "Mass Distribution",
            100,
            60,
            120
        );

    t->Draw(
        Form(
            "mass>>%s",
            h->GetName()
        ),
        "",
        "goff"
    );

    h->SetDirectory(0);

    file->Close();

    // =====================================================
    // DEFINE ROOFIT VARIABLE
    // =====================================================

    RooRealVar mass(
        "mass",
        "Mass [GeV]",
        60,
        120
    );

    // =====================================================
    // IMPORT DATA
    // =====================================================

    RooDataHist data(
        "data"+uid,
        "dataset",
        mass,
        Import(*h)
    );

    // =====================================================
    // SIGNAL PARAMETERS
    // =====================================================

    RooRealVar mean(
        "mean"+uid,
        "Z Mass",
        91.2,
        85.0,
        95.0
    );

    RooRealVar width(
        "width"+uid,
        "Width",
        2.5,
        0.1,
        10.0
    );

    // =====================================================
    // SIGNAL PDF
    // =====================================================

    RooBreitWigner signal(
        "signal"+uid,
        "Signal PDF",
        mass,
        mean,
        width
    );

    // =====================================================
    // BACKGROUND PARAMETER
    // =====================================================

    RooRealVar tau(
        "tau"+uid,
        "Background Slope",
        -0.05,
        -1.0,
        0.0
    );

    // =====================================================
    // BACKGROUND PDF
    // =====================================================

    RooExponential background(
        "background"+uid,
        "Background PDF",
        mass,
        tau
    );

    // =====================================================
    // EVENT YIELDS
    // =====================================================

    RooRealVar nsig(
        "nsig"+uid,
        "Signal Events",
        10000,
        0,
        100000
    );

    RooRealVar nbkg(
        "nbkg"+uid,
        "Background Events",
        40000,
        0,
        100000
    );

    // =====================================================
    // COMBINED MODEL
    // =====================================================

    RooAddPdf model(
        "model"+uid,
        "Signal + Background",
        RooArgList(signal, background),
        RooArgList(nsig, nbkg)
    );

    // =====================================================
    // FIT
    // =====================================================

    model.fitTo(data);

    // =====================================================
    // CREATE CANVAS
    // =====================================================

    TCanvas *c =
        new TCanvas(
            TString("c")+uid,
            "RooFit Mass Analysis",
            1800,
            700
        );

    c->Divide(3,1);

    // =====================================================
    // PANEL 1 : RAW DATA
    // =====================================================

    c->cd(1);
    RooPlot *frame1 = mass.frame();
    data.plotOn(frame1, Name("data1"));  // give it a name
    frame1->SetTitle("Raw Mass Distribution");
    frame1->Draw();

    // add manual TLegend
    TLegend *leg1 = new TLegend(0.6, 0.85, 0.9, 0.9);
    leg1->AddEntry(frame1->findObject("data1"), "Data", "ep");
    leg1->Draw();

    // =====================================================
    // PANEL 2 : FULL FIT
    // =====================================================

    c->cd(2);
    RooPlot *frame2 = mass.frame();
    data.plotOn(frame2, Name("data2"));
    model.plotOn(frame2, Name("fullfit"));
    model.plotOn(frame2, Components(background), 
                LineColor(kGreen), Name("bkg"));
    model.plotOn(frame2, Components(signal), 
                LineColor(kRed), Name("sig"));
    frame2->SetTitle("Signal + Background Fit");
    frame2->Draw();

    TLegend *leg2 = new TLegend(0.6, 0.65, 0.9, 0.9);
    leg2->AddEntry(frame2->findObject("data2"),  "Data",       "ep");
    leg2->AddEntry(frame2->findObject("fullfit"),"Full Fit",   "l");
    leg2->AddEntry(frame2->findObject("bkg"),    "Background", "l");
    leg2->AddEntry(frame2->findObject("sig"),    "Signal",     "l");
    leg2->Draw();

    // =====================================================
    // PANEL 3 : BACKGROUND SUBTRACTED SIGNAL
    // =====================================================

    c->cd(3);

    TH1F *hsig =
        (TH1F*)h->Clone(
            "hsig"+uid
        );

    hsig->SetTitle(
        "Background Subtracted Signal;Mass [GeV];Counts"
    );

    // IMPORTANT:
    // create normalization set ONCE

    RooArgSet normSet(mass);

    // subtract fitted background

    for(int i=1; i<=hsig->GetNbinsX(); i++){

        double x =
            hsig->GetBinCenter(i);

        double bw =
            hsig->GetBinWidth(i);

        // set RooFit variable
        mass.setVal(x);

        // normalized background density
        double bkg_density =
            background.getVal(normSet);

        // expected background counts
        double bkg =
            nbkg.getVal()
            *
            bkg_density
            *
            bw;

        // original count
        double y =
            hsig->GetBinContent(i);

        // subtract
        hsig->SetBinContent(
            i,
            y - bkg
        );
    }

    hsig->SetLineColor(kBlue);

    hsig->SetMarkerStyle(20);

    //hsig->Draw("E1");

    // after hsig->Draw("E1")
    gPad->Update();
    TPaveStats *st = (TPaveStats*)hsig->FindObject("stats");
    if (st) {
        st->SetName("stats");
        // either remove it entirely
        st->SetOptStat(0);  // hide stat box
        // or rename it
        st->SetLabel("Background Subtracted");
}

    // =====================================================
    // FITTED SIGNAL SHAPE
    // =====================================================

    TF1 *fsig =
        new TF1(
            "fsig"+uid,
            "[0]*TMath::BreitWigner(x,[1],[2])",
            60,
            120
        );

    double scale =
        nsig.getVal()
        *
        hsig->GetBinWidth(1);

    fsig->SetParameters(
        scale,
        mean.getVal(),
        width.getVal()
    );

    fsig->SetLineColor(kRed);

    fsig->SetLineWidth(2);

    gStyle->SetOptStat(0);  // suppress weird stat box
    hsig->Draw("E1");
    fsig->Draw("same");

    TLegend *leg3 = new TLegend(0.6, 0.75, 0.9, 0.9);
    leg3->AddEntry(hsig, "Bkg Subtracted", "ep");
    leg3->AddEntry(fsig, "Signal Fit",     "l");
    leg3->Draw();



    // =====================================================
    // UPDATE
    // =====================================================

    c->Modified();

    c->Update();

    gPad->Update();

    ofstream out("../Txt_files/Newtestresult.txt");

    out << "Fit Results:" << endl;
    out << "Mass: " << mean.getVal() << "+/-" << mean.getError() << " GeV" << endl;
    out << "Width: " << width.getVal() << "+/-" << width.getError() << " GeV" << endl;
    out << "Signal Events: " << nsig.getVal() << "+/-" << nsig.getError() << endl;
    out << "Background Events: " << nbkg.getVal() << "+/-" << nbkg.getError() << endl;

    out.close();

    // =====================================================
    // CLEANUP
    // =====================================================

}
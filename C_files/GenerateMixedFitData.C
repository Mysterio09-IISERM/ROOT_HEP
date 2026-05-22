#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TMath.h"

void GenerateMixedFitData() {

    TFile *f = new TFile("../Root_files/mixed_fit_data.root", "RECREATE");

    TRandom3 rng(0);

    double gauss_exp;
    double bw_exp;
    double doublegauss_poly;
    double asym_flat;
    double twopeak_bkg;
    double smeared_bw;
    double tiny_signal;
    double ugly_bkg;

    TTree *t = new TTree("data", "Separate fitting datasets");

    t->Branch("gauss_exp", &gauss_exp, "gauss_exp/D");
    t->Branch("bw_exp", &bw_exp, "bw_exp/D");
    t->Branch("doublegauss_poly", &doublegauss_poly, "doublegauss_poly/D");
    t->Branch("asym_flat", &asym_flat, "asym_flat/D");
    t->Branch("twopeak_bkg", &twopeak_bkg, "twopeak_bkg/D");
    t->Branch("smeared_bw", &smeared_bw, "smeared_bw/D");
    t->Branch("tiny_signal", &tiny_signal, "tiny_signal/D");
    t->Branch("ugly_bkg", &ugly_bkg, "ugly_bkg/D");

    for(int i=0; i<150000; i++) {

        // =====================================================
        // 1. Gaussian + Exponential
        // =====================================================

        if(rng.Rndm() < 0.35)
            gauss_exp = rng.Gaus(50,4);
        else
            gauss_exp = rng.Exp(35);

        // =====================================================
        // 2. Breit-Wigner + Exponential
        // =====================================================

        if(rng.Rndm() < 0.30)
            bw_exp = rng.BreitWigner(91,2.5);
        else
            bw_exp = rng.Exp(40);

        // =====================================================
        // 3. Double Gaussian + Polynomial Background
        // =====================================================

        if(rng.Rndm() < 0.40) {

            if(rng.Rndm() < 0.6)
                doublegauss_poly = rng.Gaus(40,3);
            else
                doublegauss_poly = rng.Gaus(48,7);

        } else {

            double r = rng.Rndm();

            doublegauss_poly = 100 * sqrt(r);
        }

        // =====================================================
        // 4. Asymmetric Signal + Flat Background
        // =====================================================

        if(rng.Rndm() < 0.30) {

            double g = rng.Gaus(60,4);

            if(rng.Rndm() < 0.15)
                g -= rng.Exp(6);

            asym_flat = g;

        } else {

            asym_flat = rng.Uniform(0,100);
        }

        // =====================================================
        // 5. Two Resonances + Background
        // =====================================================

        double r5 = rng.Rndm();

        if(r5 < 0.25)
            twopeak_bkg = rng.Gaus(35,2);

        else if(r5 < 0.45)
            twopeak_bkg = rng.BreitWigner(75,4);

        else
            twopeak_bkg = rng.Exp(50);

        // =====================================================
        // 6. Smeared Breit-Wigner
        // =====================================================

        if(rng.Rndm() < 0.35) {

            double true_mass = rng.BreitWigner(91,2);

            smeared_bw = rng.Gaus(true_mass,2);

        } else {

            double r = rng.Rndm();

            smeared_bw = 120 * pow(r,0.7);
        }

        // =====================================================
        // 7. Tiny Signal in Huge Background
        // =====================================================

        if(rng.Rndm() < 0.05)
            tiny_signal = rng.Gaus(60,2);
        else
            tiny_signal = rng.Exp(45);

        // =====================================================
        // 8. Broad Signal + Ugly Background
        // =====================================================

        if(rng.Rndm() < 0.30) {

            ugly_bkg = rng.Gaus(55,10);

        } else {

            if(rng.Rndm() < 0.5)
                ugly_bkg = rng.Exp(30);
            else
                ugly_bkg = rng.Uniform(20,120);
        }

        t->Fill();
    }

    t->Write();

    f->Close();

    cout << "separate_datasets.root created!" << endl;
}
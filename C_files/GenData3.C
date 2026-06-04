#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TMath.h"

void GenData3() {

    gSystem->Unlink("../Root_files/GenData3.root");

    TFile *f = new TFile("../Root_files/GenData3.root", "NEW");
    TTree *t = new TTree("data", "Intermediate complexity fitting datasets");

    TRandom3 rng(2025);

    double var1, var2, var3, var4, var5;

    t->Branch("var1", &var1, "var1/D");
    t->Branch("var2", &var2, "var2/D");
    t->Branch("var3", &var3, "var3/D");
    t->Branch("var4", &var4, "var4/D");
    t->Branch("var5", &var5, "var5/D");

    const int N = 500000;

    for (int i = 0; i < N; i++) {

        // Uncomment for debugging progress
        // if(i % 10000 == 0) printf("%d\n", i);

        // =====================================================
        // VAR1
        // =====================================================

        double r1 = rng.Rndm();

        if (r1 < 0.25) {
            do { var1 = rng.Gaus(50.0, 2.0); }
            while (var1 < 0 || var1 > 100);
        }
        else if (r1 < 0.70) {
            do { var1 = rng.Exp(30.0); }
            while (var1 > 100);
        }
        else {
            var1 = rng.Uniform(0, 100);
        }

        // =====================================================
        // VAR2
        // =====================================================

        double r2 = rng.Rndm();

        if (r2 < 0.20) {
            do { var2 = rng.Gaus(40.0, 3.0); }
            while (var2 < 0 || var2 > 100);
        }
        else if (r2 < 0.35) {
            do { var2 = rng.Gaus(52.0, 4.0); }
            while (var2 < 0 || var2 > 100);
        }
        else {
            do { var2 = rng.Exp(50.0); }
            while (var2 > 100);
        }

        // =====================================================
        // VAR3
        // =====================================================

        double r3 = rng.Rndm();

        if (r3 < 0.30) {

            do {
                double core = rng.Gaus(80.0, 3.0);

                if (rng.Rndm() < 0.18)
                    core -= rng.Exp(5.0);

                var3 = core;

            } while (var3 < 20 || var3 > 140);
        }
        else {

            do {
                double x = rng.Uniform(20, 140);

                double fx = (x - 20) * exp(-(x - 20) / 40.0);

                if (rng.Rndm() < fx / 14.7) {
                    var3 = x;
                    break;
                }

            } while (true);
        }

        // =====================================================
        // VAR4
        // =====================================================

        double r4 = rng.Rndm();

        if (r4 < 0.35) {

            do {
                double true_m = rng.BreitWigner(91.0, 2.5);
                var4 = rng.Gaus(true_m, 1.5);

            } while (var4 < 75 || var4 > 110);
        }
        else {

            do {
                double x = rng.Uniform(75, 110);

                double fx = pow(75.0 / x, 2.5);

                if (rng.Rndm() < fx) {
                    var4 = x;
                    break;
                }

            } while (true);
        }

        // =====================================================
        // VAR5
        // =====================================================

        double r5 = rng.Rndm();

        if (r5 < 0.20) {
            do { var5 = rng.Gaus(200.0, 5.0); }
            while (var5 < 100 || var5 > 350);
        }
        else if (r5 < 0.35) {
            do { var5 = rng.Gaus(240.0, 6.0); }
            while (var5 < 100 || var5 > 350);
        }
        else if (r5 < 0.45) {
            do { var5 = rng.Gaus(220.0, 25.0); }
            while (var5 < 100 || var5 > 350);
        }
        else {
            do { var5 = rng.Exp(80.0) + 100.0; }
            while (var5 > 350);
        }

        t->Fill();
    }

    f->cd();
    t->Write("", TObject::kOverwrite);
    f->Close();
    delete f;

    printf("\nGenData3.root created!\n");
    printf("Events: %d\n\n", N);
    printf("var1: Gaus + expo + flat bkg          range: 0-100\n");
    printf("var2: Two Gaus + expo bkg              range: 0-100\n");
    printf("var3: Crystal Ball + curved bkg        range: 20-140\n");
    printf("var4: Voigtian + power law bkg         range: 75-110 GeV\n");
    printf("var5: Two signals + peaking bkg + expo range: 100-350\n");
}
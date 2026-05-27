#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"

void GenNewdata(){

        // Remove old in-memory tree if it exists
    TObject *oldtree = gROOT->FindObject("data");
    if(oldtree) delete oldtree;

    // =====================================================
    // OUTPUT FILE
    // =====================================================

    TFile *f = new TFile(
        "../Root_files/Newdata.root",
        "RECREATE"
    );
    

    TTree *t = new TTree("data","Data");
    t->SetDirectory(f);

    // =====================================================
    // RANDOM GENERATOR
    // =====================================================

    // Fixed seed for reproducibility

    TRandom3 rng(15875);

    // =====================================================
    // TREE BRANCH VARIABLES
    // =====================================================

    double narrow_signal;
    double overlap_signal;
    double smeared_bw;
    double crystal_tail;
    double tiny_signal;
    double ugly_bkg;
    double reflection_peak;
    double threshold_bkg;
    double double_peak;
    double weighted_data;
    double weight;

    // =====================================================
    // CREATE BRANCHES
    // =====================================================

    t->Branch("narrow_signal",
              &narrow_signal,
              "narrow_signal/D");

    t->Branch("overlap_signal",
              &overlap_signal,
              "overlap_signal/D");

    t->Branch("smeared_bw",
              &smeared_bw,
              "smeared_bw/D");

    t->Branch("crystal_tail",
              &crystal_tail,
              "crystal_tail/D");

    t->Branch("tiny_signal",
              &tiny_signal,
              "tiny_signal/D");

    t->Branch("ugly_bkg",
              &ugly_bkg,
              "ugly_bkg/D");

    t->Branch("reflection_peak",
              &reflection_peak,
              "reflection_peak/D");

    t->Branch("threshold_bkg",
              &threshold_bkg,
              "threshold_bkg/D");

    t->Branch("double_peak",
              &double_peak,
              "double_peak/D");

    t->Branch("weighted_data",
              &weighted_data,
              "weighted_data/D");

    t->Branch("weight",
              &weight,
              "weight/D");

    // =====================================================
    // NUMBER OF EVENTS
    // =====================================================

    const int N = 500000;

    // =====================================================
    // EVENT LOOP
    // =====================================================

    for(int i=0; i<N; i++){

        // =================================================
        // 1. Narrow Signal + Curved Background
        // =================================================

        if(rng.Rndm() < 0.08){

            narrow_signal = rng.Gaus(91,1.5);

        }else{

            double r = rng.Rndm();

            narrow_signal = 120 - 100*sqrt(r);
        }

        // =================================================
        // 2. Two Overlapping Signals
        // =================================================

        double r2 = rng.Rndm();

        if(r2 < 0.20){

            overlap_signal = rng.Gaus(50,3);

        }else if(r2 < 0.35){

            overlap_signal = rng.Gaus(58,4);

        }else{

            overlap_signal = rng.Exp(40);
        }

        // =================================================
        // 3. Smeared Breit-Wigner
        // =================================================

        if(rng.Rndm() < 0.25){

            double true_mass =
                rng.BreitWigner(125,2);

            smeared_bw =
                rng.Gaus(true_mass,1.5);

        }else{

            smeared_bw =
                150*pow(rng.Rndm(),0.8);
        }

        // =================================================
        // 4. Crystal Ball Style Tail
        // =================================================

        if(rng.Rndm() < 0.30){

            double core =
                rng.Gaus(80,2);

            if(rng.Rndm() < 0.15)
                crystal_tail =
                    core - rng.Exp(4);
            else
                crystal_tail = core;

        }else{

            crystal_tail =
                rng.Uniform(40,120);
        }

        // =================================================
        // 5. Tiny Signal
        // =================================================

        if(rng.Rndm() < 0.015){

            tiny_signal =
                rng.Gaus(60,1.2);

        }else{

            tiny_signal =
                rng.Exp(35);
        }

        // =================================================
        // 6. Broad Signal + Ugly Background
        // =================================================

        double r6 = rng.Rndm();

        if(r6 < 0.30){

            ugly_bkg =
                rng.Gaus(55,10);

        }else if(r6 < 0.65){

            ugly_bkg =
                rng.Exp(30);

        }else{

            ugly_bkg =
                rng.Uniform(20,120);
        }

        // =================================================
        // 7. Reflection Peak
        // =================================================

        double r7 = rng.Rndm();

        if(r7 < 0.20){

            reflection_peak =
                rng.Gaus(90,2);

        }else if(r7 < 0.28){

            reflection_peak =
                rng.Gaus(78,5);

        }else{

            reflection_peak =
                rng.Exp(45);
        }

        // =================================================
        // 8. Threshold Background
        // =================================================

        if(rng.Rndm() < 0.15){

            threshold_bkg =
                rng.Gaus(70,2);

        }else{

            threshold_bkg =
                30 + 90*pow(rng.Rndm(),2.5);
        }

        // =================================================
        // 9. Double Peak
        // =================================================

        double r9 = rng.Rndm();

        if(r9 < 0.18){

            double_peak =
                rng.Gaus(65,2);

        }else if(r9 < 0.32){

            double_peak =
                rng.Gaus(70,2.5);

        }else{

            double_peak =
                rng.Exp(50);
        }

        // =================================================
        // 10. Weighted Dataset
        // =================================================

        if(rng.Rndm() < 0.20){

            weighted_data =
                rng.Gaus(50,3);

            weight = 2.0;

        }else{

            weighted_data =
                rng.Exp(40);

            weight = 0.5;
        }

        // =================================================
        // FILL TREE
        // =================================================

        t->Fill();
    }

    // =====================================================
    // WRITE FILE
    // =====================================================

    f->cd();

    // overwrite safely
    t->Write("", TObject::kOverwrite);

    delete t;
    f->Close();
    delete f;

    printf("\nROOT file created successfully.\n");
    printf("Seed used: 15875\n");
    printf("Events generated: %d\n\n", N);
}
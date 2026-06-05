#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TSystem.h"

void ComplexFitSamples()
{
    gSystem->Unlink("../Root_files/ComplexFitSamples.root");

    TFile *f = new TFile("../Root_files/ComplexFitSamples.root","RECREATE");
    TTree *t = new TTree("data","Advanced RooFit training datasets");

    double var1,var2,var3,var4,var5;

    t->Branch("var1",&var1,"var1/D");
    t->Branch("var2",&var2,"var2/D");
    t->Branch("var3",&var3,"var3/D");
    t->Branch("var4",&var4,"var4/D");
    t->Branch("var5",&var5,"var5/D");

    TRandom3 rng(2026);

    const int N = 500000;

    for(int i=0;i<N;i++)
    {
        //----------------------------------------------------
        // VAR1
        // CrystalBall-like + Gaussian + BreitWigner bkg
        //----------------------------------------------------

        double r1 = rng.Rndm();

        if(r1 < 0.25)
        {
            do{
                var1 = rng.Gaus(125.,2.);
                if(rng.Rndm()<0.15)
                    var1 -= rng.Exp(2.5);
            } while(var1<100 || var1>180);
        }
        else if(r1 < 0.35)
        {
            do{
                var1 = rng.Gaus(125.,4.);
            } while(var1<100 || var1>180);
        }
        else if(r1 < 0.50)
        {
            do{
                var1 = rng.BreitWigner(140.,12.);
            } while(var1<100 || var1>180);
        }
        else
        {
            while(true)
            {
                double x = rng.Uniform(100.,180.);
                double fx =
                    pow(x-100.,1.5) *
                    exp(-(x-100.)/20.);

                if(rng.Rndm() < fx/98.)
                {
                    var1 = x;
                    break;
                }
            }
        }

        //----------------------------------------------------
        // VAR2
        // Three overlapping peaks
        //----------------------------------------------------

        double r2 = rng.Rndm();

        if(r2 < 0.15)
        {
            do{
                var2 = rng.Gaus(190.,4.);
            }while(var2<150 || var2>250);
        }
        else if(r2 < 0.25)
        {
            do{
                var2 = rng.Gaus(200.,5.);
            }while(var2<150 || var2>250);
        }
        else if(r2 < 0.33)
        {
            do{
                var2 = rng.Gaus(212.,6.);
            }while(var2<150 || var2>250);
        }
        else if(r2 < 0.50)
        {
            do{
                var2 = rng.Gaus(205.,18.);
            }while(var2<150 || var2>250);
        }
        else
        {
            do{
                var2 = rng.Exp(40.) + 150.;
            }while(var2>250);
        }

        //----------------------------------------------------
        // VAR3
        // Signal on resonance shoulder
        //----------------------------------------------------

        double r3 = rng.Rndm();

        if(r3 < 0.18)
        {
            do{
                var3 = rng.BreitWigner(91.,2.5);
            }while(var3<60 || var3>140);
        }
        else if(r3 < 0.55)
        {
            do{
                var3 = rng.BreitWigner(80.,15.);
            }while(var3<60 || var3>140);
        }
        else
        {
            while(true)
            {
                double x = rng.Uniform(60.,140.);
                double fx = pow(60./x,2.3);

                if(rng.Rndm()<fx)
                {
                    var3=x;
                    break;
                }
            }
        }

        //----------------------------------------------------
        // VAR4
        // Double-sided CB-like signal
        //----------------------------------------------------

        double r4 = rng.Rndm();

        if(r4 < 0.22)
        {
            do{
                var4 = rng.Gaus(300.,4.);

                if(rng.Rndm()<0.20)
                    var4 += rng.Exp(3.);

                if(rng.Rndm()<0.15)
                    var4 -= rng.Exp(2.);

            }while(var4<250 || var4>360);
        }
        else if(r4 < 0.42)
        {
            do{
                var4 = rng.Gaus(290.,20.);
            }while(var4<250 || var4>360);
        }
        else
        {
            do{
                var4 = rng.Exp(60.) + 250.;
            }while(var4>360);
        }

        //----------------------------------------------------
        // VAR5
        // Nightmare spectrum
        //----------------------------------------------------

        double r5 = rng.Rndm();

        if(r5 < 0.08)
        {
            do{
                var5 = rng.Gaus(150.,2.);
            }while(var5<120 || var5>260);
        }
        else if(r5 < 0.15)
        {
            do{
                var5 = rng.Gaus(165.,3.);
            }while(var5<120 || var5>260);
        }
        else if(r5 < 0.21)
        {
            do{
                var5 = rng.Gaus(180.,4.);
            }while(var5<120 || var5>260);
        }
        else if(r5 < 0.27)
        {
            do{
                var5 = rng.Gaus(220.,8.);
            }while(var5<120 || var5>260);
        }
        else if(r5 < 0.42)
        {
            do{
                var5 = rng.Gaus(170.,15.);
            }while(var5<120 || var5>260);
        }
        else
        {
            while(true)
            {
                double x = rng.Uniform(120.,260.);

                double fx =
                    pow(x-120.,1.2) *
                    exp(-(x-120.)/45.);

                if(rng.Rndm() < fx/180.)
                {
                    var5 = x;
                    break;
                }
            }
        }

        t->Fill();
    }

    t->Write();
    f->Close();

    printf("\nComplexFitSamples.root created\n");
    printf("Events per variable: %d\n",N);
}
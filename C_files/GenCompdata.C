#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TMath.h"

void GenCompdata()
{
    gSystem->Unlink("../Root_files/Compdata.root");
    const Long64_t NEvents = 1000000;

    TFile *outfile = new TFile("../Root_files/Compdata.root","RECREATE");
    TTree *tree = new TTree("Events","5-variable toy dataset");

    double var1, var2, var3, var4, var5;

    tree->Branch("var1",&var1,"var1/D");
    tree->Branch("var2",&var2,"var2/D");
    tree->Branch("var3",&var3,"var3/D");
    tree->Branch("var4",&var4,"var4/D");
    tree->Branch("var5",&var5,"var5/D");

    TRandom3 rng(0);

    for(Long64_t i=0;i<NEvents;i++)
    {
        double r = rng.Rndm();

        //----------------------------------
        // VAR1
        //----------------------------------
        if(r < 0.25)
        {
            do{
                var1 = rng.Gaus(150.,3.);

                if(rng.Rndm() < 0.12)
                    var1 += rng.Exp(4.);

            }while(var1 < 100. || var1 > 180.);
        }
        else
        {
            if(rng.Rndm() < 0.4)
                var1 = 100. + 80.*pow(rng.Rndm(),1.8);
            else
                var1 = 100. + 80.*rng.Rndm();
        }

        //----------------------------------
        // VAR2
        //----------------------------------
        double r2 = rng.Rndm();

        if(r2 < 0.15)
        {
            do{
                var2 = rng.BreitWigner(91.,2.5);
            }while(var2 < 60. || var2 > 140.);
        }
        else if(r2 < 0.45)
        {
            do{
                var2 = rng.BreitWigner(80.,10.);
            }while(var2 < 60. || var2 > 140.);
        }
        else
        {
            var2 = 60. + 80.*rng.Rndm();
        }

        //----------------------------------
        // VAR3
        //----------------------------------
        double r3 = rng.Rndm();

        if(r3 < 0.20)
        {
            do{
                var3 = rng.Gaus(250.,8.);
            }while(var3 < 200. || var3 > 320.);
        }
        else
        {
            var3 = 200. + 120.*pow(rng.Rndm(),0.7);
        }

        //----------------------------------
        // VAR4
        //----------------------------------
        double r4 = rng.Rndm();

        if(r4 < 0.25)
        {
            do{
                var4 = 300.
                     + 0.25*(var1 - 150.)
                     + rng.Gaus(0.,4.);
            }while(var4 < 250. || var4 > 360.);
        }
        else
        {
            var4 = 250. + 110.*rng.Rndm();
        }

        //----------------------------------
        // VAR5
        //----------------------------------
        double r5 = rng.Rndm();

        if(r5 < 0.10)
        {
            do{
                var5 = rng.Gaus(500.,15.);
            }while(var5 < 400. || var5 > 650.);
        }
        else
        {
            var5 = 400. + rng.Exp(40.);

            if(var5 > 650.)
                var5 = 650.;
        }

        tree->Fill();
    }

    tree->Write("", TObject::kOverwrite);
    std::cout << tree->GetAutoSave() << std::endl;
    outfile->ls();
    outfile->Close();

    printf("Generated %lld events in Compdata.root\n",NEvents);
}
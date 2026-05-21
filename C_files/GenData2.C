void GenData2() {

    // =====================================================
    // OUTPUT FILE
    // =====================================================

    TFile *fout = new TFile("../Root_files/mixdata.root", "RECREATE");
    TTree *t    = new TTree("data", "Mixed Signal and Background Data");

    // =====================================================
    // BRANCHES
    // =====================================================

    Double_t mass, energy;
    t->Branch("mass",   &mass,   "mass/D");    // invariant mass -- Z peak + expo bkg
    t->Branch("energy", &energy, "energy/D");  // energy -- Gaussian signal + poly bkg

    // =====================================================
    // RANDOM NUMBER GENERATOR
    // =====================================================

    TRandom3 rng(123);  // fixed seed -- reproducible

    // =====================================================
    // GENERATE EVENTS
    // =====================================================

    Int_t nEvents = 100000;

    for (Int_t i = 0; i < nEvents; i++) {

        // -------------------------------------------------
        // MASS: 30% signal (Breit-Wigner Z peak at 91.2 GeV)
        //       70% background (exponential)
        // -------------------------------------------------

        if (rng.Uniform() < 0.3) {
            // signal -- keep regenerating until in range
            do { mass = rng.BreitWigner(91.2, 2.5); }
            while (mass < 60 || mass > 120);
        } else {
            // background -- exponential falling from 60 GeV
            do { mass = rng.Exp(15) + 60; }
            while (mass < 60 || mass > 120);
        }

        // -------------------------------------------------
        // ENERGY: 40% signal (Gaussian peak at 50 GeV)
        //         60% background (quadratic polynomial shape)
        // -------------------------------------------------

        if (rng.Uniform() < 0.4) {
            // signal -- Gaussian peak
            do { energy = rng.Gaus(50.0, 3.0); }
            while (energy < 0 || energy > 100);
        } else {
            // background -- quadratic shape
            // rejection sampling to get polynomial distribution
            Double_t x, y;
            do {
                x = rng.Uniform(0, 100);
                // background shape: peaks at edges, dips in middle
                // f(x) = 0.5 + (x-50)^2 / 5000
                y = rng.Uniform(0, 1.5);
            } while (y > 0.5 + (x-50)*(x-50)/5000.0);
            energy = x;
        }

        t->Fill();
    }

    // =====================================================
    // WRITE AND CLOSE
    // =====================================================

    t->Write();
    fout->Close();

    cout << "Generated " << nEvents << " events" << endl;
    cout << "File saved: mixdata.root" << endl;
    cout << "Branches: mass (60-120 GeV), energy (0-100 GeV)" << endl;
    cout << "Mass:   30% BW signal + 70% expo background" << endl;
    cout << "Energy: 40% Gaussian signal + 60% polynomial background" << endl;
}
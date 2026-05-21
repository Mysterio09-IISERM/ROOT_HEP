void Gendata() {
    TFile *fout = new TFile("../Root_files/generatedata.root", "RECREATE");
    TTree *t = new TTree("data", "Simulated HEP data");

    Double_t mass, pt, eta, phi;
    t->Branch("mass", &mass, "mass/D");  // invariant mass
    t->Branch("pt",   &pt,   "pt/D");    // transverse momentum
    t->Branch("eta",  &eta,  "eta/D");   // pseudorapidity
    t->Branch("phi",  &phi,  "phi/D");   // azimuthal angle

    TRandom3 rng(42);  // fixed seed for reproducibility

    for (int i = 0; i < 50000; i++) {
        // Z boson signal peak at 91.2 GeV with width 2.5 GeV
        // plus exponential background
        if (rng.Uniform() < 0.3) {
            mass = rng.BreitWigner(91.2, 2.5);  // signal
        } else {
            mass = rng.Exp(20) + 60;             // background
        }
        pt  = rng.Exp(20);
        eta = rng.Gaus(0, 2.5);
        phi = rng.Uniform(-TMath::Pi(), TMath::Pi());
        t->Fill();
    }

    t->Write();
    fout->Close();
    cout << "File created: generatedata.root with 50000 events" << endl;
}
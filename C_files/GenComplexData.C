#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TMath.h"

void GenComplexData() {

    gSystem->Unlink("../Root_files/ComplexData.root");

    TFile *f = new TFile("../Root_files/ComplexData.root", "NEW");
    TTree *t = new TTree("data", "Complex HEP-like datasets");

    TRandom3 rng(42);

    double var1, var2, var3, var4, var5;
    double weight1, weight2, weight3, weight4, weight5;

    t->Branch("var1",    &var1,    "var1/D");
    t->Branch("var2",    &var2,    "var2/D");
    t->Branch("var3",    &var3,    "var3/D");
    t->Branch("var4",    &var4,    "var4/D");
    t->Branch("var5",    &var5,    "var5/D");
    t->Branch("weight1", &weight1, "weight1/D");
    t->Branch("weight2", &weight2, "weight2/D");
    t->Branch("weight3", &weight3, "weight3/D");
    t->Branch("weight4", &weight4, "weight4/D");
    t->Branch("weight5", &weight5, "weight5/D");

    const int N = 500000;

    for (int i = 0; i < N; i++) {

        // =====================================================
        // VAR1: B meson decay -- tiny signal buried in huge
        //       combinatorial background with threshold turn-on
        //
        // Signal:     1.5% -- narrow Gaussian at 5279 MeV (B0 mass)
        //             sigma = 3 MeV (typical Belle II resolution)
        // Background: 98.5% -- ARGUS function (threshold background)
        //             rises from threshold at 5289 MeV, falls sharply
        //             ARGUS: f(x) ∝ x*sqrt(1-(x/m0)^2)*exp(c*(1-(x/m0)^2))
        //
        // Physics: This is exactly the B0 → D*-π+ invariant mass
        //          spectrum at Belle II. The ARGUS background comes
        //          from the beam energy constraint endpoint.
        // =====================================================

        double r1 = rng.Rndm();

        if (r1 < 0.015) {
            // signal: B0 mass peak
            do { var1 = rng.Gaus(5279.0, 3.0); }
            while (var1 < 5200 || var1 > 5290);
            weight1 = 1.0;
        } else {
            // ARGUS background -- rejection sampling
            double m0 = 5289.0;  // kinematic endpoint
            double c  = -15.0;   // ARGUS shape parameter
            double x, y;
            do {
                x = rng.Uniform(5200, 5289);
                double xi = 1.0 - (x/m0)*(x/m0);
                double fmax = 0.5 * m0;
                y = rng.Uniform(0, fmax);
                double fx = x * sqrt(xi) * exp(c * xi);
                if (y < fx) { var1 = x; break; }
            } while (true);
            weight1 = 1.0;
        }

        // =====================================================
        // VAR2: Dalitz plot projection -- D0 → π+π-π0
        //       Invariant mass squared m²(π+π0)
        //
        // Signal:     40% -- ρ(770) resonance (BW, mean=0.59 GeV²,
        //                    width=0.18 GeV²) -- dominates Dalitz plot
        //             15% -- f0(500)/σ (broad BW, mean=0.25, width=0.30)
        //                    -- wide scalar resonance
        // Background: 45% -- phase space (polynomial in m²)
        //
        // Physics: Directly relevant to your IDC352 project on
        //          D0 → π+π-π0 CP asymmetry measurement at Belle II
        //          The ρ(770) and f0(500) are the dominant isobars
        // =====================================================

        double r2 = rng.Rndm();

        if (r2 < 0.40) {
            // rho(770): BW in m^2 space
            // m_rho = 0.775 GeV → m^2 = 0.60 GeV^2
            // Gamma_rho = 0.149 GeV → Gamma in m^2 ≈ 2*m*Gamma = 0.23
            double m2;
            do { m2 = rng.BreitWigner(0.60, 0.115); }
            while (m2 < 0.07 || m2 > 1.80);
            var2 = m2;
            weight2 = 1.0;
        } else if (r2 < 0.55) {
            // f0(500): very broad BW
            double m2;
            do { m2 = rng.BreitWigner(0.25, 0.20); }
            while (m2 < 0.07 || m2 > 1.80);
            var2 = m2;
            weight2 = 1.0;
        } else {
            // phase space: polynomial distribution in m^2
            // for D0 → 3pi, phase space ~ linear in m^2
            double m2;
            do {
                m2 = rng.Uniform(0.07, 1.80);
                double ps = 1.0 - 0.3*(m2 - 0.93);  // slight slope
                if (rng.Rndm() < ps/1.3) { var2 = m2; break; }
            } while (true);
            weight2 = 1.0;
        }

        // =====================================================
        // VAR3: Three overlapping resonances + rising background
        //       Simulates J/ψ region with ψ(2S) and X(3872)
        //
        // Signal 1:   12% -- J/ψ: narrow BW at 3097 MeV, width=0.09 MeV
        //                    (so narrow it looks like a spike -- detector
        //                     resolution of 10 MeV dominates → Voigtian)
        // Signal 2:    8% -- ψ(2S): Voigtian at 3686 MeV, BW width=0.3 MeV
        //                    detector sigma=8 MeV
        // Signal 3:    3% -- X(3872): extremely narrow, at 3872 MeV
        //                    width < 1.2 MeV → detector dominated
        // Background: 77% -- rising power law (non-resonant e+e- → μ+μ-)
        //
        // Physics: This is the dimuon invariant mass spectrum
        //          seen in B → K μ+μ- decays -- key for LFU tests
        // =====================================================

        double r3 = rng.Rndm();

        if (r3 < 0.12) {
            // J/psi: resolution dominated
            double true_m = rng.BreitWigner(3096.9, 0.09);
            var3 = rng.Gaus(true_m, 10.0);
            weight3 = 1.0;
        } else if (r3 < 0.20) {
            // psi(2S)
            double true_m = rng.BreitWigner(3686.1, 0.3);
            var3 = rng.Gaus(true_m, 8.0);
            weight3 = 1.0;
        } else if (r3 < 0.23) {
            // X(3872): tiny signal
            double true_m = rng.BreitWigner(3871.7, 1.0);
            var3 = rng.Gaus(true_m, 8.0);
            weight3 = 1.0;
        } else {
            // rising power law background
            double x;
            do {
                x = rng.Uniform(3000, 4000);
                // power law: f(x) ∝ (x/3000)^1.5
                double fx = pow(x/3000.0, 1.5);
                if (rng.Rndm() < fx/pow(4000.0/3000.0, 1.5)) {
                    var3 = x; break;
                }
            } while (true);
            weight3 = 1.0;
        }

        // =====================================================
        // VAR4: Signal with per-event mass resolution + weights
        //       Simulates missing energy decay: B → D* l ν
        //
        // The key challenge: each event has a DIFFERENT resolution
        // because the missing neutrino momentum is estimated
        // differently per event. This creates a non-Gaussian signal.
        //
        // Signal:     25% -- Gaussian with per-event sigma
        //                    mean=0 (beam constraint variable)
        //                    sigma drawn from Gamma distribution
        //                    → results in broader-than-Gaussian tails
        // Background: 75% -- combination of flat + peaking background
        //                    (peaking bkg from B → D* π misID)
        //
        // Weights: events weighted by sin²(θ) where θ is decay angle
        //          → non-uniform weights, physically motivated
        // =====================================================

        double r4 = rng.Rndm();

        if (r4 < 0.25) {
            // signal with variable resolution
            // sigma per event from Gamma(alpha=4, beta=0.5) → mean=2
            double sigma_evt = rng.Gaus(2.0, 0.5);
            if (sigma_evt < 0.3) sigma_evt = 0.3;
            var4 = rng.Gaus(0.0, sigma_evt);
            // weight by decay angle distribution
            double costheta = rng.Uniform(-1, 1);
            weight4 = 1.0 - costheta*costheta;  // sin^2(theta)
        } else if (r4 < 0.85) {
            // combinatorial background: exponential
            var4 = rng.Exp(2.0) * (rng.Rndm() < 0.5 ? 1 : -1);
            double costheta = rng.Uniform(-1, 1);
            weight4 = 0.5 + 0.5*costheta*costheta;  // different angular dist
        } else {
            // peaking background: misidentified B→D*π
            // appears near signal but shifted
            var4 = rng.Gaus(1.5, 0.8);
            weight4 = rng.Uniform(0.3, 0.7);
        }

        // =====================================================
        // VAR5: CP asymmetry measurement
        //       Simulates time-integrated CP asymmetry in
        //       B0 → K+π- decay
        //
        // The variable is the charge asymmetry proxy:
        // q = charge of kaon (+1 or -1)
        // We measure the asymmetry between q=+1 and q=-1 events
        //
        // Signal+:    18% -- B0bar → K+π- (favored)
        //                    reconstructed mass peak at 5279 MeV
        // Signal-:    12% -- B0 → K-π+ (suppressed by CP violation)
        //                    A_CP ≈ -0.083 (PDG value)
        // Background: 70% -- charge-symmetric combinatorial + 
        //                    charge-asymmetric detector efficiency
        //
        // The tricky part: background has a small charge asymmetry
        // (~2%) from detector effects that must be corrected
        //
        // Physics: This is the measurement that first showed
        //          direct CP violation in B decays (BaBar/Belle)
        // =====================================================

        double r5 = rng.Rndm();

        if (r5 < 0.18) {
            // B0bar → K+π- signal (favored mode)
            do { var5 = rng.Gaus(5279.0, 25.0); }
            while (var5 < 5100 || var5 > 5400);
            weight5 = 1.0;
        } else if (r5 < 0.30) {
            // B0 → K-π+ signal (CP suppressed)
            // A_CP = (N+ - N-)/(N+ + N-) = -0.083
            // N-/N+ = (1 + A_CP)/(1 - A_CP) ≈ 0.84
            do { var5 = rng.Gaus(5279.0, 25.0); }
            while (var5 < 5100 || var5 > 5400);
            weight5 = 0.84;  // CP suppression via weight
        } else {
            // background: ARGUS-like combinatorial
            double m0 = 5289.0;
            double c  = -10.0;
            double x, y;
            do {
                x = rng.Uniform(5100, 5289);
                double xi = 1.0 - (x/m0)*(x/m0);
                if (xi < 0) continue;
                double fx = x * sqrt(xi) * exp(c * xi);
                y = rng.Uniform(0, 800);
                if (y < fx) {
                    var5 = x; break;
                }
            } while (true);
            // small detector charge asymmetry in background
            weight5 = 1.0 + 0.02 * rng.Gaus(0, 1);
            if (weight5 < 0) weight5 = 0.01;
        }

        t->Fill();
    }

    f->cd();
    t->Write("", TObject::kOverwrite);
    f->Close();
    delete f;

    printf("ComplexData.root created!\n");
    printf("Events: %d\n", N);
    printf("\nBranch summary:\n");
    printf("var1: B0 mass (MeV) -- tiny signal + ARGUS background\n");
    printf("var2: m^2(pi+pi0) (GeV^2) -- rho+f0+phase space (Dalitz)\n");
    printf("var3: dimuon mass (MeV) -- J/psi+psi2S+X3872+power law\n");
    printf("var4: beam constraint variable -- variable resolution + weights\n");
    printf("var5: B→Kπ mass (MeV) -- CP asymmetry measurement\n");
}
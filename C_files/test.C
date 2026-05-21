void test(){

    TCanvas *oldc1 = (TCanvas*)gROOT->FindObject("c1");
    if (oldc1) {
        delete oldc1;
    }

    TFile *F1 = new TFile("../Root_files/generatedata.root", "READ");
    TTree *t1 = (TTree*)F1->Get("data");

    TH1F *hmass = new TH1F("hmass", "Mass Distribution; Mass [GeV]; Counts", 100, 60, 120);
    t1->Draw("mass >> hmass", "", "goff");
    hmass->SetDirectory(0);
    F1->Close();
    delete F1;

    TCanvas *c1 = new TCanvas("c1", "Mass Distribution", 800, 600);
    c1->Divide(3, 1);

    TF1 *ffit = new TF1("ffit", "breitwigner(0) + expo(3)", 60, 120);
    ffit->SetParameters(
        hmass->GetMaximum(),  // [0] amplitude
        91.2,                 // [1] mass initial guess
        2.5,                  // [2] width initial guess
        5.0,                  // [3] expo constant
        -0.05);               // [4] expo slope

    ffit->SetParLimits(0, 0, 100000);  // amplitude > 0
    ffit->SetParLimits(1, 85, 95);     // mass
    ffit->SetParLimits(2, 0.1, 10);    // width > 0

    c1->cd(1);
    TH1F *hraw = (TH1F*)hmass->Clone("hraw");
    hraw->SetLineColor(kBlue);
    hraw->Draw("E1");

    c1->cd(2);
    hmass->Fit(ffit, "RS0");
    hmass->SetTitle("Combined Fit; Mass [GeV]; Counts");
    hmass->SetLineColor(kBlue);
    hmass->Draw("E1");
    ffit->SetLineColor(kRed);
    ffit->Draw("same");

    TF1 *fbkg = new TF1("fbkg", "expo", 60, 120);
    fbkg->SetParameters(ffit->GetParameter(3), ffit->GetParameter(4));
    fbkg->SetLineColor(kGreen);
    fbkg->Draw("same");

    c1->cd(3);
    TH1F *hsig = (TH1F*)hmass->Clone("hsig");
    hsig->SetTitle("Background Subtracted Signal; Mass [GeV]; Counts");

    for(int i=1; i<=hsig->GetNbinsX(); i++){
        Double_t x = hsig->GetBinCenter(i);
        Double_t bkg = fbkg->Eval(x); //// * hsig->GetBinWidth(i);
        hsig->SetBinContent(i, hsig->GetBinContent(i) - bkg);
        Double_t err = hsig->GetBinError(i);
        hsig->SetBinError(i, err); 
    }

    hsig->SetLineColor(kBlue);
    hsig->Draw("E1");

    TF1 *fsig = new TF1("fsig", "breitwigner", 60, 120);
    fsig->SetParameters(
        ffit->GetParameter(0),  // amplitude
        ffit->GetParameter(1),  // mass
        ffit->GetParameter(2)); // width
    fsig->SetLineColor(kRed);
    fsig->Draw("same");

    c1->Update();
}
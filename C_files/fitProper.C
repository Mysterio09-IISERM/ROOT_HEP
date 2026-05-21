void fitProper() {
    if (gROOT->FindObject("c1")) delete gROOT->FindObject("c1");
    if (gROOT->FindObject("h1")) delete gROOT->FindObject("h1");
    if (gROOT->FindObject("fgaus")) delete gROOT->FindObject("fgaus");
    if (gROOT->FindObject("c2")) delete gROOT->FindObject("c2");
    if (gROOT->FindObject("h2")) delete gROOT->FindObject("h2");
    if (gROOT->FindObject("flan1")) delete gROOT->FindObject("flan1");

    // step 1: get data
    TFile *f = new TFile("../Root_files/tmva_class_example.root", "READ");
    //TFile *f = TFile::Open("https://root.cern/files/tmva_class_example.root");
    TTree *t1 = (TTree*)f->Get("TreeS");
    TTree *t2 = (TTree*)f->Get("TreeB");

    // step 2: fill histogram
    TH1F *h1 = new TH1F("h1", "Signal var1; var1; Counts", 50, -5, 5);
    t1->Draw("var1 >> h1");
    h1->SetDirectory(0);  // detach from file to keep it alive after file is closed
    TH1F *h2 = new TH1F("h2", "Background var1; var1; Counts", 50, -5, 5);
    t2->Draw("var1 >> h2","","weight");
    h2->SetDirectory(0);  // detach from file to keep it alive after file is closed
    TH1F *h3 = new TH1F("h3", "Signal var2; var2; Counts", 50, -5, 5);
    t1->Draw("var2 >> h3");
    h3->SetDirectory(0);  // detach from file to keep it alive after file is closed
    TH1F *h4 = new TH1F("h4", "Background var2; var2; Counts", 50, -5, 5);
    t2->Draw("var2 >> h4","","weight");
     h4->SetDirectory(0);  // detach from file to keep it alive after file is closed
    TH1F *h5 = new TH1F("h5", "Signal var3; var3; Counts", 50, -5, 5);
     t1->Draw("var3 >> h5");
     h5->SetDirectory(0);  // detach from file to keep it alive after file is closed
    TH1F *h6 = new TH1F("h6", "Background var3; var3; Counts", 50, -5, 5);
     t2->Draw("var3 >> h6","","weight");
     h6->SetDirectory(0);  // detach from file to keep it alive after file is closed
     TH1F *h7 = new TH1F("h7", "Signal var4; var4; Counts", 50, -5, 5);
     t1->Draw("var4 >> h7");
     h7->SetDirectory(0);  // detach from file to keep it alive after file is closed
    TH1F *h8 = new TH1F("h8", "Background var4; var4; Counts", 50, -5, 5);
     t2->Draw("var4 >> h8","","weight");
     h8->SetDirectory(0);  // detach from file to keep it alive after file is closed

    f->Close();  // close file now that histograms are filled and detached
    
    // step 3: look at the data BEFORE fitting
    TCanvas *c1 = new TCanvas("c1", "Signal Fit", 800, 600);
    TCanvas *c2 = new TCanvas("c2", "Background Fit", 800, 600);
    c1->Divide(2,2);
    c2->Divide(2,2);
    c1->cd(1);
    h1->Draw();
    c1->cd(2);
    h3->Draw();
     c1->cd(3);
    h5->Draw();
     c1->cd(4);
    h7->Draw();
    c1->Update();
    c2->cd(1);
    h2->Draw();
     c2->cd(2);
    h4->Draw();
     c2->cd(3);
    h6->Draw();
     c2->cd(4);
     h8->Draw();
    c2->Update();



    // step 4: define function with range
    TF1 *fgaus = new TF1("fgaus", "gaus", -5, 5);

    // step 5: set initial guesses from histogram properties
    fgaus->SetParameter(0, h1->GetMaximum());  // amplitude ~ peak height
    fgaus->SetParameter(1, h1->GetMean());     // mean ~ center of distribution
    fgaus->SetParameter(2, h1->GetRMS());      // sigma ~ spread

    cout << "Initial guesses:" << endl;
    cout << "Amplitude = " << h1->GetMaximum() << endl;
    cout << "Mean      = " << h1->GetMean() << endl;
    cout << "Sigma     = " << h1->GetRMS() << endl;

    // step 6: fit
    // "R" = use function range
    // "S" = store result
    // "Q" = quiet (remove for verbose output while learning)
    h1->Fit(fgaus, "RS");

    // step 7: draw fit on top
    c1->cd(1);
    fgaus->SetLineColor(kRed);
    fgaus->Draw("same");

    // step 8: read and understand results
    cout << "\nFit Results:" << endl;
    cout << "Amplitude = " << fgaus->GetParameter(0) 
         << " +/- " << fgaus->GetParError(0) << endl;
    cout << "Mean      = " << fgaus->GetParameter(1) 
         << " +/- " << fgaus->GetParError(1) << endl;
    cout << "Sigma     = " << fgaus->GetParameter(2) 
         << " +/- " << fgaus->GetParError(2) << endl;
    cout << "Chi2/NDF  = " << fgaus->GetChisquare() << " / " << fgaus->GetNDF() << endl;
    cout << "Fit Prob  = " << fgaus->GetProb() << endl;

/*
    // step 4: define function with range
    TF1 *flan1 = new TF1("flan1", "landau", -5, 5);

    // step 5: set initial guesses from histogram properties
    // CORRECT -- only 2 parameters for expo
    flan1->SetParameter(0, h2->GetMaximum());  // amplitude
    flan1->SetParameter(1, h2->GetMean());     // mean
    flan1->SetParameter(2, h2->GetRMS());      // sigma

    cout << "\nInitial guesses:" << endl;
    cout << "p0 = " << h2->GetMaximum() << endl;
    cout << "p1 = " << h2->GetMean() << endl;
    cout << "p2 = " << h2->GetRMS() << endl;

    // step 6: fit
    // "R" = use function range
    // "S" = store result
    // "Q" = quiet (remove for verbose output while learning)
    h2->Fit(flan1, "RS");

    // step 7: draw fit on top
    flan1->SetLineColor(kBlue);
    flan1->Draw("same");

    // step 8: read and understand results
    cout << "\nFit Results:" << endl;
    cout << "p0 = " << flan1->GetParameter(0) 
         << " +/- " << flan1->GetParError(0) << endl;
    cout << "p1 = " << flan1->GetParameter(1) 
         << " +/- " << flan1->GetParError(1) << endl;
    cout << "p2 = " << flan1->GetParameter(2) 
         << " +/- " << flan1->GetParError(2) << endl;
    cout << "p3 = " << flan1->GetParameter(3) 
         << " +/- " << flan1->GetParError(3) << endl;
    cout << "Chi2/NDF  = " << flan1->GetChisquare() 
         << " / " << flan1->GetNDF() << endl;
    cout << "Fit Prob  = " << flan1->GetProb() << endl;*/
}
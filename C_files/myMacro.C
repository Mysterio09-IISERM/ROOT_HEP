void myMacro() {
    if (gROOT->FindObject("c1")) delete gROOT->FindObject("c1");
    if (gROOT->FindObject("f1")) delete gROOT->FindObject("f1");

    TCanvas *c1 = new TCanvas("c1", "My Plot", 800, 600);  // explicit creation

    TF1 *f1 = new TF1("f1", "sin(x)**2 + cos(x)**3", -10.0, 10.0);
    f1->SetLineColor(kRed);
    f1->SetNpx(100000);
    f1->Draw();

    c1->Update();
}
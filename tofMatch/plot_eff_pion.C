#include "../myFunction.h"
void plot_eff_pion() {
    globalSetting();
    TCanvas* c1 = new TCanvas("c1", "A Canvas",10,10,800,800);
    setPad(c1);
    char buf[250];
    char dir[250];
    char name[250];
    
    //read eff from txt file
    string str;
    float n_tmp;
    ifstream in("data/pionPid_eff.txt");
    getline(in,str);
    cout << str << endl;
    in >> n_tmp;
    const int npt = n_tmp;
    float pt_mean[npt], pt_err[npt], tpcEff[npt], tpcEffErr[npt], tofEff[npt], tofEffErr[npt];
    for(int ipt=0; ipt<npt; ipt++) {
        if(in.eof()) { cout << "read txt file error, check the file!!!" << endl; exit(1);}
        in >> pt_mean[ipt] >> pt_err[ipt];
        in >> tpcEff[ipt] >> tpcEffErr[ipt];
        in >> tofEff[ipt] >> tofEffErr[ipt];
    }
    in.close();
    
    // PID eff.
    const float markersize = 1.2;
    const float linewidth = 2;
    TGraphErrors* gTofEff = new TGraphErrors(npt, pt_mean, tofEff, pt_err, tofEffErr);
    TGraphErrors* gTpcEff = new TGraphErrors(npt, pt_mean, tpcEff, pt_err, tpcEffErr);
    gTofEff->SetMarkerStyle(kFullCircle);
    gTofEff->SetMarkerSize(markersize);
    gTofEff->SetMarkerColor(kBlack);
    gTofEff->SetLineWidth(linewidth);
    gTofEff->SetLineColor(kBlack);
    gTpcEff->SetMarkerStyle(kFullCircle);
    gTpcEff->SetMarkerSize(markersize);
    gTpcEff->SetMarkerColor(kBlack);
    gTpcEff->SetLineWidth(linewidth);
    gTpcEff->SetLineColor(kBlack);
    
    // function for fit
    TF1* fTpcEff = new TF1("fTpcEff","pol0",0,4);
    fTpcEff->SetParameter(0,0.95);
    fTpcEff->SetLineColor(kRed);
    fTpcEff->SetLineWidth(linewidth);
    TF1* fTofEff = new TF1("fTofEff","[0]+[1]/x-[2]/x/x",0,4);
    fTofEff->SetParameters(1,0.1,0.1);
    fTofEff->SetLineColor(kRed);
    fTofEff->SetLineWidth(linewidth);
    
    // fit
    const char fitOpt[10] = "NOR";
    const float fitR_lw = 0.2;
    const float fitR_up = 4.;
    const float fitR_lw1 = 0.2;
    const float fitR_up1 = 0.35;
    gTpcEff->Fit(fTpcEff,fitOpt,"",fitR_lw,fitR_up);
    gTpcEff->Fit(fTpcEff,fitOpt,"",fitR_lw,fitR_up);
    gTpcEff->Fit(fTpcEff,fitOpt,"",fitR_lw,fitR_up);
    gTofEff->Fit(fTofEff,fitOpt,"",fitR_lw1,fitR_up1);
    gTofEff->Fit(fTofEff,fitOpt,"",fitR_lw1,fitR_up1);
    gTofEff->Fit(fTofEff,fitOpt,"",fitR_lw1,fitR_up1);
    
    // plot
    TH1F* h0 = new TH1F("","",1,0,4);
    setHisto(h0, "", "p (GeV/c)", "n#sigma_{#pi} cut eff.");
    h0->GetYaxis()->SetRangeUser(0.7,1.02);
    h0->Draw();
    gTpcEff->Draw("psame");
    fTpcEff->Draw("same");
    sprintf(name,"Pion TPC PID eff. = %.3f",fTpcEff->GetParameter(0));
    drawLatex(0.2,0.2,name,132,0.05,1);
    c1->SaveAs("pic/TpcPidEff_pi.gif");
    
    TH1F* h1 = new TH1F("","",1,0,4);
    setHisto(h1, "", "p (GeV/c)", "1/#beta_{#pi} cut eff.");
    h1->GetYaxis()->SetRangeUser(0.7,1.02);
    h1->Draw();
    gTofEff->Draw("psame");
    //fTofEff->Draw("same");
    sprintf(name,"Pion TOF PID eff.");
    drawLatex(0.2,0.2,name,132,0.05,1);
    c1->SaveAs("pic/TofPidEff_pi.gif");
    
    // write to root file
    TFile* fout = new TFile("data/PIDeff_pi.root","RECREATE");
    gTpcEff->Write("gTpcEff_pi");
    gTofEff->Write("gTofEff_pi");
    fout->Close();
}

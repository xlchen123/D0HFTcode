#include "../anaCuts.h"
#include "../myConst.h"
#include "../myFunction.h"
void doPtShift1() {
    globalSetting();
    char name[250];
    char CMD[250];
    char dir[250];
    TLegend* legend;
    TH1F* h0;
    TCanvas* c1 = new TCanvas("c1", "A Canvas",10,10,800,800);
    setPad(c1);
    
    sprintf(dir,"pic");
    sprintf(CMD,"[ -d %s ] || mkdir -p %s",dir,dir);
    gSystem->Exec(CMD);
    
    //read spectra without pt shift
    TFile* fin = new TFile("../sys/D0_Spectra_Run14_HFT_beforePtShift.root");
    TGraphErrors* gD0_err[ncent];
    TGraphErrors* gD0_sys[ncent];
    for(int icent=0; icent<ncent; icent++) {
        gD0_err[icent] = (TGraphErrors*)fin->Get(Form("gD0_err_%s",nameCent1[icent]));
        gD0_sys[icent] = (TGraphErrors*)fin->Get(Form("gD0_sys_%s",nameCent1[icent]));
    }
    fin->Close();
    
    //read published
    TFile* fin2 = new TFile("../../../physics/D0_AuAu_data_publish.root");
    gD0_shifted = (TGraphErrors*)fin2->Get("cen_0_10_err");
    fin2->Close();
    
    //define fit function, now use levy function
    char funcString[200];
    char funcString_time_pt[200];
    double m0 = 1.8645;//D0-1.8645, D+/- - 1.8693;
    sprintf(funcString,"(1/(2*TMath::Pi()))*[0]*([2]-1)*([2]-2)/([2]*[1]*([2]*[1]+%f*([2]-2)))*TMath::Power(([2]*[1]+TMath::Sqrt(x[0]*x[0]+%f*%f)-%f)/([2]*[1]),-[2])",m0,m0,m0,m0); // dN/pTdpTdy
    sprintf(funcString_time_pt,"(1/(2*TMath::Pi()))*[0]*([2]-1)*([2]-2)/([2]*[1]*([2]*[1]+%f*([2]-2)))*TMath::Power(([2]*[1]+TMath::Sqrt(x[0]*x[0]+%f*%f)-%f)/([2]*[1]),-[2])*x[0]",m0,m0,m0,m0); // dN/dpTdy
    TF1 *flevy = new TF1("flevy",funcString,0,20);
    flevy->SetParameters(6.25342e-01,3.34550e-01,1.86042e+01);
    flevy->SetLineColor(kGreen+2);
    // flevy->SetLineStyle(9);
    TF1 *flevy_time_pt = new TF1("flevy_time_pt",funcString_time_pt,0,20);
    
    // ==== Create a TGraphErrors to hold the confidence intervals
    const Int_t NCL = 200;
    TGraphErrors *grint;
    grint = new TGraphErrors(NCL);
    grint->SetTitle("Fitted line with .68 conf. band");
    for (int i=0; i<NCL; i++) grint->SetPoint(i, 0 + (float)i/20., 0);
    
    //shift as published
    
    
    // ==== Fit and Write in centrality bin
    TFile* fout = new TFile("D0_Spectra_Run14HFT.root","RECREATE");
    for(icent=0; icent<ncent; icent++) {
        Int_t NPoint = gD0_err[0]->GetN();
        for(int ip=0;ip<NPoint;ip++)
        {
            gD0_err[icent]->GetX()[ip] = gD0_shifted->GetX()[ip];
            gD0_sys[icent]->GetX()[ip] = gD0_shifted->GetX()[ip];
        }
        
        flevy->SetParameters(6.25342e-01,3.34550e-01,1.86042e+01);
        for(int ifit=0;ifit<3;ifit++) gD0_err[icent]->Fit(flevy,"INOR","",0,10);
        grint = new TGraphErrors(NCL);
        grint->SetTitle("Fitted line with .68 conf. band");
        for (int i=0; i<NCL; i++) grint->SetPoint(i, 0 + (float)i/20., 0);
        (TVirtualFitter::GetFitter())->GetConfidenceIntervals(grint,0.68);
        grint->SetName(Form("flevy_err_band_%s",nameCent1[icent]));
        
        fout->cd();
        gD0_err[icent]->Write(Form("gD0_err_%s",nameCent1[icent]));
        gD0_sys[icent]->Write(Form("gD0_sys_%s",nameCent1[icent]));
        flevy->Write(Form("flevy_%s",nameCent1[icent]));
        grint->Write();
    }
    
    fout->Close();
}

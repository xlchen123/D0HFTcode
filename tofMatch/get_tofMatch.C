#include "../myFunction.h"
#include "../myConst.h"
#include "../anaCuts.h"
TH1F* hPi;
TH1F* hK;
double fpi(double* x, double* par) {
    if(x[0]>0.4 && x[0]<1.8) {
        int bin = hPi->FindBin(x[0]);
        return hPi->GetBinContent(bin);
    }
    else {
        return par[0]/(pow(x[0]-par[1],2)+par[2])-par[4]/(exp(x[0]-par[3])+par[5])+par[6];
    }
}
double fk(double* x, double* par) {
    if(x[0]>0.45 && x[0]<2.1) {
        int bin = hK->FindBin(x[0]);
        return hK->GetBinContent(bin);
    }
    else {
        return par[0]/(pow(x[0]-par[1],2)+par[2])-par[4]/(exp(x[0]-par[3])+par[5])+par[6];
    }
    
}
void get_tofMatch() {
    globalSetting();
    char buf[250];
    char dir[250];
    char name[250];
    TLegend* legend;
    
    TCanvas *c1 = new TCanvas("c1", "c1",10,10,800,800);
    setPad(c1);
    
    //read hist
    TFile* fin = new TFile("data/tofMatch_Run14.root");
    TH2F* h2pT_tpcPi = (TH2F*)fin->Get("h2pTCent_tpcPi");
    h2pT_tpcPi->SetDirectory(0);
    TH2F* h2pT_tofPi = (TH2F*)fin->Get("h2pTCent_tofPi");
    h2pT_tofPi->SetDirectory(0);
    TH2F* h2pT_tpcK = (TH2F*)fin->Get("h2pTCent_tpcK");
    h2pT_tpcK->SetDirectory(0);
    TH2F* h2pT_tofK = (TH2F*)fin->Get("h2pTCent_tofK");
    h2pT_tofK->SetDirectory(0);
    fin->Close();
    
    //calculate tofMatch eff 2D
    TH2F* h2tofEff_pi = new TH2F(*h2pT_tpcPi);
    h2tofEff_pi->Divide(h2pT_tofPi,h2pT_tpcPi,1,1);
    TH2F* h2tofEff_k = new TH2F(*h2pT_tpcK);
    h2tofEff_k->Divide(h2pT_tofK,h2pT_tpcK,1,1);
    
    
    //caculate eff in cent
    TH1F* hefftmp_pi;
    TH1F* hefftmp_k;
    TH1F* heff_pi[ncent_tof];
    TH1F* heff_k[ncent_tof];
    for(int i=0; i<9; i++) {
        for(int icent=0; icent<ncent_tof; icent++) {
            if(i>=centLw_tof[icent]-1 && i<centUp_tof[icent]-1)
                NbinSum_tof[icent] += Nbin[i];
        }
    }
    for(int icent=0; icent<ncent_tof; icent++) {
        sprintf(name,"htofMatchCent%i_pi",icent);
        heff_pi[icent] = new TH1F(name,name,npt_tof,ptLw_tof,ptUp_tof);
        heff_pi[icent]->Sumw2();
        sprintf(name,"htofMatchCent%i_kOri",icent);
        heff_k[icent] = new TH1F(name,name,npt_tof,ptLw_tof,ptUp_tof);
        heff_k[icent]->Sumw2();
        
    }
    for(int i=0; i<9; i++) {
        hefftmp_pi = (TH1F*)h2tofEff_pi->ProjectionX(Form("htofMatchCenttmp%i_pi",i),i+2,i+2);
        hefftmp_k = (TH1F*)h2tofEff_k->ProjectionX(Form("htofMatchCenttmp%i_kOri",i),i+2,i+2);
        for(int icent=0; icent<ncent_tof; icent++) {
            if(i>=centLw_tof[icent]-1 && i<centUp_tof[icent]-1) {
                heff_pi[icent]->Add(hefftmp_pi,Nbin[i]/NbinSum_tof[icent]);
                heff_k[icent]->Add(hefftmp_k,Nbin[i]/NbinSum_tof[icent]);
            }
        }
    }
    
    //fit and draw and correct kaon tof match
    char dir[250];
    char CMD[250];
    sprintf(dir,"pic");
    sprintf(CMD,"[ -d %s ] || mkdir -p %s", dir, dir);
    gSystem->Exec(CMD);
    
    TF1* funk1 = new TF1("funk1",fk,0.2,10,7);
    funk1->SetParameters(-9.10686e-02,4.30255e+00,1.35733e+01,1.48858e+00,-4.18231e-02,-2.15280e-01,6.79806e-01);
    TF1* fK = new TF1("funk","[0]/(pow(x-[1],2)+[2])-[4]/(exp(x-[3])+[5])+[6]",0.2,10);
    TH1F* heffCorr_k[ncent_tof];
    for(int icent=0; icent<ncent_tof; icent++) {
        hK = (TH1F*)heff_k[icent]->Clone(Form("htofMatchCent%i_kClone",icent));
        hK->Fit(funk1,"NOR","",0.3,5);
        //hK->Fit(funk1,"NOR","",0.3,5);
        //hK->Fit(funk1,"NOR","",0.3,5);
        for(int ipar=0; ipar<7; ipar++) {
            fK->SetParameter(ipar,funk1->GetParameter(ipar));
        }
        heff_pi[icent]->SetLineColor(COLOR[0]);
        heff_pi[icent]->SetLineWidth(2);
        setHisto(heff_pi[icent],"","p_{T} (GeV/c)", "TOF Matching Efficiency");
        heff_pi[icent]->GetXaxis()->SetRangeUser(0,4);
        heff_pi[icent]->GetYaxis()->SetRangeUser(0,1);
        heff_k[icent]->SetLineColor(COLOR[1]);
        heff_k[icent]->SetLineWidth(2);
        setHisto(heff_k[icent],"","p_{T} (GeV/c)", "TOF Matching Efficiency");
        heff_k[icent]->GetXaxis()->SetRangeUser(0,4);
        heff_k[icent]->GetYaxis()->SetRangeUser(0,1);
        fK->SetRange(0.3,4);
        fK->SetLineColor(COLOR[2]);
        fK->SetLineStyle(7);
        fK->SetLineWidth(2);
        
        heff_pi[icent]->Draw("HIST");
        heff_k[icent]->Draw("HISTSAME");
        fK->Draw("same");
        TLegend* lg1 = new TLegend(0.65,0.75,0.9,0.9);
        lg1->SetFillStyle(0);
        lg1->SetFillColor(10);
        lg1->SetBorderSize(0);
        lg1->SetTextSize(0.04);
        lg1->AddEntry(heff_pi[icent],"#pi","l");
        lg1->AddEntry(heff_k[icent],"k","l");
        lg1->AddEntry(fK,"k Fit","l");
        lg1->Draw();
        drawLatex(0.2,0.85,nameCent_tof[icent],132,0.05,1);
        sprintf(name,"%s/tofMatch_%i.gif",dir,icent);
        c1->SaveAs(name);
        
        //correct
        sprintf(name,"htofMatchCent%i_k",icent);
        heffCorr_k[icent] = new TH1F(name,name,npt_tof,ptLw_tof,ptUp_tof);
        setHisto(heffCorr_k[icent],"","p_{T} (GeV/c)", "TOF Matching Efficiency");
        heffCorr_k[icent]->GetXaxis()->SetRangeUser(0,4);
        heffCorr_k[icent]->GetYaxis()->SetRangeUser(0,1);
        heffCorr_k[icent]->SetLineColor(COLOR[2]);
        heffCorr_k[icent]->SetLineWidth(2);
        float binWidth = heffCorr_k[icent]->GetBinWidth(2);
        for(int i=0; i<heffCorr_k[icent]->GetNbinsX(); i++) {
            float pt = (i+0.5)*binWidth;
            if(pt>0.4 && pt<3.0) {
                heffCorr_k[icent]->SetBinContent(i+1,fK->Eval(pt));
                heffCorr_k[icent]->SetBinError(i+1,0);
            }
            else {
                heffCorr_k[icent]->SetBinContent(i+1,heff_k[icent]->GetBinContent(i+1));
                heffCorr_k[icent]->SetBinError(i+1,0);
            }
        }
    }
    
    //write in centrality
    TFile* fout = new TFile(Form("TOF_Match_kpi.root"),"RECREATE");
    for(int icent=0; icent<ncent_tof; icent++) {
        heff_pi[icent]->Write();
        heff_k[icent]->Write();
        heffCorr_k[icent]->Write();
    }
    fout->Close();
    
 
    /*TF1* funk1 = new TF1("funk1",fk,0.2,10,7);
    funk1->SetParameters(-9.10686e-02,4.30255e+00,1.35733e+01,1.48858e+00,-4.18231e-02,-2.15280e-01,6.79806e-01);
    TF1* fK = new TF1("funk","[0]/(pow(x-[1],2)+[2])-[4]/(exp(x-[3])+[5])+[6]",0.2,3);
    for(int icent=0; icent<9; icent++) {
        for(int icent=0; icent<ncent_tof; icent++) {
            
        }
        int cent_lw = h2tofEff_pi
        heff_pi = (TH1F*)h2tofEff_pi->ProjectionX(Form("htofMatchCent%i_pi",icent),icent+2,icent+2);
        heff_pi->SetTitle(nameCent[icent]);
        heff_pi->GetXaxis()->SetRangeUser(1,4);
        heff_pi->GetYaxis()->SetRangeUser(0,1);
        heff_pi->Write();
        heff_k = (TH1F*)h2tofEff_k->ProjectionX(Form("htofMatchCent%i_kOri",icent),icent+2,icent+2);
        heff_k->SetTitle(nameCent[icent]);
        heff_k->GetXaxis()->SetRangeUser(1,4);
        heff_k->GetYaxis()->SetRangeUser(0,1);
        heff_k->Write();
        
        
    }
    
    
    */
}

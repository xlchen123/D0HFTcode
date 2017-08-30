#include "../myFunction.h"
#include "../myConst.h"
#include "../anaCuts.h"
void plot_eff()
{
    globalSetting();
    char name[250];
    char CMD[250];
    char dir[250];
    TLegend* legend;
    TH1F* h0;
    TCanvas* c1 = new TCanvas("c1", "A Canvas",10,10,800,800);
    setPad(c1);
    
    sprintf(dir,"pic/eff");
    sprintf(CMD,"[ -d %s ] || mkdir -p %s",dir,dir);
    gSystem->Exec(CMD);
    
    //read eff
    TFile* feff = new TFile("data/eff.root");
    TH1F* heffD0[ncent];
    TH1F* heffBinD0[ncent];
    for(int icent=0; icent<ncent; icent++) {
        sprintf(name,"heffD0_%i",icent);
        heffD0[icent] = (TH1F*)feff->Get(name);
        heffD0[icent]->SetDirectory(0);
        sprintf(name,"heffBinD0_%i",icent);
        heffBinD0[icent] = (TH1F*)feff->Get(name);
        heffBinD0[icent]->SetDirectory(0);
    }
    feff->Close();
    
    //plot eff
    h0 = new TH1F("","",10,0,8);
    setHisto(h0, "", "p_{T} (GeV/c)", "efficiency");
    h0->GetYaxis()->SetRangeUser(0, 0.08);
    h0->GetXaxis()->SetRangeUser(0,8);
    h0->Draw();
    float leg_up = 0.93;
    float leg_lw = leg_up - ncent*0.05;
    TLegend* legend1 = new TLegend(.2,leg_lw,.4,leg_up);
    legend1->SetTextFont(132);
    legend1->SetTextSize(0.042);
    legend1->SetBorderSize(0);
    legend1->SetFillStyle(0);
    //legend1->SetHeader("     D^{0}_{prompt}");
    for(int icent=0; icent<ncent; icent++) {
        heffD0[icent]->Draw("esame");
        heffD0[icent]->SetLineColor(COLOR[icent]);
        heffD0[icent]->SetLineWidth(2);
        //heffBinD0[icent]->Draw("esame");
        heffBinD0[icent]->SetLineColor(COLOR[icent]);
        heffBinD0[icent]->SetLineWidth(2);
        
        legend1->AddEntry(heffD0[icent],nameCent[icent],"lp");
        //legend2->AddEntry(heffBinD0[icent],nameCent[icent],"p");
    }
    legend1->Draw();
    //legend2->Draw();
    sprintf(name,"%s/effAll.gif",dir);
    c1->SaveAs(name);
}

#include "../myFunction.h"
#include "../myConst.h"
void plot_spectra() {
    globalSetting();
    char dir[250];
    char name[250];
    char title[250];
    //TString CMD = 0;
    char CMD[250];
    TLegend* legend1;
    TLegend* legend2;
    TH1F* h0;
    
    sprintf(dir,"pic");
    sprintf(CMD,"[ -d %s ] || mkdir -p %s",dir,dir);
    gSystem->Exec(CMD);
    
    const int ncent = 6;
    const char nameCent[ncent][250] = {"0-80%","0-10%", "10-20%", "20-40%", "40-60%", "60-80%"};
    const char nameCent1[ncent][250] = {"0-80% [#times8]","0-10% [#times1]", "10-20% [/3]", "20-40% [/6]", "40-60% [/12]", "60-80% [/24]"};
    const char nameCentXL[ncent][250] = {"0_80","0_10", "10_20", "20_40", "40_60", "60_80"};
    const float scale[ncent] = {8., 1., 1/3., 1/6., 1./12, 1./24};
    
    // const int ncent = 5;
    // const char nameCent[ncent][250] = {"0-10%", "10-20%", "20-40%", "40-60%", "60-80%"};
    // const char nameCent1[ncent][250] = {"0-10%", "10-20% [/3]", "20-40% [/6]", "40-60% [/12]", "60-80% [/24]"};
    // const char nameCentXL[ncent][250] = {"0_10", "10_20", "20_40", "40_60", "60_80"};
    // const float scale[ncent] = {1., 1/3., 1/6., 1./12, 1./24};
    

    // const int ncent = 4;
    // const char nameCent[ncent][250] = {"0-10%", "10-40%", "40-80%", "0-80%"};
    // const char nameCent1[ncent][250] = {"0-10%", "10-40% [/3]", "40-80% [/12]", "0-80% [/24]"};
    // const char nameCentXL[ncent][250] = {"0_10", "10_40", "40_80", "0_80"};
    // const float scale[ncent] = {1., 1/3., 1./12, 1./24};

    //Read spectra
    //1. from xiaolong
    TGraphErrors* gD0err_xl[ncent];
    TGraphErrors* gD0sys_xl[ncent];
    TF1* fLevy[ncent];
    TFile* fin1 = new TFile("D0_Spectra_Run14HFT.root");
    for(int icent=0; icent<ncent; icent++) {
        gD0err_xl[icent] = (TGraphErrors*)fin1->Get(Form("gD0_err_%s",nameCentXL[icent]));
        gD0sys_xl[icent] = (TGraphErrors*)fin1->Get(Form("gD0_sys_%s",nameCentXL[icent]));
        fLevy[icent] = (TF1*)fin1->Get(Form("flevy_%s",nameCentXL[icent]));
    }
    fin1->Close();
    
    //scale
    for(int icent=0; icent<ncent; icent++) {
        ScaleGraph(gD0err_xl[icent], scale[icent]);
        ScaleGraph(gD0sys_xl[icent], scale[icent]);
        fLevy[icent]->SetParameter(0,fLevy[icent]->GetParameter(0)*scale[icent]);
    }
    
    //Set for Draw
    float markerSize = 2.0;
    float lineWidth = 2;
    for(int icent=0; icent<ncent; icent++) {
        gD0err_xl[icent]->SetMarkerStyle(kFullCircle);
        gD0err_xl[icent]->SetMarkerSize(markerSize);
        gD0err_xl[icent]->SetMarkerColor(COLOR[icent]);
        gD0err_xl[icent]->SetLineWidth(lineWidth);
        gD0err_xl[icent]->SetLineColor(COLOR[icent]);
        
        fLevy[icent]->SetLineColor(COLOR[icent]);
        fLevy[icent]->SetLineWidth(lineWidth);
    }
    
    //plot
    TCanvas* c1 = new TCanvas("c1", "A Canvas",10,10,800,800);
    setPad(c1);
    float ymin = 5e-10;
    float ymax = 0.99e0;
    // h0= new TH1F("","",1,0,8);
    h0= new TH1F("","",1,0,10);
    h0->Draw();
    h0->SetMinimum(ymin),
    h0->SetMaximum(ymax);
    setHisto(h0,"","p_{T} (GeV/c)", "d^{2}N/(N_{ev}2#pip_{T}dp_{T}dy) (GeV/c)^{-2}");
    legend1 = new TLegend(0.68,0.63,0.9,0.93);
    legend1->SetFillStyle(0);
    legend1->SetFillColor(10);
    legend1->SetBorderSize(0);
    legend1->SetTextSize(0.04);
    legend1->SetTextFont(132);
    //legend1->SetHeader("Xiaolong");
    for(int icent=0; icent<ncent; icent++) {
        //legend1->AddEntry(gD0err_xl[icent],Form("%s [#times%.1f]", nameCent[icent],scale[icent]),"p");
        legend1->AddEntry(gD0err_xl[icent],nameCent1[icent],"p");
        if( icent !=2 ) fLevy[icent]->Draw("same");
        gD0err_xl[icent]->Draw("psame");
        //draw systematic error
        const float sysw = 0.15;
        for(int i=0; i<gD0sys_xl[icent]->GetN(); i++) {
            const float sysl = gD0sys_xl[icent]->GetY()[i] * 0.05;
            TLine *llw = new TLine(gD0sys_xl[icent]->GetX()[i]-sysw,gD0sys_xl[icent]->GetY()[i]-gD0sys_xl[icent]->GetEY()[i],gD0sys_xl[icent]->GetX()[i]+sysw,gD0sys_xl[icent]->GetY()[i]-gD0sys_xl[icent]->GetEY()[i]);
            llw->SetLineWidth(2);
            llw->SetLineColor(COLOR[icent]);
            llw->Draw("same");
            TLine *lhi = new TLine(gD0sys_xl[icent]->GetX()[i]-sysw,gD0sys_xl[icent]->GetY()[i]+gD0sys_xl[icent]->GetEY()[i],gD0sys_xl[icent]->GetX()[i]+sysw,gD0sys_xl[icent]->GetY()[i]+gD0sys_xl[icent]->GetEY()[i]);
            lhi->SetLineWidth(2);
            lhi->SetLineColor(COLOR[icent]);
            lhi->Draw("same");
            TLine *lv = new TLine(gD0sys_xl[icent]->GetX()[i]-sysw,gD0sys_xl[icent]->GetY()[i]-gD0sys_xl[icent]->GetEY()[i],gD0sys_xl[icent]->GetX()[i]-sysw,gD0sys_xl[icent]->GetY()[i]-gD0sys_xl[icent]->GetEY()[i]+sysl);
            lv->SetLineWidth(2);
            lv->SetLineColor(COLOR[icent]);
            lv->Draw("same");
            TLine *lv = new TLine(gD0sys_xl[icent]->GetX()[i]+sysw,gD0sys_xl[icent]->GetY()[i]-gD0sys_xl[icent]->GetEY()[i],gD0sys_xl[icent]->GetX()[i]+sysw,gD0sys_xl[icent]->GetY()[i]-gD0sys_xl[icent]->GetEY()[i]+sysl);
            lv->SetLineWidth(2);
            lv->SetLineColor(COLOR[icent]);
            lv->Draw("same");
            TLine *lv = new TLine(gD0sys_xl[icent]->GetX()[i]-sysw,gD0sys_xl[icent]->GetY()[i]+gD0sys_xl[icent]->GetEY()[i],gD0sys_xl[icent]->GetX()[i]-sysw,gD0sys_xl[icent]->GetY()[i]+gD0sys_xl[icent]->GetEY()[i]-sysl);
            lv->SetLineWidth(2);
            lv->SetLineColor(COLOR[icent]);
            lv->Draw("same");
            TLine *lv = new TLine(gD0sys_xl[icent]->GetX()[i]+sysw,gD0sys_xl[icent]->GetY()[i]+gD0sys_xl[icent]->GetEY()[i],gD0sys_xl[icent]->GetX()[i]+sysw,gD0sys_xl[icent]->GetY()[i]+gD0sys_xl[icent]->GetEY()[i]-sysl);
            lv->SetLineWidth(2);
            lv->SetLineColor(COLOR[icent]);
            lv->Draw("same");
        }
    }
    legend1->Draw();
    gPad->SetLogy();
    sprintf(name,"%s/D0_spectra.pdf",dir);
    c1->SaveAs(name);
}

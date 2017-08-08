#include "../myFunction.h"
#include "../myConst.h"
#include "../anaCuts.h"
void plot_RAA_1() {
    globalSetting();
    char dir[250];
    char name[250];
    char title[250];
    //TString CMD = 0;
    char CMD[250];
    TLegend* legend;
    TH1F* h0;
    
    sprintf(dir,"pic");
    sprintf(CMD,"[ -d %s ] || mkdir -p %s",dir,dir);
    gSystem->Exec(CMD);
    
    //const int ncent = 3;
    //const char nameCent[ncent][250] = {"0-10%", "10-20%", "20-40%", "40-60%", "60-80%"};
    //const char nameCentXL[ncent][250] = {"0_10", "10_20", "20_40", "40_60", "60_80"};
    //const char nameCent[ncent][250] = {"0-20%", "20-40%", "40-80%"};
    //const char nameCentXL[ncent][250] = {"0_20", "20_40", "40_80"};
    
    //Read spectra
    //1. from xiaolong
    TGraphErrors* gD0err_xl[ncent];
    TGraphErrors* gD0sys_xl[ncent];
    TFile* fin1 = new TFile("D0RAA_Run14HFT.root");
    for(int icent=0; icent<ncent; icent++) {
        gD0err_xl[icent] = (TGraphErrors*)fin1->Get(Form("D0_RAA_err_%s",nameCent1[icent]));
        gD0sys_xl[icent] = (TGraphErrors*)fin1->Get(Form("D0_RAA_sys_%s",nameCent1[icent]));
    }
    fin1->Close();
    
    //Set for Draw
    float markerSize = 2.0;
    float lineWidth = 2;
    for(int icent=0; icent<ncent; icent++) {
        gD0err_xl[icent]->SetMarkerStyle(kFullCircle);
        gD0err_xl[icent]->SetMarkerSize(markerSize);
        gD0err_xl[icent]->SetMarkerColor(COLOR[icent]);
        gD0err_xl[icent]->SetLineWidth(lineWidth);
        gD0err_xl[icent]->SetLineColor(COLOR[icent]);
    }
    
    //plot
    TCanvas* c1 = new TCanvas("c1", "A Canvas",10,10,800,800);
    setPad(c1);
    float ymin = 0;
    float ymax = 1;
    h0= new TH1F("","",1,0,10);
    h0->SetMinimum(ymin),
    h0->SetMaximum(ymax);
    setHisto(h0,"","p_{T} (GeV/c)", "R_{AA}");
    for(int icent=0; icent<ncent; icent++) {
        const float legy_lw = 0.93-1*0.05;
        legend = new TLegend(0.7,legy_lw,0.92,0.93);
        legend->SetFillStyle(0);
        legend->SetFillColor(10);
        legend->SetBorderSize(0);
        legend->SetTextSize(0.04);
        legend->SetTextFont(132);
        //legend->SetHeader(nameCent[icent]);
        h0->Draw();
        legend->AddEntry(gD0err_xl[icent],nameCent[icent],"p");
        gD0err_xl[icent]->Draw("psame");
        //draw systematic error
        const float sysw = 0.15;
        for(int i=0; i<gD0sys_xl[icent]->GetN(); i++) {
            const float sysl = gD0sys_xl[icent]->GetY()[i] * 0.03;
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
        legend->Draw();
        sprintf(name,"%s/D0RAA_%s.pdf",dir,nameCent1[icent]);
        c1->SaveAs(name);
    }
}

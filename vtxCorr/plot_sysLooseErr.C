#include "../anaCuts.h"
#include "../myFunction.h"
#include "../myConst.h"
void plot_sysLooseErr() {
    globalSetting();
    char dir[250];
    char name[250];
    char CMD[250];
    char lname[16][250];
    //TPaveStats* ptstates;
    TLegend* legend;
    TH1F* h0;
    
    sprintf(dir,"pic/sysLoose");
    sprintf(CMD,"[ -d %s ] || mkdir -p %s",dir,dir);
    gSystem->Exec(CMD);
    
    const int nType = 2;
    const char nameDir[nType][100] = {"./out", "./out"};
    const char nameType[nType][100] = {"loose", "sys_loose"};
    const char nameType1[nType][100] = {"gRef cut 1", "gRef cut 2"};
    
    // read and cal.
    TH1F* hRatio[nType][ncent_vtx];
    for(int itype; itype<nType; itype++) {
        TFile* fin = new TFile(Form("%s/vtxCorr_%s.root",nameDir[itype],nameType[itype]));
        for(int icent=0; icent<ncent_vtx; icent++) {
            TH1F* htmpHjCuts = (TH1F*)fin->Get(Form("hptCutsHj_cent%s_%s",nameCent_vtx1[icent],nameType[itype]));
            TH1F* htmpHjNoCuts = (TH1F*)fin->Get(Form("hptNoCutsHj_cent%s_%s",nameCent_vtx1[icent],nameType[itype]));
            TH1F* htmpFsCuts = (TH1F*)fin->Get(Form("hptCutsFs_cent%s_%s",nameCent_vtx1[icent],nameType[itype]));
            TH1F* htmpFsNoCuts = (TH1F*)fin->Get(Form("hptNoCutsFs_cent%s_%s",nameCent_vtx1[icent],nameType[itype]));
            
            // rebin
            htmpHjCuts = (TH1F*)htmpHjCuts->Rebin(npt,Form("htmpHjCuts_%i_%i",itype,icent),nptbin);
            htmpHjNoCuts = (TH1F*)htmpHjNoCuts->Rebin(npt,Form("htmpHjNoCuts_%i_%i",itype,icent),nptbin);
            htmpFsCuts = (TH1F*)htmpFsCuts->Rebin(npt,Form("htmpFsCuts_%i_%i",itype,icent),nptbin);
            htmpFsNoCuts = (TH1F*)htmpFsNoCuts->Rebin(npt,Form("htmpFsNoCuts_%i_%i",itype,icent),nptbin);
            
            htmpHjCuts->Divide(htmpHjNoCuts);
            htmpFsCuts->Divide(htmpFsNoCuts);
            
            htmpHjCuts->Divide(htmpFsCuts);
            hRatio[itype][icent] = (TH1F*)htmpHjCuts->Clone(Form("hRatio_%i_%i",itype,icent));
            hRatio[itype][icent]->SetDirectory(0);
        }
        fin->Close();
    }
    
    // cal. error
    TH1F* hmean[ncent_vtx];
    TH1F* herr[ncent_vtx];
    for(int icent=0; icent<ncent_vtx; icent++) {
        hmean[icent] = (TH1F*)hRatio[0][icent]->Clone(Form("mean_%s",nameCent_vtx1[icent]));
        herr[icent] = (TH1F*)hRatio[0][icent]->Clone(Form("sysError_%s",nameCent_vtx1[icent]));
        for(int ibin=1; ibin<=hRatio[0][icent]->GetNbinsX(); ibin++) {
            float num1 = hRatio[0][icent]->GetBinContent(ibin);
            float num2 = hRatio[1][icent]->GetBinContent(ibin);
            float error;
            if((num1+num2)==0) error = 0;
            else error = fabs(num1 - num2)/(num1 + num2);
            herr[icent]->SetBinContent(ibin,error);
            hmean[icent]->SetBinContent(ibin,(num1+num2)/2.);
        }
    }
    
    // plot
    TCanvas* c1 = new TCanvas("c1", "A Canvas",10,10,800,800);
    setPad(c1);
    const float linewidth = 2;
    const float markersize =  1.5;
    for(int icent=0; icent<ncent_vtx; icent++) {
        for(int itype=0; itype<nType; itype++) {
            hRatio[itype][icent]->SetMarkerStyle(MARKERSTYLE[itype]);
            hRatio[itype][icent]->SetMarkerColor(COLOR[itype]);
            hRatio[itype][icent]->SetMarkerSize(markersize);
            hRatio[itype][icent]->SetLineColor(COLOR[itype]);
            hRatio[itype][icent]->SetLineWidth(linewidth);
        }
        herr[icent]->SetLineWidth(linewidth);
        herr[icent]->SetLineColor(COLOR[nType]);
        hmean[icent]->SetLineWidth(linewidth);
        hmean[icent]->SetLineColor(COLOR[nType+1]);
        
        for(int itype=0; itype<nType; itype++) {
            if(itype==0) {
                hRatio[itype][icent]->Draw("e");
                setHisto(hRatio[itype][icent],"","p_{T} (GeV/c)", "Eff. Ratio (pure HIJING / fast Sim)");
                hRatio[itype][icent]->GetYaxis()->SetRangeUser(0,1.8);
                //hRatio[itype][icent]->GetYaxis()->SetRangeUser(1e-3,9);
            }
            else hRatio[itype][icent]->Draw("esame");
        }
        
        hmean[icent]->Draw("histsame");
        herr[icent]->Draw("histsame");
        
        if(!legend) {
            const float legy_up = 0.88;
            const float legy_lw = legy_up-(nType+2)*0.05;
            legend = new TLegend(0.2,legy_lw,0.6,legy_up);
            legend->SetFillStyle(0);
            legend->SetFillColor(10);
            legend->SetBorderSize(0);
            legend->SetTextSize(0.05);
            legend->SetTextFont(132);
            for(int itype=0; itype<nType; itype++) legend->AddEntry(hRatio[itype][icent],nameType1[itype],"p");
            legend->AddEntry(hmean[icent],"mean","l");
            legend->AddEntry(herr[icent],"error","l");
        }
        legend->Draw();
        
        sprintf(name,"AuAu @200GeV %s", nameCent_vtx[icent]);
        drawLatex(0.2,0.89,name,132,0.05,1);
        
        //gPad->SetLogy();
        sprintf(name,"%s/SysLoose_%s.gif",dir,nameCent_vtx1[icent]);
        c1->SaveAs(name);
    }
}

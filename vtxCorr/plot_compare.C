#include "../anaCuts.h"
#include "../myFunction.h"
#include "../myConst.h"
void plot_compare() {
    globalSetting();
    char dir[250];
    char name[250];
    char CMD[250];
    char lname[16][250];
    //TPaveStats* ptstates;
    TLegend* legend;
    TH1F* h0;
    
    sprintf(dir,"pic/compare");
    sprintf(CMD,"[ -d %s ] || mkdir -p %s",dir,dir);
    gSystem->Exec(CMD);
    
    const int nType = 4;
    const char nameType[nType][100] = {"default", "tight", "loose", "pt1"};
    const char nameType1[nType][100] = {"default", "tight topology cuts", "loose topology cuts", "p_{T} > 0.3"};
    
    //read eff from pure hijing
    TH2F* h2Base_hj;
    TH2F* h2Eff_hj[nType];
    TFile* fin1 = new TFile("../out/effPureHijing.root");
    h2Base_hj = (TH2F*)fin1->Get("h2McPtCent");
    h2Base_hj->SetDirectory(0);
    for(int i=0; i<nType; i++) {
        h2Eff_hj[i] = (TH2F*)fin1->Get(Form("h2PtCent_%s",nameType[i]));
        h2Eff_hj[i]->SetDirectory(0);
    }
    fin1->Close();
    
    //read eff from fast simu based on hijing
    TH2F* h2Base_fs;
    TH2F* h2Eff_fs[nType];
    TFile* fin2 = new TFile("../out/eff.FastSimu.Hijing.root");
    h2Base_fs = (TH2F*)fin2->Get("h2McPtCent");
    h2Base_fs->SetDirectory(0);
    for(int i=0; i<nType; i++) {
        h2Eff_fs[i] = (TH2F*)fin2->Get(Form("h2PtCent_%s",nameType[i]));
        h2Eff_fs[i]->SetDirectory(0);
    }
    fin2->Close();
    
    //calculate eff
    TH1F* heffHj[ncent_vtx][nType];
    TH1F* heffFs[ncent_vtx][nType];
    for(int icent=0; icent<ncent_vtx; icent++) {
        TH1F* hbase_hj = (TH1F*)h2Base_hj->ProjectionX(Form("hbaseHj_%i",icent),icent+1,icent+1)->Rebin(npt,Form("hbaseHjRebin_%i",icent),nptbin);
        TH1F* hbase_fs = (TH1F*)h2Base_fs->ProjectionX(Form("hbaseFs_%i",icent),icent+1,icent+1)->Rebin(npt,Form("hbaseFsRebin_%i",icent),nptbin);
        for(int itype=0; itype<nType; itype++) {
            heffHj[icent][itype] = (TH1F*)h2Eff_hj[itype]->ProjectionX(Form("heffHj_%i_%i",icent,itype),icent+1,icent+1)->Rebin(npt,Form("heffHjRebin_%i_%i",icent,itype),nptbin);
            heffHj[icent][itype]->Divide(hbase_hj);
            
            heffFs[icent][itype] = (TH1F*)h2Eff_fs[itype]->ProjectionX(Form("heffFs_%i_%i",icent,itype),icent+1,icent+1)->Rebin(npt,Form("heffFsRebin_%i_%i",icent,itype),nptbin);
            heffFs[icent][itype]->Divide(hbase_fs);
        }
    }
    
    TH1F* hRatio[ncent_vtx][nType];
    for(int icent=0; icent<ncent_vtx; icent++) {
        for(int itype=0; itype<nType; itype++) {
            hRatio[icent][itype] = (TH1F*)heffHj[icent][itype]->Clone(Form("hRatio_%i_%i",itype,icent));
            hRatio[icent][itype]->Divide(heffFs[icent][itype]);
        }
    }
    
    // plot
    TCanvas* c1 = new TCanvas("c1", "A Canvas",10,10,800,800);
    setPad(c1);
    const float linewidth = 2;
    const float markersize =  1.5;
    for(int icent=0; icent<ncent_vtx; icent++) {
        for(int itype=0; itype<nType; itype++) {
            hRatio[icent][itype]->SetMarkerStyle(MARKERSTYLE[itype]);
            hRatio[icent][itype]->SetMarkerColor(COLOR[itype]);
            hRatio[icent][itype]->SetMarkerSize(markersize);
            hRatio[icent][itype]->SetLineColor(COLOR[itype]);
            hRatio[icent][itype]->SetLineWidth(linewidth);
        }
        
        for(int itype=0; itype<nType; itype++) {
            if(itype==0) {
                hRatio[icent][itype]->Draw("e");
                setHisto(hRatio[icent][itype],"","p_{T} (GeV/c)", "Eff. Ratio (pure HIJING / fast Sim)");
                hRatio[icent][itype]->GetYaxis()->SetRangeUser(0,1.8);
                //hRatio[icent][itype]->GetYaxis()->SetRangeUser(1e-3,9);
            }
            else hRatio[icent][itype]->Draw("esame");
        }
        
        //herr[icent]->Draw("histsame");
        
        if(!legend) {
            const float legy_up = 0.88;
            const float legy_lw = legy_up-nType*0.05;
            legend = new TLegend(0.2,legy_lw,0.6,legy_up);
            legend->SetFillStyle(0);
            legend->SetFillColor(10);
            legend->SetBorderSize(0);
            legend->SetTextSize(0.05);
            legend->SetTextFont(132);
            for(int itype=0; itype<nType; itype++) legend->AddEntry(hRatio[icent][itype],nameType1[itype],"p");
            //legend->AddEntry(herr[icent],"error","l");
        }
        legend->Draw();
        
        sprintf(name,"AuAu @200GeV %s", nameCent_vtx[icent]);
        drawLatex(0.2,0.89,name,132,0.05,1);
        
        //gPad->SetLogy();
        sprintf(name,"%s/compare_%s.gif",dir,nameCent_vtx1[icent]);
        c1->SaveAs(name);
    }
}

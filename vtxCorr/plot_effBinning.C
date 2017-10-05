#include "../myFunction.h"
#include "../anaCuts.h"
void plot_effBinning()
{
    globalSetting();
    char dir[250];
    char name[250];
    char CMD[250];
    char lname[16][250];
    //TPaveStats* ptstates;
    TLegend* legend;
    TH1F* h0;
    
    sprintf(dir,"pic/eff");
    sprintf(CMD,"[ -d %s ] || mkdir -p %s",dir,dir);
    gSystem->Exec(CMD);
    
    const int neff = 4;
    const char effName[neff][100] = {"h2PtCentAccCut","h2PtCentTpcCut","h2PtCentTpcHftCut","h2PtCentTpcHftTopoCut"};
    const char effLegend[neff][100] = {"Acceptence Only", "TPC", "TPC + HFT", "TPC + HFT + TopoCuts"};
    const char effName1[neff][100] = {"AcceptOnly", "TpcOnly", "TpcHft", "TpcHftTopo"};
    
    //read eff from pure hijing
    TH2F* h2Base_hj;
    TH2F* h2Eff_hj[neff];
    TFile* fin1 = new TFile("./out/effPureHijing.root");
    h2Base_hj = (TH2F*)fin1->Get("h2McPtCent");
    h2Base_hj->SetDirectory(0);
    for(int i=0; i<neff; i++) {
        h2Eff_hj[i] = (TH2F*)fin1->Get(effName[i]);
        h2Eff_hj[i]->SetDirectory(0);
    }
    fin1->Close();
    
    //read eff from fast simu based on hijing
    TH2F* h2Base_fs;
    TH2F* h2Eff_fs[neff];
    TFile* fin2 = new TFile("./out/eff.FastSimu.Hijing.root");
    h2Base_fs = (TH2F*)fin2->Get("h2McPtCent");
    h2Base_fs->SetDirectory(0);
    for(int i=0; i<neff; i++) {
        h2Eff_fs[i] = (TH2F*)fin2->Get(effName[i]);
        h2Eff_fs[i]->SetDirectory(0);
    }
    fin2->Close();
    
    //calculate eff
    TH1F* heffHj[ncent_vtx][neff];
    TH1F* heffFs[ncent_vtx][neff];
    for(int icent=0; icent<ncent_vtx; icent++) {
        TH1F* hbase_hj = (TH1F*)h2Base_hj->ProjectionX(Form("hbaseHj_%i",icent),icent+1,icent+1)->Rebin(npt,Form("hbaseHjRebin_%i",icent),nptbin);
        TH1F* hbase_fs = (TH1F*)h2Base_fs->ProjectionX(Form("hbaseFs_%i",icent),icent+1,icent+1)->Rebin(npt,Form("hbaseFsRebin_%i",icent),nptbin);
        for(int ieff=0; ieff<neff; ieff++) {
            heffHj[icent][ieff] = (TH1F*)h2Eff_hj[ieff]->ProjectionX(Form("heffHj_%i_%i",icent,ieff),icent+1,icent+1)->Rebin(npt,Form("heffHjRebin_%i_%i",icent,ieff),nptbin);
            heffHj[icent][ieff]->Divide(hbase_hj);
            
            heffFs[icent][ieff] = (TH1F*)h2Eff_fs[ieff]->ProjectionX(Form("heffFs_%i_%i",icent,ieff),icent+1,icent+1)->Rebin(npt,Form("heffFsRebin_%i_%i",icent,ieff),nptbin);
            heffFs[icent][ieff]->Divide(hbase_fs);
        }
    }
    
    //calculate double ratio (Pure hijing / Fast Simu)
    TH1F* hratio[ncent_vtx][neff];
    for(int icent=0; icent<ncent_vtx; icent++) {
        for(int ieff=0; ieff<neff; ieff++) {
            hratio[icent][ieff] = (TH1F*)heffHj[icent][ieff]->Clone(Form("hratio%i%i",ieff,icent));
            hratio[icent][ieff]->Divide(heffFs[icent][ieff]);
        }
    }
    
    //fit function
    TF1* fratio = new TF1("fratio","pol0",0,8);
    fratio->SetParName(0,"R");
    fratio->SetLineColor(2);
    fratio->SetLineWidth(2);
    fratio->SetParameter(0,0.5);
    fratio->SetParLimits(0,0,2);
    
    //cout << heffHj[0]->GetNbinsX() << "\t" << heffHj[0]->GetBinWidth(2) << endl;
    //cout << heffFs[0]->GetNbinsX() << "\t" << heffFs[0]->GetBinWidth(2) << endl;
    
    //plot
    TCanvas *c1 = new TCanvas("c1", "c1",10,10,600,800);
    gStyle->SetOptFit(0);
    gStyle->SetOptStat(0);
    gStyle->SetEndErrorSize(0.01);
    c1->SetFillColor(10);
    c1->SetBorderMode(0);
    c1->SetBorderSize(2);
    c1->SetFrameFillColor(0);
    c1->SetFrameBorderMode(0);
    c1->SetLogy();
    //  c1->SetGridx();
    //  c1->SetGridy();
    c1->SetLeftMargin(0.18);
    c1->SetBottomMargin(0.15);
    c1->SetTopMargin(0.025);
    c1->SetRightMargin(0.025);
    
    TPad *pad1 = new TPad("pad1","",0.,0.5,1.,1.);
    pad1->SetFillStyle(4000);
    pad1->SetFillColor(10);
    pad1->SetBorderMode(0);
    pad1->SetBorderSize(0);
    pad1->SetFrameFillColor(0);
    pad1->SetFrameBorderMode(0);
    pad1->SetFrameBorderSize(0);
    pad1->SetFrameLineWidth(1);
    pad1->SetBottomMargin(0.00);
    pad1->SetLeftMargin(0.15);
    pad1->SetRightMargin(0.02);
    pad1->SetTopMargin(0.01);
    pad1->SetTickx();
    pad1->SetTicky();
    
    TPad *pad2 = new TPad("pad2","",0.,0.,1.,0.5);
    pad2->SetFillColor(10);
    pad2->SetBorderMode(0);
    pad2->SetBorderSize(0);
    pad2->SetFrameFillColor(0);
    pad2->SetFrameBorderMode(0);
    pad2->SetFrameBorderSize(0);
    pad2->SetFrameLineWidth(1);
    pad2->SetTopMargin(0.00);
    pad2->SetLeftMargin(0.15);
    pad2->SetRightMargin(0.02);
    pad2->SetBottomMargin(0.15);
    pad2->SetTickx();
    pad2->SetTicky();
    
    pad1->Draw();
    pad2->Draw();
    
    for(int icent=0; icent<ncent_vtx; icent++) {
        for(int ieff=0; ieff<neff; ieff++) {
            const float ymin = 0;
            const float ymax = 1.2*heffFs[icent][ieff]->GetMaximum();
            pad1->cd();
            //pad1->SetLogy();
            heffHj[icent][ieff]->Draw("e");
            heffHj[icent][ieff]->GetXaxis()->SetRangeUser(0,8.);
            heffHj[icent][ieff]->GetYaxis()->SetRangeUser(ymin,ymax);
            //cout << ymax << endl;
            setHisto(heffHj[icent][ieff],"","D^{0} p_{T} (GeV/c)", "efficiency");
            heffHj[icent][ieff]->GetXaxis()->SetLabelOffset(999);
            heffHj[icent][ieff]->GetXaxis()->SetLabelSize(0.0);
            heffHj[icent][ieff]->GetXaxis()->SetLabelFont(22);
            heffHj[icent][ieff]->GetXaxis()->SetTitleFont(22);
            heffHj[icent][ieff]->GetYaxis()->SetTitleOffset(1.2);
            heffHj[icent][ieff]->GetYaxis()->SetTitleSize(0.06);
            //heffHj[icent][ieff]->GetYaxis()->SetLabelOffset(0.12);
            heffHj[icent][ieff]->GetYaxis()->SetLabelSize(0.06);
            heffHj[icent][ieff]->GetYaxis()->SetLabelFont(22);
            heffHj[icent][ieff]->GetYaxis()->SetTitleFont(22);
            heffHj[icent][ieff]->SetMarkerStyle(kFullCircle);
            heffHj[icent][ieff]->SetMarkerSize(1.2);
            heffHj[icent][ieff]->SetMarkerColor(1);
            heffHj[icent][ieff]->SetLineColor(1);
            heffHj[icent][ieff]->SetLineWidth(2);
            heffFs[icent][ieff]->Draw("esame");
            heffFs[icent][ieff]->SetMarkerStyle(kOpenCircle);
            heffFs[icent][ieff]->SetMarkerSize(1.2);
            heffFs[icent][ieff]->SetMarkerColor(2);
            heffFs[icent][ieff]->SetLineColor(2);
            heffFs[icent][ieff]->SetLineWidth(2);
            TLegend* lg1 = new TLegend(0.58,0.83,0.9,0.95);
            lg1->SetFillStyle(0);
            lg1->SetFillColor(10);
            lg1->SetBorderSize(0);
            lg1->SetTextSize(0.05);
            lg1->SetTextFont(132);
            //lg1->SetHeader(Form("%s",nameCent[i]));
            lg1->AddEntry(heffHj[icent][ieff],"Pure hijing","p");
            lg1->AddEntry(heffFs[icent][ieff],"Fast simu.","p");
            lg1->Draw();
            sprintf(name,"AuAu @200GeV  %s",nameCent_vtx[icent]);
            drawLatex(0.2,0.9,name,132,0.05,1);
            sprintf(name,"%s",effLegend[ieff]);
            drawLatex(0.2,0.82,name,132,0.05,1);
            
            pad2->cd();
            hratio[icent][ieff]->Fit(fratio,"INOR","",0,5);
            hratio[icent][ieff]->Fit(fratio,"INOR","",0,5);
            hratio[icent][ieff]->Fit(fratio,"INOR","",0,5);
            hratio[icent][ieff]->SetMarkerStyle(kFullCircle);
            hratio[icent][ieff]->SetMarkerSize(1.2);
            hratio[icent][ieff]->SetMarkerColor(1);
            hratio[icent][ieff]->SetLineColor(1);
            hratio[icent][ieff]->SetLineWidth(2);
            hratio[icent][ieff]->Draw("e");
            fratio->Draw("same");
            hratio[icent][ieff]->GetXaxis()->SetRangeUser(0,8.);
            hratio[icent][ieff]->GetYaxis()->SetRangeUser(0,1.5);
            setHisto(hratio[icent][ieff],"","D^{0} p_{T} (GeV/c)", "Ratio (PureHijing/FastSimu)");
            hratio[icent][ieff]->GetXaxis()->SetLabelSize(0.06);
            hratio[icent][ieff]->GetXaxis()->SetLabelFont(22);
            hratio[icent][ieff]->GetXaxis()->SetTitleFont(22);
            hratio[icent][ieff]->GetYaxis()->SetTitleOffset(1.2);
            hratio[icent][ieff]->GetYaxis()->SetTitleSize(0.06);
            hratio[icent][ieff]->GetXaxis()->SetTitleOffset(1.1);
            hratio[icent][ieff]->GetYaxis()->SetLabelSize(0.06);
            hratio[icent][ieff]->GetYaxis()->SetLabelFont(22);
            hratio[icent][ieff]->GetYaxis()->SetTitleFont(22);
            hratio[icent][ieff]->SetMarkerStyle(20);
            hratio[icent][ieff]->SetMarkerSize(1.2);
            hratio[icent][ieff]->SetMarkerColor(1);
            sprintf(name,"Ratio = %.2f #pm %.2f",fratio->GetParameter(0),fratio->GetParError(0));
            drawLatex(0.2,0.9,name,132,0.05,1);
            sprintf(name,"#chi^{2}/ndf = %.1f/%i",fratio->GetChisquare(),fratio->GetNDF());
            drawLatex(0.2,0.82,name,132,0.05,1);
            float xmin = 0;
            float xmax = 8;
            TLine* l1 = new TLine(xmin,1,xmax,1);
            l1->SetLineStyle(7);
            //l1->Draw("same");
            TLine* l2 = new TLine(xmin,0.8,xmax,0.8);
            l2->SetLineStyle(7);
            //l2->Draw("same");
            TLine* l3 = new TLine(xmin,1.2,xmax,1.2);
            l3->SetLineStyle(7);
            //l3->Draw("same");
            sprintf(name,"%s/%s_%s.gif",dir,effName1[ieff],nameCent_vtx1[icent]);
            c1->SaveAs(name);
        }
    }

    TFile* fOut = new TFile("Mc_d0Eff_peripher_newCuts.root", "RECREATE");
    for(int icent=0; icent<ncent_vtx; icent++) {
        for(int ieff=0; ieff<neff; ieff++) {
            heffHj[icent][ieff]->Write();
            heffFs[icent][ieff]->Write();
            hratio[icent][ieff]->Write();
        }
    }
    fOut->Close();
    
    
    //delete
    for(int icent=0; icent<ncent_vtx; icent++) {
        for(int ieff=0; ieff<neff; ieff++) {
            if(heffHj[icent][ieff]) delete heffHj[icent][ieff];
            if(heffFs[icent][ieff]) delete heffFs[icent][ieff];
            if(hratio[icent][ieff]) delete hratio[icent][ieff];
        }
    }
}

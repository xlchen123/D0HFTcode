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
    
    const int nType = 4;
    const char nameDir[nType][100] = {"./out", "./out", "./out", "./out"};
    const char nameType[nType][100] = {"default", "loose", "tight", "pt0_3"};
    const char nameType1[nType][100] = {"default", "looseTopo", "tightTopo", "ptCut300MeV"};
    
    // read and cal.
    TH1F* heffHj[nType][ncent_vtx];
    TH1F* heffFs[nType][ncent_vtx];
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
            heffHj[itype][icent] = (TH1F*)htmpHjCuts->Clone(Form("heffHj_%i_%i",itype,icent));
            heffHj[itype][icent]->SetDirectory(0);
            heffFs[itype][icent] = (TH1F*)htmpFsCuts->Clone(Form("heffFs_%i_%i",itype,icent));
            heffFs[itype][icent]->SetDirectory(0);
            
            htmpHjCuts->Divide(htmpFsCuts);
            hRatio[itype][icent] = (TH1F*)htmpHjCuts->Clone(Form("hRatio_%i_%i",itype,icent));
            hRatio[itype][icent]->SetDirectory(0);
        }
        fin->Close();
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
    pad1->SetTopMargin(0.025);
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
    
    for(int itype=0; itype<nType; itype++) {
        for(int icent=0; icent<ncent_vtx; icent++) {
            const float ymin = 0;
            const float ymax = 1.4*heffFs[itype][icent]->GetMaximum();
            pad1->cd();
            //pad1->SetLogy();
            heffHj[itype][icent]->Draw("e");
            heffHj[itype][icent]->GetXaxis()->SetRangeUser(0,8.);
            heffHj[itype][icent]->GetYaxis()->SetRangeUser(0,0.1);
            //cout << ymax << endl;
            setHisto(heffHj[itype][icent],"","D^{0} p_{T} (GeV/c)", "efficiency");
            heffHj[itype][icent]->GetXaxis()->SetLabelOffset(999);
            heffHj[itype][icent]->GetXaxis()->SetLabelSize(0.0);
            heffHj[itype][icent]->GetXaxis()->SetLabelFont(22);
            heffHj[itype][icent]->GetXaxis()->SetTitleFont(22);
            heffHj[itype][icent]->GetYaxis()->SetTitleOffset(1.2);
            heffHj[itype][icent]->GetYaxis()->SetTitleSize(0.06);
            //heffHj[itype][icent]->GetYaxis()->SetLabelOffset(0.12);
            heffHj[itype][icent]->GetYaxis()->SetLabelSize(0.06);
            heffHj[itype][icent]->GetYaxis()->SetLabelFont(22);
            heffHj[itype][icent]->GetYaxis()->SetTitleFont(22);
            heffHj[itype][icent]->SetMarkerStyle(kFullCircle);
            heffHj[itype][icent]->SetMarkerSize(1.2);
            heffHj[itype][icent]->SetMarkerColor(1);
            heffHj[itype][icent]->SetLineColor(1);
            heffHj[itype][icent]->SetLineWidth(2);
            heffFs[itype][icent]->Draw("esame");
            heffFs[itype][icent]->SetMarkerStyle(kOpenCircle);
            heffFs[itype][icent]->SetMarkerSize(1.2);
            heffFs[itype][icent]->SetMarkerColor(2);
            heffFs[itype][icent]->SetLineColor(2);
            heffFs[itype][icent]->SetLineWidth(2);
            TLegend* lg1 = new TLegend(0.65,0.75,0.9,0.95);
            lg1->SetFillStyle(0);
            lg1->SetFillColor(10);
            lg1->SetBorderSize(0);
            lg1->SetTextSize(0.07);
            lg1->SetTextFont(132);
            //lg1->SetHeader(Form("%s",nameCent[i]));
            lg1->AddEntry(heffHj[itype][icent],"Pure hijing","p");
            lg1->AddEntry(heffFs[itype][icent],"Fast simu.","p");
            lg1->Draw();
            sprintf(name,"AuAu @200GeV  %s",nameCent_vtx1[icent]);
            drawLatex(0.18,0.88,name,132,0.07,1);
            //sprintf(name,"%s",effLegend[ieff]);
            //drawLatex(0.2,0.82,name,132,0.05,1);
            
            pad2->cd();
            hRatio[itype][icent]->Fit(fratio,"INOR","",0,5);
            hRatio[itype][icent]->Fit(fratio,"INOR","",0,5);
            hRatio[itype][icent]->Fit(fratio,"INOR","",0,5);
            hRatio[itype][icent]->SetMarkerStyle(kFullCircle);
            hRatio[itype][icent]->SetMarkerSize(1.2);
            hRatio[itype][icent]->SetMarkerColor(1);
            hRatio[itype][icent]->SetLineColor(1);
            hRatio[itype][icent]->SetLineWidth(2);
            hRatio[itype][icent]->Draw("e");
            fratio->Draw("same");
            hRatio[itype][icent]->GetXaxis()->SetRangeUser(0,8.);
            hRatio[itype][icent]->GetYaxis()->SetRangeUser(0,1.6);
            setHisto(hRatio[itype][icent],"","D^{0} p_{T} (GeV/c)", "Ratio (PureHijing/FastSimu)");
            hRatio[itype][icent]->GetXaxis()->SetLabelSize(0.06);
            hRatio[itype][icent]->GetXaxis()->SetLabelFont(22);
            hRatio[itype][icent]->GetXaxis()->SetTitleFont(22);
            hRatio[itype][icent]->GetYaxis()->SetTitleOffset(1.2);
            hRatio[itype][icent]->GetYaxis()->SetTitleSize(0.06);
            hRatio[itype][icent]->GetXaxis()->SetTitleOffset(1.1);
            hRatio[itype][icent]->GetYaxis()->SetLabelSize(0.06);
            hRatio[itype][icent]->GetYaxis()->SetLabelFont(22);
            hRatio[itype][icent]->GetYaxis()->SetTitleFont(22);
            hRatio[itype][icent]->SetMarkerStyle(20);
            hRatio[itype][icent]->SetMarkerSize(1.2);
            hRatio[itype][icent]->SetMarkerColor(1);
            sprintf(name,"Ratio = %.2f #pm %.2f",fratio->GetParameter(0),fratio->GetParError(0));
            drawLatex(0.18,0.88,name,132,0.07,1);
            sprintf(name,"#chi^{2}/ndf = %.1f/%i",fratio->GetChisquare(),fratio->GetNDF());
            drawLatex(0.18,0.78,name,132,0.07,1);
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
            sprintf(name,"%s/eff_%s_cent%s.gif",dir,nameType1[itype],nameCent_vtx1[icent]);
            c1->SaveAs(name);
        }
    }
    
    //delete
    for(int itype=0; itype<nType; itype++) {
        for(int icent=0; icent<ncent_vtx; icent++) {
            if(heffHj[itype][icent]) delete heffHj[itype][icent];
            if(heffFs[itype][icent]) delete heffFs[itype][icent];
            if(hRatio[itype][icent]) delete hRatio[itype][icent];
        }
    }
}

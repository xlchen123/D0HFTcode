#include "../myFunction.h"
#include "../anaCuts.h"

// std::pair<float, float> gMassExclusionRange(1.80, 1.92);
std::pair<float, float> gMassExclusionRange(1.82, 1.90);
//std::pair<float, float> gMassYieldCountRange(1.82, 1.91);
//std::pair<float, float> gMassYieldSideBand1(1.71, 1.80);
//std::pair<float, float> gMassYieldSideBand2(1.93, 2.02);
Double_t funResidualBg(Double_t *x, Double_t *par)
{
    if (x[0] > gMassExclusionRange.first && x[0] < gMassExclusionRange.second)
    {
        TF1::RejectPoint();
        return 0;
    }
    
    return par[0] + par[1] * x[0];
}

void plot_rawy()
{
    globalSetting();
    char buf[250];
    char dir[250];
    char name[250];
    char title[250];
    //TString CMD = 0;
    char CMD[250];
    TLegend* legend;
    char lname[16][250];
    TPaveStats* ptstates;
    TH1F* htemp;
    TH1F* h0;
    TH1F* h1;
    const double PI = TMath::Pi();
    float pt_lw,    pt_up,    m_lw,    m_up,    cent_lw,    cent_up ;
    int   ptbin_lw, ptbin_up, mbin_lw, mbin_up, centbin_lw, centbin_up;
    
    //fit and plot
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
    pad1->SetLeftMargin(0.13);
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
    pad2->SetLeftMargin(0.13);
    pad2->SetRightMargin(0.02);
    pad2->SetBottomMargin(0.15);
    pad2->SetTickx();
    pad2->SetTicky();
    
    pad1->Draw();
    pad2->Draw();
    
    TF1 *fD0 = new TF1("fD0","gausn(0)*[3]+pol1(4)",1.6,2.1);
    TF1 *gausfun = new TF1("gausfun","gausn(0)*[3]",1.6,2.1);
    TF1 *resifun = new TF1("resifun","pol1(0)",1.6,2.1);
    TF1 *resifunFit = new TF1("resifunFit",funResidualBg,1.6,2.1,2);
    resifunFit->SetNpx(501);
    resifunFit->SetLineColor(kRed);
    //resifunFit->SetLineStyle(7);
    resifunFit->SetLineWidth(2);
    fD0->SetParNames("N","#mu","#sigma","BinWidth","A0","A1");
    // float fitRange_lw = 1.72;
    float fitRange_lw = 1.76;
    // float fitRange_up = 2.1;
    float fitRange_up = 2.0;
    
    for(int icent=0; icent<ncent; icent++) {
// if(icent != 0) continue;
        sprintf(dir,Form("pic/%s",nameCent1[icent]));
        sprintf(CMD,"[ -d %s ] || mkdir -p %s", dir,dir);
        gSystem->Exec(CMD);
        sprintf(CMD,"[ -d data ] || mkdir -p data");
        gSystem->Exec(CMD);
        
        float nptmean[npt], npterr[npt], mumean[npt], muerr[npt], sigmamean[npt], sigmaerr[npt];
        float rawY[npt], rawYerr[npt];
        for(int ipt=0; ipt<npt; ipt++) {
            nptmean[ipt] = 0.5*(nptbin[ipt]+nptbin[ipt+1]);
            npterr[ipt] = 0.5*(-nptbin[ipt]+nptbin[ipt+1]);
        }
        
        TH1F* hInvMass[npt];
        TH1F* hInvMassUL[npt];
        TH1F* hInvMassLS[npt];
        TH1F* hInvMassMixUL[npt];
        TFile* fin = new TFile("data/D0_HFT_signal.root");
        fin->cd();
        for(int ipt=0; ipt<npt; ipt++) {
            sprintf(name,"%s/hInvMassUL_%i_%i",nameCent1[icent],icent,ipt);
            hInvMassUL[ipt] = (TH1F*)fin->Get(name);
            hInvMassUL[ipt]->SetDirectory(0);
            sprintf(name,"%s/hInvMassLS_%i_%i",nameCent1[icent],icent,ipt);
            hInvMassLS[ipt] = (TH1F*)fin->Get(name);
            hInvMassLS[ipt]->SetDirectory(0);
            sprintf(name,"%s/hInvMassMixUL_%i_%i",nameCent1[icent],icent,ipt);
            hInvMassMixUL[ipt] = (TH1F*)fin->Get(name);
            hInvMassMixUL[ipt]->SetDirectory(0);
            sprintf(name,"%s/hInvMassMix_%i_%i",nameCent1[icent],icent,ipt);
            hInvMass[ipt] = (TH1F*)fin->Get(name);
            hInvMass[ipt]->SetDirectory(0);
        }
        fin->Close();
        
        Float_t mean, sigma, lowx, higx;
        double S, B;
        for(int ipt=0; ipt<npt; ipt++) {
            pad1->cd();
            pad1->SetLogy(0);
            hInvMassUL[ipt]->Draw("E");
            setHisto(hInvMassUL[ipt],"","M_{k#pi} (GeV/c^{2})", "Counts (per 10 MeV/c^{2})");
            hInvMassUL[ipt]->SetMaximum(1.5*hInvMassUL[ipt]->GetMaximum());
            hInvMassUL[ipt]->GetXaxis()->SetTitleSize(0.0);
            hInvMassUL[ipt]->GetXaxis()->SetLabelOffset(999);
            hInvMassUL[ipt]->GetXaxis()->SetLabelSize(0.0);
            hInvMassUL[ipt]->GetXaxis()->SetLabelFont(22);
            hInvMassUL[ipt]->GetXaxis()->SetTitleFont(22);
            hInvMassUL[ipt]->GetYaxis()->SetTitleOffset(1.0);
            hInvMassUL[ipt]->GetYaxis()->SetTitleSize(0.06);
            //hInvMassUL[ipt]->GetYaxis()->SetLabelOffset(0.12);
            hInvMassUL[ipt]->GetYaxis()->SetLabelSize(0.06);
            hInvMassUL[ipt]->GetYaxis()->SetLabelFont(22);
            hInvMassUL[ipt]->GetYaxis()->SetTitleFont(22);
            hInvMassUL[ipt]->SetMarkerStyle(20);
            hInvMassUL[ipt]->SetMarkerSize(1.2);
            hInvMassUL[ipt]->SetMarkerColor(1);
            hInvMassMixUL[ipt]->Draw("HISTSAME");
            hInvMassMixUL[ipt]->SetLineColor(2);
            hInvMassMixUL[ipt]->SetLineWidth(2);
            hInvMassLS[ipt]->Draw("HISTSAME");
            hInvMassLS[ipt]->SetLineColor(4);
            hInvMassLS[ipt]->SetLineWidth(2);
            //drawLatex(0.15,0.9,"UnLike Sign",22,0.06,1);
            //drawLatex(0.15,0.84,"Like Sign",22,0.06,2);
            drawLatex(0.55,0.9,"STAR Au+Au @ 200 GeV",22,0.05,1);
            drawLatex(0.55,0.81,Form("%s",nameCent[icent]),22,0.05,4);
            drawLatex(0.55,0.73,Form("%.1f < p_{T} < %.1f GeV/c",nptbin[ipt],nptbin[ipt+1]),22,0.05,4);
            //drawLatex(0.55,0.64,Form("%.3f < DCA_{D^{0}} < %.3f",ndcabin[idca],ndcabin[idca+1]),22,0.06,4);
            //drawLatex(0.55,0.5,"STAR Preliminary",92,0.06,1);
            TLegend* legend2 = new TLegend(0.15,.78,.5,.95);
            legend2->SetTextFont(22);
            legend2->SetTextSize(0.05);
            legend2->SetBorderSize(0);
            legend2->SetFillStyle(0);
            legend2->AddEntry(hInvMassUL[ipt]," US (Same-Event)","lep");
            legend2->AddEntry(hInvMassLS[ipt]," LS (Same-Event)","l");
            legend2->AddEntry(hInvMassMixUL[ipt]," US (Mix-Event)","l");
            //legend2->AddEntry(hInvMass[ipt]," Unlike - Like","lep");
            legend2->Draw();
            
            pad2->cd();
            setHisto(hInvMass[ipt],"","M_{k#pi} (GeV/c^{2})", "Counts (per 10 MeV/c^{2})");
            hInvMass[ipt]->GetXaxis()->SetTitleSize(0.06);
            //hInvMass[ipt]->GetXaxis()->SetLabelOffset(1.0);
            hInvMass[ipt]->GetXaxis()->SetLabelSize(0.06);
            hInvMass[ipt]->GetXaxis()->SetLabelFont(22);
            hInvMass[ipt]->GetXaxis()->SetTitleFont(22);
            hInvMass[ipt]->GetYaxis()->SetTitleOffset(1.0);
            hInvMass[ipt]->GetYaxis()->SetTitleSize(0.06);
            //hInvMass[ipt]->GetYaxis()->SetLabelOffset(0.12);
            hInvMass[ipt]->GetYaxis()->SetLabelSize(0.06);
            hInvMass[ipt]->GetYaxis()->SetLabelFont(22);
            hInvMass[ipt]->GetYaxis()->SetTitleFont(22);
            hInvMass[ipt]->SetMarkerStyle(20);
            hInvMass[ipt]->SetMarkerSize(1.2);
            hInvMass[ipt]->SetMarkerColor(1);
            
            // Fit residual bg and signal (gauss + pol1)
            float N = hInvMass[ipt]->GetBinContent(hInvMass[ipt]->FindBin(1.86));
            resifunFit->SetParameters(0.5*N, -0.5*N);
            hInvMass[ipt]->Fit(resifunFit,"INOR","",fitRange_lw, fitRange_up);
            hInvMass[ipt]->Fit(resifunFit,"INOR","",fitRange_lw, fitRange_up);
            // hInvMass[ipt]->Fit(resifunFit,"INOR","",fitRange_lw, fitRange_up);
            //cout << resifunFit->GetParameter(0) << "\t" << resifunFit->GetParameter(1) << endl;
            // fD0->SetParameters(N, 1.866, 0.014, hInvMass[ipt]->GetBinWidth(4),resifunFit->GetParameter(0),resifunFit->GetParameter(1));
            fD0->SetParameters(N, 1.866, 0.014, hInvMass[ipt]->GetBinWidth(4),1,1);
            //fD0->SetParLimits(0,0,1.e9);
            // fD0->SetParLimits(1,1.85,1.88);
            // fD0->SetParLimits(2,0.003,0.03);
            fD0->FixParameter(3,hInvMass[ipt]->GetBinWidth(4));
            // fD0->FixParameter(4,resifunFit->GetParameter(0));
            // fD0->FixParameter(5,resifunFit->GetParameter(1));
            hInvMass[ipt]->Fit(fD0,"INOR","",fitRange_lw, fitRange_up);
            hInvMass[ipt]->Fit(fD0,"INOR","",fitRange_lw, fitRange_up);
            //hInvMass[ipt]->Fit(fD0,"INOR","",fitRange_lw, fitRange_up);
            
            // caculate significance
            double* par = fD0->GetParameters();
            mean = par[1];
            sigma=par[2];
            const float nsigma = 3*sigma;
            const float shift = 0.01;
            lowx = mean - nsigma;
            higx = mean + nsigma;
            float lowx1 = lowx - 2*nsigma - shift;
            float higx1 = lowx - shift;
            float lowx2 = higx + shift;
            float higx2 = higx + 2*nsigma + shift;
            float signalAndBg, signalAndBg_err, sideBand1, sideBand1_err, sideBand2, sideBand2_err;
            hIntegralAndError(hInvMass[ipt],lowx,higx,signalAndBg,signalAndBg_err);
            hIntegralAndError(hInvMass[ipt],lowx1,higx1,sideBand1,sideBand1_err);
            hIntegralAndError(hInvMass[ipt],lowx2,higx2,sideBand2,sideBand2_err);
            cout << hInvMass[ipt]->FindBin(lowx) << "\t" << hInvMass[ipt]->FindBin(higx) << endl;
            float bg = 0.5*(sideBand1+sideBand2);
            float bg_err = 0.5*sqrt(pow(sideBand1_err,2)+pow(sideBand2_err,2));
            //float yield = signalAndBg - bg;
            //float yielderr = sqrt(pow(signalAndBg_err,2)+pow(bg_err,2));
            float yield = fD0->GetParameter(0);
            float yielderr = fD0->GetParError(0);
            
            rawY[ipt] = yield; rawYerr[ipt] = yielderr;
            mumean[ipt] = fD0->GetParameter(1); muerr[ipt] = fD0->GetParError(1);
            sigmamean[ipt] = 1000.*fD0->GetParameter(2); sigmaerr[ipt] = 1000.*fD0->GetParError(2);
            fD0->SetLineColor(2);
            
            //plot
            float y_up = 1.8*hInvMass[ipt]->GetMaximum();
            float y_lw = 0. - 0.2*hInvMass[ipt]->GetMaximum();
            cout << y_lw << "\t" << y_up << endl;
            hInvMass[ipt]->GetYaxis()->SetRangeUser(y_lw,y_up);
            hInvMass[ipt]->Draw("E");
            resifunFit->SetRange(fitRange_lw,fitRange_up);
            resifunFit->DrawCopy("same");
            fD0->SetRange(fitRange_lw,fitRange_up);
            fD0->DrawCopy("same");
            S = yield;
            float mixbg, mixbg_err;
            hIntegralAndError(hInvMassMixUL[ipt],lowx,higx,mixbg,mixbg_err);
            B = bg + mixbg;
            cout << "bin width: " << hInvMass[ipt]->GetBinWidth(2) << endl;
            cout << "S = " << S << endl;
            cout << "B = " << B << endl;
            int N_y = round(yield);
            int N_error = round(yielderr);
            int Nwr = N_y - N_y%10;
            int Nerrwr = N_error - N_error%10;
            drawLatex(0.55,0.9,Form("RawYield = %i #pm %i",Nwr, Nerrwr),22,0.05,4);
            //drawLatex(0.55,0.8,Form("RawYield = 810 #pm 130"),22,0.06,4);
            //drawLatex(0.55,0.78,"STAR Preliminary",92,0.05,1);
            //drawLatex(0.55,0.855,Form("S/#sqrt{S+2B}=%.1f",S/sqrt(S+2*B)),22,0.06,4);
            TLegend* legend = new TLegend(0.15,.8,.5,.95);
            //TLegend* legend = new TLegend(0.55,.83,.8,.96);
            legend->SetTextFont(22);
            legend->SetTextSize(0.05);
            legend->SetBorderSize(0);
            legend->SetFillStyle(0);
            legend->AddEntry(hInvMass[ipt]," US (SE) - US (ME)","lep");
            legend->AddEntry(fD0," Fit","l");
            //legend2->AddEntry(hInvMass[ipt]," Unlike - Like","lep");
            legend->Draw();
            
            sprintf(name,"%s/cent%s_pt_%.1f_%.1f.gif",dir,nameCent1[icent],nptbin[ipt],nptbin[ipt+1]);
            c1->SaveAs(name);
        }
        
        pad1->cd();
        if(h0) delete h0;
        h0= new TH1F("","",100,0,8);
        // h0= new TH1F("","",100,0,10);
        h0->Draw();
        h0->SetMinimum(1.78),
        h0->SetMaximum(1.96);
        setHisto(h0,"","p_{T} (GeV/c)", "D^{0} mass mean (GeV/c^{2})");
        h0->GetXaxis()->SetTitleSize(0.0);
        h0->GetXaxis()->SetLabelOffset(999);
        h0->GetXaxis()->SetLabelSize(0.0);
        h0->GetXaxis()->SetLabelFont(22);
        h0->GetXaxis()->SetTitleFont(22);
        h0->GetYaxis()->SetTitleOffset(1.0);
        h0->GetYaxis()->SetTitleSize(0.06);
        //h0->GetYaxis()->SetLabelOffset(0.12);
        h0->GetYaxis()->SetLabelSize(0.06);
        h0->GetYaxis()->SetLabelFont(22);
        h0->GetYaxis()->SetTitleFont(22);
        TGraphErrors* gr1 = new TGraphErrors(npt,nptmean,mumean,npterr,muerr);
        gr1->SetMarkerStyle(20);
        gr1->SetMarkerColor(1);
        gr1->SetMarkerSize(1.8);
        gr1->SetLineWidth(3);
        gr1->Draw("psame");
        sprintf(name,"%s",nameCent[icent]);
        drawLatex(0.8,0.85,name,22,0.07,1);
        drawLatex(0.15,0.85,"Mean",22,0.07,1);
        //sprintf(name,"%s/D0_mean_%s.eps",dir,nameCent1[icent]);
        //c1->SaveAs(name);
        
        pad2->cd();
        if(h1) delete h1;
        h1= new TH1F("","",100,0,8);
        // h1= new TH1F("","",100,0,10);
        h1->Draw();
        h1->SetMinimum(0),
        h1->SetMaximum(32);
        setHisto(h1,"","p_{T} (GeV/c)", "D^{0} mass sigma (MeV/c^{2})");
        //h1->SetMaximum(1.8*hInvMass[i]->GetMaximum());
        h1->GetXaxis()->SetTitleSize(0.06);
        //h1->GetXaxis()->SetLabelOffset(1.0);
        h1->GetXaxis()->SetLabelSize(0.06);
        h1->GetXaxis()->SetLabelFont(22);
        h1->GetXaxis()->SetTitleFont(22);
        h1->GetYaxis()->SetTitleOffset(1.0);
        h1->GetYaxis()->SetTitleSize(0.06);
        //h1->GetYaxis()->SetLabelOffset(0.12);
        h1->GetYaxis()->SetLabelSize(0.06);
        h1->GetYaxis()->SetLabelFont(22);
        h1->GetYaxis()->SetTitleFont(22);
        TGraphErrors* gr2 = new TGraphErrors(npt,nptmean,sigmamean,npterr,sigmaerr);
        gr2->SetMarkerStyle(20);
        gr2->SetMarkerColor(1);
        gr2->SetMarkerSize(1.8);
        gr2->SetLineWidth(3);
        gr2->Draw("psame");
        sprintf(name,"%s",nameCent[icent]);
        drawLatex(0.8,0.85,name,22,0.07,1);
        drawLatex(0.15,0.85,"Sigma",22,0.07,1);
        sprintf(name,"%s/D0_mean_sigma_%s.gif",dir,nameCent1[icent]);
        c1->SaveAs(name);
        
        ofstream frawy(Form("data/rawY_%s.txt",nameCent1[icent]));
        for(int i=0; i<npt; i++) frawy << rawY[i] << "\t" << rawYerr[i] << endl;
        frawy.close();
    }
}

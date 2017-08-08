#include "../myFunction.h"

void plot_Nsig_kaon() {
    globalSetting();
    TCanvas* c1 = new TCanvas("c1", "A Canvas",10,10,800,800);
    setPad(c1);
    char buf[250];
    char dir[250];
    char name[250];
    char lname[16][250];
    TLegend* legend;
    TH1F* htemp;
    TPaveStats* ptstates;
    const double PI = TMath::Pi();
    const double twoPI = 2*PI;
    float massPi = 0.13957;
    float massK  = 0.493677;
    
    const int np = 22;
    float pbin[np+1];
    pbin[0] = 0.3;
    pbin[1] = 0.32;
    pbin[2] = 0.35;
    //pbin[2] = 0.37;
    for(int i=3; i<20; i++) pbin[i] = 0.4+0.1*(i-3);
    pbin[20] = 2.2;
    pbin[21] = 2.5;
    pbin[22] = 3.0;
    float p_mean[np], p_err[np];
    for(int i=0; i<np; i++) {
        p_mean[i] = 0.5*(pbin[i]+pbin[i+1]);
        p_err[i] = 0.5*(-pbin[i]+pbin[i+1]);
    }
    int p_low, p_up;
    float width;
    float mean1[np], sigma1[np], N1[np], meanErr1[np], sigmaErr1[np], NErr1[np];
    float mean2[np], sigma2[np], N2[np], meanErr2[np], sigmaErr2[np], NErr2[np];
    float eff1[np], efferr1[np], eff2[np], efferr2[np];
    
    TFile* fin = new TFile("data/kaon_noHFT.root");
    fin->cd();
    TH2F* hNsigUL=(TH2F*)fin->Get("hNsigUL");
    hNsigUL->SetDirectory(0);
    TH2F* hNsigLS=(TH2F*)fin->Get("hNsigLS");
    hNsigLS->SetDirectory(0);
    TH2F* hNsig = new TH2F("hNsig","hNsig",300,0,3,80,-4,4);
    hNsig->SetDirectory(0);
    hNsig->Sumw2();
    hNsig->Add(hNsigUL,hNsigLS,1,-1);
    TH2F* hNsigTofUL=(TH2F*)fin->Get("hNsigTofUL");
    hNsigTofUL->SetDirectory(0);
    TH2F* hNsigTofLS=(TH2F*)fin->Get("hNsigTofLS");
    hNsigTofLS->SetDirectory(0);
    TH2F* hNsigTof = new TH2F("hNsigTof","hNsigTof",300,0,3,80,-4,4);
    hNsigTof->SetDirectory(0);
    hNsigTof->Sumw2();
    hNsigTof->Add(hNsigTofUL,hNsigTofLS,1,-1);
    fin->Close();
    
    //TF1* f_lw = new TF1("f1_lw",fpi_lw,0.2,20,0); //1
    //TF1* f_up = new TF1("f1_up",fpi_up,0.2,20,2); //1
    //f_up->SetParameters(5.429,-2.143);
    TF1* f_lw = new TF1("f2_lw",fk_lw,0.2,20,3); //2
    f_lw->SetParameters(-7.536934,5.828048,-1.307178);
    TF1* f_up = new TF1("f2_up",fk_up,0.2,20,3);  //2
    f_up->SetParameters(8.690754,-6.023949,1.316979);
    
    pad = new TPad("pad","",0.00,0.00,1.00,1.00);
    pad->Draw();
    TPDF *pdf = new TPDF("PDF/KaonPid.pdf");
    pad->cd();
    setPad(gPad,0.1,0.1,0.05,0.12, 0, 0);
    TBox *bLabela = new TBox(0.01, 0.88, 0.99, 0.99);
    bLabela->SetFillColor(kBlue);
    bLabela->Draw();
    TBox *bLabelb = new TBox(0.01, 0.01, 0.99, 0.12);
    bLabelb->SetFillColor(kBlue);
    bLabelb->Draw();
    TLatex tTitle;
    tTitle.SetNDC();
    tTitle.SetTextColor(kBlack);
    tTitle.SetTextSize(0.05);
    tTitle.SetTextFont(22);
    sprintf(name, "Kaon PID and efficiency");
    tTitle.DrawLatex(0.3, 0.55, name);
    
    tTitle.SetTextSize(0.035);
    tTitle.DrawLatex(0.15, 0.15, (new TDatime())->AsSQLString());
    tTitle.SetTextColor(kBlack);
    tTitle.SetTextSize(0.035);
    tTitle.DrawLatex(0.15, 0.2, "By Xiaolong Chen");
    pad->Update();
    
    const int nh = 2;
    TH1F* h1[nh];     TH1F* h2[nh];
    TF1* f1[nh];      TF1* f2[nh];
    int ii;
    int nRebin = 4;
    char fitOp[50] = "INOR";
    //char fitOp[50] = "NOR";
    for(int i=0; i<np; i++) {
        if(i%nh==0) {
            pdf->NewPage();
            pad->cd();
            pad->Clear();
            //c1->cd();
            pad->Divide(2,2,0.001,0.001);
            pad->cd(1);
        }
        p_low = (int)hNsig->GetXaxis()->FindBin(pbin[i]+1e-6);
        p_up = (int)hNsig->GetXaxis()->FindBin(pbin[i+1]-1e-6);
        float n, m;
        
        if(i%nh==0) ii = 0;  else ii = 1;
        if(i!=1&&h1[ii]) delete h1[ii];
        if(i!=1&&f1[ii]) delete f1[ii];
        if(i!=1&&h2[ii]) delete h2[ii];
        if(i!=1&&f2[ii]) delete f2[ii];
        
        f1[ii] = new TF1("f1", oneGaus,  -10,10,4);
        h1[ii] = new TH1F(Form("h1_%i",i),    Form(name,"h1_%i",i),    80,-4,4);
        h1[ii]->Sumw2();
        //cout << Form("h1_%i",i) << endl;
        f2[ii] = new TF1("f2", oneGaus,  -10,10,4);
        h2[ii] = new TH1F(Form(name,"h2_%i",i),    Form(name,"h2_%i",i),    80,-4,4);
        h2[ii]->Sumw2();
        
        if(ii==0) pad->cd(1);  else pad->cd(3);
        setPad(gPad,0.15,0.15,0.15,0.15,0);
        htemp = (TH1F*)hNsig->ProjectionY("_py",p_low,p_up);
        h1[ii]->Add(htemp);
        delete htemp;
        h1[ii]->Rebin(nRebin);
        //h1[ii]->SetMinimum(minimum);
        setHisto(h1[ii], "", "n#sigma_{k}", "Counts");
        width    = h1[ii]->GetBinWidth(2);
        mean1[i]  = h1[ii]->GetMean();
        sigma1[i] = h1[ii]->GetRMS();
        N1[i]     = h1[ii]->GetEntries();
        f1[ii]->SetParameters(N1[i],mean1[i],sigma1[i],width);
        //f1[ii]->SetParameters(N1[i],-0.1,0.99,width);
        f1[ii]->SetParLimits(2,0.5,1.5);
        f1[ii]->SetParNames("N","#mu","#sigma","width");
        f1[ii]->FixParameter(3,width);
        f1[ii]->SetLineColor(1);
        f1[ii]->SetLineWidth(4);
        //f1[ii]->SetLineStyle(7);
        h1[ii]->Fit(f1[ii],fitOp,"",-4,4);
        h1[ii]->Fit(f1[ii],fitOp,"",-4,4);
        h1[ii]->Fit(f1[ii],fitOp,"",-4,4);
        //h1[ii]->Fit(f1[ii],fitOp,"",f1[ii]->GetParameter(1)-2,f1[ii]->GetParameter(1)+2);
        h1[ii]->SetMarkerStyle(20);
        h1[ii]->SetMarkerSize(1.2);
        h1[ii]->SetMarkerColor(1);
        h1[ii]->Draw("e");
        f1[ii]->Draw("same");
        mean1[i] = f1[ii]->GetParameter(1);  meanErr1[i] = f1[ii]->GetParError(1);
        sigma1[i] = f1[ii]->GetParameter(2);  sigmaErr1[i] = f1[ii]->GetParError(2);
        m = f1[ii]->Integral(-2,2)/width;
        n = f1[ii]->GetParameter(0);
        eff1[i] = (m+1)/(n+2);
        efferr1[i] = eff1[i]*(n-m+1)/(n+2)/(n+3);
        sprintf(lname[0],"nSigma_k");
        sprintf(lname[1],"#chi^{2}/ndf = %.0f / %.0f",f1[ii]->GetChisquare(),f1[ii]->GetNDF());
        sprintf(lname[2],"mean = %.3f #pm %.3f",mean1[i],meanErr1[i]);
        sprintf(lname[3],"sigma = %.3f #pm %.3f",sigma1[i],sigmaErr1[i]);
        if(ptstates) delete ptstats;
        ptstats = setPaveStats(h1[ii], (float)0.55, (float)0.7, (float)0.9, (float)0.95, (int)4, lname);
        sprintf(name,"%.2f #leq p < %.2f",pbin[i],pbin[i+1]);
        drawLatex(0.18,0.7,name,22,0.05,4);
        //gPad->Update();
        
        if(ii==0) pad->cd(2);  else pad->cd(4);
        setPad(gPad,0.15,0.15,0.15,0.15,0);
        htemp = (TH1F*)hNsigTof->ProjectionY("_py",p_low,p_up);
        h2[ii]->Add(htemp);
        delete htemp;
        h2[ii]->Rebin(nRebin);
        //h1[ii]->SetMinimum(minimum);
        setHisto(h2[ii], "", "n#sigma_{k}^{TOF}", "Counts");
        width    = h2[ii]->GetBinWidth(2);
        mean2[i]  = h2[ii]->GetMean();
        sigma2[i] = h2[ii]->GetRMS();
        N2[i]     = h2[ii]->GetEntries();
        f2[ii]->SetParameters(N2[i],mean2[i],sigma2[i],width);
        f2[ii]->SetParNames("N","#mu","#sigma","width");
        f2[ii]->FixParameter(3,width);
        f2[ii]->SetLineColor(2);
        f2[ii]->SetLineWidth(4);
        //f2[ii]->SetLineStyle(7);
        if(i==13) {h2[ii]->Fit(f2[ii],fitOp,"",-1.7,4);h2[ii]->Fit(f2[ii],fitOp,"",-1.7,4);h2[ii]->Fit(f2[ii],fitOp,"",-1.7,4);}
        else {h2[ii]->Fit(f2[ii],fitOp,"",-4,4);h2[ii]->Fit(f2[ii],fitOp,"",-4,4);h2[ii]->Fit(f2[ii],fitOp,"",-4,4);}
        //h2[ii]->Fit(f2[ii],fitOp,"",f2[ii]->GetParameter(1)-2,f2[ii]->GetParameter(1)+2);
        h2[ii]->SetMarkerStyle(20);
        h2[ii]->SetMarkerSize(1.2);
        h2[ii]->SetMarkerColor(2);
        h2[ii]->Draw("e");
        f2[ii]->Draw("same");
        mean2[i] = f2[ii]->GetParameter(1);  meanErr2[i] = f2[ii]->GetParError(1);
        sigma2[i] = f2[ii]->GetParameter(2);  sigmaErr2[i] = f2[ii]->GetParError(2);
        //sigma2[i] *= 2.;
        //f2[ii]->SetParameter(2,sigma2[i]);
        m = f2[ii]->Integral(-0.03/0.013,0.03/0.013)/width;
        n = f2[ii]->GetParameter(0);
        eff2[i] = (m+1)/(n+2);
        efferr2[i] = eff2[i]*(n-m+1)/(n+2)/(n+3);
        sprintf(lname[0],"nSigmaTof_k");
        sprintf(lname[1],"#chi^{2}/ndf = %.0f / %.0f",f2[ii]->GetChisquare(),f2[ii]->GetNDF());
        sprintf(lname[2],"mean = %.3f #pm %.3f",mean2[i],meanErr2[i]);
        sprintf(lname[3],"sigma = %.3f #pm %.3f",sigma2[i],sigmaErr2[i]);
        if(ptstates) delete ptstats;
        ptstats = setPaveStats(h2[ii], (float)0.55, (float)0.7, (float)0.9, (float)0.95, (int)4, lname);
        sprintf(name,"%.2f #leq p < %.2f",pbin[i],pbin[i+1]);
        drawLatex(0.18,0.7,name,22,0.05,4);
        gPad->Update();
        
        if(i%nh==1){
            sprintf(name, "%i", i/2+1);
            //pad->SaveAs(name);
            tTitle.DrawLatex(0.02, 0.02, name);
            pad->cd();
            pad->Update();
        }
    }
    if(i%4!=0){
        pad->cd(4);
        sprintf(name, "%i", i/2+1);
        //pad->SaveAs(name);
        tTitle.DrawLatex(0.02, 0.02, name);
        pad->cd();
        pad->Update();  //if i%4 !=0, this is necessary
    }
    
    pdf->Off();
    
    pdf->On();
    pdf->NewPage();
    pad->cd();
    pad->Clear();
    pad->Divide(2,2,0.001,0.001);
    pad->cd(1);
    
    pad->cd(1);
    setPad(gPad,0.15,0.15,0.15,0.15,0);
    h0= new TH1F("","",100,0,3.2);
    h0->SetMinimum(-1),
    h0->SetMaximum(1);
    setHisto(h0,"","p (GeV/c)", "mean");
    //setHisto(h0,"","p (GeV/c)", "ratio",1,0.05,0.05, 1,0.9,1.1, 1,0,3.2,0.,2);
    TGraphErrors* gr1 = new TGraphErrors(np,p_mean,mean1,p_err,meanErr1);
    gr1->SetMarkerStyle(20);
    gr1->SetMarkerColor(1);
    gr1->SetMarkerSize(1.);
    h0->Draw();
    gr1->Draw("psame");
    sprintf(name,"n#sigma_{k} mean");
    drawLatex(0.55,0.75,name,22,0.05,4);
    gPad->Update();
    
    pad->cd(2);
    setPad(gPad,0.15,0.15,0.15,0.15,0);
    h0= new TH1F("","",100,0,3.2);
    h0->SetMinimum(0.),
    h0->SetMaximum(2.);
    setHisto(h0,"","p (GeV/c)", "sigma");
    //setHisto(h0,"","p (GeV/c)", "ratio",1,0.05,0.05, 1,0.9,1.1, 1,0,3.2,0.,2);
    TGraphErrors* gr2 = new TGraphErrors(np,p_mean,sigma1,p_err,sigmaErr1);
    gr2->SetMarkerStyle(20);
    gr2->SetMarkerColor(1);
    gr2->SetMarkerSize(1.);
    h0->Draw();
    gr2->Draw("psame");
    sprintf(name,"n#sigma_{k} sigma");
    drawLatex(0.55,0.75,name,22,0.05,4);
    gPad->Update();
    
    pad->cd(3);
    setPad(gPad,0.15,0.15,0.15,0.15,0);
    h0= new TH1F("","",100,0,3.2);
    h0->SetMinimum(-1),
    h0->SetMaximum(1);
    setHisto(h0,"","p (GeV/c)", "mean");
    //setHisto(h0,"","p (GeV/c)", "ratio",1,0.05,0.05, 1,0.9,1.1, 1,0,3.2,0.,2);
    TGraphErrors* gr3 = new TGraphErrors(np,p_mean,mean2,p_err,meanErr2);
    gr3->SetMarkerStyle(20);
    gr3->SetMarkerColor(2);
    gr3->SetMarkerSize(1.);
    h0->Draw();
    gr3->Draw("psame");
    sprintf(name,"n#sigma_{k}^{TOF} mean");
    drawLatex(0.55,0.75,name,22,0.05,4);
    gPad->Update();
    
    pad->cd(4);
    setPad(gPad,0.15,0.15,0.15,0.15,0);
    h0= new TH1F("","",100,0,3.2);
    h0->SetMinimum(0.),
    h0->SetMaximum(2.);
    setHisto(h0,"","p (GeV/c)", "sigma");
    //setHisto(h0,"","p (GeV/c)", "ratio",1,0.05,0.05, 1,0.9,1.1, 1,0,3.2,0.,2);
    TGraphErrors* gr4 = new TGraphErrors(np,p_mean,sigma2,p_err,sigmaErr2);
    gr4->SetMarkerStyle(20);
    gr4->SetMarkerColor(2);
    gr4->SetMarkerSize(1.);
    h0->Draw();
    gr4->Draw("psame");
    sprintf(name,"n#sigma_{k}^{TOF} sigma");
    drawLatex(0.55,0.75,name,22,0.05,4);
    gPad->Update();
    
    pad->cd();
    pad->Update();
    
    pdf->On();
    pdf->NewPage();
    pad->cd();
    pad->Clear();
    pad->SetPad(0,0.25,1,0.75);
    pad->Divide(2,1);
    
    pad->cd(1);
    setPad(gPad,0.15,0.15,0.15,0.15,0);
    h0= new TH1F("","",100,0,3.2);
    h0->SetMinimum(0.7),
    h0->SetMaximum(1.02);
    setHisto(h0,"","p (GeV/c)", "n#sigma_{k} cut eff");
    //setHisto(h0,"","p (GeV/c)", "ratio",1,0.05,0.05, 1,0.9,1.1, 1,0,3.2,0.,2);
    TGraphErrors* gr5 = new TGraphErrors(np,p_mean,eff1,p_err,efferr1);
    gr5->SetMarkerStyle(20);
    gr5->SetMarkerColor(1);
    gr5->SetMarkerSize(1.);
    h0->Draw();
    gr5->Draw("psame");
    sprintf(name,"n#sigma_{k} eff");
    drawLatex(0.55,0.75,name,22,0.05,4);
    TLine* ll = new TLine(0.3,0.7,0.3,1.0);
    ll->SetLineColor(4);
    ll->SetLineStyle(7);
    ll->SetLineWidth(2);
    ll->Draw("same");
    gPad->Update();
    
    pad->cd(2);
    setPad(gPad,0.15,0.15,0.15,0.15,0);
    h0= new TH1F("","",100,0,3.2);
    h0->SetMinimum(0.7),
    h0->SetMaximum(1.02);
    setHisto(h0,"","p (GeV/c)", "n#sigma_{k}^{TOF} cut eff");
    //setHisto(h0,"","p (GeV/c)", "ratio",1,0.05,0.05, 1,0.9,1.1, 1,0,3.2,0.,2);
    TGraphErrors* gr6 = new TGraphErrors(np,p_mean,eff2,p_err,efferr2);
    gr6->SetMarkerStyle(20);
    gr6->SetMarkerColor(2);
    gr6->SetMarkerSize(1.);
    h0->Draw();
    gr6->Draw("psame");
    sprintf(name,"n#sigma_{k}^{TOF} eff");
    drawLatex(0.55,0.75,name,22,0.05,4);
    ll->Draw("same");
    gPad->Update();
    
    cout << "p\tTOF PID eff.\tTPC PID eff." << endl;
    for(int i=0; i<np; i++) cout << p_mean[i] << "\t" << eff2[i] << "\t" << eff1[i] << endl;
    ofstream out1("data/kaonPid_eff.txt");
    out1 << "p" << "\t" << "p_err" << "\t" << "TPC PID eff." << "\t" << "err" << "\t" << "TOF PID eff." << "\t" << "err" << endl;
    out1<<np<<endl;
    for(int i=0; i<np; i++) {
        out1<<p_mean[i]<<"\t"<<p_err[i]<<"\t";
        out1<<eff1[i]<<"\t"<<efferr1[i]<<"\t";
        out1<<eff2[i]<<"\t"<<efferr2[i]<<endl;
    }
    out1.close();
    
    pad->cd();
    pad->Update();
    
    pdf->On();
    pdf->Close();
    delete pdf;
}

#include "../myFunction.h"
#include "../myConst.h"
#include "../anaCuts.h"
TF1* ftemp;
double funMultiplyPt(double* x, double* par) {
    return ftemp->Eval(x[0]) * 2.* TMath::Pi() * x[0];
}
void write_RAA_Integral() {
    globalSetting();
    const float PI = TMath::Pi();
    
    //define Nbin
    
    //read D0 yield in AuAu
    TGraphErrors* gD0_Run14HFT_err[ncent];
    TGraphErrors* gD0_Run14HFT_sys[ncent];
    TFile* fin = new TFile("../sys/D0_Spectra_Run14_HFT_beforePtShift.root");
    for(int icent=0; icent<ncent; icent++) {
        gD0_Run14HFT_err[icent] = (TGraphErrors*)fin->Get(Form("gD0_err_%s",nameCent1[icent]));
        gD0_Run14HFT_sys[icent] = (TGraphErrors*)fin->Get(Form("gD0_sys_%s",nameCent1[icent]));
    }
    fin->Close();
    
    //integral
    float pt_meanX[npt];
    for(int ipt=0; ipt<npt; ipt++) {
        pt_meanX[ipt] = 0.5*(nptbin[ipt]+nptbin[ipt+1]);
        float ptW = nptbin[ipt+1] - nptbin[ipt];
        for(int icent=0; icent<ncent; icent++) {
            gD0_Run14HFT_err[icent]->GetY()[ipt] *= (ptW*pt_meanX[ipt]*2*PI);
            gD0_Run14HFT_err[icent]->GetEY()[ipt] *= (ptW*pt_meanX[ipt]*2*PI);
            gD0_Run14HFT_sys[icent]->GetY()[ipt] *= (ptW*pt_meanX[ipt]*2*PI);
            gD0_Run14HFT_sys[icent]->GetEY()[ipt] *= (ptW*pt_meanX[ipt]*2*PI);
        }
    }
    
    //read D0 yield in pp
    TFile* fin1 = new TFile("data/B_and_D0_ptSpectra.root");
    TF1* fppbase = (TF1*)fin1->Get("fD0_pp_Run09");
    fppbase->SetName("fppbase");
    TF1* fppbase_lw = (TF1*)fin1->Get("fD0_pp_lw_Run09");
    fppbase_lw->SetName("fppbase_lw");
    TF1* fppbase_up = (TF1*)fin1->Get("fD0_pp_up_Run09");
    fppbase_up->SetName("fppbase_up");
    fin1->Close();
    
    //calculate pp baseline
    float ppbase[npt], ppbase_lw[npt], ppbase_up[npt];
    const float factor = 0.565/42.;
    /*for(int i=0; i<npt; i++) {
        float pt = gD0_Run14HFT_err[0]->GetX()[i];
        ppbase[i] = fppbase->Eval(pt)*factor;
        float ppbaseErr_lw = fabs(fppbase->Eval(pt)*factor - fppbase_lw->Eval(pt)*factor);
        ppbaseErr_lw = sqrt(pow(ppbaseErr_lw,2) - pow(0.06*ppbase[i],2));
        ppbase_lw[i] = ppbase[i] - ppbaseErr_lw; //cancel TPC track sys error
        float ppbaseErr_up = fabs(fppbase->Eval(pt)*factor - fppbase_lw->Eval(pt)*factor);
        ppbaseErr_up = sqrt(pow(ppbaseErr_up,2) - pow(0.06*ppbase[i],2));
        ppbase_up[i] = ppbase[i] + ppbaseErr_up; //cancel TPC track sys error
    }*/
    ftemp = (TF1*)fppbase->Clone("fppbaseClone");
    TF1* fppbaseI = new TF1("fppbaseI",funMultiplyPt,0,10,0);
    for(int ipt=0; ipt<npt; ipt++) {
        float ptW = nptbin[ipt+1] - nptbin[ipt];
        ppbase[ipt] = fppbaseI->Integral(nptbin[ipt],nptbin[ipt+1])*factor;
    }
    ftemp = (TF1*)fppbase_lw->Clone("fppbase_lwClone");
    TF1* fppbaseI_lw = new TF1("fppbaseI_lw",funMultiplyPt,0,10,0);
    for(int ipt=0; ipt<npt; ipt++) {
        float ptW = nptbin[ipt+1] - nptbin[ipt];
        ppbase_lw[ipt] = fppbaseI_lw->Integral(nptbin[ipt],nptbin[ipt+1])*factor;
        float ppbaseErr_lw = fabs(ppbase[ipt] - ppbase_lw[ipt]);
        ppbaseErr_lw = sqrt(pow(ppbaseErr_lw,2) - pow(0.06*ppbase[ipt],2));
        ppbase_lw[ipt] = ppbase[ipt] - ppbaseErr_lw; //cancel TPC track sys error
    }
    ftemp = (TF1*)fppbase_up->Clone("fppbase_upClone");
    TF1* fppbaseI_up = new TF1("fppbaseI_up",funMultiplyPt,0,10,0);
    for(int ipt=0; ipt<npt; ipt++) {
        float ptW = nptbin[ipt+1] - nptbin[ipt];
        ppbase_up[ipt] = fppbaseI_up->Integral(nptbin[ipt],nptbin[ipt+1])*factor;
        float ppbaseErr_up = fabs(ppbase[ipt] - ppbase_up[ipt]);
        ppbaseErr_up = sqrt(pow(ppbaseErr_up,2) - pow(0.06*ppbase[ipt],2));
        ppbase_up[ipt] = ppbase[ipt] + ppbaseErr_up; //cancel TPC track sys error
    }
    
    //caculate pp baseline error band
    float ppErr_up[npt], ppErr_lw[npt];
    for(int i=0; i<npt; i++) {
        ppErr_lw[i] = (ppbase[i] - ppbase_lw[i]);// /ppbase[i];
        ppErr_up[i] = (ppbase_up[i] - ppbase[i]);// /ppbase[i];
        //cout << ppErr_lw[i] << "\t" << ppErr_up[i] << endl;
    }
    
    //RAA and its errors
    float pt_mean[ncent][npt], y[ncent][npt], yerr[ncent][npt], ysys[ncent][npt];
    float raa[ncent][npt], raaErr[ncent][npt], raaSys[ncent][npt]; //ratio err
    float raa_pp_up[ncent][npt], raa_pp_lw[ncent][npt]; //pp err
    ofstream out("D0RAA_Run14HFT.txt");
    for(int icent=0; icent<ncent; icent++) {
        out << nameCent[icent] << endl;
        out << "pT \t RAA \t Statistical error \t Systematic error \t pp error (low) \t pp error (up)" << endl;
        for(int i=0; i<npt; i++) {
            pt_mean[icent][i] = gD0_Run14HFT_err[icent]->GetX()[i];
            y[icent][i] = gD0_Run14HFT_err[icent]->GetY()[i];
            yerr[icent][i] = gD0_Run14HFT_err[icent]->GetEY()[i];
            ysys[icent][i] = sqrt(pow(gD0_Run14HFT_sys[icent]->GetEY()[i],2) - pow(y[icent][i]*0.06,2));  //cancel TPC track sys error
            //cout << y[icent][i] << "\t" << ppbase[i] << endl;
            raa[icent][i] = y[icent][i]/ppbase[i] /NbinMean[icent];
            raaErr[icent][i] = raa[icent][i] * yerr[icent][i]/y[icent][i];
            raaSys[icent][i] = raa[icent][i] * ysys[icent][i]/y[icent][i];
            
            raa_pp_lw[icent][i] = raa[icent][i] * ppErr_up[i]/ppbase[i]; //pp is in denominator
            raa_pp_up[icent][i] = raa[icent][i] * ppErr_lw[i]/ppbase[i];
            
            out << pt_mean[icent][i] << "\t" << raa[icent][i] << "\t" << raaErr[icent][i] << "\t" << raaSys[icent][i] << "\t" << raa_pp_lw[icent][i] << "\t" << raa_pp_up[icent][i] << endl;
        }
        out << endl;
    }
    out.close();
    TGraphErrors* gRAA[ncent];
    TGraphErrors* gRAA_sys[ncent];
    TGraphAsymmErrors* gRAA_pp[ncent];
    for(int icent=0; icent<ncent; icent++) {
        gRAA[icent] = new TGraphErrors(npt,pt_mean[icent],raa[icent],0,raaErr[icent]);
        gRAA_sys[icent] = new TGraphErrors(npt,pt_mean[icent],raa[icent],0,raaSys[icent]);
        gRAA_pp[icent] = new TGraphAsymmErrors(npt,pt_mean[icent],raa[icent],0,0,raa_pp_lw[icent],raa_pp_up[icent]);
    }
    
    //Write
    TFile* fout = new TFile("D0RAA_Run14HFT.root","RECREATE");
    fout->cd();
    for(int icent=0; icent<ncent; icent++) {
        gRAA[icent]->SetMarkerStyle(MARKERSTYLE[0]);
        gRAA[icent]->SetMarkerColor(COLOR[0]);
        gRAA[icent]->SetMarkerSize(2.5);
        gRAA[icent]->Write(Form("D0_RAA_err_%s",nameCent1[icent]));
        gRAA_sys[icent]->SetMarkerStyle(MARKERSTYLE[0]);
        gRAA_sys[icent]->SetMarkerColor(COLOR[0]);
        gRAA_sys[icent]->SetMarkerSize(2.5);
        gRAA_sys[icent]->Write(Form("D0_RAA_sys_%s",nameCent1[icent]));
        gRAA_pp[icent]->SetMarkerStyle(MARKERSTYLE[0]);
        gRAA_pp[icent]->SetMarkerColor(COLOR[0]);
        gRAA_pp[icent]->SetMarkerSize(2.5);
        gRAA_pp[icent]->Write(Form("D0_RAA_pperr_%s",nameCent1[icent]));
    }
    fout->Close();
}

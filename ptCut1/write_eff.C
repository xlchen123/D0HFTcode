#include "../anaCuts.h"
void write_eff() {
    char name[250];
    
    for(int icent=0; icent<9; icent++) {
        for(int i=0; i<ncent; i++) {
            if(icent>=centLw[i]&&icent<centUp[i]) Nbin_Sum[i] += Nbin[icent];
        }
    }
    
    // read vtx Corr
    TH1D* hmean[ncentAll_vtx];
    TH1D* herr[ncentAll_vtx];
    TFile* fin = new TFile("../vtxCorr/data/vtxCorr_pt0_3.root");
    for(int icent=0; icent<ncentAll_vtx; icent++) {
        hmean[icent] = (TH1D*)fin->Get(Form("mean_%s",nameCentAll_vtx1[icent]));
        hmean[icent]->SetDirectory(0);
        herr[icent] = (TH1D*)fin->Get(Form("sysError_%s",nameCentAll_vtx1[icent]));
        herr[icent]->SetDirectory(0);
    }
    fin->Close();
    
    //read hist
    // TFile* fsimuD0 = new TFile("../D0_eff.root"); //D0
    TFile* fsimuD0 = new TFile("../D0_eff_secondTrack.root"); //D0
    TH2F* h2Pt_D0 = (TH2F*)fsimuD0->Get("h2Pt");
    h2Pt_D0->SetDirectory(0);
    TH2F* h2PtCut_D0 = (TH2F*)fsimuD0->Get("h2PtCut_pt1");
    h2PtCut_D0->SetDirectory(0);
    fsimuD0->Close();
    
    //Rebin
    int nRebin = 1;
    h2Pt_D0->RebinX(nRebin); h2PtCut_D0->RebinX(nRebin);
    
    //eff in 9 cent bin and cent combined
    TH1F* hpt[9];
    TH1F* hptCut[9];
    TH1D* hptBin[9];
    TH1D* hptBinCut[9];
    TH1F* heffD0_inCent[9];
    TH1D* heffBinD0_inCent[9];
    TH1F* heffD0[ncent];
    TH1D* heffBinD0[ncent];
    for(int icent=0; icent<ncent; icent++) {
        sprintf(name,"heffD0_%i",icent);
        heffD0[icent] = new TH1F(name,name,npt_eff, ptLw_eff, ptUp_eff);
        heffD0[icent]->Sumw2();
        sprintf(name,"heffBinD0_%i",icent);
        heffBinD0[icent] = new TH1D(name,name,npt,nptbin);
        heffBinD0[icent]->Sumw2();
    }
    
    for(int i=0; i<9; i++) {
        //
        hpt[i] = (TH1F*)h2Pt_D0->ProjectionX(Form("hpt_%i",i),i+1,i+1);
        hptCut[i] = (TH1F*)h2PtCut_D0->ProjectionX(Form("hptCut_%i",i),i+1,i+1);
        // bin binning--rebin
        hptBin[i] = (TH1D*)hpt[i]->Rebin(npt,Form("hptBin_%i",i),nptbin);
        hptBinCut[i] = (TH1D*)hptCut[i]->Rebin(npt,Form("hptBinCut_%i",i),nptbin);
        
        //caculate eff in 9 cent bin
        heffD0_inCent[i] = new TH1F(*hpt[i]);
        heffD0_inCent[i]->Divide(hptCut[i],hpt[i],1,1);
        heffBinD0_inCent[i] = new TH1D(*hptBin[i]);
        heffBinD0_inCent[i]->Divide(hptBinCut[i],hptBin[i],1,1);
    }
    
    // deal with the drop at > 5GeV in cent 70-80%
    // int icent=0;
    // const float ptLW = 5.;
    // int bin_lw = heffD0_inCent[icent]->FindBin(ptLW+1.e-6);
    // int bin_up = heffD0_inCent[icent]->GetNbinsX();
    // for(int ibin=bin_lw; ibin<=bin_up; ibin++) {
    //     heffD0_inCent[icent]->SetBinContent(ibin,heffD0_inCent[icent+1]->GetBinContent(ibin));
    //     heffD0_inCent[icent]->SetBinError(ibin,heffD0_inCent[icent+1]->GetBinError(ibin));
    // }
    // bin_lw = heffBinD0_inCent[icent]->FindBin(ptLW+1.e-6);
    // bin_up = heffBinD0_inCent[icent]->GetNbinsX();
    // for(int ibin=bin_lw; ibin<=bin_up; ibin++) {
    //     heffBinD0_inCent[icent]->SetBinContent(ibin,heffBinD0_inCent[icent+1]->GetBinContent(ibin));
    //     heffBinD0_inCent[icent]->SetBinError(ibin,heffBinD0_inCent[icent+1]->GetBinError(ibin));
    // }
    
    /////////////////////////  multiply vtx corr ratio  /////////////////////////
    /////////////////////////////////////////////////////////////////////////////
    // binning
    //heffBinD0_inCent[icent]->Multiply(hmean[icent]);
    for(int icent=0; icent<9; icent++) {
        for(int ibin=1; ibin<=heffBinD0_inCent[icent]->GetNbinsX(); ibin++) {
            float ptCenter = heffBinD0_inCent[icent]->GetBinCenter(ibin);
            int binFind = hmean[icent]->FindBin(ptCenter);
            float R = hmean[icent]->GetBinContent(binFind);
            float Rerr = 0;//hmean[icent]->GetBinError(binFind);
            float num = heffBinD0_inCent[icent]->GetBinContent(ibin);
            float err = heffBinD0_inCent[icent]->GetBinError(ibin);
            err = num*R*sqrt(pow(Rerr/R,2)+pow(err/num,2));
            num = num*R;
            heffBinD0_inCent[icent]->SetBinContent(ibin,num);
            heffBinD0_inCent[icent]->SetBinError(ibin,err);
        }
        // detailed bin efficiency
        for(int ibin=1; ibin<=heffD0_inCent[icent]->GetNbinsX(); ibin++) {
            float ptCenter = heffD0_inCent[icent]->GetBinCenter(ibin);
            int binFind = ptCenter<8. ? hmean[icent]->FindBin(ptCenter): hmean[icent]->GetNbinsX();
            float R = hmean[icent]->GetBinContent(binFind);
            float Rerr = 0;//hmean[icent]->GetBinError(binFind);
            float num = heffD0_inCent[icent]->GetBinContent(ibin);
            float err = heffD0_inCent[icent]->GetBinError(ibin);
            err = num*R*sqrt(pow(Rerr/R,2)+pow(err/num,2));
            num = num*R;
            heffD0_inCent[icent]->SetBinContent(ibin,num);
            heffD0_inCent[icent]->SetBinError(ibin,err);
        }
    }
    /////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////
    
    ////////////////////////////  combine centrality  ///////////////////////////
    /////////////////////////////////////////////////////////////////////////////
    for(int i=0; i<9; i++) {
        //combine cent bin
        for(int icent=0; icent<ncent; icent++) {
            if(i>=centLw[icent]&&i<centUp[icent]) {
                // no correlation
                //heffD0[icent]->Add(heffD0_inCent[i],Nbin[i]/Nbin_Sum[icent]);
                //heffBinD0[icent]->Add(heffBinD0_inCent[i],Nbin[i]/Nbin_Sum[icent]);
                
                // total correlation
                for(int ibin=1; ibin<=heffD0[icent]->GetNbinsX(); ibin++) {
                    float sum = heffD0[icent]->GetBinContent(ibin);
                    float sumErr = heffD0[icent]->GetBinError(ibin);
                    float num = heffD0_inCent[i]->GetBinContent(ibin)*Nbin[i]/Nbin_Sum[icent];
                    float numErr = heffD0_inCent[i]->GetBinError(ibin)*Nbin[i]/Nbin_Sum[icent];
                    sum += num;
                    sumErr += numErr;
                    heffD0[icent]->SetBinContent(ibin,sum);
                    heffD0[icent]->SetBinError(ibin,sumErr);
                }
                for(int ibin=1; ibin<=heffBinD0[icent]->GetNbinsX(); ibin++) {
                    float sum = heffBinD0[icent]->GetBinContent(ibin);
                    float sumErr = heffBinD0[icent]->GetBinError(ibin);
                    float num = heffBinD0_inCent[i]->GetBinContent(ibin)*Nbin[i]/Nbin_Sum[icent];
                    float numErr = heffBinD0_inCent[i]->GetBinError(ibin)*Nbin[i]/Nbin_Sum[icent];
                    sum += num;
                    sumErr += numErr;
                    heffBinD0[icent]->SetBinContent(ibin,sum);
                    heffBinD0[icent]->SetBinError(ibin,sumErr);
                }
            }
        }
    }
    // read efficiency  ---end: see default/write_eff.C
    
    //write eff
    TFile* fout = new TFile("data/eff.root","RECREATE");
    for(int icent=0; icent<ncent; icent++) {
        heffD0[icent]->SetTitle(nameCent[icent]);
        heffD0[icent]->Write();
        heffBinD0[icent]->SetTitle(nameCent[icent]);
        heffBinD0[icent]->Write();
    }
    fout->Close();
}

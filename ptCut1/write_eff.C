#include "../anaCuts.h"
void write_eff() {
    char name[250];
    
    for(int icent=0; icent<9; icent++) {
        for(int i=0; i<ncent; i++) {
            if(icent>=centLw[i]&&icent<centUp[i]) Nbin_Sum[i] += Nbin[icent];
        }
    }
    
    //read hist
    // TFile* fsimuD0 = new TFile("../D0_eff.root"); //D0
    // TFile* fsimuD0 = new TFile("../D0_eff_secondTrack.root"); //D0
    // TFile* fsimuD0 = new TFile("../D0_eff_secondTrack_v3.root"); //D0 Y
    // TFile* fsimuD0 = new TFile("../D0_eff_pos.root"); //D0 pos
    TFile* fsimuD0 = new TFile("../D0_eff_neg.root"); //D0 pos
    TH3F* h3Pt_D0 = (TH3F*)fsimuD0->Get("h3PtCentY");
    h3Pt_D0->SetDirectory(0);
    TH3F* h3PtCut_D0 = (TH3F*)fsimuD0->Get("h3PtCentYCut_pt1");
    h3PtCut_D0->SetDirectory(0);
    fsimuD0->Close();
    
    // const float YStart = -0.5;
    // const float YEnd = 0.5;
    const float YStart = -1.0;
    const float YEnd = 1.0;

    //Rebin
    int nRebin = 1;
    h3Pt_D0->RebinX(nRebin); h3PtCut_D0->RebinX(nRebin);
    
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
        int ybin_lw = h3Pt_D0->GetZaxis()->FindBin(YStart+1e-6);
        int ybin_up = h3Pt_D0->GetZaxis()->FindBin(YEnd-1e-6);

        hpt[i] = (TH1F*)h3Pt_D0->ProjectionX(Form("hpt_%i",i),i+1,i+1,ybin_lw,ybin_up);
        hptCut[i] = (TH1F*)h3PtCut_D0->ProjectionX(Form("hptCut_%i",i),i+1,i+1,ybin_lw,ybin_up);
        // bin binning--rebin
        hptBin[i] = (TH1D*)hpt[i]->Rebin(npt,Form("hptBin_%i",i),nptbin);
        hptBinCut[i] = (TH1D*)hptCut[i]->Rebin(npt,Form("hptBinCut_%i",i),nptbin);
        
        //caculate eff in 9 cent bin
        heffD0_inCent[i] = new TH1F(*hpt[i]);
        heffD0_inCent[i]->Divide(hptCut[i],hpt[i],vtxWeight[i],1);
        heffBinD0_inCent[i] = new TH1D(*hptBin[i]);
        heffBinD0_inCent[i]->Divide(hptBinCut[i],hptBin[i],vtxWeight[i],1);
        
        //combine cent bin
        for(int icent=0; icent<ncent; icent++) {
            if(i>=centLw[icent]&&i<centUp[icent]) {
                heffD0[icent]->Add(heffD0_inCent[i],Nbin[i]/Nbin_Sum[icent]);
                heffBinD0[icent]->Add(heffBinD0_inCent[i],Nbin[i]/Nbin_Sum[icent]);
            }
        }
    }
    
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

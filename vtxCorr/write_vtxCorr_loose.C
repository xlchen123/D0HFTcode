#include "../anaCuts.h"
void write_vtxCorr_loose() {
    char dir[250];
    char CMD[250];
    char name[250];
    char title[250];
    
    sprintf(dir,"data");
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
    
    // cal. correction ratio and its sys error: 30-80%
    TH1F* hmean[ncent];
    TH1F* herr[ncent];
    for(int icent=0; icent<ncent_vtx-1; icent++) {
        hmean[icent] = (TH1F*)hRatio[0][icent]->Clone(Form("mean_%s",nameCentAll_vtx1[icent]));
        herr[icent] = (TH1F*)hRatio[0][icent]->Clone(Form("sysError_%s",nameCentAll_vtx1[icent]));
        for(int ibin=1; ibin<=hRatio[0][icent]->GetNbinsX(); ibin++) {
            float num1 = hRatio[0][icent]->GetBinContent(ibin);
            float num2 = hRatio[1][icent]->GetBinContent(ibin);
            float error;
            if((num1+num2)==0) error = 0;
            else error = fabs(num1 - num2)/(num1 + num2);
            herr[icent]->SetBinContent(ibin,error);
            hmean[icent]->SetBinContent(ibin,(num1+num2)/2.);
            herr[icent]->SetBinError(ibin,0);
            hmean[icent]->SetBinError(ibin,fabs(num1 - num2)/2.);
        }
    }
    
    // cal. correction ratio and its sys error: 20-30%
    int icent = ncent_vtx-1;
    hmean[icent] = (TH1F*)hRatio[0][icent]->Clone(Form("mean_%s",nameCentAll_vtx1[icent]));
    herr[icent] = (TH1F*)hRatio[0][icent]->Clone(Form("sysError_%s",nameCentAll_vtx1[icent]));
    for(int ibin=1; ibin<=hmean[icent]->GetNbinsX(); ibin++) {
        float num1 = hRatio[0][icent]->GetBinContent(ibin);
        float num2 = 1;
        float error;
        if((num1+num2)==0) error = 0;
        else error = fabs(num1 - num2)/(num1 + num2);
        herr[icent]->SetBinContent(ibin,error);
        hmean[icent]->SetBinContent(ibin,(num1+num2)/2.);
        herr[icent]->SetBinError(ibin,0);
        hmean[icent]->SetBinError(ibin,fabs(num1 - num2)/2.);
    }
    
    // cal. correction ratio and its sys error: 0-20%
    for(int icent=ncent_vtx; icent<ncentAll_vtx; icent++) {
        hmean[icent] = (TH1F*)hRatio[0][icent]->Clone(Form("mean_%s",nameCentAll_vtx1[icent]));
        herr[icent] = (TH1F*)hRatio[0][icent]->Clone(Form("sysError_%s",nameCentAll_vtx1[icent]));
        for(int ibin=1; ibin<=hmean[icent]->GetNbinsX(); ibin++) {
            //float num1 = 1;
            //float num2 = 1;
            //float error;
            //if((num1+num2)==0) error = 0;
            //else error = fabs(num1 - num2)/(num1 + num2);
            //herr[icent]->SetBinContent(ibin,error);
            //hmean[icent]->SetBinContent(ibin,(num1+num2)/2.);
            float num = 1;
            float error = 0.05;
            herr[icent]->SetBinContent(ibin,error);
            hmean[icent]->SetBinContent(ibin,num);
            herr[icent]->SetBinError(ibin,0);
            hmean[icent]->SetBinError(ibin,num*error);
        }
    }
    
    // write
    sprintf(name,"%s/vtxCorr_%s.root",dir,nameType[0]);
    TFile* fout = new TFile(name,"RECREATE");
    for(int icent=0; icent<ncentAll_vtx; icent++) {
        hmean[icent]->Write();
        herr[icent]->Write();
    }
    fout->Close();
}

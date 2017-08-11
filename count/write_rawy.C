#include "../anaCuts.h"
void write_rawy()
{
    char buf[250];
    char dir1[250];
    char name[250];
    char title[250];
    char CMD[250];
    const double PI = TMath::Pi();
    float pt_lw,    pt_up,    m_lw,    m_up,    cent_lw,    cent_up ;
    int   ptbin_lw, ptbin_up, mbin_lw, mbin_up, centbin_lw, centbin_up;
    
    const float dauPtCut = 0.6;
    const float mR_lw = 1.6;
    const float mR_up = 2.1;
    //read
    TFile* fin = new TFile("../D0_data_mix.root");
    THnF* hnmassUL = (THnF*)fin->Get("hD0CentPtEtaMDphiDaug_standard");
    //hnmassUL->SetDirectory(0);
    THnF* hnmassLS = (THnF*)fin->Get("hD0CentPtEtaMDphiDaugLikeSign_standard");
    //hnmassLS->SetDirectory(0);
    THnF* hnmassMixUL = (THnF*)fin->Get("hD0CentPtEtaMDphiDaugMixed_standard");
    //hnmassMixUL->SetDirectory(0);
    THnF* hnmassMixLS = (THnF*)fin->Get("hD0CentPtEtaMDphiDaugLikeSignMixed_standard");
    //hnmassMixLS->SetDirectory(0);
    fin->Close();
    
    //cout << hnmassMixLS->GetNbins() << endl;
    
    //write
    sprintf(dir1,"data");
    sprintf(CMD,"[ -d %s ] || mkdir -p %s",dir1,dir1);
    gSystem->Exec(CMD);
    
    TH1F* hInvMass;
    TH1F* hInvMassUL;
    TH1F* hInvMassLS;
    TH1F* hInvMassMix;
    TH1F* hInvMassMixUL;
    TH1F* hInvMassMixLS;
    TFile* fout = new TFile("data/D0_HFT_signal.root","RECREATE");
    TDirectory* dir;
    for(int icent=0; icent<ncent; icent++) {
        fout->cd();
        TDirectory* dir = NULL;
        if(!(dir = (TDirectory*)fout->Get(Form("%s",nameCent1[icent]))))
            dir = (TDirectory*)fout->mkdir(Form("%s",nameCent1[icent]));
        for(int ipt=0; ipt<npt; ipt++) {
            centbin_lw = hnmassUL->GetAxis(0)->FindBin(centLw[icent]+1e-6);
            centbin_up = hnmassUL->GetAxis(0)->FindBin(centUp[icent]-1e-6);
            ptbin_lw = hnmassUL->GetAxis(1)->FindBin(nptbin[ipt]+1e-6);
            ptbin_up = hnmassUL->GetAxis(1)->FindBin(nptbin[ipt+1]-1e-6);
            int dau1PtBin = hnmassUL->GetAxis(2)->FindBin(dauPtCut+1e-6);
            int dau2PtBin = hnmassUL->GetAxis(4)->FindBin(dauPtCut+1e-6);
            int dau1NPtBin = 11;  //11
            int dau2NPtBin = 11;  //11
            
            cout << "pt bin: " << ptbin_lw << "-" << ptbin_up << ", pt range: " << nptbin[ipt] << "-" << nptbin[ipt+1] << endl;
            
            // projection
            // 1 - unlike sign
            hnmassUL->GetAxis(0)->SetRange(centbin_lw,centbin_up);
            hnmassUL->GetAxis(1)->SetRange(ptbin_lw,ptbin_up);
            hnmassUL->GetAxis(2)->SetRange(dau1PtBin,dau1NPtBin);
            hnmassUL->GetAxis(4)->SetRange(dau2PtBin,dau2NPtBin);
            hInvMassUL = (TH1F*)hnmassUL->Projection(3,"E");
            sprintf(name,"hInvMassUL_%i_%i",icent,ipt);
            hInvMassUL->SetName(name);
            sprintf(title,"unlike sign %s %.1f<p_{T}<%.1f",nameCent[icent],nptbin[ipt],nptbin[ipt+1]);
            hInvMassUL->SetTitle(title);
            dir->cd();
            hInvMassUL->GetXaxis()->SetRangeUser(1.6,2.1);
            hInvMassUL->Write();
            dir->Save();
            //cout << "dauPtCut = " << dauPtCut << endl;
            //cout << "dauPtCut bin = " << dau1PtBin << "\t" << dau2PtBin << endl;
            //cout << "centbin = " << centbin_lw << "-" << centbin_up << endl;
            //cout << hInvMassUL->Integral(161,210) << endl;
            
            // 2 - like sign
            hnmassLS->GetAxis(0)->SetRange(centbin_lw,centbin_up);
            hnmassLS->GetAxis(1)->SetRange(ptbin_lw,ptbin_up);
            hnmassLS->GetAxis(2)->SetRange(dau1PtBin,dau1NPtBin);
            hnmassLS->GetAxis(4)->SetRange(dau2PtBin,dau2NPtBin);
            hInvMassLS = (TH1F*)hnmassLS->Projection(3,"E");
            sprintf(name,"hInvMassLS_%i_%i",icent,ipt);
            hInvMassLS->SetName(name);
            sprintf(title,"like sign %s %.1f<p_{T}<%.1f",nameCent[icent],nptbin[ipt],nptbin[ipt+1]);
            hInvMassLS->SetTitle(title);
            dir->cd();
            hInvMassLS->GetXaxis()->SetRangeUser(1.6,2.1);
            hInvMassLS->Write();
            dir->Save();
            
            // 3 - unlike-sign - like-sign
            hInvMass = (TH1F*)hInvMassUL->Clone(Form("hInvMassUL_clone1_%i%i",icent,ipt));
            sprintf(name,"hInvMass_%i_%i",icent,ipt);
            sprintf(title,"signal(UL-LS) %s %.1f<p_{T}<%.1f",nameCent[icent],nptbin[ipt],nptbin[ipt+1]);
            hInvMass->SetName(name);
            hInvMass->SetTitle(title);
            hInvMass->Add(hInvMassUL,hInvMassLS,1,-1);
            dir->cd();
            hInvMass->GetXaxis()->SetRangeUser(1.6,2.1);
            hInvMass->Write();
            dir->Save();
            
            // 4 - mixEvent like sign
            hnmassMixLS->GetAxis(0)->SetRange(centbin_lw,centbin_up);
            hnmassMixLS->GetAxis(1)->SetRange(ptbin_lw,ptbin_up);
            hnmassMixLS->GetAxis(2)->SetRange(dau1PtBin,dau1NPtBin);
            hnmassMixLS->GetAxis(4)->SetRange(dau2PtBin,dau2NPtBin);
            hInvMassMixLS = (TH1F*)hnmassMixLS->Projection(3,"E");
            sprintf(name,"hInvMassMixLS_%i_%i",icent,ipt);
            hInvMassMixLS->SetName(name);
            sprintf(title,"mix-event like sign %s %.1f<p_{T}<%.1f",nameCent[icent],nptbin[ipt],nptbin[ipt+1]);
            hInvMassMixLS->SetTitle(title);
            
            // 5 - mixEvent unlike sign
            hnmassMixUL->GetAxis(0)->SetRange(centbin_lw,centbin_up);
            hnmassMixUL->GetAxis(1)->SetRange(ptbin_lw,ptbin_up);
            hnmassMixUL->GetAxis(2)->SetRange(dau1PtBin,dau1NPtBin);
            hnmassMixUL->GetAxis(4)->SetRange(dau2PtBin,dau2NPtBin);
            hInvMassMixUL = (TH1F*)hnmassMixUL->Projection(3,"E");
            sprintf(name,"hInvMassMixUL_%i_%i",icent,ipt);
            hInvMassMixUL->SetName(name);
            sprintf(title,"mix-event unlike sign %s %.1f<p_{T}<%.1f",nameCent[icent],nptbin[ipt],nptbin[ipt+1]);
            hInvMassMixUL->SetTitle(title);

            cout << " cent : " << nameCent[icent] << endl;
            cout << "axis 0 : " << centbin_lw << " " << centbin_up << endl;
            cout << "axis 1 : " << ptbin_lw << " " <<ptbin_up << endl;
            cout << "axis 2 : " << dau1PtBin << " " << dau1NPtBin << endl;
            cout << "axis 4 : " << dau2PtBin << " " <<dau2NPtBin << endl;
            
            // scale mix-event unlike
            mbin_lw = hInvMassLS->FindBin(mR_lw+1.e-6);
            mbin_up = hInvMassLS->FindBin(mR_up-1.e-6);
            //cout << "N bin number: " << hInvMassLS->GetNbinsX() << endl;
            //cout << "mass range: " << mR_lw << "\t" << mR_up << endl;
            //cout << "mass bin range:" << mbin_lw << "\t" << mbin_up << endl;
            //cout << hInvMassMixLS->Integral(mbin_lw,mbin_up) << "\t" << hInvMassMixUL->Integral(mbin_lw,mbin_up) << endl;
            const float norm = hInvMassMixLS->Integral(mbin_lw,mbin_up) ? hInvMassLS->Integral(mbin_lw,mbin_up) / hInvMassMixLS->Integral(mbin_lw,mbin_up) : 0;
            cout << "norm.=" << norm << endl;
            hInvMassMixUL->Scale(norm);
            hInvMassMixLS->Scale(norm);
            //cout << hInvMassMixLS->Integral(mbin_lw,mbin_up) << "\t" << hInvMassMixUL->Integral(mbin_lw,mbin_up) << "\t" << hInvMassLS->Integral(mbin_lw,mbin_up) << endl;
            
            // write mix-event
            dir->cd();
            hInvMassMixUL->GetXaxis()->SetRangeUser(1.6,2.1);
            hInvMassMixUL->Write();
            dir->Save();
            
            hInvMassMix = (TH1F*)hInvMassUL->Clone(Form("hInvMassUL_clone2_%i%i",icent,ipt));
            sprintf(name,"hInvMassMix_%i_%i",icent,ipt);
            sprintf(title,"signal(UL-MixUL) %s %.1f<p_{T}<%.1f",nameCent[icent],nptbin[ipt],nptbin[ipt+1]);
            hInvMassMix->SetName(name);
            hInvMassMix->SetTitle(title);
            hInvMassMix->Add(hInvMassUL,hInvMassMixUL,1,-1);
            dir->cd();
            hInvMassMix->GetXaxis()->SetRangeUser(1.6,2.1);
            hInvMassMix->Write();
            dir->Save();
        }
        dir->Close();
        cout << nameCent[icent] << ", cent bin: " << centbin_lw << "-" << centbin_up << endl;
    }
    fout->Close();
    
    if(hnmassUL) {delete hnmassUL; hnmassUL=NULL;}
    if(hnmassLS) {delete hnmassLS; hnmassLS=NULL;}
    if(hnmassMixUL) {delete hnmassMixUL; hnmassMixUL=NULL;}
    if(hnmassMixLS) {delete hnmassMixLS; hnmassMixLS=NULL;}
    /*if(hInvMassUL) {delete hInvMassUL; hInvMassUL=NULL;}
    if(hInvMassLS) {delete hInvMassLS; hInvMassLS=NULL;}
    if(hInvMass) {delete hInvMass; hInvMass=NULL;}
    if(hInvMassMixUL) {delete hInvMassMixUL; hInvMassMixUL=NULL;}
    if(hInvMassMix) {delete hInvMassMix; hInvMassMix=NULL;}*/
}

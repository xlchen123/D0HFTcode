#include "../anaCuts.h"
void fitRawCounts() {
    ifstream in;
    float ptbin[npt], rawy[npt], rawyErr[npt];
    float eff, efferr, y, yerr;
    for(int icent=0; icent<1; icent++) {
        in.open(Form("data/rawY_%s.txt",nameCent1[icent]));
        for(int ipt=0; ipt<npt; ipt++) in >> rawy[ipt] >> rawyErr[ipt];
        in.close();
        for(int ipt=0; ipt<npt; ipt++) {
            ptbin[ipt] = 0.5*(nptbin[ipt]+nptbin[ipt+1]);
            cout << nameCent[icent] << ", pt: " << nptbin[ipt] << "-" << nptbin[ipt+1] <<  endl;
        }
    }

    TGraphErrors *g1 = new TGraphErrors(npt,ptbin,rawy,0, rawyErr);
    g1->SetName("g1");
    g1->SetTitle("0-10\% signal counts");
    g1->Draw("APL");

    TF1 *f1= new TF1("f1","[0]+[1]*x*[2]*x*x+[3]*x*x*x",0,8);
    g1->Fit("f1","R");

    TFile* fout = new TFile("rawSignalCounts.root","RECREATE");
    g1->Write();
}

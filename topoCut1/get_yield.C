#include "../anaCuts.h"
void get_yield() {
    //read eff
    TH1D* heff[ncent];
    TFile* fin = new TFile("data/eff.root");
    for(int icent=0; icent<ncent; icent++) {
        heff[icent] = (TH1D*)fin->Get(Form("heffBinD0_%i",icent));
        heff[icent]->SetDirectory(0);
    }
    
    //do the eff correction
    ifstream in;
    ofstream out;
    float rawy[npt], rawyErr[npt];
    float eff, efferr, y, yerr;
    for(int icent=0; icent<ncent; icent++) {
        in.open(Form("data/rawY_%s.txt",nameCent1[icent]));
        for(int ipt=0; ipt<npt; ipt++) in >> rawy[ipt] >> rawyErr[ipt];
        in.close();
        out.open(Form("data/yield_%s.txt",nameCent1[icent]));
        for(int ipt=0; ipt<npt; ipt++) {
            float pt = 0.5*(nptbin[ipt]+nptbin[ipt+1]);
            int bin = heff[icent]->FindBin(pt);
            eff = heff[icent]->GetBinContent(bin);
            efferr = 0;//heff[icent]->GetBinError(bin);
            cout << nameCent[icent] << ", pt: " << nptbin[ipt] << "-" << nptbin[ipt+1] << ", eff error: " << efferr/eff << endl;
            y = rawy[ipt]/eff;
            yerr = y*sqrt(pow(rawyErr[ipt]/rawy[ipt],2)+pow(efferr/eff,2));
            out << y << "\t" << yerr << endl;
        }
        out.close();
    }
}

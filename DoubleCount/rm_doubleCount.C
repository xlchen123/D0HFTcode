#include "../anaCuts.h"
void rm_doubleCount() {
    //read Base double count TGraph -- mean and sys
    TFile* fin = new TFile("MisPID_SB_Final.root");
    TGraph* gdcRatio_mea[ncentBase_dc];
    TGraph* gdcRatio_sys[ncentBase_dc];
    for(int icent=0; icent<ncentBase_dc; icent++) {
        gdcRatio_mea[icent] = (TGraph*)fin->Get(Form("DoubleCounting_Cen_%i_SB",icent));
        gdcRatio_sys[icent] = (TGraph*)fin->Get(Form("DoubleCounting_Cen_%i_Full",icent));
    }
    fin->Close();
    
    //subtract double count and calculate its sys error
    ifstream in;
    ofstream out;
    float y[npt], yerr[npt], ybase[npt], yerrbase[npt];
    float sys;
    gSystem->Exec("[ -d data ] || mkdir -p data");
    for(int icent=0; icent<ncent; icent++) {
        // in.open(Form("../default/data/yield_%s.txt",nameCent1[icent]));
        in.open(Form("../default/data/re_yield_%s.txt",nameCent1[icent]));
        for(int ipt=0; ipt<npt; ipt++)  {
            if(in.eof()) break;
            in >> ybase[ipt] >> yerrbase[ipt];
        }
        in.close();
        
        cout << " TEST HERE !! " <<  " icent = " << icent << endl;

        out.open(Form("data/yieldSys_%s.txt",nameCent1[icent]));
        for(int ipt=0; ipt<npt; ipt++) {
            float dcR_mea = 0;
            float dcR_sys = 0;
            float pt_mean = 0.5*(nptbin[ipt]+nptbin[ipt+1]);
            int NbinSum_dc = 0;
            for(int i=centLw_dc[icent]; i<centUp_dc[icent]; i++) {
                NbinSum_dc += NbinMeanBase_dc[i];
                dcR_mea += gdcRatio_mea[i]->Eval(pt_mean) * NbinMeanBase_dc[i];
                dcR_sys += gdcRatio_sys[i]->Eval(pt_mean) * NbinMeanBase_dc[i];
            }
            dcR_mea /= NbinSum_dc;
            dcR_sys /= NbinSum_dc;
            
            y[ipt] = ybase[ipt]*(1.-dcR_mea);
            yerr[ipt] = yerrbase[ipt]*(1.-dcR_mea);
            sys = fabs(dcR_mea-dcR_sys)/(1.-dcR_mea);
            
            //cout << dcR_mea << "\t" << dcR_sys << endl;
            out << y[ipt] << "\t" << yerr[ipt] << "\t" << sys << "\t" << dcR_mea << endl;
        }
        out.close();
    }
}

#include "../anaCuts.h"
void get_Rcp() {
    
    ifstream in;
    ofstream out;
    float y[ncent][npt], yerr[ncent][npt];
    float Rcp[ncent][npt], Rcperr[ncent][npt];
    for(int icent=0; icent<ncent; icent++) {
        in.open(Form("data/yield_%s.txt",nameCent1[icent]));
        for(int ipt=0; ipt<npt; ipt++) in >> y[icent][ipt] >> yerr[icent][ipt];
        in.close();
    }

    for(int icent=0; icent<ncent; icent++) {
        out.open(Form("data/Rcp1_%s.txt",nameCent1[icent]));
        for(int ipt=0; ipt<npt; ipt++) {
            Rcp[icent][ipt] = y[icent][ipt]/y[4][ipt];
            Rcperr[icent][ipt] = Rcp[icent][ipt] * sqrt( pow(yerr[icent][ipt]/y[icent][ipt],2) +  pow(yerr[4][ipt]/y[4][ipt],2) );
            // out << y[icent][ipt] << "\t" << yerr[icent][ipt] << endl;
            out << Rcp[icent][ipt] << "\t" << Rcperr[icent][ipt] << endl;
        }
        out.close();
    }

    for(int icent=0; icent<ncent; icent++) {
        out.open(Form("data/Rcp2_%s.txt",nameCent1[icent]));
        for(int ipt=0; ipt<npt; ipt++) {
            Rcp[icent][ipt] = y[icent][ipt]/y[6][ipt];
            Rcperr[icent][ipt] = Rcp[icent][ipt] * sqrt( pow(yerr[icent][ipt]/y[icent][ipt],2) +  pow(yerr[6][ipt]/y[6][ipt],2) );
            // out << y[icent][ipt] << "\t" << yerr[icent][ipt] << endl;
            out << Rcp[icent][ipt] << "\t" << Rcperr[icent][ipt] << endl;
        }
        out.close();
    }
}

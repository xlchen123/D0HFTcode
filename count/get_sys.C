#include "../anaCuts.h"
void get_sys() {
    ifstream in;
    ofstream out;
    float y[npt], yerr[npt], ybase[npt], yerrbase[npt];
    float sys;
    for(int icent=0; icent<ncent; icent++) {
        in.open(Form("../default/data/re_yield_%s.txt",nameCent1[icent]));
        for(int ipt=0; ipt<npt; ipt++) in >> ybase[ipt] >> yerrbase[ipt];
        in.close();
        
        in.open(Form("data/yield_%s.txt",nameCent1[icent]));
        for(int ipt=0; ipt<npt; ipt++) in >> y[ipt] >> yerr[ipt];
        in.close();
        
        out.open(Form("data/yieldSys_%s.txt",nameCent1[icent]));
        for(int ipt=0; ipt<npt; ipt++) {
            sys = fabs(y[ipt]-ybase[ipt])/ybase[ipt];
            out << ybase[ipt] << "\t" << sys << endl;
            //cout << y[ipt] << "\t" << ybase[ipt] << endl;
        }
        out.close();
    }
}

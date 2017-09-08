#include "../anaCuts.h"
#include "../myConst.h"
void write_Rcp_err() {
    char CMD[250];
    char dir[250];
    sprintf(dir,"data");
    sprintf(CMD,"[ -d %s ] || mkdir -p %s",dir,dir);
    // gSystem->Exec(CMD);
    
    float ptmean[npt], pterr[npt];
    for(int ipt=0; ipt<npt; ipt++) {
        ptmean[ipt] = 0.5*(nptbin[ipt]+nptbin[ipt+1]);
        pterr[ipt] = 0.5*(-nptbin[ipt]+nptbin[ipt+1]);
    }
    
    TFile* fPionEmd = new TFile("../pid_sys/piplus_sys.root");
    TF1* fpiTpc[ncent];
    TF1* fpiTpcRcp1[ncent];
    TF1* fpiTpcRcp2[ncent];
    for(int icent=0; icent<ncent; icent++) {
      fpiTpc[icent] = (TF1*)fPionEmd->Get(Form("fCombine3_%s",nameCent1[icent]));
      fpiTpc[icent] ->SetName(Form("fpiTpc_%s",nameCent1[icent]));
      fpiTpcRcp1[icent] = (TF1*)fPionEmd->Get(Form("fCombine3Rcp1_%s",nameCent1[icent]));
      fpiTpcRcp1[icent] ->SetName(Form("fpiTpcRcp1_%s",nameCent1[icent]));
      fpiTpcRcp2[icent] = (TF1*)fPionEmd->Get(Form("fCombine3Rcp2_%s",nameCent1[icent]));
      fpiTpcRcp2[icent] ->SetName(Form("fpiTpcRcp2_%s",nameCent1[icent]));
    }
    fPionEmd->Close();

    TFile* fKaonEmd = new TFile("../pid_sys/kaonminus_sys.root");
    TF1* fkTpc[ncent];
    TF1* fkTpcRcp1[ncent];
    TF1* fkTpcRcp2[ncent];
    for(int icent=0; icent<ncent; icent++) {
      fkTpc[icent] = (TF1*)fKaonEmd->Get(Form("fCombine3_%s",nameCent1[icent]));
      fkTpc[icent] ->SetName(Form("fkTpc_%s",nameCent1[icent]));
      fkTpcRcp1[icent] = (TF1*)fKaonEmd->Get(Form("fCombine3Rcp1_%s",nameCent1[icent]));
      fkTpcRcp1[icent] ->SetName(Form("fkTpcRcp1_%s",nameCent1[icent]));
      fkTpcRcp2[icent] = (TF1*)fKaonEmd->Get(Form("fCombine3Rcp2_%s",nameCent1[icent]));
      fkTpcRcp2[icent] ->SetName(Form("fkTpcRcp2_%s",nameCent1[icent]));
    }
    fKaonEmd->Close();
    
    TFile* filepipid = new TFile("../pid_sys/pion_PidEff_Ks_170822_use_forSys.root");
    TF1* fpipid;
    fpipid = (TF1*)filepipid->Get("fpionNsigTof_effratio");
    fpipid ->SetName(Form("fpipid"));//9.87102992642156396e-01
    filepipid->Close();

    TFile* filekaonpid = new TFile("../pid_sys/kaon_PidEff_phi_170822_use_forSys.root");
    TF1* fkpid;
    fkpid = (TF1*)filekaonpid->Get("fkaonNsigTof_effratio");
    fkpid ->SetName(Form("fkpid")); // 9.87039958588842303e-01
    filekaonpid->Close();

    //caculate err
    float y[ncent][npt], yerr[ncent][npt], ysys[ncent][npt];
    float y0sys[ncent][npt], y1sys[ncent][npt];
    float yRcp[ncent][npt], yRcp1sys[ncent][npt], yRcp2sys[ncent][npt];
    float vtxSys[ncent][npt];
    float tmp[ncent][npt];
    float tmp1[ncent][npt];
    float tmp2[ncent][npt];
    float tmp3[ncent][npt];
    ifstream in;
    ofstream out;
    TFile* fout = new TFile("D0_Rcp_Sys.root","RECREATE");

    float unit[npt];
    for(int ipt=0; ipt<npt; ipt++)  unit[ipt] = 1.;

    // for yield extraction and double count and Vtxsys
    for(int icent=0; icent<ncent; icent++) {
        //sys 2 
        //-- 2.1 count with side band
        in.open(Form("../count/data/yieldSys_%s.txt",nameCent1[icent]));
        if(in.eof()) { cout << "No cout sys error file!!!" << endl; exit(1);}
        for(int ipt=0; ipt<npt; ipt++) {
            if(in.eof()) break;
            in >> y[icent][ipt] >> tmp1[icent][ipt];
        }
        in.close();
        
        //-- 2.2 count with fit change fit range
        in.open(Form("../fitRange/data/yieldSys_%s.txt",nameCent1[icent]));
        if(in.eof()) { cout << "No change fit range sys error file!!!" << endl; exit(1);}
        for(int ipt=0; ipt<npt; ipt++) {
            if(in.eof()) break;
            in >> y[icent][ipt] >> tmp2[icent][ipt];
        }
        in.close();

        //-- 2.3 count with likesign bkg subtraction
        in.open(Form("../likeSign/data/yieldSys_%s.txt",nameCent1[icent]));
        if(in.eof()) { cout << "No like-sign sys error file!!!" << endl; exit(1);}
        for(int ipt=0; ipt<npt; ipt++) {
            if(in.eof()) break;
            in >> y[icent][ipt] >> tmp3[icent][ipt];
        }
        in.close();

        for(int ipt=0; ipt<npt; ipt++) {
          tmp[icent][ipt] = tmp1[icent][ipt] > tmp2[icent][ipt] ? tmp1[icent][ipt] : tmp2[icent][ipt];
          tmp[icent][ipt] = tmp[icent][ipt] > tmp3[icent][ipt] ? tmp[icent][ipt] : tmp3[icent][ipt];
          y0sys[icent][ipt] = tmp[icent][ipt];
        }

        //sys 5 -- double count
        in.open(Form("../DoubleCount/data/yieldSys_%s.txt",nameCent1[icent]));
        if(in.eof()) { cout << "No double count sys error file!!!" << endl; exit(1);}
        for(int ipt=0; ipt<npt; ipt++) {
            if(in.eof()) break;
            in >> y[icent][ipt] >> yerr[icent][ipt] >> tmp[icent][ipt] >> tmp1[icent][ipt];
            y1sys[icent][ipt] = tmp[icent][ipt];
        }
        in.close();

        // do the vertex reso. correction
        // vtx sys. error
        in.open(Form("../vtxCorr/data/vtxSys_%s.txt",nameCent1[icent]));
        if(in.eof()) { cout << "No vtx sys error file!!!" << endl; exit(1);}
        for(int ipt=0; ipt<npt; ipt++) {
            if(in.eof()) break;
            in >> vtxSys[icent][ipt];
            // yRcp1sys[icent][ipt] = sqrt(pow(tmp[icent][ipt],2)+pow(yRcp1sys[icent][ipt],2));
        }
        in.close();
    }
        
    
    //for Rcp1
    for(int icent=0; icent<ncent; icent++) {
        //init
        for(int ipt=0; ipt<npt; ipt++) {
            yRcp[icent][ipt] = 1;
            yRcp1sys[icent][ipt] = 0;
        }

        //statistics error
        in.open(Form("../default/data/Rcp1_%s.txt",nameCent1[icent]));
        if(in.eof()) { cout << "No default yield error file!!!" << endl; exit(1);}
        for(int ipt=0; ipt<npt; ipt++) {
            if(in.eof()) cout << {
                break;
            }
            in >> tmp1[icent][ipt] >> tmp[icent][ipt];// the second is useless
            // yerr[ipt] /= y[ipt];
        }
        in.close();

        double TmpCom =  sqrt(pow(fpiTpcRcp1[icent]->Eval(1),2) + pow(1.-fpipid->Eval(1),2)) + sqrt(pow(fkTpcRcp1[icent]->Eval(1),2) + pow(1.-fkpid->Eval(1),2));
        cout<< "TmpCom = " << TmpCom << endl;

        //sys 1 -- tpc track
        for(int ipt=0; ipt<npt; ipt++) {
            yRcp1sys[icent][ipt] = sqrt(pow(TmpCom,2)+pow(yRcp1sys[icent][ipt],2));
        }

        //sys 2 -- yield extaction 
        for(int ipt=0; ipt<npt; ipt++) {
            yRcp1sys[icent][ipt] = sqrt(pow(y0sys[icent][ipt],2) + pow(y0sys[4][ipt],2) +pow(yRcp1sys[icent][ipt],2));
        }

        //sys 3, Daughter pt Cut scan // choose the maximum difference
        //sys 3.1 -- daughter pt cut1
        in.open(Form("../ptCut1/data/Rcp1_%s.txt",nameCent1[icent]));
        if(in.eof()) { cout << "No pt = 0.3 sys error file!!!" << endl; exit(1);}
        for(int ipt=0; ipt<npt; ipt++) {
            if(in.eof()) break;
            in >> tmp2[icent][ipt] >> tmp[icent][ipt];
        }
        in.close();

        //sys 3.2 -- daughter pt cut2
        in.open(Form("../ptCut2/data/Rcp1_%s.txt",nameCent1[icent]));
        if(in.eof()) { cout << "No pt = 0.5 sys error file!!!" << endl; exit(1);}
        for(int ipt=0; ipt<npt; ipt++) {
            if(in.eof()) break;
            in >> tmp3[icent][ipt] >> tmp[icent][ipt];
        }
        in.close();

        for(int ipt=0; ipt<npt; ipt++) {
          tmp[icent][ipt] = fabs(tmp1[icent][ipt] - tmp2[icent][ipt])>fabs(tmp1[icent][ipt] - tmp3[icent][ipt])? tmp2[icent][ipt] : tmp3[icent][ipt];
          yRcp1sys[icent][ipt] = sqrt(pow(tmp[icent][ipt]/tmp1[icent][ipt]-1.,2)+pow(yRcp1sys[icent][ipt],2));
        }

        //sys 4 topological cut
        //4.1 -- tight topo cuts
        in.open(Form("../topoCut1/data/Rcp1_%s.txt",nameCent1[icent]));
        if(in.eof()) { cout << "No tight topo sys error file!!!" << endl; exit(1);}
        for(int ipt=0; ipt<npt; ipt++) {
            if(in.eof()) break;
            in >> tmp2[icent][ipt] >> tmp[icent][ipt];
            // if(icent == 4 && ipt >= npt-3) tmp2[icent][ipt] =0; // for the tight cut, the signal is not good, no need to include this sys for some certain pt bin
        }
        in.close();

        //4.2 -- loose topo cuts
        in.open(Form("../topoCut2/data/Rcp1_%s.txt",nameCent1[icent]));
        if(in.eof()) { cout << "No loose topo sys error file!!!" << endl; exit(1);}
        for(int ipt=0; ipt<npt; ipt++) {
            if(in.eof()) break;
            in >> tmp3[icent][ipt] >> tmp[icent][ipt];
        }
        in.close();

        for(int ipt=0; ipt<npt; ipt++) {

          // tmp[icent][ipt] = fabs(tmp1[icent][ipt] - tmp2[icent][ipt])>fabs(tmp1[icent][ipt] - tmp3[icent][ipt])? tmp2[icent][ipt] : tmp3[icent][ipt];
          tmp[icent][ipt] = (tmp2[icent][ipt] + tmp3[icent][ipt])/2.;
          yRcp1sys[icent][ipt] = sqrt(pow(tmp[icent][ipt]/tmp1[icent][ipt]-1.,2)+pow(yRcp1sys[icent][ipt],2));
        }

        //sys 5 -- double count
        for(int ipt=0; ipt<npt; ipt++) {
            yRcp1sys[icent][ipt] = sqrt(pow(y1sys[icent][ipt],2) + pow(y1sys[4][ipt],2) +pow(yRcp1sys[icent][ipt],2));
        }

        // do the vertex reso. correction
        // vtx sys. error
        for(int ipt=0; ipt<npt; ipt++) {
            // yRcp1sys[icent][ipt] = sqrt(pow(vtxSys[icent][ipt],2)+pow(vtxSys[4][ipt],2)+pow(yRcp1sys[icent][ipt],2));
            yRcp1sys[icent][ipt] = sqrt(pow(vtxSys[icent][ipt],2)+pow(yRcp1sys[icent][ipt],2));//seperate 60-80% vtx contribution as band
        }

        // vtx stat. error
        // in.open(Form("../vtxCorr/data/vtxStat_%s.txt",nameCent1[icent]));
        // if(in.eof()) { cout << "No vtx sys error file!!!" << endl; exit(1);}
        // for(int ipt=0; ipt<npt; ipt++) {
        //     if(in.eof()) break;
        //     in >> tmp[ipt];
        //     ysys[ipt] = sqrt(pow(tmp[ipt],2)+pow(ysys[ipt],2));
        // }
        // in.close();


        //out put
        out.open(Form("data/Rcp1_sys_%s.txt",nameCent1[icent]));
        out << "yield" << "\t\t" << "sys.relative" << endl;
        for(int ipt=0; ipt<npt; ipt++) {
            out << 1 << "\t" << yRcp1sys[icent][ipt] << endl;
        }
        out.close();

        fout->cd();
        TGraphErrors* gD0Rcpsys = new TGraphErrors(npt, ptmean, unit, 0, yRcp1sys[icent]);
        gD0Rcpsys->SetMarkerStyle(MARKERSTYLE[0]);
        gD0Rcpsys->SetMarkerSize(2.);
        gD0Rcpsys->SetMarkerColor(COLOR[icent]);
        gD0Rcpsys->SetTitle(nameCent[icent]);
        gD0Rcpsys->Write(Form("gD0Rcp1_sys_%s",nameCent1[icent]));

     }

        TGraphErrors* gD0Rcpsys_vtx = new TGraphErrors(npt, ptmean, unit, 0, vtxSys[4]);
        gD0Rcpsys_vtx->SetMarkerStyle(MARKERSTYLE[0]);
        gD0Rcpsys_vtx->SetMarkerSize(2.);
        gD0Rcpsys_vtx->SetMarkerColor(COLOR[4]);
        gD0Rcpsys_vtx->SetTitle(nameCent[4]);
        gD0Rcpsys_vtx->Write(Form("gD0Rcp1_sys_vtx_%s",nameCent1[4]));

    //for Rcp2
    for(int icent=0; icent<ncent; icent++) {
        //init
        for(int ipt=0; ipt<npt; ipt++) {
            yRcp[icent][ipt] = 1;
            yRcp2sys[icent][ipt] = 0;
        }

        //statistics error
        in.open(Form("../default/data/Rcp2_%s.txt",nameCent1[icent]));
        if(in.eof()) { cout << "No default yield error file!!!" << endl; exit(1);}
        for(int ipt=0; ipt<npt; ipt++) {
            if(in.eof()) cout << {
                break;
            }
            in >> tmp1[icent][ipt] >> tmp[icent][ipt];// the second is useless
            // yerr[ipt] /= y[ipt];
        }
        in.close();

        double TmpCom =  sqrt(pow(fpiTpcRcp2[icent]->Eval(1),2) + pow(1.-fpipid->Eval(1),2)) + sqrt(pow(fkTpcRcp2[icent]->Eval(1),2) + pow(1.-fkpid->Eval(1),2));
        cout<< "TmpCom = " << TmpCom << endl;

        //sys 1 -- tpc track
        for(int ipt=0; ipt<npt; ipt++) {
            yRcp2sys[icent][ipt] = sqrt(pow(TmpCom,2)+pow(yRcp2sys[icent][ipt],2));
        }

        //sys 2 -- yield extaction 
        for(int ipt=0; ipt<npt; ipt++) {
            yRcp2sys[icent][ipt] = sqrt(pow(y0sys[icent][ipt],2) + pow(y0sys[6][ipt],2) +pow(yRcp2sys[icent][ipt],2));
        }

        //sys 3, Daughter pt Cut scan // choose the maximum difference
        //sys 3.1 -- daughter pt cut1
        in.open(Form("../ptCut1/data/Rcp2_%s.txt",nameCent1[icent]));
        if(in.eof()) { cout << "No pt = 0.3 sys error file!!!" << endl; exit(1);}
        for(int ipt=0; ipt<npt; ipt++) {
            if(in.eof()) break;
            in >> tmp2[icent][ipt] >> tmp[icent][ipt];
        }
        in.close();

        //sys 3.2 -- daughter pt cut2
        in.open(Form("../ptCut2/data/Rcp2_%s.txt",nameCent1[icent]));
        if(in.eof()) { cout << "No pt = 0.5 sys error file!!!" << endl; exit(1);}
        for(int ipt=0; ipt<npt; ipt++) {
            if(in.eof()) break;
            in >> tmp3[icent][ipt] >> tmp[icent][ipt];
        }
        in.close();

        for(int ipt=0; ipt<npt; ipt++) {
          tmp[icent][ipt] = fabs(tmp1[icent][ipt] - tmp2[icent][ipt])>fabs(tmp1[icent][ipt] - tmp3[icent][ipt])? tmp2[icent][ipt] : tmp3[icent][ipt];
          yRcp2sys[icent][ipt] = sqrt(pow(tmp[icent][ipt]/tmp1[icent][ipt]-1.,2)+pow(yRcp2sys[icent][ipt],2));
        }

        //sys 4 topological cut
        //4.1 -- tight topo cuts
        in.open(Form("../topoCut1/data/Rcp2_%s.txt",nameCent1[icent]));
        if(in.eof()) { cout << "No tight topo sys error file!!!" << endl; exit(1);}
        for(int ipt=0; ipt<npt; ipt++) {
            if(in.eof()) break;
            in >> tmp2[icent][ipt] >> tmp[icent][ipt];
            // if(icent == 4 && ipt >= npt-3) tmp2[icent][ipt] =0; // for the tight cut, the signal is not good, no need to include this sys for some certain pt bin
        }
        in.close();

        //4.2 -- loose topo cuts
        in.open(Form("../topoCut2/data/Rcp2_%s.txt",nameCent1[icent]));
        if(in.eof()) { cout << "No loose topo sys error file!!!" << endl; exit(1);}
        for(int ipt=0; ipt<npt; ipt++) {
            if(in.eof()) break;
            in >> tmp3[icent][ipt] >> tmp[icent][ipt];
        }
        in.close();

        for(int ipt=0; ipt<npt; ipt++) {

          // tmp[icent][ipt] = fabs(tmp1[icent][ipt] - tmp2[icent][ipt])>fabs(tmp1[icent][ipt] - tmp3[icent][ipt])? tmp2[icent][ipt] : tmp3[icent][ipt];
          tmp[icent][ipt] = (tmp2[icent][ipt] + tmp3[icent][ipt])/2.;
          yRcp2sys[icent][ipt] = sqrt(pow(tmp[icent][ipt]/tmp1[icent][ipt]-1.,2)+pow(yRcp2sys[icent][ipt],2));
        }

        //sys 5 -- double count
        for(int ipt=0; ipt<npt; ipt++) {
            yRcp2sys[icent][ipt] = sqrt(pow(y1sys[icent][ipt],2) + pow(y1sys[6][ipt],2) +pow(yRcp2sys[icent][ipt],2));
        }

        // do the vertex reso. correction
        // vtx sys. error
        for(int ipt=0; ipt<npt; ipt++) {
            // yRcp2sys[icent][ipt] = sqrt(pow(vtxSys[icent][ipt],2)+pow(vtxSys[6][ipt],2)+pow(yRcp2sys[icent][ipt],2));
            yRcp2sys[icent][ipt] = sqrt(pow(vtxSys[icent][ipt],2)+pow(yRcp2sys[icent][ipt],2));//seperate 60-80% vtx contribution as band
        }
        // vtx stat. error
        // in.open(Form("../vtxCorr/data/vtxStat_%s.txt",nameCent1[icent]));
        // if(in.eof()) { cout << "No vtx sys error file!!!" << endl; exit(1);}
        // for(int ipt=0; ipt<npt; ipt++) {
        //     if(in.eof()) break;
        //     in >> tmp[ipt];
        //     ysys[ipt] = sqrt(pow(tmp[ipt],2)+pow(ysys[ipt],2));
        // }
        // in.close();


        //out put
        out.open(Form("data/Rcp2_sys_%s.txt",nameCent1[icent]));
        out << "yield" << "\t\t" << "sys.relative" << endl;
        for(int ipt=0; ipt<npt; ipt++) {
            out << 1 << "\t" << yRcp2sys[icent][ipt] << endl;
        }
        out.close();

        fout->cd();
        TGraphErrors* gD0Rcpsys = new TGraphErrors(npt, ptmean, unit, 0, yRcp2sys[icent]);
        gD0Rcpsys->SetMarkerStyle(MARKERSTYLE[0]);
        gD0Rcpsys->SetMarkerSize(2.);
        gD0Rcpsys->SetMarkerColor(COLOR[icent]);
        gD0Rcpsys->SetTitle(nameCent[icent]);
        gD0Rcpsys->Write(Form("gD0Rcp2_sys_%s",nameCent1[icent]));
     }

        TGraphErrors* gD0Rcpsys_vtx = new TGraphErrors(npt, ptmean, unit, 0, vtxSys[6]);
        gD0Rcpsys_vtx->SetMarkerStyle(MARKERSTYLE[0]);
        gD0Rcpsys_vtx->SetMarkerSize(2.);
        gD0Rcpsys_vtx->SetMarkerColor(COLOR[6]);
        gD0Rcpsys_vtx->SetTitle(nameCent[6]);
        gD0Rcpsys_vtx->Write(Form("gD0Rcp2_sys_vtx_%s",nameCent1[6]));

    fout->Close();
    
}

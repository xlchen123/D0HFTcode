#include "../myFunction.h"
#include "../anaCuts.h"
#include "../myConst.h"
void plot_Rcp_err()
{
   globalSetting();
   char CMD[250];
   char dir[250];
   char title[250];
    char name[250];

   sprintf(dir, "pic");
   sprintf(CMD, "[ -d %s ] || mkdir -p %s", dir, dir);
   gSystem->Exec(CMD);

   TLegend* legend;
   TH1F* h0;

   const int nerr = 6;
   const char nameErr[nerr][250] = {"TPC track", "raw yield extract", "single track p_{T}", "50%/150% topo. eff.", "double count", "vertex correction"};

   float ptmean[npt], pterr[npt];
   for (int ipt = 0; ipt < npt; ipt++)
   {
      ptmean[ipt] = 0.5 * (nptbin[ipt] + nptbin[ipt + 1]);
      pterr[ipt] = 0.5 * (-nptbin[ipt] + nptbin[ipt + 1]);
   }

   TFile* fPionEmd = new TFile("../pid_sys/piplus_sys.root");
   TF1* fpiTpc[ncent];
   TF1* fpiTpcRcp1[ncent];
   TF1* fpiTpcRcp2[ncent];
   for (int icent = 0; icent < ncent; icent++)
   {
      fpiTpc[icent] = (TF1*)fPionEmd->Get(Form("fCombine3_%s", nameCent1[icent]));
      fpiTpc[icent] ->SetName(Form("fpiTpc_%s", nameCent1[icent]));
      fpiTpcRcp1[icent] = (TF1*)fPionEmd->Get(Form("fCombine3Rcp1_%s", nameCent1[icent]));
      fpiTpcRcp1[icent] ->SetName(Form("fpiTpcRcp1_%s", nameCent1[icent]));
      fpiTpcRcp2[icent] = (TF1*)fPionEmd->Get(Form("fCombine3Rcp2_%s", nameCent1[icent]));
      fpiTpcRcp2[icent] ->SetName(Form("fpiTpcRcp2_%s", nameCent1[icent]));
   }
   fPionEmd->Close();

   TFile* fKaonEmd = new TFile("../pid_sys/kaonminus_sys.root");
   TF1* fkTpc[ncent];
   TF1* fkTpcRcp1[ncent];
   TF1* fkTpcRcp2[ncent];
   for (int icent = 0; icent < ncent; icent++)
   {
      fkTpc[icent] = (TF1*)fKaonEmd->Get(Form("fCombine3_%s", nameCent1[icent]));
      fkTpc[icent] ->SetName(Form("fkTpc_%s", nameCent1[icent]));
      fkTpcRcp1[icent] = (TF1*)fKaonEmd->Get(Form("fCombine3Rcp1_%s", nameCent1[icent]));
      fkTpcRcp1[icent] ->SetName(Form("fkTpcRcp1_%s", nameCent1[icent]));
      fkTpcRcp2[icent] = (TF1*)fKaonEmd->Get(Form("fCombine3Rcp2_%s", nameCent1[icent]));
      fkTpcRcp2[icent] ->SetName(Form("fkTpcRcp2_%s", nameCent1[icent]));
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
   float yRcp[ncent][npt];
   float yRcp1sys[ncent][nerr][npt];
   float yRcp2sys[ncent][nerr][npt];
   float vtxSys[ncent][npt];
   float tmp[ncent][npt];
   float tmp1[ncent][npt];
   float tmp2[ncent][npt];
   float tmp3[ncent][npt];
   ifstream in;
   ofstream out;

   float unit[npt];
   for (int ipt = 0; ipt < npt; ipt++)  unit[ipt] = 1.;

   // for yield extraction and double count
   for (int icent = 0; icent < ncent; icent++)
   {
      //sys 2
      //-- 2.1 count with side band
      in.open(Form("../count/data/yieldSys_%s.txt", nameCent1[icent]));
      if (in.eof())
      {
         cout << "No cout sys error file!!!" << endl;
         exit(1);
      }
      for (int ipt = 0; ipt < npt; ipt++)
      {
         if (in.eof()) break;
         in >> y[icent][ipt] >> tmp1[icent][ipt];
      }
      in.close();

      //-- 2.2 count with fit change fit range
      in.open(Form("../fitRange/data/yieldSys_%s.txt", nameCent1[icent]));
      if (in.eof())
      {
         cout << "No change fit range sys error file!!!" << endl;
         exit(1);
      }
      for (int ipt = 0; ipt < npt; ipt++)
      {
         if (in.eof()) break;
         in >> y[icent][ipt] >> tmp2[icent][ipt];
      }
      in.close();

      //-- 2.3 count with likesign bkg subtraction
      in.open(Form("../likeSign/data/yieldSys_%s.txt", nameCent1[icent]));
      if (in.eof())
      {
         cout << "No like-sign sys error file!!!" << endl;
         exit(1);
      }
      for (int ipt = 0; ipt < npt; ipt++)
      {
         if (in.eof()) break;
         in >> y[icent][ipt] >> tmp3[icent][ipt];
      }
      in.close();

      for (int ipt = 0; ipt < npt; ipt++)
      {
         tmp[icent][ipt] = tmp1[icent][ipt] > tmp2[icent][ipt] ? tmp1[icent][ipt] : tmp2[icent][ipt];
         tmp[icent][ipt] = tmp[icent][ipt] > tmp3[icent][ipt] ? tmp[icent][ipt] : tmp3[icent][ipt];
         y0sys[icent][ipt] = tmp[icent][ipt];
      }

      //sys 5 -- double count
      in.open(Form("../DoubleCount/data/yieldSys_%s.txt", nameCent1[icent]));
      if (in.eof())
      {
         cout << "No double count sys error file!!!" << endl;
         exit(1);
      }
      for (int ipt = 0; ipt < npt; ipt++)
      {
         if (in.eof()) break;
         in >> y[icent][ipt] >> yerr[icent][ipt] >> tmp[icent][ipt] >> tmp1[icent][ipt];
         y1sys[icent][ipt] = tmp[icent][ipt];
      }
      in.close();
      // do the vertex reso. correction
      // vtx sys. error
      in.open(Form("../vtxCorr/data/vtxSys_%s.txt", nameCent1[icent]));
      if (in.eof())
      {
         cout << "No vtx sys error file!!!" << endl;
         exit(1);
      }
      for (int ipt = 0; ipt < npt; ipt++)
      {
         if (in.eof()) break;
         in >> vtxSys[icent][ipt];
      }
      in.close();
   }


   //for Rcp1
   for (int icent = 0; icent < ncent; icent++)
   {
      int ierr = 0;
      //init
      for (int ipt = 0; ipt < npt; ipt++)
      {
         yRcp[icent][ipt] = 1;
      }

      //statistics error
      in.open(Form("../default/data/Rcp1_%s.txt", nameCent1[icent]));
      if (in.eof())
      {
         cout << "No default yield error file!!!" << endl;
         exit(1);
      }
      for (int ipt = 0; ipt < npt; ipt++)
      {
         if (in.eof()) cout <<
         {
            break;
         }
         in >> tmp1[icent][ipt] >> tmp[icent][ipt];// the second is useless
         // yerr[ipt] /= y[ipt];
      }
      in.close();

      double TmpCom =  sqrt(pow(fpiTpcRcp1[icent]->Eval(1), 2) + pow(1. - fpipid->Eval(1), 2)) + sqrt(pow(fkTpcRcp1[icent]->Eval(1), 2) + pow(1. - fkpid->Eval(1), 2));
      cout << "TmpCom = " << TmpCom << endl;

      //sys 1 -- tpc track
      for (int ipt = 0; ipt < npt; ipt++)
      {
         yRcp1sys[icent][ierr][ipt] = TmpCom;
      }
      ierr++;

      //sys 2 -- yield extaction
      for (int ipt = 0; ipt < npt; ipt++)
      {
         yRcp1sys[icent][ierr][ipt] = sqrt(pow(y0sys[icent][ipt], 2) + pow(y0sys[4][ipt], 2));
      }
      ierr++;

      //sys 3, Daughter pt Cut scan // choose the maximum difference
      //sys 3.1 -- daughter pt cut1
      in.open(Form("../ptCut1/data/Rcp1_%s.txt", nameCent1[icent]));
      if (in.eof())
      {
         cout << "No pt = 0.3 sys error file!!!" << endl;
         exit(1);
      }
      for (int ipt = 0; ipt < npt; ipt++)
      {
         if (in.eof()) break;
         in >> tmp2[icent][ipt] >> tmp[icent][ipt];
      }
      in.close();

      //sys 3.2 -- daughter pt cut2
      in.open(Form("../ptCut2/data/Rcp1_%s.txt", nameCent1[icent]));
      if (in.eof())
      {
         cout << "No pt = 0.5 sys error file!!!" << endl;
         exit(1);
      }
      for (int ipt = 0; ipt < npt; ipt++)
      {
         if (in.eof()) break;
         in >> tmp3[icent][ipt] >> tmp[icent][ipt];
      }
      in.close();

      for (int ipt = 0; ipt < npt; ipt++)
      {
         tmp[icent][ipt] = fabs(tmp1[icent][ipt] - tmp2[icent][ipt]) > fabs(tmp1[icent][ipt] - tmp3[icent][ipt]) ? tmp2[icent][ipt] : tmp3[icent][ipt];
         yRcp1sys[icent][ierr][ipt] = sqrt(pow(tmp[icent][ipt] / tmp1[icent][ipt] - 1., 2));
      }
      ierr++;

      //sys 4 topological cut
      //4.1 -- tight topo cuts
      in.open(Form("../topoCut1/data/Rcp1_%s.txt", nameCent1[icent]));
      if (in.eof())
      {
         cout << "No tight topo sys error file!!!" << endl;
         exit(1);
      }
      for (int ipt = 0; ipt < npt; ipt++)
      {
         if (in.eof()) break;
         in >> tmp2[icent][ipt] >> tmp[icent][ipt];
         // if(icent == 4 && ipt >= npt-3) tmp2[icent][ipt] =0; // for the tight cut, the signal is not good, no need to include this sys for some certain pt bin
      }
      in.close();

      //4.2 -- loose topo cuts
      in.open(Form("../topoCut2/data/Rcp1_%s.txt", nameCent1[icent]));
      if (in.eof())
      {
         cout << "No loose topo sys error file!!!" << endl;
         exit(1);
      }
      for (int ipt = 0; ipt < npt; ipt++)
      {
         if (in.eof()) break;
         in >> tmp3[icent][ipt] >> tmp[icent][ipt];
      }
      in.close();

      for (int ipt = 0; ipt < npt; ipt++)
      {
         // if(icent == 4 && ipt >= npt-3) // for the tight cut, the signal is not good, no need to include this sys for some certain pt bin
         // tmp[icent][ipt] = tmp3[icent][ipt];
         // else
         // tmp[icent][ipt] = fabs(tmp1[icent][ipt] - tmp2[icent][ipt])>fabs(tmp1[icent][ipt] - tmp3[icent][ipt])? tmp2[icent][ipt] : tmp3[icent][ipt];
         tmp[icent][ipt] = (tmp2[icent][ipt] + tmp3[icent][ipt]) / 2.;
         yRcp1sys[icent][ierr][ipt] = sqrt(pow(tmp[icent][ipt] / tmp1[icent][ipt] - 1., 2));
      }
      ierr++;

      //sys 5 -- double count
      for (int ipt = 0; ipt < npt; ipt++)
      {
         yRcp1sys[icent][ierr][ipt] = sqrt(pow(y1sys[icent][ipt], 2) + pow(y1sys[4][ipt], 2));
      }
      ierr++;

      // do the vertex reso. correction
      // vtx sys. error
      for (int ipt = 0; ipt < npt; ipt++)
      {
         yRcp1sys[icent][ierr][ipt] = sqrt(pow(vtxSys[icent][ipt], 2)+pow(vtxSys[4][ipt], 2));
      }
      ierr++;

      cout << "ierr = " << ierr << endl;

   }

   TGraph* gSys[ncent][nerr];
   for (int icent = 0; icent < ncent; icent++)
   {
      for (int ierr = 0; ierr < nerr; ierr++)
      {
         gSys[icent][ierr] = new TGraph(npt, ptmean, yRcp1sys[icent][ierr]);
      }
   }

   
   //plot
   TCanvas* c1 = new TCanvas("c1", "A Canvas",10,10,600,800);
   setPad(c1);
   gPad->SetLogy();
   for(int icent=0; icent<ncent; icent++) {
       float ymin = 1.01e-3;
       float ymax = 9.5e0;
       h0= new TH1F("","",1,0,nptbin[npt]);
       h0->Draw();
       h0->SetMinimum(ymin),
       h0->SetMaximum(ymax);
       setHisto(h0,"","p_{T} (GeV/c)", "Sys. error (#times 100%)");
       const float legy_lw = 0.93 - 0.04*nerr;
       legend = new TLegend(0.6,legy_lw,0.9,0.93);
       legend->SetFillStyle(0);
       legend->SetFillColor(10);
       legend->SetBorderSize(0);
       legend->SetTextSize(0.04);
       legend->SetTextFont(132);
       for(int ierr=0; ierr<nerr; ierr++) {
           legend->AddEntry(gSys[icent][ierr],nameErr[ierr],"p");
           gSys[icent][ierr]->Draw("psame");
           gSys[icent][ierr]->SetMarkerStyle(MARKERSTYLE[ierr]);
           gSys[icent][ierr]->SetMarkerColor(COLOR[ierr]);
           gSys[icent][ierr]->SetMarkerSize(1.5);
       }
       legend->Draw();
       drawLatex(0.18,0.9,Form("AuAu @200GeV  %s",nameCent[icent]),132,0.04,1);
       gPad->SetLogy();
       sprintf(name,"%s/sysErr_%s_Rcp1.pdf",dir,nameCent1[icent]);
       c1->SaveAs(name);
   }
   

   //for Rcp2
   for (int icent = 0; icent < ncent; icent++)
   {
      int ierr = 0;
      //init
      for (int ipt = 0; ipt < npt; ipt++)
      {
         yRcp[icent][ipt] = 1;
      }

      //statistics error
      in.open(Form("../default/data/Rcp2_%s.txt", nameCent1[icent]));
      if (in.eof())
      {
         cout << "No default yield error file!!!" << endl;
         exit(1);
      }
      for (int ipt = 0; ipt < npt; ipt++)
      {
         if (in.eof()) cout <<
         {
            break;
         }
         in >> tmp1[icent][ipt] >> tmp[icent][ipt];// the second is useless
         // yerr[ipt] /= y[ipt];
      }
      in.close();

      double TmpCom =  sqrt(pow(fpiTpcRcp2[icent]->Eval(1), 2) + pow(1. - fpipid->Eval(1), 2)) + sqrt(pow(fkTpcRcp2[icent]->Eval(1), 2) + pow(1. - fkpid->Eval(1), 2));
      cout << "TmpCom = " << TmpCom << endl;

      //sys 1 -- tpc track
      for (int ipt = 0; ipt < npt; ipt++)
      {
         yRcp2sys[icent][ierr][ipt] = TmpCom;
      }
      ierr++;

      //sys 2 -- yield extaction
      for (int ipt = 0; ipt < npt; ipt++)
      {
         yRcp2sys[icent][ierr][ipt] = sqrt(pow(y0sys[icent][ipt], 2) + pow(y0sys[6][ipt], 2));
      }
      ierr++;

      //sys 3, Daughter pt Cut scan // choose the maximum difference
      //sys 3.1 -- daughter pt cut1
      in.open(Form("../ptCut1/data/Rcp2_%s.txt", nameCent1[icent]));
      if (in.eof())
      {
         cout << "No pt = 0.3 sys error file!!!" << endl;
         exit(1);
      }
      for (int ipt = 0; ipt < npt; ipt++)
      {
         if (in.eof()) break;
         in >> tmp2[icent][ipt] >> tmp[icent][ipt];
      }
      in.close();

      //sys 3.2 -- daughter pt cut2
      in.open(Form("../ptCut2/data/Rcp2_%s.txt", nameCent1[icent]));
      if (in.eof())
      {
         cout << "No pt = 0.5 sys error file!!!" << endl;
         exit(1);
      }
      for (int ipt = 0; ipt < npt; ipt++)
      {
         if (in.eof()) break;
         in >> tmp3[icent][ipt] >> tmp[icent][ipt];
      }
      in.close();

      for (int ipt = 0; ipt < npt; ipt++)
      {
         tmp[icent][ipt] = fabs(tmp1[icent][ipt] - tmp2[icent][ipt]) > fabs(tmp1[icent][ipt] - tmp3[icent][ipt]) ? tmp2[icent][ipt] : tmp3[icent][ipt];
         yRcp2sys[icent][ierr][ipt] = sqrt(pow(tmp[icent][ipt] / tmp1[icent][ipt] - 1., 2));
      }
      ierr++;

      //sys 4 topological cut
      //4.1 -- tight topo cuts
      in.open(Form("../topoCut1/data/Rcp2_%s.txt", nameCent1[icent]));
      if (in.eof())
      {
         cout << "No tight topo sys error file!!!" << endl;
         exit(1);
      }
      for (int ipt = 0; ipt < npt; ipt++)
      {
         if (in.eof()) break;
         in >> tmp2[icent][ipt] >> tmp[icent][ipt];
         // if(icent == 4 && ipt >= npt-3) tmp2[icent][ipt] =0; // for the tight cut, the signal is not good, no need to include this sys for some certain pt bin
      }
      in.close();

      //4.2 -- loose topo cuts
      in.open(Form("../topoCut2/data/Rcp2_%s.txt", nameCent1[icent]));
      if (in.eof())
      {
         cout << "No loose topo sys error file!!!" << endl;
         exit(1);
      }
      for (int ipt = 0; ipt < npt; ipt++)
      {
         if (in.eof()) break;
         in >> tmp3[icent][ipt] >> tmp[icent][ipt];
      }
      in.close();

      for (int ipt = 0; ipt < npt; ipt++)
      {
         // tmp[icent][ipt] = fabs(tmp1[icent][ipt] - tmp2[icent][ipt])>fabs(tmp1[icent][ipt] - tmp3[icent][ipt])? tmp2[icent][ipt] : tmp3[icent][ipt];
         tmp[icent][ipt] = (tmp2[icent][ipt] + tmp3[icent][ipt]) / 2.;
         yRcp2sys[icent][ierr][ipt] = sqrt(pow(tmp[icent][ipt] / tmp1[icent][ipt] - 1., 2));
      }
      ierr++;

      //sys 5 -- double count
      for (int ipt = 0; ipt < npt; ipt++)
      {
         yRcp2sys[icent][ierr][ipt] = sqrt(pow(y1sys[icent][ipt], 2) + pow(y1sys[6][ipt], 2));
      }
      ierr++;

      // do the vertex reso. correction
      // vtx sys. error
      for (int ipt = 0; ipt < npt; ipt++)
      {
         yRcp2sys[icent][ierr][ipt] = sqrt(pow(vtxSys[icent][ipt], 2)+pow(vtxSys[6][ipt], 2));
      }
      in.close();
      ierr++;

      cout << "ierr = " << ierr << endl;

   }

   TGraph* gSys[ncent][nerr];
   for (int icent = 0; icent < ncent; icent++)
   {
      for (int ierr = 0; ierr < nerr; ierr++)
      {
         gSys[icent][ierr] = new TGraph(npt, ptmean, yRcp2sys[icent][ierr]);
      }
   }

   
   //plot
   TCanvas* c1 = new TCanvas("c1", "A Canvas",10,10,600,800);
   setPad(c1);
   gPad->SetLogy();
   for(int icent=0; icent<ncent; icent++) {
       float ymin = 1.01e-3;
       float ymax = 9.5e0;
       h0= new TH1F("","",1,0,nptbin[npt]);
       h0->Draw();
       h0->SetMinimum(ymin),
       h0->SetMaximum(ymax);
       setHisto(h0,"","p_{T} (GeV/c)", "Sys. error (#times 100%)");
       const float legy_lw = 0.93 - 0.04*nerr;
       legend = new TLegend(0.6,legy_lw,0.9,0.93);
       legend->SetFillStyle(0);
       legend->SetFillColor(10);
       legend->SetBorderSize(0);
       legend->SetTextSize(0.04);
       legend->SetTextFont(132);
       for(int ierr=0; ierr<nerr; ierr++) {
           legend->AddEntry(gSys[icent][ierr],nameErr[ierr],"p");
           gSys[icent][ierr]->Draw("psame");
           gSys[icent][ierr]->SetMarkerStyle(MARKERSTYLE[ierr]);
           gSys[icent][ierr]->SetMarkerColor(COLOR[ierr]);
           gSys[icent][ierr]->SetMarkerSize(1.5);
       }
       legend->Draw();
       drawLatex(0.18,0.9,Form("AuAu @200GeV  %s",nameCent[icent]),132,0.04,1);
       gPad->SetLogy();
       sprintf(name,"%s/sysErr_%s_Rcp2.pdf",dir,nameCent1[icent]);
       c1->SaveAs(name);
   }

}

#include "../myFunction.h"
#include "../myConst.h"
void plot_Rcp1_pTshift()
{
   globalSetting();
   char dir[250];
   char name[250];
   char title[250];
   //TString CMD = 0;
   char CMD[250];
   TLegend* legend;
   TH1F* h0;

   sprintf(dir, "pic");
   sprintf(CMD, "[ -d %s ] || mkdir -p %s", dir, dir);
   gSystem->Exec(CMD);

   //change the base line--last cent bin
   const int ncent = 5;
   const char nameCent[ncent][250] = {"0-10%", "10-20%", "20-40%", "40-60%", "60-80%"};
   const char nameCentXL[ncent][250] = {"0_10", "10_20", "20_40", "40_60", "60_80"};
   float NbinMean[ncent] = {938.80170, 579.89409, 288.35051, 91.37100, 21.37396};

   //Read spectra
   TGraphErrors* gD0err_xl[ncent];
   TGraphErrors* gD0sys_xl[ncent];
   TF1* fLevy[ncent];
   TFile* fin1 = new TFile("../ptShift/D0_Spectra_Run14HFT.root");
   for (int icent = 0; icent < ncent; icent++)
   {
      gD0err_xl[icent] = (TGraphErrors*)fin1->Get(Form("gD0_err_%s", nameCentXL[icent]));
      gD0sys_xl[icent] = (TGraphErrors*)fin1->Get(Form("gD0_sys_%s", nameCentXL[icent]));
      //fLevy[icent] = (TF1*)fin1->Get(Form("flevy_%s",nameCentXL[icent]));
   }
   fin1->Close();

   TFile* fin2 = new TFile("../sys/D0_Rcp_Sys.root");
   TGraphErrors* gD0Rcp1_sys[ncent];
   TGraphErrors* gD0Rcp2_sys[ncent];
   for (int icent = 0; icent < ncent; icent++)
   {
      gD0Rcp1_sys[icent] = (TGraphErrors*)fin2->Get(Form("gD0Rcp1_sys_%s", nameCentXL[icent]));
      gD0Rcp2_sys[icent] = (TGraphErrors*)fin2->Get(Form("gD0Rcp2_sys_%s", nameCentXL[icent]));
   }
   TGraphErrors* gD0Rcp1_sys_vtx = (TGraphErrors*)fin2->Get(Form("gD0Rcp1_sys_vtx_60_80"));
   TGraphErrors* gD0Rcp2_sys_vtx = (TGraphErrors*)fin2->Get(Form("gD0Rcp2_sys_vtx_40_80"));
   fin2->Close();

   //calculate Rcp
   for (int icent = 0; icent < ncent; icent++)
   {
      cout << nameCent[icent] << endl;
      for (int ipt = 0; ipt < gD0err_xl[0]->GetN(); ipt++)
      {
         float y = gD0err_xl[icent]->GetY()[ipt] / NbinMean[icent];
         float yErr = gD0err_xl[icent]->GetEY()[ipt] / NbinMean[icent];
         float ySys = gD0sys_xl[icent]->GetEY()[ipt] / NbinMean[icent];
         float base = gD0err_xl[ncent - 1]->GetY()[ipt] / NbinMean[ncent - 1]; //60-80%
         float baseErr = gD0err_xl[ncent - 1]->GetEY()[ipt] / NbinMean[ncent - 1];
         float baseSys = gD0sys_xl[ncent - 1]->GetEY()[ipt] / NbinMean[ncent - 1];
         float Rcp = y / base;
         float RcpErr = Rcp * sqrt(pow(yErr / y, 2) + pow(baseErr / base, 2));
         // float RcpSys = Rcp * sqrt(pow(ySys / y, 2) + pow(baseSys / base, 2));
         float RcpSys = Rcp * gD0Rcp1_sys[icent]->GetEY()[ipt];
         gD0err_xl[icent]->GetY()[ipt] = Rcp;
         gD0err_xl[icent]->GetEY()[ipt] = RcpErr;
         gD0sys_xl[icent]->GetY()[ipt] = Rcp;
         gD0sys_xl[icent]->GetEY()[ipt] = RcpSys;
         cout << gD0err_xl[icent]->GetX()[ipt] << "\t" << Rcp << "\t" << RcpErr << "\t" << RcpSys << endl;
      }
      cout << endl;
   }

   //set for plot
   float markerSize = 2.0;
   float lineWidth = 2;
   for (int icent = 0; icent < ncent; icent++)
   {
      gD0err_xl[icent]->SetMarkerStyle(kFullCircle);
      gD0err_xl[icent]->SetMarkerSize(markerSize);
      gD0err_xl[icent]->SetMarkerColor(COLOR[icent]);
      gD0err_xl[icent]->SetLineWidth(lineWidth);
      gD0err_xl[icent]->SetLineColor(COLOR[icent]);

      gD0sys_xl[icent]->SetMarkerStyle(kFullCircle);
      gD0sys_xl[icent]->SetMarkerSize(markerSize);
      gD0sys_xl[icent]->SetMarkerColor(COLOR[icent]);
      gD0sys_xl[icent]->SetLineWidth(lineWidth);
      gD0sys_xl[icent]->SetLineColor(COLOR[icent]);

      //fLevy[icent]->SetLineColor(COLOR[icent]);
      //fLevy[icent]->SetLineWidth(lineWidth);
   }

   //plot Rcp
   TCanvas* c1 = new TCanvas("c1", "A Canvas", 10, 10, 800, 800);
   setPad(c1);
   float ymin = 0;
   float ymax = 1.5;
   h0 = new TH1F("", "", 1, 0, 8);
   h0->Draw();
   h0->SetMinimum(ymin),
   h0->SetMaximum(ymax);
   setHisto(h0, "", "p_{T} (GeV/c)", Form("R_{cp} (/%s)", nameCent[ncent - 1]));
   const float legy_lw = 0.93 - 0.05 * (ncent - 1);
   legend = new TLegend(0.7, legy_lw, 0.9, 0.93);
   legend->SetFillStyle(0);
   legend->SetFillColor(10);
   legend->SetBorderSize(0);
   legend->SetTextSize(0.04);
   legend->SetTextFont(132);
   //legend->SetHeader("Guannan");
   gD0Rcp1_sys_vtx->Draw("e3");
   gD0Rcp1_sys_vtx->SetLineColor(16);
   gD0Rcp1_sys_vtx->SetFillColor(16);
   const float sysw = 0.15;
   for (int icent = 0; icent < ncent - 1; icent++)
   {
      legend->AddEntry(gD0err_xl[icent], nameCent[icent], "p");
      gD0err_xl[icent]->Draw("psame");
      // gD0sys_xl[icent]->Draw("psame");
      for (int i = 0; i < gD0sys_xl[icent]->GetN(); i++)
      {
         const float sysl = gD0sys_xl[icent]->GetY()[i] * 0.05;
         TLine *llw = new TLine(gD0sys_xl[icent]->GetX()[i] - sysw, gD0sys_xl[icent]->GetY()[i] - gD0sys_xl[icent]->GetEY()[i], gD0sys_xl[icent]->GetX()[i] + sysw, gD0sys_xl[icent]->GetY()[i] - gD0sys_xl[icent]->GetEY()[i]);
         llw->SetLineWidth(2);
         llw->SetLineColor(COLOR[icent]);
         llw->Draw("same");
         TLine *lhi = new TLine(gD0sys_xl[icent]->GetX()[i] - sysw, gD0sys_xl[icent]->GetY()[i] + gD0sys_xl[icent]->GetEY()[i], gD0sys_xl[icent]->GetX()[i] + sysw, gD0sys_xl[icent]->GetY()[i] + gD0sys_xl[icent]->GetEY()[i]);
         lhi->SetLineWidth(2);
         lhi->SetLineColor(COLOR[icent]);
         lhi->Draw("same");
         TLine *lv = new TLine(gD0sys_xl[icent]->GetX()[i] - sysw, gD0sys_xl[icent]->GetY()[i] - gD0sys_xl[icent]->GetEY()[i], gD0sys_xl[icent]->GetX()[i] - sysw, gD0sys_xl[icent]->GetY()[i] - gD0sys_xl[icent]->GetEY()[i] + sysl);
         lv->SetLineWidth(2);
         lv->SetLineColor(COLOR[icent]);
         lv->Draw("same");
         TLine *lv = new TLine(gD0sys_xl[icent]->GetX()[i] + sysw, gD0sys_xl[icent]->GetY()[i] - gD0sys_xl[icent]->GetEY()[i], gD0sys_xl[icent]->GetX()[i] + sysw, gD0sys_xl[icent]->GetY()[i] - gD0sys_xl[icent]->GetEY()[i] + sysl);
         lv->SetLineWidth(2);
         lv->SetLineColor(COLOR[icent]);
         lv->Draw("same");
         TLine *lv = new TLine(gD0sys_xl[icent]->GetX()[i] - sysw, gD0sys_xl[icent]->GetY()[i] + gD0sys_xl[icent]->GetEY()[i], gD0sys_xl[icent]->GetX()[i] - sysw, gD0sys_xl[icent]->GetY()[i] + gD0sys_xl[icent]->GetEY()[i] - sysl);
         lv->SetLineWidth(2);
         lv->SetLineColor(COLOR[icent]);
         lv->Draw("same");
         TLine *lv = new TLine(gD0sys_xl[icent]->GetX()[i] + sysw, gD0sys_xl[icent]->GetY()[i] + gD0sys_xl[icent]->GetEY()[i], gD0sys_xl[icent]->GetX()[i] + sysw, gD0sys_xl[icent]->GetY()[i] + gD0sys_xl[icent]->GetEY()[i] - sysl);
         lv->SetLineWidth(2);
         lv->SetLineColor(COLOR[icent]);
         lv->Draw("same");
      }
   }
   legend->Draw("same");


   sprintf(name, "%s/D0_Rcp1_pTshift.pdf", dir);
   c1->SaveAs(name);

   // ==== Write
   TFile* fout = new TFile("D0_Rcp1_pTshift.root", "RECREATE");
   fout->cd();
   for (int icent = 0; icent < ncent; icent++)
   {
      gD0err_xl[icent]->Write(Form("gD0_Rcp_err_%s", nameCentXL[icent]));
      gD0sys_xl[icent]->Write(Form("gD0_Rcp_sys_%s", nameCentXL[icent]));
   }
   gD0Rcp1_sys_vtx->Write("gD0Rcp1_sys_vtx_60_80");
}

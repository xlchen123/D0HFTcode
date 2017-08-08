//===============================
//Version: 1.0 
//Time: Sun Nov  6 18:45:20 PST 2016 
//Author: Long Zhou 
//Discribe: Compare D0 spectra from each mthod  

#include <stdio.h>
#include <iostream>
#include <TH2.h>
#include <TStyle.h>
#include "TDatime.h"
#include "TSystem.h"
#include "TPad.h"
#include "TVirtualPad.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"
#include "TString.h"
#include "TLatex.h"
#include "TPDF.h"
#include "TGaxis.h"
float findGraphError(TGraphErrors *gr, float pos);
Int_t decomposeTGraphErrors(TGraphErrors *grAll, TGraph *gr_m, TGraph *gr_u,TGraph *gr_d);
TLatex *drawLatex(Double_t x, Double_t y, char *text, Int_t textFont, Double_t textSize, Int_t colorIndex);
TH1D* graph2hist(TGraphAsymmErrors *gr); // when the asymmerrors is a symmerrors
TGraphErrors *scaleGraph(TGraphErrors *gr1,TGraphErrors *gr2,float scale1,float scale2);
void scaleGraph(TGraphErrors *gr, float scale);
TGraphErrors *scaleGraph(TGraphErrors *gr, TF1 *f, float x1=1, float x2=1);
TLine* drawLine(Double_t xlow,Double_t ylow, Double_t xup, Double_t yup, Double_t lineWidth, Int_t lineColor,Int_t linestyle);
void ptShift()
{
  gStyle->SetOptFit(00000);
  TFile *ds_file_run14_lz  = new TFile("Ds_spectra.root");
  TFile *fout = new TFile("Ds_spectra_shifted.root","RECREATE");

  const Int_t ncenbins = 4;
  const Int_t nloops   = 10;
  const Int_t nptbins  = 4;

  TString centrality[ncenbins] = {"0_80","0_10","10_40","40_80"};
  TString cenname = "";
  TString outname = "Ds_spectra_shifted";
  Double_t ptbins[nptbins][2]={{1.5,2.5},{2.5,3.5},{3.5,5.0},{5.0,8.0}};

  TGraphErrors *ds_spectra_shifted[ncenbins];
  TGraphErrors *gtemp;
  char funcString[200];
  char funcString_time_pt[200];
  double m0 = 1.9682;//D0-1.8645, D+/- - 1.8693;
  sprintf(funcString,"(1/(2*TMath::Pi()))*[0]*([2]-1)*([2]-2)/([2]*[1]*([2]*[1]+%f*([2]-2)))*TMath::Power(([2]*[1]+TMath::Sqrt(x[0]*x[0]+%f*%f)-%f)/([2]*[1]),-[2])",m0,m0,m0,m0); // [0] - dsig/dy [2] - n, [1] - C,
  sprintf(funcString_time_pt,"(1/(2*TMath::Pi()))*[0]*([2]-1)*([2]-2)/([2]*[1]*([2]*[1]+%f*([2]-2)))*TMath::Power(([2]*[1]+TMath::Sqrt(x[0]*x[0]+%f*%f)-%f)/([2]*[1]),-[2])*x[0]",m0,m0,m0,m0); // [0] - dsig/dy [2] - n, [1] - C,
  
  TF1 *flevy = new TF1("flevy",funcString,0,20);
  flevy->SetParameters(6.25342e-01,3.34550e-01,1.86042e+01);
  flevy->SetLineColor(kGreen+2);
  // flevy->SetLineStyle(9);

  TF1 *flevy_time_pt = new TF1("flevy_time_pt",funcString_time_pt,0,20);

  //canvas only needed for this documentation
  TCanvas* c1 = new TCanvas("example","",800,800);
  c1->SetFillStyle(1001);
  c1->SetFillColor(kWhite);
  c1->SetGrid(0,0);

  cenname = centrality[0];
  TGraphErrors *ds_spectra  = (TGraphErrors *)ds_file_run14_lz->Get(Form("Ds_spectra_%s",cenname.Data()));
  ds_spectra->SetName(Form("Ds_spectra_%s",cenname.Data()));
  // ds_spectra->SetMarkerStyle(33);
  // ds_spectra_shifted->SetName(Form("Ds_spectra_shift_%s",cenname.Data()));
  flevy->SetName(Form("Levy_%s",cenname.Data()));
  
  // Plots 1 : Draw the spectra without any shift 
  // ------------------------------------------------
  TH1F *frame = (TH1F *)gPad->DrawFrame(0,1e-7,10,1);
  frame->SetTitle("Ds spectra;p_{T} [GeV/c];d^{2}N/(N_{ev}2#pip_{T}dp_{T}dy)");
  gPad->SetLogy(1);
  Int_t styleindex[6] = {24,25,26,29,27,28};
  Int_t colorindex[10] = {kRed-2,kOrange-2,kSpring-2,kGreen-2,kTeal-2,kCyan-2,kBlue-2,kViolet-2,kMagenta-2,kBlack};
  
  ds_spectra->SetLineColor(2);
  ds_spectra->SetMarkerSize(2);
  ds_spectra->SetMarkerStyle(20);
  ds_spectra->SetMarkerColor(2);
  ds_spectra->DrawClone("psame");
  ds_spectra->Fit(flevy,"NQ");
  flevy->DrawClone("same");
  SaveToPDF(c1,outname,0);

  flevy->SetName(Form("Levy_%s",cenname.Data()));
  fout->WriteTObject(flevy);	  // 0-80%, the final fit function
  // flevy->SetNpx();
  
  // plots 2 : Draw all shifts
  // --------------------------
  gPad->DrawFrame(0,1e-7,10,1,"Ds spectra;p_{T} [GeV/c];d^{2}N/(N_{ev}2#pip_{T}dp_{T}dy)");
  gtemp = (TGraphErrors *)ds_spectra->Clone();
  gtemp->SetMarkerStyle(33);
  ds_spectra->Draw("psame");
  Int_t NPoint = ds_spectra->GetN();
  for(int n=0;n<nloops;n++)
    {
      flevy->SetParameters(6.25342e-01,3.34550e-01,1.86042e+01);
      gtemp->Fit(flevy,"NQ");
      // flevy->Print();
      flevy_time_pt->SetParameters(flevy->GetParameters());
      // flevy_time_pt->Print();

      for(int j=0;j<NPoint;j++)
	{
	  Float_t x = gtemp->GetX()[j];
	  Float_t y = gtemp->GetY()[j];
	  Float_t ey = gtemp->GetEY()[j];
	  Float_t x_corr = flevy_time_pt->Integral(ptbins[j][0],ptbins[j][1])/flevy->Integral(ptbins[j][0],ptbins[j][1]);
	  gtemp->SetPoint(j,x_corr,y);
	  gtemp->SetPointError(j,0,ey);
	}
      gtemp->SetMarkerColor(colorindex[n]);
      gtemp->SetFillColor(colorindex[n]);
      flevy->SetLineColor(colorindex[n]);
      gtemp->DrawClone("psame");
      flevy->DrawClone("same");
    }
  SaveToPDF(c1,outname);

  // ==== Create a TGraphErrors to hold the confidence intervals
  const Int_t NCL = 200;
  TGraphErrors *grint;
  grint = new TGraphErrors(NCL);
  grint->SetTitle("Fitted line with .68 conf. band");
  for (int i=0; i<NCL; i++) grint->SetPoint(i, 0 + (float)i/20., 0);

  (TVirtualFitter::GetFitter())->GetConfidenceIntervals(grint,0.68);
  grint->SetName(Form("flevy_err_band_%s",cenname.Data()));
  
  // Plots 3 : Draw before and after: 
  // --------------------------------
  flevy->SetName(Form("Levy_shifted_%s",cenname.Data()));
  fout->WriteTObject(flevy);	  // 0-80%, the final fit function
  gPad->DrawFrame(0,1e-7,10,1,"Ds spectra;p_{T} [GeV/c];d^{2}N/(N_{ev}2#pip_{T}dp_{T}dy)");
  gtemp->DrawClone("psame");
  flevy->DrawClone("same");
  ds_spectra->DrawClone("psame");
  flevy->SetParameters(6.25342e-01,3.34550e-01,1.86042e+01);
  ds_spectra->Fit(flevy,"NQ");
  flevy->SetLineColor(kGray);
  flevy->DrawClone("same");
  
  TLegend *leg = new TLegend(0.60,0.70,0.80,0.80);
  leg->SetFillColor(10);
  leg->SetTextFont(62);
  leg->SetTextSize(0.04);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(ds_spectra,"w/o shift","pl");
  leg->AddEntry(gtemp,"w shift","pl");
  leg->DrawClone();
  SaveToPDF(c1,outname,2);
  gtemp->SetName(Form("Ds_spectra_shifted_%s",cenname.Data()));

  fout->WriteTObject(ds_spectra); // 0-80%-raw
  fout->WriteTObject(gtemp);	  // 0-80%-shifted
  fout->WriteTObject(grint);	  // 0-80%-fit error band 
 
  // Plots 4 : Shift the other 3 centrality bins
  // --------------------------------------------
  for(int ic=1;ic<ncenbins;ic++)
    {
      cenname = centrality[ic];
      ds_spectra_shifted[ic] = (TGraphErrors *)ds_file_run14_lz->Get(Form("Ds_spectra_%s",cenname.Data()));
      ds_spectra_shifted[ic]->SetName(Form("Ds_spectra_%s",cenname.Data()));
      flevy->SetName(Form("Levy_shifted_%s",cenname.Data()));
      fout->WriteTObject(ds_spectra_shifted[ic]); // raw
 
      Int_t NPoint = ds_spectra->GetN();
      for(int ip=0;ip<NPoint;ip++)
	{
	  Double_t pt        = gtemp->GetX()[ip];
	  Double_t yield     = ds_spectra_shifted[ic]->GetY()[ip];
	  Double_t err_yield = ds_spectra_shifted[ic]->GetEY()[ip];
	  ds_spectra_shifted[ic]->SetPoint(ip,pt,yield);
	  ds_spectra_shifted[ic]->SetPointError(ip,0,err_yield);
	}
      ds_spectra_shifted[ic]->SetName(Form("Ds_spectra_shifted_%s",cenname.Data()));

      for(int ifit=0;ifit<3;ifit++) ds_spectra_shifted[ic]->Fit(flevy,"NQ");
      grint = new TGraphErrors(NCL);
      grint->SetTitle("Fitted line with .68 conf. band");
      for (int i=0; i<NCL; i++) grint->SetPoint(i, 0 + (float)i/20., 0);
      (TVirtualFitter::GetFitter())->GetConfidenceIntervals(grint,0.68);
      grint->SetName(Form("flevy_err_band_%s",cenname.Data()));

      fout->WriteTObject(ds_spectra_shifted[ic]); // shifted 
      fout->WriteTObject(flevy);
      fout->WriteTObject(grint);
    }
  return;
}

TGraphErrors* scaleGraph(TGraphErrors *gr1,TGraphErrors *gr2,float scale1,float scale2)
{
  Int_t NPoint = gr1->GetN();
  TGraphErrors *gr = new TGraphErrors(NPoint);
  gr->SetMarkerStyle(gr1->GetMarkerStyle());
  gr->SetMarkerSize(gr1->GetMarkerSize());
  gr->SetMarkerColor(gr1->GetMarkerColor());
  gr->SetLineColor(gr1->GetLineColor());
  for(int i=0;i<NPoint;i++)
    {
      float x = gr1->GetX()[i];
      float y = scale1*gr1->GetY()[i];      
      float ex = gr1->GetEX()[i];
      float ey = scale1*gr1->GetEY()[i];
      float y2 = scale2*gr2->Eval(x);      
      float ey2 = findGraphError(gr2,x);
      float err = y/y2 * sqrt(pow(ey/y,2)+pow(ey2/y2,2));

      gr->SetPoint(i,x,y/y2);
      gr->SetPointError(i,ex,err);
    }
  return gr;
}

float findGraphError(TGraphErrors *gr, float pos)
{
  Int_t NPoint = gr->GetN();
  TGraph *gr_m = new TGraph(NPoint);
  TGraph *gr_u = new TGraph(NPoint);
  TGraph *gr_d = new TGraph(NPoint);
  decomposeTGraphErrors(gr, gr_m, gr_u, gr_d);
  return (gr_u->Eval(pos) - gr_d->Eval(pos))/2.0;
}

TH1D* graph2hist(TGraphAsymmErrors *gr)
{
  Int_t NPoint = gr->GetN();
  float xmin = gr->GetX()[0];
  float ex   = gr->GetEXlow()[0];
  float xmax = gr->GetX()[NPoint-1];

  TH1D* hgraph = new TH1D("hgraph","",NPoint,xmin-ex,xmax+ex);
  for (int i = 0; i < gr->GetN(); i++) {
    float x = gr->GetX()[i];
    float y = gr->GetY()[i];
    float ey = gr->GetEYhigh()[i];
    int bin = hgraph->GetXaxis()->FindBin(x);
    hgraph->SetBinContent(bin,y);
    hgraph->SetBinError(bin,ey);
  }
  
  return hgraph;
}

TLatex *drawLatex(Double_t x, Double_t y, char *text, Int_t textFont, Double_t textSize, Int_t colorIndex)
{
  TLatex *latex = new TLatex(x, y, text);
  latex->SetNDC();
  latex->SetTextFont(textFont);
  latex->SetTextSize(textSize);
  latex->SetTextColor(colorIndex);
  latex->Draw("same");
  return latex;
}

void scaleGraph(TGraphErrors *gr, float scale)
{
  Int_t NPoint = gr->GetN();
  for (int i = 0; i < NPoint; i++) {
    float x = gr->GetX()[i];
    float y = gr->GetY()[i];
    float ex = gr->GetEX()[i];
    float ey = gr->GetEY()[i];
    y = scale*y;
    gr->SetPoint(i, x, y);
    gr->SetPointError(i, ex, ey*scale);
  }
}

TGraphErrors* scaleGraph(TGraphErrors *gr1,TF1 *gr2,float x1,float x2)
{
  Int_t NPoint = gr1->GetN();
  TGraphErrors *gr = new TGraphErrors(NPoint);
  gr->SetMarkerStyle(gr1->GetMarkerStyle());
  gr->SetMarkerSize(gr1->GetMarkerSize());
  gr->SetMarkerColor(gr1->GetMarkerColor());
  gr->SetLineColor(gr1->GetLineColor());
  for(int i=0;i<NPoint;i++)
    {
      float x = gr1->GetX()[i];
      float y = x1*gr1->GetY()[i];      
      float ex = gr1->GetEX()[i];
      float ey = x1*gr1->GetEY()[i];
      float y2 = x2*gr2->Eval(x);      
      float err = x1/y2*ey;
      // float ey2 = findGraphError(gr2,x);
      // float err = y/y2 * sqrt(pow(ey/y,2)+pow(ey2/y2,2));

      gr->SetPoint(i,x,y/y2);
      gr->SetPointError(i,ex,err);
    }
  return gr;
}

Int_t decomposeTGraphErrors(TGraphErrors *grAll, TGraph *gr_m, TGraph *gr_u,TGraph *gr_d)
{
  // cout<< " Decompose TGraphErrors to 3 TGraphs "<<endl;
  Int_t N1 = gr_m->GetN();
  Int_t N2 = gr_u->GetN();
  Int_t N3 = gr_d->GetN();
  Int_t NPoint = grAll->GetN();

  if (((N1-NPoint) !=0) || ((N1-N2) !=0) || ((N2-N3) != 0) ) {cout << "The number of point in final graph is not correct , please check it !" << endl; return 0;}

  for (int i = 0; i < NPoint; i++) {
    float x   = grAll->GetX()[i];
    float y   = grAll->GetY()[i];
    float ex  = grAll->GetEX()[i];
    float ey  = grAll->GetEY()[i];
    gr_m->SetPoint(i, x, y);
    gr_u->SetPoint(i, x, y + ey);
    gr_d->SetPoint(i, x, y - ey);
  }
  return 1;
}

TLine* drawLine(Double_t xlow,Double_t ylow, Double_t xup, Double_t yup, Double_t lineWidth, Int_t lineColor,Int_t linestyle){
  TLine *l1 = new TLine(xlow,ylow,xup,yup);
  l1->SetLineWidth(lineWidth);
  l1->SetLineColor(lineColor);
  l1->SetLineStyle(linestyle);
  l1->Draw("same");
  return l1;
}


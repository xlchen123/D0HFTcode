#include <stdio>
#include <stdlib>
#include <iostream>
#include <fstream>
#include <math.h>
#include "iomanip.h"

#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "THn.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TStyle.h"

TH1D* histo(char *name, Double_t xlow, Double_t xup, Double_t ylow, Double_t yup, char* xTitle, char* yTitle);
TLatex* drawLatex(Double_t x, Double_t y, char* text, Int_t textFont, Double_t textSize, Int_t colorIndex);
TLine* drawLine(Double_t xlow, Double_t ylow, Double_t xup, Double_t yup, Double_t lineWidth, Int_t lineColor);
void drawLines(Double_t xlow, Double_t ylow, Double_t xup, Double_t yup, Double_t lineWidth, Int_t lineColor);
void setpad(TPad *pad, float left, float right, float top, float bottom);

Double_t line(Double_t *x, Double_t *par)
{
   return par[0] +  par[1] * x[0];
}
Double_t line2(Double_t *x, Double_t *par)
{
   return par[0] + par[1] * x[0]  + par[2] * x[0] * x[0];
}

Double_t expol(Double_t *x, Double_t *par)
{
   return exp(par[0] + x[0] * par[1]);
}

std::pair<float, float> gMassCorrectionFitRange(1.71, 2.02);
std::pair<float, float> gMassExclusionRange(1.80, 1.92);
std::pair<float, float> gMassYieldCountRange(1.82, 1.91);

std::pair<float, float> gMassYieldSideBand1(1.71, 1.80);
std::pair<float, float> gMassYieldSideBand2(1.93, 2.02);

Double_t residualBg(Double_t *x, Double_t *par)
{
   if (x[0] > gMassExclusionRange.first && x[0] < gMassExclusionRange.second)
   {
      TF1::RejectPoint();
      return 0;
   }

   return par[0] + par[1] * x[0];
}
//http://paste.ubuntu.com/12347407/
// TF1 fitExclusion("tmpFit",residualBg, gMassCorrectionFitRange.first, gMassCorrectionFitRange.second, 2);
// fitExclusion.SetParameters(1,1);
// hist->Fit(&fitExclusion,"QRNL");

Double_t LineBg(Double_t *x, Double_t *par)
{
   return par[0] + par[1] * x[0];
}


void plotInvmMB(char* fileName = "../../D0_data_mix")
{

   // const Int_t nCent = 10;
   // const Int_t iCent[nCent + 1] = {2, 6, 11, 4, 11, 6, 9, 4, 9, 2, 9 };
   // const char cent[nCent][20] = {"40-80%", "0-40%", "Null", "0-60%", "NULL", "10-40%", "NULL", "10-60%", "NULL", "10-80%"};
   // Double_t Cen[nCent] = { 0.6 , 0.2 , 0.0,  0.3 , 0.0, 0.25, 0.0,  0.35, 0.0, 0.45 };
   // Double_t Cenerr[nCent] = {0.2,  0.2 , 0.0, 0.3, 0.0, 0.15, 0.0,  0.25, 0.0, 0.35 };
   //
   // Double_t nNcollCent[nCent] = {59.01580, 401.44672 , 1,  392.44598 , 1, 401.44672 , 1, 278.421 , 1, 205.446};
   // Double_t ndNcollCent[nCent] = {1};

   const Int_t nCent = 5;
   const Int_t iCent[nCent + 1] = {2,  6, 9, 11, 2, 11};
   const char cent[nCent][20] = {"40-80%", "10-40%", "0-10%", "Null", "0-80%"};
   Double_t Cen[nCent] = { 0.6 , 0.25 , 0.05, 0.0,  0.4  };
   Double_t Cenerr[nCent] = {0.2,  0.15 , 0.05, 0.0, 0.4};

   Double_t nNcollCent[nCent] = {59.01580, 401.44672 , 959.42547, 1, 303.78965};
   Double_t ndNcollCent[nCent] = {14.33639, 31.22933, 27.80131, 1, 21.22558};
   //
   // // const Int_t nCent = 6;
   // const Int_t iCent[nCent + 1] = {2,  6, 8, 9, 11, 2, 11};
   // const char cent[nCent][20] = {"40-80%", "20-40%", "10-20%", "0-10%", "Null", "0-80%"};
   // Double_t Cen[nCent] = { 0.6 , 0.3 , 0.15 , 0.05, 0.0,  0.4  };
   // Double_t Cenerr[nCent] = {0.2,  0.1 , 0.05 , 0.05, 0.0, 0.4};
   //
   // Double_t nNcollCent[nCent] = {59.01580, 299.07347, 606.93118, 959.42547, 1, 303.78965};
   // Double_t ndNcollCent[nCent] = {14.33639, 31.49844, 30.60806, 27.80131, 1, 21.22558};
   //
   const Int_t nPtBins = 13;
   const Double_t ptEdge[nPtBins + 1] =  { 0., 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 7.0, 8.0, 10 };
   // const Int_t nPtBins = 10;
   // const Double_t ptEdge[nPtBins + 1] =  { 0., 0.5, 1.0, 1.50, 2.0, 2.50, 3.0, 4.0, 6.0, 8., 10 };
   // const Int_t nPtBins = 9;
   // const Double_t ptEdge[nPtBins + 1] =  { 0., 1.0,  2.0, 2.50, 3.0, 3.50, 4.0, 5.0, 8.,10 };


   //Next parameters are fit with 0-80% data
   double MBFitMean1[nPtBins] =
   {
      1.86471 - 3.*0.00301412,
      1.86471 - 3.*0.00288264,
      1.85808 - 3.*0.00802594,
      1.86168 - 3.*0.00154701,
      1.86467 - 3.*0.00127562,
      1.8664 - 3.*0.00437981,
      1.86219 - 3.*0.00186508,
      1.86509 - 3.*0.00171032,
      1.86371 - 3.*0.00151382,
      1.8664 - 3.*0.00437981,
      1.86219 - 3.*0.00186508,
      1.86509 - 3.*0.00171032,
      1.86371 - 3.*0.00151382,
      //  1.8656 -3.*0.00172028,
   };
   double MBFitMean2[nPtBins] =
   {
      1.86471 + 3.*0.00301412,
      1.86471 + 3.*0.00288264,
      1.85808 + 3.*0.00802594,
      1.86168 + 3.*0.00154701,
      1.86467 + 3.*0.00127562,
      1.8664 + 3.*0.00437981,
      1.86219 + 3.*0.00186508,
      1.86509 + 3.*0.00171032,
      1.86371 + 3.*0.00151382,
      1.8664 + 3.*0.00437981,
      1.86219 + 3.*0.00186508,
      1.86509 + 3.*0.00171032,
      1.86371 + 3.*0.00151382,
//  1.8656 +3.*0.00172028,
   };
   double MBFitSigma1[nPtBins] =
   {
      0.0175737 - 3.*0.0129536 ,
      0.0128005 - 3.*0.0112506 ,
      0.0170953 - 3.*0.00251076,
      0.0122122 - 3.*0.00124907,
      0.0115828 - 3.*0.00124784,
      0.0186458 - 3.*0.00179435,
      0.0178611 - 3.*0.00169469,
      0.0169888 - 3.*0.00176045,
      0.0199233 - 3.*0.00170984,
      0.0186458 - 3.*0.00179435,
      0.0178611 - 3.*0.00169469,
      0.0169888 - 3.*0.00176045,
      0.0199233 - 3.*0.00170984,
//  0.0212652-3.*0.00191156,
   };
   double MBFitSigma2[nPtBins] =
   {
      0.0175737 + 3.*0.0129536 ,
      0.0128005 + 3.*0.0112506 ,
      0.0170953 + 3.*0.00251076,
      0.0122122 + 3.*0.00124907,
      0.0115828 + 3.*0.00124784,
      0.0186458 + 3.*0.00179435,
      0.0178611 + 3.*0.00169469,
      0.0169888 + 3.*0.00176045,
      0.0199233 + 3.*0.00170984,
      0.0186458 + 3.*0.00179435,
      0.0178611 + 3.*0.00169469,
      0.0169888 + 3.*0.00176045,
      0.0199233 + 3.*0.00170984,
//  0.0212652+3.*0.00191156,
   };

   // Int_t mColor[nCent] = {1, 2, 4, 6, 1, 1};
   // Int_t mStyle[nCent] = {20, 24, 28, 20, 34, 34};
   // Int_t mSize[nCent] = {1.3, 0.8, 2.5, 0.8, 2.5, 2.5};
   Int_t mColor[nCent] = {1, 2,  6, 1, 1};
   Int_t mStyle[nCent] = {20, 24,  20, 34, 34};
   Int_t mSize[nCent] = {1.3, 0.8,  0.8, 2.5, 2.5};

   for (Int_t i = 0; i < nCent; i++)
   {
      Cen[i] *= 100;
      Cenerr[i] *= 100;
   }

   Double_t pt[nPtBins];
   Double_t pterr[nPtBins];
   for (int i = 0; i < nPtBins; i++)
   {
      pt[i] = 0.5 * (ptEdge[i] + ptEdge[i + 1]);
      pterr[i] = pt[i] - ptEdge[i];
   }


   Double_t countCen[nCent][nPtBins] = {0};
   Double_t counterrCen[nCent][nPtBins] = {0};
   Double_t signifCen[nCent][nPtBins] = {0};
   Double_t signiferrCen[nCent][nPtBins] = {0};
   Double_t dNcountCen[nCent][nPtBins] = {0};
   Double_t dNcounterrCen[nCent][nPtBins] = {0};
   Double_t dN2countCen[nCent][nPtBins] = {0};
   Double_t dN2counterrCen[nCent][nPtBins] = {0};

   double   gMean[nCent][nPtBins] = {0};
   double   gMeanErr[nCent][nPtBins] = {0};
   double   gSigma[nCent][nPtBins] = {0};
   double   gSigmaErr[nCent][nPtBins] = {0};
   TGraphErrors *GraphD0Mean[nCent];
   TGraphErrors *GraphD0Sigma[nCent];

   const Int_t nRebin = 2;


//==============begin of read histrgrams
   // TH3F* mh3InvariantMassVsPtVsCent;
   // TH3F* mh3InvariantMassVsPtVsCentLike;
   // TH3F* mh3InvariantMassVsPtVsCentTof;
   // TH3F* mh3InvariantMassVsPtVsCentTofLike;
//==============end of read files

   THnF *hD0CentPtEtaMDphi, *hD0CentPtEtaMDphiLikeSign, *hD0CentPtEtaMDphiMixed, *hD0CentPtEtaMDphiLikeSignMixed;
   THnF *hNumInvMassvsPtMB, *hDenInvMassvsPtMBLike, *hDenInvMassvsPtMBMixUnLike, *hDenInvMassvsPtMBMixLike;
   TH1F *hNumInvMassMB[nCent][nPtBins], *hDenInvMassMBLike[nCent][nPtBins], *hDenInvMassMBMixUnLike[nCent][nPtBins], *hDenInvMassMBMixLike[nCent][nPtBins];

   TH1F *hNumInvMassMBCopy[nCent][nPtBins], *hDenInvMassMBLikeCopy[nCent][nPtBins], *hDenInvMassMBMixUnLikeCopy[nCent][nPtBins], *hDenInvMassMBMixLikeCopy[nCent][nPtBins];
   TH1F *mh1Cent;
   Double_t nEventsCent[nCent];

   char buf[1024];
   sprintf(buf, "%s.root", fileName);
   TFile *f = new TFile(buf);
   // Range from MixedEvent Code
   // int nBinsDaug[nDimDaug] = {9, 100, 10, 250, 10};  //******//cent, pt, daughterpt1, m, daughterpt2
   // double xMinDaug[nDimDaug] = {0, 0, 0.6, 0, 0.6};
   // double xMaxDaug[nDimDaug] = {9, 10, 1.6, 2.5, 1.6};
   // mh3InvariantMassVsPtVsCentTof = (TH3F*) f->Get("mh3InvariantMassVsPtVsCentTof");
   // mh3InvariantMassVsPtVsCentTofLike = (TH3F*) f->Get("mh3InvariantMassVsPtVsCentTofLike");
   //
   // mh3InvariantMassVsPtVsCentTof->SetName("mh3InvariantMassVsPtVsCentTof");
   // mh3InvariantMassVsPtVsCentTofLike->SetName("mh3InvariantMassVsPtVsCentTofLike");

   hNumInvMassvsPtMB = (THnF*) f->Get("hD0CentPtEtaMDphiDaug_standard");
   hDenInvMassvsPtMBLike = (THnF*) f->Get("hD0CentPtEtaMDphiDaugLikeSign_standard");
   hDenInvMassvsPtMBMixUnLike = (THnF*) f->Get("hD0CentPtEtaMDphiDaugMixed_standard");
   hDenInvMassvsPtMBMixLike = (THnF*) f->Get("hD0CentPtEtaMDphiDaugLikeSignMixed_standard");

   hNumInvMassvsPtMB->SetName("hNumInvMassvsPtMB");
   hDenInvMassvsPtMBLike->SetName("hDenInvMassvsPtMBLike");
   hDenInvMassvsPtMBMixUnLike->SetName("hDenInvMassvsPtMBMixUnLike");
   hDenInvMassvsPtMBMixLike->SetName("hDenInvMassvsPtMBMixLike");
   mh1Cent = (TH1F*)f->Get("hCentralityWeighted");
   mh1Cent->SetName("mh1Cent");

//sprintf(buf,"D0_Counts_Cent_%s.dat",fileName);
   sprintf(buf, "D0_Counts_Cent.dat");
   ofstream outdata(buf);
   outdata <<  "centrality   " << "nEvents   " << "NColl   " << "dNColl   " << "pTLow   " << "pTHigh   " << "count   " << "counterr   " << "dN2count   " << "dN2counterr   " <<  endl;
   sprintf(buf, "D0_FitParameter.dat");
   ofstream outFit(buf);


   //Hard Code, Daughter pT Cut
   //========================
   double DaugptCut = 0.3;
   //========================
   const int nDaugPtBins = 10;
   const int nMixRange = nPtBins - 4;//exclude the last pT bin [5-8], use side band method for [5-8]GeV
   // const Double_t ptEdge[nPtBins + 1] =  { 0., 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 7.0, 8.0, 10 };
   int iDaugpTbin1;
   int iDaugpTbin2;
   for (int ic = 0; ic < nCent; ic++)
   {
      if (iCent[ic + 1] <= iCent[ic]) continue;
      for (int i = 0; i < nPtBins; i++)
      {
         if (ptEdge[i + 1] <= ptEdge[i]) continue;

         //SameEvent Unlike
         sprintf(buf, "hNumInvMassMB_Cent%d_Pt%d", ic, i);
         Int_t ibin1 = hNumInvMassvsPtMB->GetAxis(1)->FindBin(ptEdge[i] + 1e-6);
         Int_t ibin2 = hNumInvMassvsPtMB->GetAxis(1)->FindBin(ptEdge[i + 1] - 1e-6);
         hNumInvMassvsPtMB->GetAxis(1)->SetRange(ibin1, ibin2);
         hNumInvMassvsPtMB->GetAxis(0)->SetRange(iCent[ic] - 1, iCent[ic + 1] - 2);
         // cout << "icen =============== " << iCent[ic] - 1 << "   ;" << iCent[ic + 1] - 2 << endl;
         iDaugpTbin1 = hNumInvMassvsPtMB->GetAxis(2)->FindBin(DaugptCut + 1e-6);
         iDaugpTbin2 = hNumInvMassvsPtMB->GetAxis(4)->FindBin(DaugptCut + 1e-6);
         hNumInvMassvsPtMB->GetAxis(2)->SetRange(iDaugpTbin1, nDaugPtBins + 1);
         hNumInvMassvsPtMB->GetAxis(4)->SetRange(iDaugpTbin2, nDaugPtBins + 1);
         hNumInvMassMB[ic][i] = (TH1F*) hNumInvMassvsPtMB->Projection(3, "E");
         hNumInvMassMB[ic][i]->SetName(buf);
         sprintf(buf, "hNumInvMassMB_Cent%s_Pt%2.1f_%2.1f", cent[ic], ptEdge[i], ptEdge[i + 1]);
         hNumInvMassMB[ic][i]->SetTitle(buf);
         // cout << "Maxium========" << hNumInvMassMB[ic][i]->GetMaximum() << endl;
         // cout << "Minum========" << hNumInvMassMB[ic][i]->GetMinimum() << endl;

         //SameEvent like
         sprintf(buf, "hDenInvMassMBLike_Cent%d_Pt%d", ic, i);
         Int_t ibin3 = hDenInvMassvsPtMBLike->GetAxis(1)->FindBin(ptEdge[i] + 1e-6);
         Int_t ibin4 = hDenInvMassvsPtMBLike->GetAxis(1)->FindBin(ptEdge[i + 1] - 1e-6);
         if (ibin3 != ibin1 || ibin4 != ibin2) cout << "Bins has some issue,need check" << endl;
         hDenInvMassvsPtMBLike->GetAxis(1)->SetRange(ibin1, ibin2);
         hDenInvMassvsPtMBLike->GetAxis(0)->SetRange(iCent[ic] - 1, iCent[ic + 1] - 2);

         iDaugpTbin1 = hDenInvMassvsPtMBLike->GetAxis(2)->FindBin(DaugptCut + 1e-6);
         iDaugpTbin2 = hDenInvMassvsPtMBLike->GetAxis(4)->FindBin(DaugptCut + 1e-6);
         hDenInvMassvsPtMBLike->GetAxis(2)->SetRange(iDaugpTbin1, nDaugPtBins + 1);
         hDenInvMassvsPtMBLike->GetAxis(4)->SetRange(iDaugpTbin2, nDaugPtBins + 1);

         hDenInvMassMBLike[ic][i] = (TH1F*) hDenInvMassvsPtMBLike->Projection(3, "E");
         hDenInvMassMBLike[ic][i]->SetName(buf);
         sprintf(buf, "hDenInvMassMBLike_Cent%s_Pt%2.1f_%2.1f", cent[ic], ptEdge[i], ptEdge[i + 1]);
         hDenInvMassMBLike[ic][i]->SetTitle(buf);

         //MixEvent Unlike
         sprintf(buf, "hDenInvMassMBMixUnLike_Cent%d_Pt%d", ic, i);
         Int_t ibin5 = hDenInvMassvsPtMBMixUnLike->GetAxis(1)->FindBin(ptEdge[i] + 1e-6);
         Int_t ibin6 = hDenInvMassvsPtMBMixUnLike->GetAxis(1)->FindBin(ptEdge[i + 1] - 1e-6);
         if (ibin5 != ibin1 || ibin6 != ibin2) cout << "Bins has some issue,need check" << endl;
         hDenInvMassvsPtMBMixUnLike->GetAxis(1)->SetRange(ibin1, ibin2);
         hDenInvMassvsPtMBMixUnLike->GetAxis(0)->SetRange(iCent[ic] - 1, iCent[ic + 1] - 2);
         iDaugpTbin1 = hDenInvMassvsPtMBMixUnLike->GetAxis(2)->FindBin(DaugptCut + 1e-6);
         iDaugpTbin2 = hDenInvMassvsPtMBMixUnLike->GetAxis(4)->FindBin(DaugptCut + 1e-6);
         hDenInvMassvsPtMBMixUnLike->GetAxis(2)->SetRange(iDaugpTbin1, nDaugPtBins + 1);
         hDenInvMassvsPtMBMixUnLike->GetAxis(4)->SetRange(iDaugpTbin2, nDaugPtBins + 1);

         hDenInvMassMBMixUnLike[ic][i] = (TH1F*) hDenInvMassvsPtMBMixUnLike->Projection(3, "E");
         hDenInvMassMBMixUnLike[ic][i]->SetName(buf);
         sprintf(buf, "hDenInvMassMBMixUnLike_Cent%s_Pt%2.1f_%2.1f", cent[ic], ptEdge[i], ptEdge[i + 1]);
         hDenInvMassMBMixUnLike[ic][i]->SetTitle(buf);

         //MixEvent like
         sprintf(buf, "hDenInvMassMBMixLike_Cent%d_Pt%d", ic, i);
         Int_t ibin7 = hDenInvMassvsPtMBMixLike->GetAxis(1)->FindBin(ptEdge[i] + 1e-6);
         Int_t ibin8 = hDenInvMassvsPtMBMixLike->GetAxis(1)->FindBin(ptEdge[i + 1] - 1e-6);
         if (ibin7 != ibin1 || ibin8 != ibin2) cout << "Bins has some issue,need check" << endl;
         hDenInvMassvsPtMBMixLike->GetAxis(1)->SetRange(ibin1, ibin2);
         hDenInvMassvsPtMBMixLike->GetAxis(0)->SetRange(iCent[ic] - 1, iCent[ic + 1] - 2);
         iDaugpTbin1 = hDenInvMassvsPtMBMixLike->GetAxis(2)->FindBin(DaugptCut + 1e-6);
         iDaugpTbin2 = hDenInvMassvsPtMBMixLike->GetAxis(4)->FindBin(DaugptCut + 1e-6);
         hDenInvMassvsPtMBMixLike->GetAxis(2)->SetRange(iDaugpTbin1, nDaugPtBins + 1);
         hDenInvMassvsPtMBMixLike->GetAxis(4)->SetRange(iDaugpTbin2, nDaugPtBins + 1);

         hDenInvMassMBMixLike[ic][i] = (TH1F*) hDenInvMassvsPtMBMixLike->Projection(3, "E");
         hDenInvMassMBMixLike[ic][i]->SetName(buf);
         sprintf(buf, "hDenInvMassMBMixLike_Cent%s_Pt%2.1f_%2.1f", cent[ic], ptEdge[i], ptEdge[i + 1]);
         hDenInvMassMBMixLike[ic][i]->SetTitle(buf);

      }
   }

//making plots
   gStyle->Reset("plain");
   TCanvas *c1 = new TCanvas("c1", "c1", 0, 0, 700, 600);
   gStyle->SetOptFit(0);
   gStyle->SetOptStat(0);
   gStyle->SetEndErrorSize(2);

   setpad(c1, 0.14, 0.02, 0.01, 0.12);

   TF1 *lineFun1 = new TF1("lineFun1", line, 1.7, 2.1, 2);

   TF1 *myGausLine1 = new TF1("myGausLine1", "[0]*exp(-0.5*pow((x-[1])/[2],2))/sqrt(2.*TMath::Pi())/[2]+[3]+[4]*x", 1.7, 2.1);
   myGausLine1->SetLineColor(2);
   myGausLine1->SetLineWidth(2);
   myGausLine1->SetParNames("#D0", "Mass", "Width", "a0", "a1");
   myGausLine1->SetParameters(10, 1.87, 0.015, 1, 1);

   TF1 *lineBg = new TF1("lineBg", LineBg, gMassCorrectionFitRange.first, gMassCorrectionFitRange.second, 2);
   TF1 *fitExclusion = new TF1("fitExclusion", residualBg, gMassCorrectionFitRange.first, gMassCorrectionFitRange.second, 2);
   lineBg->SetLineColor(3);
   lineBg->SetLineStyle(2);
   lineBg->SetLineWidth(2);
   lineBg->SetParNames("p0", "p1");
// fitExclusion.SetParameters(1,1);
// hist->Fit(&fitExclusion,"QRNL");
   fitExclusion->SetLineColor(4);

   TFile *Xin_file = new TFile("D0_Counts_Cent.root", "recreate");
   for (int ic = 0; ic < nCent; ic++)
   {
      if (iCent[ic + 1] <= iCent[ic]) continue;
      nEventsCent[ic] = mh1Cent->Integral(iCent[ic] - 1, iCent[ic + 1] - 2);
      cout << "nEventsCent[ic]=" << nEventsCent[ic] << endl;
      // Use mixed Event method to subtract Yield
      //Begin 0-5GeV// mixed event yield
      for (int i = 0; i < nMixRange; i++)//for Mixed Event
      {
         if (ptEdge[i + 1] <= ptEdge[i]) continue;

         if (ic == 0)
         {
            //      hNumInvMassMB[ic][i]->Rebin(nRebin);
            //      hDenInvMassMBMixLike[ic][i]->Rebin(nRebin);
         }

         // if (!(ic == 2 && i == 3)) //especial for pt2.0-2.5GeV //use likesign to subtract background
         {
            Double_t x1 = 1.66, x2 = 2.08;
            c1->SetLogy(0);
            Int_t ibin = hNumInvMassMB[ic][i]->FindBin(x1);
            x1 = hNumInvMassMB[ic][i]->GetBinLowEdge(ibin);
            hNumInvMassMB[ic][i]->GetXaxis()->SetRangeUser(x1, x2);
            Double_t max = hNumInvMassMB[ic][i]->GetMaximum();
            //      cout<<"max =============== "<<max<<endl;
            ibin = hNumInvMassMB[ic][i]->FindBin(x2);
            x2 = hNumInvMassMB[ic][i]->GetBinLowEdge(ibin + 1);
            Double_t min = hNumInvMassMB[ic][i]->GetBinContent(ibin - 1);
            //cout<<"ibin = "<<ibin<<endl;
            Double_t y1 = 0.7 * min;
            if (min <= 0) y1 = -5.0;
            Double_t y2 = 1.5 * max;
            cout << "ymin = " << y1 << ", ymax = " << y2 << endl;

            hNumInvMassMB[ic][i]->Sumw2();
            hNumInvMassMB[ic][i]->GetXaxis()->SetRangeUser(x1, x2);
            hNumInvMassMB[ic][i]->GetYaxis()->SetRangeUser(y1, y2);
            hNumInvMassMB[ic][i]->SetMarkerStyle(20);
            hNumInvMassMB[ic][i]->SetMarkerSize(1);
            hNumInvMassMB[ic][i]->SetLineColor(1);
            hNumInvMassMB[ic][i]->SetLineWidth(2.);
            hNumInvMassMB[ic][i]->SetMarkerColor(1);
            // hNumInvMassMB[ic][i]->SetTitle("");
            hNumInvMassMB[ic][i]->Draw();

            cout << hNumInvMassMB[ic][i]->GetMaximum() << endl;
            cout << hNumInvMassMB[ic][i]->GetYaxis()->GetBinUpEdge(1) << endl;
            sprintf(buf, "Counts/(%d MeV/c^{2})", 1000 * hNumInvMassMB[ic][i]->GetXaxis()->GetBinWidth(2));
            hNumInvMassMB[ic][i]->GetYaxis()->SetTitle(buf);
            hNumInvMassMB[ic][i]->GetYaxis()->SetTitleOffset(0.8);
            hNumInvMassMB[ic][i]->GetXaxis()->SetTitle("M_{K#pi} (GeV/c^{2})");
            hNumInvMassMB[ic][i]->GetXaxis()->SetTitleOffset(1.0);
            hNumInvMassMB[ic][i]->GetXaxis()->SetTitleSize(0.06);
            hNumInvMassMB[ic][i]->GetYaxis()->SetTitleSize(0.06);
            hNumInvMassMB[ic][i]->GetXaxis()->CenterTitle();
            hNumInvMassMB[ic][i]->GetYaxis()->CenterTitle();
            //Double_t max = hNumInvMassMB[ic][i]->GetMaximum()*1.2;
            //hNumInvMassMB[ic][i]->SetMaximum(max);
            //hNumInvMassMB[ic][i]->SetMinimum(0);
            Int_t ibin = hNumInvMassMB[ic][i]->FindBin(x1);
            x1 = hNumInvMassMB[ic][i]->GetBinLowEdge(ibin);
            ibin = hNumInvMassMB[ic][i]->FindBin(x2);
            x2 = hNumInvMassMB[ic][i]->GetBinLowEdge(ibin + 1);
            drawLines(x1, y1, x2, y2, 3, 1);

            hDenInvMassMBLike[ic][i]->SetLineColor(4);
            hDenInvMassMBLike[ic][i]->SetLineWidth(2.);
            hDenInvMassMBLike[ic][i]->Draw("same");

            //normalize mixed-event background
            Double_t norm = 0.1;
            /* Double_t m1=1.95, m2=2.10; */
            Double_t m1 = 1.6, m2 = 2.1;
            // Double_t m1 = 1.8, m2 = 2.0;
            Int_t mbin1 = hDenInvMassMBMixLike[ic][i]->FindBin(m1 + 0.0001);
            Int_t mbin2 = hDenInvMassMBMixLike[ic][i]->FindBin(m2 - 0.0001);
            norm = hDenInvMassMBLike[ic][i]->Integral(mbin1, mbin2) / hDenInvMassMBMixLike[ic][i]->Integral(mbin1, mbin2);
            cout << "norm.=" << norm << endl;

            hDenInvMassMBMixUnLike[ic][i]->Scale(norm);

            hDenInvMassMBMixUnLike[ic][i]->SetLineColor(2);
            hDenInvMassMBMixUnLike[ic][i]->SetLineWidth(2.);
            hDenInvMassMBMixUnLike[ic][i]->Draw("same");

            //Copy to Two panles
            sprintf(buf, "hNumInvMassMBCopy_Cent%d_Pt%d", ic, i);
            hNumInvMassMBCopy[ic][i] = (TH1F*) hNumInvMassMB[ic][i]->Clone(buf);
            hNumInvMassMBCopy[ic][i] ->SetName(buf);
            sprintf(buf, "hNumInvMassMBCopy_Cent%s_Pt%2.1f_%2.1f", cent[ic], ptEdge[i], ptEdge[i + 1]);
            hNumInvMassMBCopy[ic][i] ->SetTitle(buf);
            hNumInvMassMBCopy[ic][i]->GetYaxis()->SetTitleOffset(0.8);
            hNumInvMassMBCopy[ic][i]->GetXaxis()->SetTitleSize(0.09);
            hNumInvMassMBCopy[ic][i]->GetXaxis()->SetRangeUser(x1, x2);
            hNumInvMassMBCopy[ic][i]->GetYaxis()->SetRangeUser(y1, y2);
            //Copy to two panles
            sprintf(buf, "hDenInvMassMBLikeCopy_Cent%d_Pt%d", ic, i);
            hDenInvMassMBLikeCopy[ic][i] = (TH1F*) hDenInvMassMBLike[ic][i]->Clone(buf);
            hDenInvMassMBLikeCopy[ic][i]->SetName(buf);
            sprintf(buf, "hDenInvMassMBLikeCopy_Cent%s_Pt%2.1f_%2.1f", cent[ic], ptEdge[i], ptEdge[i + 1]);
            hDenInvMassMBLikeCopy[ic][i]->SetTitle(buf);
            //Copy to two panles
            sprintf(buf, "hDenInvMassMBMixUnLikeCopy_Cent%d_Pt%d", ic, i);
            hDenInvMassMBMixUnLikeCopy[ic][i] = (TH1F*) hDenInvMassMBMixUnLike[ic][i]->Clone(buf);
            hDenInvMassMBMixUnLikeCopy[ic][i]->SetName(buf);
            sprintf(buf, "hDenInvMassMBMixUnLikeCopy_Cent%s_Pt%2.1f_%2.1f", cent[ic], ptEdge[i], ptEdge[i + 1]);
            hDenInvMassMBMixUnLikeCopy[ic][i]->SetTitle(buf);
            hDenInvMassMBMixUnLikeCopy[ic][i]->Sumw2();

            /*             */

            myGausLine1->SetNpx(501);
            lineFun1->SetNpx(501);
            fitExclusion->SetNpx(501);
            lineBg->SetNpx(501);

            // Begin of Yield extraction
            c1->Clear();
            c1->Update();
            setpad(c1, 0.14, 0.02, 0.01, 0.12);
            //gStyle->SetOptFit(1111);
            gStyle->SetOptFit(0);
            hDenInvMassMBLike[ic][i]->Sumw2();
            hDenInvMassMBMixUnLike[ic][i]->Scale(1. / norm);
            hDenInvMassMBMixUnLike[ic][i]->Sumw2();
            hDenInvMassMBMixUnLike[ic][i]->Scale(norm);
            // hNumInvMassMB[ic][i]->Add(hDenInvMassMBLike[ic][i], -1);
            hNumInvMassMB[ic][i]->Add(hDenInvMassMBMixUnLike[ic][i], -1);
            Double_t max = hNumInvMassMB[ic][i]->GetMaximum() * 1.6;
            hNumInvMassMB[ic][i]->SetMaximum(max);
            Double_t min = hNumInvMassMB[ic][i]->GetMinimum();
            //manual tune the Yaxis-range
            if (min < -10 && i == 0) min = -100;
            else if (min < -10 && i == 1) min = -100;
            else if (min < -10 && i == 2) min = -100;
            else min = min * 1.5;
            hNumInvMassMB[ic][i]->SetMinimum(min);
            for (int ip = 0; ip < 5; ip++)
            {
               myGausLine1->ReleaseParameter(i);
            }
            myGausLine1->SetParameters(10, 1.87, 0.015, 1, 1);
            myGausLine1->SetParLimits(0, 0, 1e8);
            myGausLine1->SetParLimits(1, MBFitMean1[i], MBFitMean2[i]);
            myGausLine1->SetParLimits(2, MBFitSigma1[i], MBFitSigma2[i]);
            // cout<<MBFitMean1[i]<<","<<MBFitMean2[i]<<","<<MBFitSigma1[i]<<","<<MBFitSigma2[i]<<endl;
            hNumInvMassMB[ic][i]->Fit(myGausLine1, "r0");

            lineFun1->SetParameter(0, myGausLine1->GetParameter(3));
            lineFun1->SetParameter(1, myGausLine1->GetParameter(4));
            lineFun1->SetLineStyle(2);
            lineFun1->SetLineColor(2);
            lineFun1->SetLineWidth(2);

            Float_t chi2 = myGausLine1->GetChisquare();
            Float_t ndf = myGausLine1->GetNDF();
            Float_t amp = myGausLine1->GetParameter(0);
            Float_t amperror = myGausLine1->GetParError(0);
            Float_t mean = myGausLine1->GetParameter(1);
            Float_t meanerror = myGausLine1->GetParError(1);
            Float_t sigma = myGausLine1->GetParameter(2);
            Float_t sigmaerror = myGausLine1->GetParError(2);
            //print out 0-80% parameters, used to contrain others
            if (ic == 4)outFit << cent[ic] << " " << ptEdge[i] << "  " << ptEdge[i + 1] << " SubtractFit:  " << myGausLine1->GetParameter(0) << "+-" << myGausLine1->GetParError(0) << " ;" << myGausLine1->GetParameter(1) << "+2.*" << myGausLine1->GetParError(1) << " ," << myGausLine1->GetParameter(1) << "-2.*" << myGausLine1->GetParError(1) << " ," << myGausLine1->GetParameter(2) << "+2.*" << myGausLine1->GetParError(2) << " ," << myGausLine1->GetParameter(2) << "-2.*" << myGausLine1->GetParError(2) << " ," << myGausLine1->GetParameter(3) << "+-" << myGausLine1->GetParError(3) << " ;" << myGausLine1->GetParameter(4) << "+-" << myGausLine1->GetParError(4) << " ;" << endl;
            gMean[ic][i] = mean;
            gMeanErr[ic][i] = meanerror;
            gSigma[ic][i] = sigma;
            gSigmaErr[ic][i] = sigmaerror;
            sprintf(buf, "Gaus+line:");
            drawLatex(0.65, 0.84, buf, 42, 0.045, 1);
            sprintf(buf, "#chi^{2}/ndf : %3.1f/%3.1f", chi2, ndf);
            drawLatex(0.65, 0.79, buf, 42, 0.045, 1);
            sprintf(buf, "const  : %3.2f#pm%3.3f", amp, amperror);
            drawLatex(0.65, 0.73, buf, 42, 0.045, 1);
            sprintf(buf, "mean  : %3.2f#pm%3.3f", mean, meanerror);
            drawLatex(0.65, 0.68, buf, 42, 0.045, 1);
            sprintf(buf, "#sigma        : %3.3f#pm%3.3f", sigma, sigmaerror);
            drawLatex(0.65, 0.63, buf, 42, 0.045, 1);

            int ibintmp1 = hNumInvMassMB[ic][i]->GetXaxis()->FindBin(mean - 3.0 * sigma + 1.e-6);
            int ibintmp2 = hNumInvMassMB[ic][i]->GetXaxis()->FindBin(mean + 3.0 * sigma - 1.e-6);

            double ttmp1 = hNumInvMassMB[ic][i]->GetBinLowEdge(hNumInvMassMB[ic][i]->FindBin(mean - 3.0 * sigma + 1.e-6));
            double ttmp2 = hNumInvMassMB[ic][i]->GetBinLowEdge(hNumInvMassMB[ic][i]->FindBin(mean + 3.0 * sigma - 1.e-6) + 1);

            double fitSignalSub = hNumInvMassMB[ic][i]->Integral(ibintmp1, ibintmp2) - lineFun1->Integral(ttmp1, ttmp2) / hNumInvMassMB[ic][i]->GetBinWidth(2);

            // double fitSignalSub = hNumInvMassMB[ic][i]->Integral(ibintmp1, ibintmp2) - lineFun1->Integral(mean - 3.0 * sigma, mean + 3.0 * sigma) / hNumInvMassMB[ic][i]->GetBinWidth(2);
            double fitSignalSubErr = 0 ;
            for (int ijk = ibintmp1; ijk <= ibintmp2; ++ijk)
            {
               fitSignalSubErr += TMath::Power(hNumInvMassMB[ic][i]->GetBinError(ijk), 2);
            }
            // fitSignalSubErr = sqrt(fitSignalSubErr);

            fitExclusion->SetParameter(0, myGausLine1->GetParameter(3));
            fitExclusion->SetParameter(1, myGausLine1->GetParameter(4));
            // TFitResultPtr r = hNumInvMassMB[ic][i]->Fit(fitExclusion, "QRLS0");
            TFitResultPtr r = hNumInvMassMB[ic][i]->Fit(fitExclusion, "QRFS0");
            // cout << "**" << fitExclusion->GetParameter(0);
            // cout << "**" << fitExclusion->GetParameter(1);
            // cout << "**" << fitExclusion->GetParError(0);
            // cout << "**" << fitExclusion->GetParError(1) << endl;
            lineBg->SetParameter(0, fitExclusion->GetParameter(0));
            lineBg->SetParameter(1, fitExclusion->GetParameter(1));
            lineBg->SetParError(0, fitExclusion->GetParError(0));
            lineBg->SetParError(1, fitExclusion->GetParError(1));
            // double ResidualBg = lineBg->Integral(gMassYieldCountRange.first, gMassYieldCountRange.second) / hNumInvMassMB[ic][i]->GetBinWidth(2);
            // double ResidualBgErr = lineBg->IntegralError(gMassYieldCountRange.first, gMassYieldCountRange.second, fitExclusion->GetParameters(), r->GetCovarianceMatrix().GetMatrixArray()) / hNumInvMassMB[ic][i]->GetBinWidth(2);
            double ttmp3 = hNumInvMassMB[ic][i]->GetBinLowEdge(hNumInvMassMB[ic][i]->FindBin(gMassYieldCountRange.first + 1.e-6));
            double ttmp4 = hNumInvMassMB[ic][i]->GetBinLowEdge(hNumInvMassMB[ic][i]->FindBin(gMassYieldCountRange.second - 1.e-6) + 1);
            double ResidualBg = lineBg->Integral(ttmp3, ttmp4) / hNumInvMassMB[ic][i]->GetBinWidth(2);

            double ResidualBgErr = lineBg->IntegralError(ttmp3, ttmp4, fitExclusion->GetParameters(), r->GetCovarianceMatrix().GetMatrixArray()) / hNumInvMassMB[ic][i]->GetBinWidth(2);
            ibintmp1 = hNumInvMassMB[ic][i]->GetXaxis()->FindBin(gMassYieldCountRange.first + 1e-6);
            ibintmp2 = hNumInvMassMB[ic][i]->GetXaxis()->FindBin(gMassYieldCountRange.second - 1e-6);
            fitSignalSub = hNumInvMassMB[ic][i]->Integral(ibintmp1, ibintmp2) - ResidualBg;
            fitSignalSubErr += TMath::Power(ResidualBgErr, 2);
            fitSignalSubErr = sqrt(fitSignalSubErr);
            lineBg->Draw("same");

            sprintf(buf, "#color[2]{subtract fit residual:} #D0 = %d #pm %d", fitSignalSub, fitSignalSubErr);
            drawLatex(0.2, 0.56, buf, 42, 0.04, 1);

            drawLine(mean - 3.0 * sigma, min, mean - 3.0 * sigma, min + (max - min) * 0.6, 2, 4);
            drawLine(mean + 3.0 * sigma, min, mean + 3.0 * sigma, min + (max - min) * 0.6, 2, 4);

            drawLines(x1, min, x2, max, 3, 1);

            sprintf(buf, "Au+Au 200GeV, %s MB", cent[ic]);
            drawLatex(0.16, 0.94, buf, 62, 0.045, 1);
            sprintf(buf, "%3.2f<p_{T}<%3.2f GeV/c", ptEdge[i], ptEdge[i + 1]);
            drawLatex(0.16, 0.88, buf, 42, 0.045, 1);

            sprintf(buf, "pic/D0_InvmSub_cent_%d_%d_pt_%3.2f_%3.2f.gif", iCent[ic] - 1, iCent[ic + 1] - 2, ptEdge[i], ptEdge[i + 1]);
            c1->SaveAs(buf);
            // End of Yield extraction


            ///From  Next is to Make two pannel Plots
            //
            //Make Two panel plots
            gStyle->Reset("plain");
            TCanvas *c11 = new TCanvas("c11", "c11", 0, 0, 650, 850);
            gStyle->SetOptStat(0);
            gStyle->SetOptTitle(0);
            setpad(c11, 0.14, 0.02, 0.01, 0.12);
            c11->SetLogy(0);
            float small = 0;
            c11->Divide(1, 2, small, small);
            //Make tow panels
            c11->cd(1);
            c11_1->SetFillColor(-1);
            gPad->SetTickx();
            gPad->SetBottomMargin(small);
            hNumInvMassMBCopy[ic][i]->GetYaxis()->SetTitleOffset(0.8);
            hNumInvMassMBCopy[ic][i]->GetYaxis()->SetTitleSize(0.07);
            hNumInvMassMBCopy[ic][i]->Draw();
            hDenInvMassMBLikeCopy[ic][i]->Draw("same");
            hDenInvMassMBMixUnLikeCopy[ic][i]->Draw("same");
            drawLines(x1, y1, x2, y2, 3, 1);
            TLegend *leg = new TLegend(0.58, 0.65, 0.75, 0.94);
            leg->SetFillColor(10);
            leg->SetBorderSize(0);
            leg->SetTextFont(42);
            leg->SetTextSize(0.064);
            leg->AddEntry(hNumInvMassMBCopy[ic][i], "unlike-sign, same evt", "pl");
            leg->AddEntry(hDenInvMassMBLikeCopy[ic][i], "like-sign, same evt", "l");
            leg->AddEntry(hDenInvMassMBMixUnLikeCopy[ic][i], "unlike-sign, mixed evt", "l");
            leg->Draw("same");
            sprintf(buf, "Au+Au 200GeV, %s MB", cent[ic]);
            drawLatex(0.18, 0.94, buf, 62, 0.065, 1);
            sprintf(buf, "%3.2f<p_{T}<%3.2f GeV/c", ptEdge[i], ptEdge[i + 1]);
            drawLatex(0.18, 0.85, buf, 42, 0.065, 1);

            c11->cd(2);
            c11_2->SetFillColor(-1);
            gPad->SetTickx();
            gPad->SetTopMargin(small);
            hNumInvMassMB[ic][i]->GetYaxis()->SetTitleOffset(0.8);
            hNumInvMassMB[ic][i]->GetYaxis()->SetTitleSize(0.06);
            hNumInvMassMB[ic][i]->Draw();
            hNumInvMassMB[ic][i]->Write();
            myGausLine1->Draw("same");
            lineBg->Draw("same");
            fitExclusion->Draw("same");
            lineFun1->Draw("same");
            drawLines(x1, min, x2, max, 3, 1);
            sprintf(buf, "Au+Au 200GeV, %s MB", cent[ic]);
            drawLatex(0.18, 0.92, buf, 62, 0.06, 1);
            sprintf(buf, "%3.2f<p_{T}<%3.2f GeV/c", ptEdge[i], ptEdge[i + 1]);
            drawLatex(0.18, 0.83, buf, 42, 0.06, 1);
            sprintf(buf, "#color[2]{#D0 = %d #pm %d}", fitSignalSub, fitSignalSubErr);
            drawLatex(0.2, 0.56, buf, 42, 0.04, 1);

            sprintf(buf, "Gaus+pol1:");
            // drawLatex(0.65, 0.74, buf, 42, 0.05, 1);
            sprintf(buf, "#chi^{2}/ndf : %3.1f/%3.1f", chi2, ndf);
            drawLatex(0.65, 0.69, buf, 42, 0.05, 1);
            sprintf(buf, "mean  : %3.2f", mean);
            drawLatex(0.65, 0.63, buf, 42, 0.05, 1);
            sprintf(buf, "#sigma        : %3.3f", sigma);
            drawLatex(0.65, 0.58, buf, 42, 0.05, 1);

            sprintf(buf, "pic/D0_InvmSub_TwoPanel_cent_%d_%d_pt_%3.2f_%3.2f.gif", iCent[ic] - 1, iCent[ic + 1] - 2, ptEdge[i], ptEdge[i + 1]);
            // if (ic == 3 || ic == 5) c11->SaveAs(buf); //only save 0-10% and 0-80% plots
            c11->SaveAs(buf); //only save 0-10% and 0-80% plots
         }

         if (fitSignalSub >= 0 && fitSignalSubErr >= 0)
         {
            //From subtracted background
            countCen[ic][i] = fitSignalSub;
            counterrCen[ic][i] = fitSignalSubErr;
            signifCen[ic][i] = fitSignalSub / fitSignalSubErr;
            signiferrCen[ic][i] = 0;
            //calculate dN/2pipTdpT/Events
            dNcountCen[ic][i] = countCen[ic][i] / 2. / TMath::Pi() / pt[i] / (2.*pterr[i]) / nEventsCent[ic] / 2.;
            dNcounterrCen[ic][i] = counterrCen[ic][i] / 2. / TMath::Pi() / pt[i] / (2.*pterr[i]) / nEventsCent[ic] / 2.;
            dN2countCen[ic][i] = dNcountCen[ic][i] / nNcollCent[ic];
            dN2counterrCen[ic][i] = dNcounterrCen[ic][i] / nNcollCent[ic];
         }
         else
         {
            countCen[ic][i] = 0;
            counterrCen[ic][i] = 0;
            signifCen[ic][i] = 0;
            signiferrCen[ic][i] = 0;
            //calculate dN/2pipTdpT/Events
            dNcountCen[ic][i] = countCen[ic][i] / 2. / TMath::Pi() / pt[i] / (2.*pterr[i]) / nEventsCent[ic] / 2.;
            dNcounterrCen[ic][i] = counterrCen[ic][i] / 2. / TMath::Pi() / pt[i] / (2.*pterr[i]) / nEventsCent[ic] / 2.;
            dN2countCen[ic][i] = dNcountCen[ic][i] / nNcollCent[ic];
            dN2counterrCen[ic][i] = dNcounterrCen[ic][i] / nNcollCent[ic];
         }

         outdata << cent[ic] << "  " << nEventsCent[ic] << "  " << nNcollCent[ic] << "  " << ndNcollCent[ic] << "  " << ptEdge[i] << "  " << ptEdge[i + 1] << "  " << countCen[ic][i] << "  " << counterrCen[ic][i] << "  " << dN2countCen[ic][i] << "  " << dN2counterrCen[ic][i] <<  endl;


      }
      //End of 0-5GeV// mixed event yield


      // Use same Event method to subtract Yield  //Begin 5-8GeV//
      for (int i = nMixRange; i < nPtBins; i++)// For LikeSign Subtract // pt Large than 5 Gev/ donot use mixedEvent to subtract background, use sideband
      {
         if (ptEdge[i + 1] <= ptEdge[i]) continue;

         if (ic == 0)
         {
            //      hNumInvMassMB[ic][i]->Rebin(nRebin);
            //      hDenInvMassMBMixLike[ic][i]->Rebin(nRebin);
         }

         Double_t x1 = 1.66, x2 = 2.08;
         c1->SetLogy(0);
         Int_t ibin = hNumInvMassMB[ic][i]->FindBin(x1);
         x1 = hNumInvMassMB[ic][i]->GetBinLowEdge(ibin);
         hNumInvMassMB[ic][i]->GetXaxis()->SetRangeUser(x1, x2);
         Double_t max = hNumInvMassMB[ic][i]->GetMaximum();
         ibin = hNumInvMassMB[ic][i]->FindBin(x2);
         x2 = hNumInvMassMB[ic][i]->GetBinLowEdge(ibin + 1);
         Double_t min = hNumInvMassMB[ic][i]->GetBinContent(ibin - 1);
         Double_t y1 = 0.7 * min;
         if (min <= 0 && i == 5) y1 = -5.0;
         Double_t y2 = 1.5 * max;

         hNumInvMassMB[ic][i]->Sumw2();
         hDenInvMassMBLike[ic][i]->Sumw2();
         hNumInvMassMB[ic][i]->GetXaxis()->SetRangeUser(x1, x2);
         hNumInvMassMB[ic][i]->GetYaxis()->SetRangeUser(y1, y2);
         hNumInvMassMB[ic][i]->SetMarkerStyle(20);
         hNumInvMassMB[ic][i]->SetMarkerSize(1);
         hNumInvMassMB[ic][i]->SetLineColor(1);
         hNumInvMassMB[ic][i]->SetLineWidth(2.);
         hNumInvMassMB[ic][i]->SetMarkerColor(1);
         hNumInvMassMB[ic][i]->Draw();
         cout << hNumInvMassMB[ic][i]->GetMaximum() << endl;
         cout << hNumInvMassMB[ic][i]->GetYaxis()->GetBinUpEdge(1) << endl;
         sprintf(buf, "Counts/(%d MeV/c^{2})", 1000 * hNumInvMassMB[ic][i]->GetXaxis()->GetBinWidth(2));
         hNumInvMassMB[ic][i]->GetYaxis()->SetTitle(buf);
         hNumInvMassMB[ic][i]->GetYaxis()->SetTitleOffset(1.1);
         hNumInvMassMB[ic][i]->GetXaxis()->SetTitleSize(0.05);
         hNumInvMassMB[ic][i]->GetYaxis()->SetTitleSize(0.06);
         hNumInvMassMB[ic][i]->GetXaxis()->SetTitle("M_{K#pi} (GeV/c^{2})");
         hNumInvMassMB[ic][i]->GetXaxis()->CenterTitle();
         hNumInvMassMB[ic][i]->GetYaxis()->CenterTitle();
         hNumInvMassMB[ic][i]->GetXaxis()->SetTitleOffset(1.2);
         Int_t ibin = hNumInvMassMB[ic][i]->FindBin(x1);
         x1 = hNumInvMassMB[ic][i]->GetBinLowEdge(ibin);
         ibin = hNumInvMassMB[ic][i]->FindBin(x2);
         x2 = hNumInvMassMB[ic][i]->GetBinLowEdge(ibin + 1);
         drawLines(x1, y1, x2, y2, 3, 1);

         hDenInvMassMBLike[ic][i]->SetLineColor(4);
         hDenInvMassMBLike[ic][i]->SetLineWidth(2.);
         hDenInvMassMBLike[ic][i]->Draw("same");

         //normalize mixed-event background
         Double_t norm = 0.1;
         Double_t m1 = 1.6, m2 = 2.1;
         Int_t mbin1 = hDenInvMassMBMixLike[ic][i]->FindBin(m1 + 0.0001);
         Int_t mbin2 = hDenInvMassMBMixLike[ic][i]->FindBin(m2 - 0.0001);
         if (hDenInvMassMBMixLike[ic][i]->Integral(mbin1, mbin2) != 0) norm = hDenInvMassMBLike[ic][i]->Integral(mbin1, mbin2) / hDenInvMassMBMixLike[ic][i]->Integral(mbin1, mbin2);
         cout << "norm.=" << norm << endl;

         hDenInvMassMBMixUnLike[ic][i]->Scale(norm);

         hDenInvMassMBMixUnLike[ic][i]->SetLineColor(2);
         hDenInvMassMBMixUnLike[ic][i]->SetLineWidth(2.);
         hDenInvMassMBMixUnLike[ic][i]->Draw("same");

         //Copy to Two panles
         sprintf(buf, "hNumInvMassMBCopy_Cent%d_Pt%d", ic, i);
         hNumInvMassMBCopy[ic][i] = (TH1F*) hNumInvMassMB[ic][i]->Clone(buf);
         hNumInvMassMBCopy[ic][i] ->SetName(buf);
         sprintf(buf, "hNumInvMassMBCopy_Cent%s_Pt%2.1f_%2.1f", cent[ic], ptEdge[i], ptEdge[i + 1]);
         hNumInvMassMBCopy[ic][i] ->SetTitle(buf);
         hNumInvMassMBCopy[ic][i]->GetYaxis()->SetTitleOffset(1.1);
         hNumInvMassMBCopy[ic][i]->GetXaxis()->SetTitleSize(0.05);
         hNumInvMassMBCopy[ic][i]->GetXaxis()->SetRangeUser(x1, x2);
         hNumInvMassMBCopy[ic][i]->GetYaxis()->SetRangeUser(y1, y2);
         //Copy to two panles
         sprintf(buf, "hDenInvMassMBLikeCopy_Cent%d_Pt%d", ic, i);
         hDenInvMassMBLikeCopy[ic][i] = (TH1F*) hDenInvMassMBLike[ic][i]->Clone(buf);
         hDenInvMassMBLikeCopy[ic][i]->SetName(buf);
         sprintf(buf, "hDenInvMassMBLikeCopy_Cent%s_Pt%2.1f_%2.1f", cent[ic], ptEdge[i], ptEdge[i + 1]);
         hDenInvMassMBLikeCopy[ic][i]->SetTitle(buf);
         //Copy to two panles
         sprintf(buf, "hDenInvMassMBMixUnLikeCopy_Cent%d_Pt%d", ic, i);
         hDenInvMassMBMixUnLikeCopy[ic][i] = (TH1F*) hDenInvMassMBMixUnLike[ic][i]->Clone(buf);
         hDenInvMassMBMixUnLikeCopy[ic][i]->SetName(buf);
         sprintf(buf, "hDenInvMassMBMixUnLikeCopy_Cent%s_Pt%2.1f_%2.1f", cent[ic], ptEdge[i], ptEdge[i + 1]);
         hDenInvMassMBMixUnLikeCopy[ic][i]->SetTitle(buf);
         hDenInvMassMBMixUnLikeCopy[ic][i]->Sumw2();
         /*
             */

         c1->Clear();
         c1->Update();
         setpad(c1, 0.14, 0.02, 0.01, 0.12);
         gStyle->SetOptFit(0);
         hNumInvMassMB[ic][i]->Add(hDenInvMassMBLike[ic][i], -1);
         Double_t max = hNumInvMassMB[ic][i]->GetMaximum() * 1.6;
         hNumInvMassMB[ic][i]->SetMaximum(max);
         Double_t min = hNumInvMassMB[ic][i]->GetMinimum();
         if (min < -10) min = -10;
         else min = min * 1.5;
         hNumInvMassMB[ic][i]->SetMinimum(min);

         for (int ip = 0; ip < 5; ip++)
         {
            myGausLine1->ReleaseParameter(i);
         }
         myGausLine1->SetParameters(10, 1.87, 0.015, 1, 1);
         myGausLine1->SetParLimits(0, 0, 1e8);
         // myGausLine1->SetParLimits(1, 1.84, 1.89);
         // myGausLine1->SetParLimits(2, 0, 0.035);

         myGausLine1->SetParLimits(1, MBFitMean1[i], MBFitMean2[i]);
         myGausLine1->SetParLimits(2, MBFitSigma1[i], MBFitSigma2[i]);

         hNumInvMassMB[ic][i]->Fit(myGausLine1, "r0");

         // hNumInvMassMBCopy[ic][i]->Fit(myGausLine1, "r0");

         lineFun1->SetParameter(0, myGausLine1->GetParameter(3));
         lineFun1->SetParameter(1, myGausLine1->GetParameter(4));
         lineFun1->SetLineStyle(2);
         lineFun1->SetLineColor(2);
         lineFun1->SetLineWidth(2);
         // lineFun1->Draw("same");

         Float_t chi2 = myGausLine1->GetChisquare();
         Float_t ndf = myGausLine1->GetNDF();
         Float_t mean = myGausLine1->GetParameter(1);
         Float_t meanerror = myGausLine1->GetParError(1);
         Float_t sigma = myGausLine1->GetParameter(2);
         Float_t sigmaerror = myGausLine1->GetParError(2);
         if (ic == 4)outFit << cent[ic] << " " << ptEdge[i] << "  " << ptEdge[i + 1] << " SubtractFit:  " << myGausLine1->GetParameter(0) << "+-" << myGausLine1->GetParError(0) << " ;" << myGausLine1->GetParameter(1) << "+2.*" << myGausLine1->GetParError(1) << " ," << myGausLine1->GetParameter(1) << "-2.*" << myGausLine1->GetParError(1) << " ," << myGausLine1->GetParameter(2) << "+2.*" << myGausLine1->GetParError(2) << " ," << myGausLine1->GetParameter(2) << "-2.*" << myGausLine1->GetParError(2) << " ," << myGausLine1->GetParameter(3) << "+-" << myGausLine1->GetParError(3) << " ;" << myGausLine1->GetParameter(4) << "+-" << myGausLine1->GetParError(4) << " ;" << myGausLine1->GetParameter(5) << "+-" << myGausLine1->GetParError(5) << endl;
         gMean[ic][i] = mean;
         gMeanErr[ic][i] = meanerror;
         gSigma[ic][i] = sigma;
         gSigmaErr[ic][i] = sigmaerror;
         sprintf(buf, "Gaus+line:");
         drawLatex(0.65, 0.84, buf, 42, 0.045, 1);
         sprintf(buf, "#chi^{2}/ndf : %3.1f/%3.1f", chi2, ndf);
         drawLatex(0.65, 0.79, buf, 42, 0.045, 1);
         sprintf(buf, "const  : %3.2f#pm%3.3f", amp, amperror);
         drawLatex(0.65, 0.73, buf, 42, 0.045, 1);
         sprintf(buf, "mean  : %3.2f#pm%3.3f", mean, meanerror);
         drawLatex(0.65, 0.68, buf, 42, 0.045, 1);
         sprintf(buf, "#sigma        : %3.3f#pm%3.3f", sigma, sigmaerror);
         drawLatex(0.65, 0.63, buf, 42, 0.045, 1);

         int ibintmp1 = hNumInvMassMB[ic][i]->GetXaxis()->FindBin(mean - 3.0 * sigma + 1.e-6);
         int ibintmp2 = hNumInvMassMB[ic][i]->GetXaxis()->FindBin(mean + 3.0 * sigma - 1.e-6);

         double ttmp1 = hNumInvMassMB[ic][i]->GetBinLowEdge(hNumInvMassMB[ic][i]->FindBin(mean - 3.0 * sigma + 1.e-6));
         double ttmp2 = hNumInvMassMB[ic][i]->GetBinLowEdge(hNumInvMassMB[ic][i]->FindBin(mean + 3.0 * sigma - 1.e-6) + 1);

         double fitSignalSub = hNumInvMassMB[ic][i]->Integral(ibintmp1, ibintmp2) - lineFun1->Integral(ttmp1, ttmp2) / hNumInvMassMB[ic][i]->GetBinWidth(2);
         // double fitSignalSub = hNumInvMassMB[ic][i]->Integral(ibintmp1, ibintmp2) - lineFun1->Integral(mean - 3.0 * sigma, mean + 3.0 * sigma) / hNumInvMassMB[ic][i]->GetBinWidth(2);

         double fitSignalSubErr = 0 ;
         for (int ijk = ibintmp1; ijk <= ibintmp2; ++ijk)
         {
            fitSignalSubErr += TMath::Power(hNumInvMassMB[ic][i]->GetBinError(ijk), 2);
         }

         fitExclusion->SetParameter(0, myGausLine1->GetParameter(3));
         fitExclusion->SetParameter(1, myGausLine1->GetParameter(4));
         // TFitResultPtr r = hNumInvMassMB[ic][i]->Fit(fitExclusion, "QRLS0");
         TFitResultPtr r = hNumInvMassMB[ic][i]->Fit(fitExclusion, "QRFS0");
         lineBg->SetParameter(0, fitExclusion->GetParameter(0));
         lineBg->SetParameter(1, fitExclusion->GetParameter(1));
         lineBg->SetParError(0, fitExclusion->GetParError(0));
         lineBg->SetParError(1, fitExclusion->GetParError(1));

         double ttmp3 = hNumInvMassMB[ic][i]->GetBinLowEdge(hNumInvMassMB[ic][i]->FindBin(gMassYieldCountRange.first + 1.e-6));
         double ttmp4 = hNumInvMassMB[ic][i]->GetBinLowEdge(hNumInvMassMB[ic][i]->FindBin(gMassYieldCountRange.second - 1.e-6) + 1);
         double ResidualBg = lineBg->Integral(ttmp3, ttmp4) / hNumInvMassMB[ic][i]->GetBinWidth(2);

         double ResidualBgErr = lineBg->IntegralError(ttmp3, ttmp4, fitExclusion->GetParameters(), r->GetCovarianceMatrix().GetMatrixArray()) / hNumInvMassMB[ic][i]->GetBinWidth(2);
         // double ResidualBg = lineBg->Integral(gMassYieldCountRange.first, gMassYieldCountRange.second) / hNumInvMassMB[ic][i]->GetBinWidth(2);
         // double ResidualBgErr = lineBg->IntegralError(gMassYieldCountRange.first, gMassYieldCountRange.second, fitExclusion->GetParameters(), r->GetCovarianceMatrix().GetMatrixArray()) / hNumInvMassMB[ic][i]->GetBinWidth(2);
         ibintmp1 = hNumInvMassMB[ic][i]->GetXaxis()->FindBin(gMassYieldCountRange.first + 1e-6);
         ibintmp2 = hNumInvMassMB[ic][i]->GetXaxis()->FindBin(gMassYieldCountRange.second - 1e-6);
         fitSignalSub = hNumInvMassMB[ic][i]->Integral(ibintmp1, ibintmp2) - ResidualBg;
         int ibinSideBandBegin1 = hNumInvMassMB[ic][i]->GetXaxis()->FindBin(gMassYieldSideBand1.first + 1e-6);
         int ibinSideBandEnd1 = hNumInvMassMB[ic][i]->GetXaxis()->FindBin(gMassYieldSideBand1.second - 1e-6);
         int ibinSideBandBegin2 = hNumInvMassMB[ic][i]->GetXaxis()->FindBin(gMassYieldSideBand2.first + 1e-6);
         int ibinSideBandEnd2 = hNumInvMassMB[ic][i]->GetXaxis()->FindBin(gMassYieldSideBand2.second - 1e-6);
         double SideBandRatio = (gMassYieldSideBand1.second - gMassYieldSideBand1.first + gMassYieldSideBand2.second - gMassYieldSideBand2.first) * 1.0 / (gMassYieldCountRange.second - gMassYieldCountRange.first);
         double SideBandBg = (hNumInvMassMB[ic][i]->Integral(ibinSideBandBegin1, ibinSideBandEnd1) + hNumInvMassMB[ic][i]->Integral(ibinSideBandBegin2, ibinSideBandEnd2)) / SideBandRatio;
         double SideBandBgErr = 0;
         for (int ijk = ibinSideBandBegin1; ijk <= ibinSideBandEnd1; ++ijk)
         {
            SideBandBgErr += TMath::Power(hNumInvMassMB[ic][i]->GetBinError(ijk), 2);
         }
         for (int ijk = ibinSideBandBegin2; ijk <= ibinSideBandEnd2; ++ijk)
         {
            SideBandBgErr += TMath::Power(hNumInvMassMB[ic][i]->GetBinError(ijk), 2);
         }

         fitSignalSub = hNumInvMassMB[ic][i]->Integral(ibintmp1, ibintmp2) - SideBandBg;
         SideBandBgErr *= TMath::Power(1.0 / SideBandRatio, 2);
         fitSignalSubErr += SideBandBgErr;
         fitSignalSubErr = sqrt(fitSignalSubErr);

         sprintf(buf, "#color[2]{subtract fit residual:} #D0 = %d #pm %d", fitSignalSub, fitSignalSubErr);
         drawLatex(0.2, 0.56, buf, 42, 0.04, 1);

         drawLine(mean - 3.0 * sigma, min, mean - 3.0 * sigma, min + (max - min) * 0.6, 2, 4);
         drawLine(mean + 3.0 * sigma, min, mean + 3.0 * sigma, min + (max - min) * 0.6, 2, 4);

         drawLines(x1, min, x2, max, 3, 1);

         sprintf(buf, "Au+Au 200GeV, %s MB", cent[ic]);
         drawLatex(0.16, 0.94, buf, 62, 0.045, 1);
         sprintf(buf, "%3.2f<p_{T}<%3.2f GeV/c", ptEdge[i], ptEdge[i + 1]);
         drawLatex(0.16, 0.88, buf, 42, 0.045, 1);

         sprintf(buf, "pic/D0_InvmSub_%s_cent_%d_%d_pt_%3.2f_%3.2f.gif", fileName, iCent[ic] - 1, iCent[ic + 1] - 2, ptEdge[i], ptEdge[i + 1]);
         // c1->SaveAs(buf);
         // End of Yield extraction


         ///From  Next is to Make two pannel Plots
         //
         gStyle->Reset("plain");
         TCanvas *c11 = new TCanvas("c11", "c11", 0, 0, 650, 850);
         gStyle->SetOptStat(0);
         gStyle->SetOptTitle(0);
         setpad(c11, 0.14, 0.02, 0.01, 0.12);
         c11->SetLogy(0);
         float small = 0;
         c11->Divide(1, 2, small, small);
         //Make tow panels
         c11->cd(1);
         c11_1->SetFillColor(-1);
         gPad->SetTickx();
         gPad->SetBottomMargin(small);
         hNumInvMassMBCopy[ic][i]->GetYaxis()->SetTitleOffset(0.8);
         hNumInvMassMBCopy[ic][i]->GetYaxis()->SetTitleSize(0.07);
         hNumInvMassMBCopy[ic][i]->Draw();
         hDenInvMassMBLikeCopy[ic][i]->Draw("same");
         hDenInvMassMBMixUnLikeCopy[ic][i]->Draw("same");
         drawLines(x1, y1, x2, y2, 3, 1);
         TLegend *leg = new TLegend(0.58, 0.65, 0.75, 0.94);
         leg->SetFillColor(10);
         leg->SetBorderSize(0);
         leg->SetTextFont(42);
         leg->SetTextSize(0.064);
         leg->AddEntry(hNumInvMassMBCopy[ic][i], "unlike-sign, same evt", "pl");
         leg->AddEntry(hDenInvMassMBLikeCopy[ic][i], "like-sign, same evt", "l");
         leg->AddEntry(hDenInvMassMBMixUnLikeCopy[ic][i], "unlike-sign, mixed evt", "l");
         leg->Draw("same");
         sprintf(buf, "Au+Au 200GeV, %s MB", cent[ic]);
         drawLatex(0.18, 0.94, buf, 62, 0.065, 1);
         sprintf(buf, "%3.2f<p_{T}<%3.2f GeV/c", ptEdge[i], ptEdge[i + 1]);
         drawLatex(0.18, 0.85, buf, 42, 0.065, 1);

         c11->cd(2);
         c11_2->SetFillColor(-1);
         gPad->SetTickx();
         gPad->SetTopMargin(small);
         // hNumInvMassMBCopy[ic][i]->GetYaxis()->SetTitleOffset(0.8);
         // hNumInvMassMBCopy[ic][i]->GetYaxis()->SetTitleSize(0.06);
         // hNumInvMassMBCopy[ic][i]->Draw();
         // hNumInvMassMBCopy[ic][i]->Write();
         hNumInvMassMB[ic][i]->GetYaxis()->SetTitleOffset(0.8);
         hNumInvMassMB[ic][i]->GetYaxis()->SetTitleSize(0.06);
         hNumInvMassMB[ic][i]->Draw();
         hNumInvMassMB[ic][i]->Write();
         myGausLine1->Draw("same");
         lineBg->Draw("same");
         fitExclusion->Draw("same");
         lineFun1->Draw("same");

         drawLines(x1, min, x2, max, 3, 1);
         sprintf(buf, "Au+Au 200GeV, %s MB", cent[ic]);
         drawLatex(0.18, 0.92, buf, 62, 0.06, 1);
         sprintf(buf, "%3.2f<p_{T}<%3.2f GeV/c", ptEdge[i], ptEdge[i + 1]);
         drawLatex(0.18, 0.83, buf, 42, 0.06, 1);
         sprintf(buf, "#color[2]{#D0 = %d #pm %d}", fitSignalSub, fitSignalSubErr);
         drawLatex(0.2, 0.56, buf, 42, 0.04, 1);

         sprintf(buf, "Gaus+pol1:");
         // drawLatex(0.65, 0.84, buf, 42, 0.05, 1);
         sprintf(buf, "#chi^{2}/ndf : %3.1f/%3.1f", chi2, ndf);
         drawLatex(0.65, 0.69, buf, 42, 0.05, 1);
         sprintf(buf, "mean  : %3.2f", mean);
         drawLatex(0.65, 0.63, buf, 42, 0.05, 1);
         sprintf(buf, "#sigma        : %3.3f", sigma);
         drawLatex(0.65, 0.58, buf, 42, 0.05, 1);

         sprintf(buf, "pic/D0_InvmSub_TwoPanel_%s_cent_%d_%d_pt_%3.2f_%3.2f.gif", fileName, iCent[ic] - 1, iCent[ic + 1] - 2, ptEdge[i], ptEdge[i + 1]);
         // if (ic == 3 || ic == 5) c11->SaveAs(buf); //only save 0-10% and 0-80% plots
         c11->SaveAs(buf); //only save 0-10% and 0-80% plots


         if (fitSignalSub >= 0 && fitSignalSubErr >= 0)
         {
            countCen[ic][i] = fitSignalSub;
            counterrCen[ic][i] = fitSignalSubErr;
            signifCen[ic][i] = fitSignalSub / fitSignalSubErr;
            signiferrCen[ic][i] = 0;

            //calculate dN/2pipTdpT/Events
            dNcountCen[ic][i] = countCen[ic][i] / 2. / TMath::Pi() / pt[i] / (2.*pterr[i]) / nEventsCent[ic] / 2.;
            dNcounterrCen[ic][i] = counterrCen[ic][i] / 2. / TMath::Pi() / pt[i] / (2.*pterr[i]) / nEventsCent[ic] / 2.;
            dN2countCen[ic][i] = dNcountCen[ic][i] / nNcollCent[ic];
            dN2counterrCen[ic][i] = dNcounterrCen[ic][i] / nNcollCent[ic];
         }
         else
         {
            countCen[ic][i] = 0;
            counterrCen[ic][i] = 0;
            signifCen[ic][i] = 0;
            signiferrCen[ic][i] = 0;
            //calculate dN/2pipTdpT/Events
            dNcountCen[ic][i] = countCen[ic][i] / 2. / TMath::Pi() / pt[i] / (2.*pterr[i]) / nEventsCent[ic] / 2.;
            dNcounterrCen[ic][i] = counterrCen[ic][i] / 2. / TMath::Pi() / pt[i] / (2.*pterr[i]) / nEventsCent[ic] / 2.;
            dN2countCen[ic][i] = dNcountCen[ic][i] / nNcollCent[ic];
            dN2counterrCen[ic][i] = dNcounterrCen[ic][i] / nNcollCent[ic];
         }

         outdata << cent[ic] << "  " << nEventsCent[ic] << "  " << nNcollCent[ic] << "  " << ndNcollCent[ic] << "  " << ptEdge[i] << "  " << ptEdge[i + 1] << "  " << countCen[ic][i] << "  " << counterrCen[ic][i] << "  " << dN2countCen[ic][i] << "  " << dN2counterrCen[ic][i] << endl;

      }


   }
   ////End of all the Yield extraction

   outdata.close();
   outFit.close();


//From next Next are draw D0 counts
//gStyle->Reset("plain");
   TCanvas *c2 = new TCanvas("c2", "c2", 0, 0, 800, 650);
   gStyle->SetOptStat(0);
   gStyle->SetEndErrorSize(2);
   gStyle->SetOptTitle(0);
   setpad(c2, 0.14, 0.02, 0.01, 0.12);
   c2->SetLogy(1);
//gStyle->SetOptFit(0);
   TH1F *hist = new TH1F("hist", ";p_{T}(Gev/c);D0 Counts", 12, 0, 12);
   hist->Draw();
   double xx1 = 0;
   double xx2 = 12;
   double yy1 = 1;
   double yy2 = 1e5;
   TAxis *a1 = hist->GetXaxis();
   a1->SetRangeUser(xx1, xx2);
   TAxis *a2 = hist->GetYaxis();
   a2->SetRangeUser(yy1, yy2);
   hist->GetXaxis()->SetTitleSize(0.05);
   hist->GetYaxis()->SetTitleSize(0.05);
   TGraphErrors *Graphcount[nCent];

   TLegend *leg = new TLegend(0.68, 0.53, 0.90, 0.84);
   leg->SetFillColor(10);
   leg->SetBorderSize(0);
   leg->SetTextFont(42);
   leg->SetTextSize(0.04);
   for (int ic = 0; ic < nCent; ic++)
   {
      // if (ic == 4) continue;
      if (iCent[ic + 1] <= iCent[ic]) continue;
      Graphcount[ic] = new TGraphErrors(nPtBins, pt, countCen[ic] , pterr, counterrCen[ic]);
      Graphcount[ic]->SetMarkerStyle(mStyle[ic]);
      Graphcount[ic]->SetMarkerColor(mColor[ic]);
      Graphcount[ic]->SetLineColor(mColor[ic]);
      Graphcount[ic]->SetMarkerSize(1.5);
      Graphcount[ic]->SetLineWidth(1.8);
      sprintf(buf, "count %s", cent[ic]);
      Graphcount[ic]->SetTitle(buf);
      sprintf(buf, "Graphcount_%d", ic);
      Graphcount[ic]->SetName(buf);
      Graphcount[ic]->GetXaxis()->SetTitle("p_{T}(Gev/c)");
      Graphcount[ic]->GetYaxis()->SetTitle("D0 Counts");
      if (ic == 0) Graphcount[ic]->Draw("P");
      else Graphcount[ic]->Draw("Psame");
      leg->AddEntry(Graphcount[ic], cent[ic], "pl");
   }
   leg->Draw("same");
   drawLines(xx1 , yy1, xx2, yy2, 3, 1);
   drawLatex(0.16, 0.92, "Au+Au 200GeV, D0 Raw Counts ", 62, 0.05, 1);

   sprintf(buf, "STAR Preliminary");
//      drawLatex(0.2,0.2,buf,42,0.05,1);

   sprintf(buf, "pic/Counts_vs_Pt_Diff_Centrality.gif");
   c2->SaveAs(buf);

//Next is draw D0 raw spectra
//Divide Ncoll
   TH1F *histd2 = new TH1F("histd2", ";p_{T}(Gev/c);dN/(N_{evt}*2#pip_{T}dp_{T}*dy*N_{coll}) [(GeV/c)^{-2}]", 12, 0, 12);
   histd2->Draw();
   double yy1d2 = 1e-14;
   double yy2d2 = 1e-6;
   TAxis *a1 = histd2->GetXaxis();
   a1->SetRangeUser(xx1, xx2);
   TAxis *a2 = histd2->GetYaxis();
   a2->SetRangeUser(yy1d2, yy2d2);
   histd2->GetXaxis()->SetTitleSize(0.05);
   histd2->GetYaxis()->SetTitleSize(0.05);
   TGraphErrors *Graphd2count[nCent];
   TLegend *leg = new TLegend(0.68, 0.53, 0.90, 0.84);
   leg->SetFillColor(10);
   leg->SetBorderSize(0);
   leg->SetTextFont(42);
   leg->SetTextSize(0.04);
   for (int ic = 0; ic < nCent; ic++)
   {
      // if (ic == 4) continue;
      if (iCent[ic + 1] <= iCent[ic]) continue;
      Graphd2count[ic] = new TGraphErrors(nPtBins, pt, dN2countCen[ic] , pterr, dN2counterrCen[ic]);
      Graphd2count[ic]->SetMarkerStyle(mStyle[ic]);
      Graphd2count[ic]->SetMarkerColor(mColor[ic]);
      Graphd2count[ic]->SetLineColor(mColor[ic]);
      Graphd2count[ic]->SetMarkerSize(1.5);
      Graphd2count[ic]->SetLineWidth(1.8);
      sprintf(buf, "Spectra %s", cent[ic]);
      Graphd2count[ic]->SetTitle(buf);
      sprintf(buf, "Graphd2Spectra_%d", ic);
      Graphd2count[ic]->SetName(buf);
      Graphd2count[ic]->GetXaxis()->SetTitle("p_{T}(Gev/c)");
      Graphd2count[ic]->GetYaxis()->SetTitle("dN/(N_{evt}*2#pip_{T}dp_{T}*dy*N_{coll}) [(GeV/c)^{-2}]");
      if (ic == 0) Graphd2count[ic]->Draw("P");
      else Graphd2count[ic]->Draw("Psame");
      sprintf(buf, "%s, N_{coll}=%3.1f", cent[ic], nNcollCent[ic]);
      leg->AddEntry(Graphd2count[ic], buf, "pl");
   }
   leg->Draw("same");
   drawLines(xx1 , yy1d2, xx2, yy2d2, 3, 1);
   drawLatex(0.16, 0.92, "Au+Au 200GeV, D0 dN/(N_{evt}*2#pip_{T}dp_{T}*N_{coll}) ", 62, 0.05, 1);

   sprintf(buf, "STAR Preliminary");
//      drawLatex(0.2,0.2,buf,42,0.05,1);
   sprintf(buf, "pic/dN_NColl_Counts_vs_Pt_Diff_Centrality.gif");
   c2->SaveAs(buf);

   TFile *out_file = new TFile("Run14_D0_Rcp_MixedEvent_OptmizedCut2_Yy_pT1.0.root", "recreate");
   for (int ic = 0; ic < nCent; ic++)
   {
      // if (ic == 4) continue;
      if (iCent[ic + 1] <= iCent[ic]) continue;
      Graphd2count[ic]->Write();
      Graphcount[ic]->Write();
   }
   out_file->Close();

}


//===========================================================================
TH1D* histo(char *name, Double_t xlow, Double_t xup, Double_t ylow, Double_t yup, char* xTitle, char* yTitle)
{
   TH1D *dd = new TH1D(name, "", 100, xlow, xup);
   dd->SetMinimum(ylow);
   dd->SetMaximum(yup);
   dd->GetXaxis()->SetTitle(xTitle);
   dd->GetYaxis()->SetTitle(yTitle);

   dd->GetXaxis()->SetTitleSize(0.055);
   dd->GetXaxis()->SetTitleOffset(0.9);
   dd->GetXaxis()->SetLabelSize(0.045);
   dd->GetYaxis()->SetTitleSize(0.055);
   dd->GetYaxis()->SetTitleOffset(1);
   dd->GetYaxis()->SetLabelSize(0.045);
   //dd->GetXaxis()->CenterTitle(kTRUE);
   dd->GetXaxis()->SetNdivisions(512);
   return dd;
}
TLatex* drawLatex(Double_t x, Double_t y, char* text, Int_t textFont, Double_t textSize, Int_t colorIndex)
{
   TLatex *latex = new TLatex(x, y, text);
   latex->SetNDC();
   latex->SetTextFont(textFont);
   latex->SetTextSize(textSize);
   latex->SetTextColor(colorIndex);
   latex->Draw("same");
   return latex;
}

TLine* drawLine(Double_t xlow, Double_t ylow, Double_t xup, Double_t yup, Double_t lineWidth, Int_t lineColor)
{
   TLine *l1 = new TLine(xlow, ylow, xup, yup);
   l1->SetLineWidth(lineWidth);
   l1->SetLineColor(lineColor);
   l1->Draw("same");
   return l1;
}

void drawLines(Double_t xlow, Double_t ylow, Double_t xup, Double_t yup, Double_t lineWidth, Int_t lineColor)
{
   drawLine(xlow, ylow, xup, ylow, lineWidth, lineColor);
   drawLine(xlow, yup, xup, yup, lineWidth, lineColor);
   drawLine(xlow, ylow, xlow, yup, lineWidth, lineColor);
   drawLine(xup, ylow, xup, yup, lineWidth, lineColor);
}

void setpad(TPad *pad, float left, float right, float top, float bottom)
{
   pad->SetFillColor(10);
   pad->SetBorderMode(0);
   pad->SetBorderSize(0);
   pad->SetFrameFillColor(10);
   pad->SetFrameBorderMode(0);
   pad->SetFrameBorderSize(0);
   pad->SetLeftMargin(left);
   pad->SetRightMargin(right);
   pad->SetTopMargin(top);
   pad->SetBottomMargin(bottom);
}







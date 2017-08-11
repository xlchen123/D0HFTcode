#include "../myFunction.h"
#include "../myConst.h"
#include "../anaCuts.h"
void plot_err() {
    globalSetting();
    char dir[250];
    char name[250];
    char title[250];
    //TString CMD = 0;
    char CMD[250];
    TLegend* legend;
    TH1F* h0;
    
    sprintf(dir,"pic");
    sprintf(CMD,"[ -d %s ] || mkdir -p %s",dir,dir);
    gSystem->Exec(CMD);
    
    // const int nerr = 6;
    // const char nameErr[nerr][250] = {"TPC track", "raw yield extract", "single track p_{T}", "50% topo. eff.", "150% topo. eff.", "double count"};
    //
    // const int nerr = 7;
    // const char nameErr[nerr][250] = {"TPC track", "raw yield extract", "single track p_{T}1","single track p_{T}2", "50% topo. eff.", "150% topo. eff.", "double count"};
    
    const int nerr = 5;
    const char nameErr[nerr][250] = {"TPC track", "raw yield extract", "single track p_{T}", "50%/150% topo. eff.", "double count"};
    
    float pt_mean[npt], pt_err[npt];
    for(int i=0; i<npt; i++) {
        pt_mean[i] = 0.5*(nptbin[i] + nptbin[i+1]);
        pt_err[i] = 0.5*(-nptbin[i] + nptbin[i+1]);
    }
    
    // read sys err
    ifstream in;
    float systmp, ytmp;
    float sys[ncent][nerr][npt];
    float y[npt], yerr[npt], ysys[npt];
    float tmp[npt];
    float tmp1[npt];
    float tmp2[npt];
    float tmp3[npt];
    for(int icent=0; icent<ncent; icent++) {
        int ierr = 0;
        
        //sys 1 -- tpc track
        for(int ipt=0; ipt<npt; ipt++) {
            sys[icent][ierr][ipt] = 0.072;
            //ysys[ipt] = sqrt(pow(0.04,2)+pow(ysys[ipt],2));
        }
        ierr++;
        
        //sys 2 
        //-- 2.1 count with fit
        in.open(Form("../count/data/yieldSys_%s.txt",nameCent1[icent]));
        for(int ipt=0; ipt<npt; ipt++) {
            if(in.eof()) break;
            in >> y[ipt] >> tmp1[ipt];
        }
        in.close();
        
        //-- 2.2 count with fit change fit range
        in.open(Form("../fitRange/data/yieldSys_%s.txt",nameCent1[icent]));
        for(int ipt=0; ipt<npt; ipt++) {
            if(in.eof()) break;
            in >> y[ipt] >> tmp2[ipt];
        }
        in.close();

        //-- 2.3 count with likesign bkg subtraction
        in.open(Form("../likeSign/data/yieldSys_%s.txt",nameCent1[icent]));
        for(int ipt=0; ipt<npt; ipt++) {
            if(in.eof()) break;
            in >> y[ipt] >> tmp3[ipt];
        }
        in.close();

        for(int ipt=0; ipt<npt; ipt++) {
          tmp[ipt] = tmp1[ipt] > tmp2[ipt] ? tmp1[ipt] : tmp2[ipt];
          tmp[ipt] = tmp[ipt] > tmp3[ipt] ? tmp[ipt] : tmp3[ipt];
          // ysys[ipt] = sqrt(pow(tmp[ipt],2)+pow(ysys[ipt],2));
        }

        for(int ipt=0; ipt<npt; ipt++) {
          sys[icent][ierr][ipt] = tmp[ipt];
        }
        ierr++;

        //sys 3, Daughter pt Cut scan // choose the maximum difference
        //sys 3.1 -- daughter pt cut1
        in.open(Form("../ptCut1/data/yieldSys_%s.txt",nameCent1[icent]));
        for(int ipt=0; ipt<npt; ipt++) {
            if(in.eof()) break;
            in >> y[ipt] >> tmp1[ipt];
        }
        in.close();
        
        //sys 3.2 -- daughter pt cut2
        in.open(Form("../ptCut2/data/yieldSys_%s.txt",nameCent1[icent]));
        for(int ipt=0; ipt<npt; ipt++) {
            if(in.eof()) break;
            in >> y[ipt] >> tmp2[ipt];
        }
        in.close();

        for(int ipt=0; ipt<npt; ipt++) {
          tmp[ipt] = tmp1[ipt] > tmp2[ipt] ? tmp1[ipt] : tmp2[ipt];
          // ysys[ipt] = sqrt(pow(tmp[ipt],2)+pow(ysys[ipt],2));
        }

        for(int ipt=0; ipt<npt; ipt++) {
          sys[icent][ierr][ipt] = tmp[ipt];
        }
        ierr++;

        
        //sys 4 topological cut
        //4.1 -- tight topo cuts
        in.open(Form("../topoCut1/data/yieldSys_%s.txt",nameCent1[icent]));
        for(int ipt=0; ipt<npt; ipt++) {
            if(in.eof()) break;
            in >> y[ipt] >> tmp1[ipt];
            if(icent == 4 && ipt >= npt-3) tmp1[ipt] =0; // for the tight cut, the signal is not good, no need to include this sys for some certain pt bin
        }
        in.close();
        
        //4.2 -- loose topo cuts
        in.open(Form("../topoCut2/data/yieldSys_%s.txt",nameCent1[icent]));
        for(int ipt=0; ipt<npt; ipt++) {
            if(in.eof()) break;
            in >> y[ipt] >> tmp2[ipt];
        }
        in.close();
        
        for(int ipt=0; ipt<npt; ipt++) {
          tmp[ipt] = tmp1[ipt] > tmp2[ipt] ? tmp1[ipt] : tmp2[ipt];
          // ysys[ipt] = sqrt(pow(tmp[ipt],2)+pow(ysys[ipt],2));
        }
        
        for(int ipt=0; ipt<npt; ipt++) {
          sys[icent][ierr][ipt] = tmp[ipt];
        }
        ierr++;

        //sys 5 -- double count
        in.open(Form("../DoubleCount/data/yieldSys_%s.txt",nameCent1[icent]));
        for(int ipt=0; ipt<npt; ipt++) {
            if(in.eof()) break;
            float errtmp;
            float tmptmp;
            in >> ytmp >> errtmp >> systmp >> tmptmp;
            sys[icent][ierr][ipt] = systmp;
        }
        in.close();
        ierr++;
    }
    
    //define sys err graph
    TGraph* gSys[ncent][nerr];
    for(int icent=0; icent<ncent; icent++) {
        for(int ierr=0; ierr<nerr; ierr++) {
            gSys[icent][ierr] = new TGraph(npt,pt_mean,sys[icent][ierr]);
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
        sprintf(name,"%s/sysErr_%s.pdf",dir,nameCent1[icent]);
        c1->SaveAs(name);
    }
}

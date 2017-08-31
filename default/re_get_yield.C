#include "../anaCuts.h"
void re_get_yield()  //re calculate the default mean value, especially for some certain centrality/pt bin
{
   ifstream in;
   ofstream out;
   const int nyield = 4;

   float y1[npt], y1err[npt];//default fit
   float y2[npt], y2err[npt];//count
   float y3[npt], y3err[npt];//fitRange
   float y4[npt], y4err[npt];//likeSign

   float y[npt], yerr[npt];//re sign

   for (int icent = 0; icent < ncent; icent++)
   {
      in.open(Form("../default/data/yield_%s.txt", nameCent1[icent]));
      for (int ipt = 0; ipt < npt; ipt++)
         in >> y1[ipt] >> y1err[ipt];
      in.close();

      in.open(Form("../count/data/yield_%s.txt", nameCent1[icent]));
      for (int ipt = 0; ipt < npt; ipt++)
         in >> y2[ipt] >> y2err[ipt];
      in.close();

      in.open(Form("../fitRange/data/yield_%s.txt", nameCent1[icent]));
      for (int ipt = 0; ipt < npt; ipt++)
         in >> y3[ipt] >> y3err[ipt];
      in.close();

      in.open(Form("../likeSign/data/yield_%s.txt", nameCent1[icent]));
      for (int ipt = 0; ipt < npt; ipt++)
         in >> y4[ipt] >> y4err[ipt];
      in.close();


      //re sign some specific data point
      for (int ipt = 0; ipt < npt; ipt++)
      {
        y[ipt] = y1[ipt];
        yerr[ipt] = y1err[ipt];
        if( icent == 1 && ipt == 0)  //10-20%, 0-0.5GeV data point // average all the data point
        {
          y[ipt] = 1./3 * (y1[ipt] + y2[ipt] + y3[ipt]);
          yerr[ipt] = y[ipt] * 1./3. * (y1err[ipt]*1./y1[ipt] + y2err[ipt]*1./y2[ipt]+ y3err[ipt]*1./y3[ipt]);
        }
        if( icent ==0  && ipt == 0)  //0-10%, 0-0.5GeV data point // average 3 of the data points, sicne the last signal from likesign is less then 3sigma
        {
          y[ipt] = 1./2 * (y1[ipt] + y2[ipt]);
          yerr[ipt] = y[ipt] * 1./2. * (y1err[ipt]*1./y1[ipt] + y2err[ipt]*1./y2[ipt]);
        }
        if(icent == 4 && ipt > npt -3 )  //60-80% 5-6-8GeV data point// 
        {
          y[ipt] = 1./3 * (y1[ipt] + y2[ipt] +  y4[ipt]);
          yerr[ipt] = y[ipt] * 1./3. * (y1err[ipt]*1./y1[ipt] + y2err[ipt]*1./y2[ipt]+ y4err[ipt]*1./y4[ipt]);
        }
        if(icent == 4 && ipt == npt -3 )  //60-80% 4-5GeV data point// 
        {
          y[ipt] = y2[ipt];// fit didn't catch the raw counts
          yerr[ipt] = y2err[ipt]*1.;
        }
      }

      out.open(Form("../default/data/re_yield_%s.txt", nameCent1[icent]));
      for (int ipt = 0; ipt < npt; ipt++)
         out << y[ipt] << "\t" << yerr[ipt] << endl;
      out.close();
   } 

}

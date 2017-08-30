//part 1: N binary
const float Nbin[9] = {12.91047, 29.32408, 62.25545, 120.65795, 215.95050, 360.58912, 579.89409, 812.73278, 1042.75372};

//part 2: vertex reso. correction factor
//const float vtxWeight[9] = {0.29, 0.55, 0.77, 0.89, 0.93, 1., 1., 1., 1.};
//const float vtxWeight[9] = {1, 1, 1, 1, 1, 1., 1., 1., 1.};
const int ncent_vtx = 6;
const char nameCent_vtx[ncent_vtx][100] = {"70-80%","60-70%","50-60%","40-50%","30-40%","20-30%"};
const char nameCent_vtx1[ncent_vtx][100] = {"70_80","60_70","50_60","40_50","30_40","20_30"};
const int ncentAll_vtx = 9;
const char nameCentAll_vtx[ncentAll_vtx][100] = {"70-80%","60-70%","50-60%","40-50%","30-40%","20-30%","10-20%","5-10%","0-5%"};
const char nameCentAll_vtx1[ncentAll_vtx][100] = {"70_80","60_70","50_60","40_50","30_40","20_30","10_20","5_0","0_5"};

//part 3: tof match
//int cent_lw, cent_up;
const int ncent_tof = 9;
float centLw_tof[ncent_tof] = {1,2,3,4,5,6,7,8,9}; // >=   0-80%: 1-9
float centUp_tof[ncent_tof] = {2,3,4,5,6,7,8,9,10}; // <
char nameCent_tof[ncent_tof][250] = {"70-80%","60-70%","50-60%","40-50%","30-40%","20-30%","10-20%","5-10%","0-5%"};
float NbinSum_tof[ncent_tof] = {0};
const int npt_tof = 200;
float ptLw_tof = 0;
float ptUp_tof = 10;

//part 4: raw spectra -- centrality and pT
//part 4a: centrality
const int ncent = 10;
char nameCent[ncent][250] = {"0-10%","10-20%","20-40%","40-60%","60-80%","10-40%","40-80%","0-80%","10-80%","30-50%"};
char nameCent1[ncent][250] = {"0_10","10_20","20_40","40_60","60_80","10_40","40_80","0_80","10_80","30_50"};
float centLw[ncent] = {7,6,4,2,0,4,0,0,0,3}; // >=    0-80%: 1-9
float centUp[ncent] = {9,7,6,4,2,7,4,9,7,5}; // <
float NbinMean[ncent] = {938.80170, 579.89409, 288.35051, 91.37100, 21.37396, 386.08527, 56.99229, 301.05848,197, 168};
float Nbin_Sum[ncent] = {0};

// const int ncent = 4;
// char nameCent[ncent][250] = {"0-10%","10-40%","40-80%","0-80%"};
// char nameCent1[ncent][250] = {"0_10","10_40","40_80","0_80"};
// float centLw[ncent] = {7,4,0,0}; // >=    0-80%: 1-9
// float centUp[ncent] = {9,7,4,9}; // <
// float NbinMean[ncent] = {938.80170, 386.08527, 56.99229, 301.05848};
// float Nbin_Sum[ncent] = {0};

//part 4b: pT
//const int npt = 2;
//double nptbin[npt+1] = {3, 5, 8};
//const int npt = 10;
//double nptbin[npt+1] = {0, 1., 1.5, 2., 2.5, 3., 3.5, 4.0, 5., 6., 8.0};
//const int npt = 5;
//double nptbin[npt+1] = {0, 1, 2, 3, 5, 8.0};
//const int npt = 7;
//double nptbin[npt+1] = {0, 0.7, 1.1, 1.6, 2.2, 3.0, 5.0, 8.0};
// const int npt = 1;
// double nptbin[npt+1] = {1.5, 5.0};
//const int npt = 8;
//double nptbin[npt+1] = {0, 1, 2, 3, 4, 5, 6, 8,  10.};
// const int npt = 10;
// double nptbin[npt+1] = {0, 0.5, 1., 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8};
// const int npt = 13;
// double nptbin[npt+1] = { 0., 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 7.0, 8.0, 10 };

//const int npt = 9;
//const double nptbin[npt+1] = { 0., 0.5, 1., 1.5, 2., 2.5, 3., 4.0, 5.5, 8.0 };//
const int npt = 10;
const double nptbin[npt+1] = { 0., 0.5, 1., 1.5, 2., 2.5, 3.0, 4.0, 5.0, 6.0, 8.0 };//

//const int npt = 4;
//double nptbin[npt+1] = {1.5, 2.5, 3.5, 5, 8};

//part 5: eff
const int npt_eff = 200;
float ptLw_eff = 0;
float ptUp_eff = 10;

//spectra: ub
//float NbinMean[ncent] = {938.80170, 768.78266, 301.05848};

//part 6: double count
const int ncentBase_dc = 5; //base cent divide
char nameCentBase_dc[ncentBase_dc][250] = {"0-10%","10-20%","20-40%","40-60%","60-80%"}; //base cent divide
float NbinMeanBase_dc[ncentBase_dc] = {938.80170, 579.89409, 288.35051, 91.37100, 21.37396};

float centLw_dc[ncent] = {0,1,2,3,4,1,3,0,1,2}; // >= //keep consistent with centrality divide in raw spectra
float centUp_dc[ncent] = {1,2,3,4,5,3,5,5,5,4}; // <

// float centLw_dc[ncent] = {0,1,3,0}; // >= //keep consistent with centrality divide in raw spectra
// float centUp_dc[ncent] = {1,3,5,5}; // <

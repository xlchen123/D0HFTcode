# $D^0$ HFT analysis code
LBNL - Nuclear Science Division, RNC Soft Hadron Physics Group

USTC - Modern Physics Department, Office 410

## Declaration
This code was original for Our STAR HFT D0 analysis

All Rights Reserved !

###Code Authors:  
[Xiaolong Chen](https://github.com/xlchen123) (xlchen@lbl.gov)  
[Guannan Xie](https://github.com/GuannanXie) (xieguannanpp@gmail.com)  
- - -
### Presentations:  
#### STAR Protected:  
first version of root file without charge and Y 
/global/homes/x/xgn1992/rnc_Global/SL16dRun14/myAnalysis_AnalysisMeeting_Nov2/MixedEvent_hybrid/hadd/Mixed_Evt.2017Mar23.root

second version of root file with charge and Y 
/global/homes/x/xgn1992/rnc_Global/SL16dRun14/myAnalysis_AnalysisMeeting_Nov2/MixedEvent_hybrid_V2/hadd/Mixed_Evt.2017Aug16.root

- - -

###How to use this code:  
```bash
NOTE: the mixdata root file "D0_data_mix.root" was not here since the limit of Github file size

1. How to run (after copy the full files), if run on PDSF/RCF: starver SL16d ==> simplify printed things...
./run.sh

2. All can be set at anaCuts.h
2.1. Vary centrality division at part4a and part6
2.2. Vary pT division at part4b
2.3. Vary TOF match eff. centrality binning at part3
2.4. Vary the double ratio efficiency between pure hijing and fast simulation at part2

3. spectra before pT shift at: sys/D0_Spectra_Run14_HFT_beforePtShift.root 
   spectra after pT shift at: ptShift/D0_Spectra_Run14HFT.root
   single sysErr in sys/pic


4. RAA need check whether the pT shift is right in ptShift
   Or just set in ./run.sh to get integral RAA at each pT bin 

5. Rcp in Rcp/plot_Rcp.C, change the baseline: nameCent[ncent-1] && nameCent1[ncent-1] && NbinMean[ncent-1]
   Rcp/plot_Rcp.C : use pT spectra without pT shift
   Rcp/plot_Rcp_pTshift.C : use pT spectra after pT shift

```

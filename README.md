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
<!-- 1. [sPHENXIX tracker (steal from Tony)](https://github.com/GuannanXie/sPHENIX_FastSimu/blob/master/Slides/2016Sept_sPHENIX_tracking_simulations_Tony.pdf), Tony Frawley, 2016&#45;09&#45;07  -->

- - -

###How to use this code:  
```bash

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

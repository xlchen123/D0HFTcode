#!/bin/bash
root -b -q plot_effBinning.C

root -b -q plot_statErr.C
root -b -q plot_sysErr.C
root -b -q plot_sysLooseErr.C
root -b -q plot_sysTightErr.C

root -b -q write_vtxCorr_default.C
root -b -q write_vtxCorr_loose.C
root -b -q write_vtxCorr_tight.C

root -b -q write_vtxStat.C

root -b -q get_sys.C
root -b -q get_stat.C

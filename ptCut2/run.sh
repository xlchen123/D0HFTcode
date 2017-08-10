#!/bin/bash
root -b -q write_rawy.C
root -b -q plot_rawy.C
root -b -q write_eff.C
root -b -q plot_eff.C
root -b -q get_yield.C
root -b -q get_sys.C

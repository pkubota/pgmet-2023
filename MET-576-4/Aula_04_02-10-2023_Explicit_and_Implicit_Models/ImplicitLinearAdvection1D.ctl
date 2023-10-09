dset ^ImplicitLinearAdvection1D.bin
title  EDO
undef  -9999.9
xdef       200 linear 0.00 0.001
ydef  1 linear  -1.27 1
tdef     400 linear  00z01jan0001 1hr
zdef  1 levels 1000 
vars 2
FA  0 99 resultado analitico da edol yc
FC  0 99 resultado numerico  da edol yc
endvars

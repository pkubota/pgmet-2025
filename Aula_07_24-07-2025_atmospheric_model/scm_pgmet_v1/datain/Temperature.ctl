dset ^Temperature.bin
*OPTIONS  big_endian          
undef  -99999.0
title GPCP MODEL 1.0  degree
xdef   161 linear   -65.0   0.25
ydef   161 linear   -50.0   0.25
tdef   168  LINEAR  00Z18APR2023  1hr
zdef    37 levels 1000 975 950 925 900 875 850 825 800 775 
                   750 700 650 600 550 500 450 400 350 300 
                   250 225 200 175 150 125 100  70  50  30 
                    20  10   7   5   3   2   1             
VARS 1
temp    37   99    Temperature                 (K       )    
ENDVARS

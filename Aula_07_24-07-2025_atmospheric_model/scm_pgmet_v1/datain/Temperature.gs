'reinit'
'set display color white'
'c'
'sdfopen Temperature.nc'
DirName='./'
FileName='Temperature'
'rm -f 'DirName''FileName''

'set fwrite  'DirName''FileName'.bin'
'set gxout fwrite'

it=409
'set t 'it''
'q time'
time0=subwrd(result,3)
while(it<=577)
  iz=1
    while(iz<=37)
      'set y 41  201'
      'set x 221 381'
      'd t(t='it',z='iz')'
      iz=iz+1 
    endwhile
  it=it+1
endwhile
'disable fwrite'

******************************************************************************************************************************
*    ctl
******************************************************************************************************************************
filename=DirName'/'FileName
filename2=FileName'.bin'
say filename
say filename2
******************************************************************************************************************************
*    ctl
******************************************************************************************************************************
     'set t 1'
     rc=write(filename'.ctl','dset ^'filename2 )
     rc=write(filename'.ctl','*OPTIONS  big_endian          '                                                ,append )
     rc=write(filename'.ctl','undef  -99999.0'                                                               ,append )
     rc=write(filename'.ctl','title GPCP MODEL 1.0  degree'                                                   ,append )
     rc=write(filename'.ctl','xdef   161 linear   -65.0   0.25'                                              ,append )
     rc=write(filename'.ctl','ydef   161 linear   -50.0   0.25'                                              ,append )
     rc=write(filename'.ctl','tdef   168  LINEAR  'time0'  1hr'                                                    ,append )
     rc=write(filename'.ctl','zdef    37 levels 1000 975 950 925 900 875 850 825 800 775 '       ,append )
     rc=write(filename'.ctl','                   750 700 650 600 550 500 450 400 350 300 '       ,append )
     rc=write(filename'.ctl','                   250 225 200 175 150 125 100  70  50  30 '       ,append )
     rc=write(filename'.ctl','                    20  10   7   5   3   2   1             '       ,append )
     rc=write(filename'.ctl','VARS 1'                                                                                 ,append )
     rc=write(filename'.ctl','temp    37   99    Temperature                 (K       )    '       ,append )
     rc=write(filename'.ctl','ENDVARS'                                                                                ,append )
******************************************************************************************************************************

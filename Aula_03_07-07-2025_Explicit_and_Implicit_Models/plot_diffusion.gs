'reinit'
'set display color white'
'c'
*'open ImplicitLinearDiffusion1D.ctl'
'open ImplicitLinearDiffusion1DMethodCrankNicolson_NWP.ctl'
it=0
while(it<=1000)
'c'
'set vrange -2 2'
'd a(t=1)'
'd a(t='it')'
if(it=999)
  'q pos'  
  it=0
endif
it=it+1
endwhile

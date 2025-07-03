'reinit'
'set display color white'
'c'
'open AdvecLinearConceitual1D.ctl '

it=1
while(it<=200)

'd phia(t='it')'
'd phic(t='it')'
'!sleep 1'
'c'
it=it+1
endwhile

function V = integrate_eLIF(x,I,V0,Treset)

for i,Iin in enumerate(I[1:]):
    Vmodel[i+1] = Vmodel[i] + h*(func(Vmodel[i],*popt2)+Iin/(C/10.0))
    if Vmodel[i+1] > 10:
        Vmodel[i+1] = popt2[2]
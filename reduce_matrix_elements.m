
Sz0 = sum(mod(abs(B0(1,:)),2))/2; % Sz in the neutral ground state
S0 = round(10*E0(1,2))/10;
n=1
for m = 1:2

for k = 2:4

%Lehmann_amplitudes(C0(:,n),CM(:,m),B0,BM,1,1)/clebsch_gordan(0.5, S0, round(10*EM(m,2))/10,  0.5, Sz0, Sz0 + 0.5)

Lehmann_amplitudes(C0(:,n),CM(:,m),B0,BM,k,-1)/clebsch_gordan(0.5, S0, round(10*EM(m,2))/10, -0.5, Sz0, Sz0 - 0.5)

Lehmann_amplitudes(Cm(:,m),C0(:,n),Bm,B0,k,1)/clebsch_gordan(0.5,round(10*Em(m,2))/10,  S0,  0.5, Sz0 - 0.5, Sz0)

Lehmann_amplitudes(Cm(:,m),C0(:,n),Bm,B0,k,-1)/clebsch_gordan(0.5,round(10*Em(m,2))/10,  S0,  0.5, Sz0 + 0.5, Sz0)

end % k

end % m

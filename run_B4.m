
G= load('B4bas.txt'); G = G(:,[2,3])

[E_st,C,E,AE,E0,E1,C0,C1,B0,B1,H0,H1,n0u,n0d,n1u,n1d,gap_oe,Cnat,Enat,Rho,Fm,FM,CM,BM,Cm,Bm,EM,Em] = Huckel_Hubbard(G,-2.8,4.3,4);

%[E_st,C,E,AE,E00,E1,C0,C1,B0,B1,H0,H1,n0u,n0d,n1u,n1d,gap_oe,Cnat,Enat,Rho,Fm,FM,CM,BM,Cm,Bm,EM0,Em0] = Huckel_Hubbard(G,-2.8,4.3,4);

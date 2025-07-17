function[D,Ni] = unpaired_electrons(E)

E=E; E=E(1:end);
D = E.*(2-E);
D = 1-abs(1-E);
D = (E.*(2-E)).^2;
Ni = sum(D);

end %end function
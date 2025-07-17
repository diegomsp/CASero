function[Sx,Sy,Sz] = Pauli_matrices_Hfile(site)

Sx = [ 0.5, site, 0, 0, -site; 0.5, -site, 0, 0, site ];

Sy = [ -0.5*1i, site, 0, 0, -site; 0.5*1i, -site, 0, 0, site ];

Sz = [ 0.5, site, 0, 0, site; -0.5, -site, 0, 0, -site ];

end %end function Pauli_matrices_Hfile
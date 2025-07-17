function[Sx,Sy,Sz] = Pauli_matrices_Hfile_nondiag(site1,site2)

Sx = [ 0.5, site1, 0, 0, -site2; 0.5, -site1, 0, 0, site2 ];

Sy = [ -0.5*1i, site1, 0, 0, -site2; 0.5*1i, -site1, 0, 0, site2 ];

Sz = [ 0.5, site1, 0, 0, site2; -0.5, -site1, 0, 0, -site2 ];

end %end function Pauli_matrices_Hfile

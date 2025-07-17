function[Cij,Sij,Sxi,Syi,Szi] = spin_correl_in_extended_basis(C,B,nstate)

%CALCULATE SPIN CORRELATION < Psi | Si Sj | Psi  > - < Psi | Si| Psi  > < Psi | Sj | Psi  >,
% the wavefucntion Psi is expanded in basis of Slater determinants given by B.

%i j run over list StatesC
Nst = numel(B(1,:)); %number of sites
Nel = sum(abs(B(1,:)));

% extend the basis
[Badd]=construct_basis2(Nst,Nel,2,[-2,2,4,-4,1,-1,3,-3],0);
Btotal = [B; Badd];
Psi = [C(:,nstate); zeros(size(Badd,1),1)];

for isite = 1:Nst
   [Sxhfile,Syhfile,Szhfile] = Pauli_matrices_Hfile(isite);
   [Sx(:,:,isite)] = full(create_hamiltonian_sparse_matrix_general_nonreduced(Sxhfile,[],Btotal,[],[],0,Nst));
   [Sy(:,:,isite)] = full(create_hamiltonian_sparse_matrix_general_nonreduced(Syhfile,[],Btotal,[],[],0,Nst));
   [Sz(:,:,isite)] = full(create_hamiltonian_sparse_matrix_general_nonreduced(Szhfile,[],Btotal,[],[],0,Nst));
end % isite

for isite1 = 1:Nst
  Sxi(isite1) = (Psi')*((Sx(:,:,isite1))*Psi);
  Syi(isite1) = (Psi')*((Sy(:,:,isite1))*Psi);
  Szi(isite1) = (Psi')*((Sz(:,:,isite1))*Psi);
  for isite2 = isite1:Nst
    Sij(isite1,isite2) = (Psi')*((Sx(:,:,isite1)*Sx(:,:,isite2) + Sy(:,:,isite1)*Sy(:,:,isite2) + Sz(:,:,isite1)*Sz(:,:,isite2))*Psi);
    Sij(isite2,isite1) =  Sij(isite1,isite2);
  end % isite2
end % isite1

for isite1 = 1:Nst
  for isite2 = isite1:Nst
    Cij(isite1,isite2) = Sij(isite1,isite2) - ( Sxi(isite1)*Sxi(isite2) + Syi(isite1)*Syi(isite2) + Szi(isite1)*Szi(isite2) );
    Cij(isite2,isite1) =  Cij(isite1,isite2);
  end % isite2
end % isite1



end %end function spin_correl


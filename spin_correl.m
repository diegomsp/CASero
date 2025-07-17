function[Cij,Sij,Sxi,Syi,Szi] = spin_correl(C,B,nstate)

%CALCULATE SPIN CORRELATION < C(:,nstate) | Si Sj | C(:,nstate)  > - < C(:,nstate) | Si| C(:,nstate)  > < C(:,nstate) | Sj | C(:,nstate)  >,
% the wavefucntion C(:,nstate) is expanded in basis of Slater determinants given by B.

%i j run over list StatesC
Nst = numel(B(1,:)); %number of sites

for isite = 1:Nst
   [Sxhfile,Syhfile,Szhfile] = Pauli_matrices_Hfile(isite);
   [Sx(:,:,isite)] = create_hamiltonian_sparse_matrix_general_nonreduced(Sxhfile,[],B,[],[],0,Nst);
   [Sy(:,:,isite)] = create_hamiltonian_sparse_matrix_general_nonreduced(Syhfile,[],B,[],[],0,Nst);
   [Sz(:,:,isite)] = create_hamiltonian_sparse_matrix_general_nonreduced(Szhfile,[],B,[],[],0,Nst);
end % isite

for isite1 = 1:Nst
  Sxi(isite1) = (C(:,nstate)')*((Sx(:,:,isite1))*C(:,nstate));
  Syi(isite1) = (C(:,nstate)')*((Sy(:,:,isite1))*C(:,nstate));
  Szi(isite1) = (C(:,nstate)')*((Sz(:,:,isite1))*C(:,nstate));
  for isite2 = isite1:Nsites
    Sij(isite1,isite2) = (C(:,nstate)')*((Sx(:,:,isite1)*Sx(:,:,isite2) + Sy(:,:,isite1)*Sy(:,:,isite2) + Sz(:,:,isite1)*Sz(:,:,isite2))*C(:,nstate));
  end % isite2
end % isite1

for isite1 = 1:Nst
  for isite2 = isite1:Nsites
    Cij(isite1,isite2) = Sij(isite1,isite2) - ( Sxi(isite1)*Sxi(isite2) + Syi(isite1)*Syi(isite2) + Szi(isite1)*Szi(isite2) );
  end % isite2
end % isite1



end %end function spin_correl


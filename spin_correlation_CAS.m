function[Smunu] = spin_correlation_CAS(C,B,nstate,AE,mu,nu)

%i j run over list StatesC

Px = 0.5*[ 0 1; 1 0];
Nst = numel(B(1,:)); %number of sites
Nel = sum(abs(B(1,:)));
Py = 0.5*[0 -1i; 1i 0];
Pz = 0.5*[1 0; 0 -1];

Smunu = 0;
for sigma = 1:2
 for sigmap = 1:2
  for tau = 1:2
   for taup = 1:2
    if ( (-2*sigma + 3) + (-2*sigmap + 3) + (-2*tau + 3) + (-2*taup + 3) == 0 )
    for j = 1:Nst
     for jp = 1:Nst
      for k = 1:Nst
       for kp = 1:Nst
         Smunu = Smunu + ...
         ( Px(sigma,sigmap)*Px(tau,taup) + Py(sigma,sigmap)*Py(tau,taup) + Pz(sigma,sigmap)*Pz(tau,taup) )*AE(mu,j)*AE(mu,jp)*AE(nu,k)*AE(nu,kp)...
         *scalar_product([C(:,nstate),B],  apply_hamiltonian_2body(...
         [C(:,nstate),B],[1, j*(-2*sigma + 3), jp*(-2*sigmap + 3), k*(-2*tau + 3), kp*(-2*taup + 3), ],Nst) );
       end % kp
      end % k
     end % Nst
    end % Nst
   end % if
   end % taup
  end % tau
 end % sigmap
end % sigma



end %end function spin_correl

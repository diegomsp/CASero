function[Jmat]=calculate_J_matrix_Kondo_orbitals(C0,B0,E0,CM,BM,EM,Cm,Bm,Em,varargin)

Nstates0 = 1;
Nstates_charged = size(BM,1);
Nst = size(B0,2);

if (length(varargin)>0)
  Nstates0 = varargin{1};
  Nstates_charged = varargin{2};
end

% Determine S0
tol = 0.001;
Sz0 = sum(mod(abs(B0(1,:)),2))/2; % Sz in the neutral ground state
S0 = round(10*E0(1,2))/10;
Jmat = zeros(size(B0,2));
for n = 1:Nstates0 % loop over states of the neutral system
    for m = 1:Nstates_charged % loop over states of the charged system
      if ( abs(EM(m,2) - E0(1,2) - 0.5) < tol   ) CoefM = -(S0 + 0.5);  end
      if ( abs(EM(m,2) - E0(1,2) + 0.5) < tol   ) CoefM =  (S0 + 0.5);  end
      if ( abs(Em(m,2) - E0(1,2) - 0.5) < tol   ) Coefm = -(S0 + 1.0);  end
      if ( abs(Em(m,2) - E0(1,2) + 0.5) < tol   ) Coefm =   S0       ;  end
      cgMu = clebsch_gordan(0.5, S0, round(10*EM(m,2))/10,  0.5, Sz0, Sz0 + 0.5);
      cgMd = clebsch_gordan(0.5, S0, round(10*EM(m,2))/10, -0.5, Sz0, Sz0 - 0.5);
      cgmu = clebsch_gordan(0.5,round(10*Em(m,2))/10,  S0,  0.5, Sz0 - 0.5, Sz0);
      cgmd = clebsch_gordan(0.5,round(10*Em(m,2))/10,  S0,  0.5, Sz0 + 0.5, Sz0);
        for j = 1:Nst   % loop over sites of active space
           for k = 1:Nst   % loop over sites of active space

            if (abs(cgMu) > 0.001) Jmat(j,k) = Jmat(j,k) + (Lehmann_amplitudes(C0(:,n),CM(:,m),B0,BM,k,1) *Lehmann_amplitudes(C0(:,n),CM(:,m),B0,BM,j,1) /cgMu^2)/( CoefM*(EM(m,1) - E0(1,1)) ); end
            if (abs(cgMd) > 0.001) Jmat(j,k) = Jmat(j,k) + (Lehmann_amplitudes(C0(:,n),CM(:,m),B0,BM,k,-1)*Lehmann_amplitudes(C0(:,n),CM(:,m),B0,BM,j,-1)/cgMd^2)/( CoefM*(EM(m,1) - E0(1,1)) ); end
            if (abs(cgmu) > 0.001) Jmat(j,k) = Jmat(j,k) + (Lehmann_amplitudes(Cm(:,m),C0(:,n),Bm,B0,k,1) *Lehmann_amplitudes(Cm(:,m),C0(:,n),Bm,B0,j,1) /cgmu^2)/( Coefm*(Em(m,1) - E0(1,1)) ); end
            if (abs(cgmd) > 0.001) Jmat(j,k) = Jmat(j,k) + (Lehmann_amplitudes(Cm(:,m),C0(:,n),Bm,B0,k,-1)*Lehmann_amplitudes(Cm(:,m),C0(:,n),Bm,B0,j,-1)/cgmd^2)/( Coefm*(Em(m,1) - E0(1,1)) ); end

            %Jmat(j,k) = Jmat(j,k) - (Lehmann_amplitudes(C0(:,n),Cm(:,m),B0,Bm,k,1) *Lehmann_amplitudes(C0(:,n),Cm(:,m),B0,Bm,j,1) /clebsch_gordan(0.5, S0, round(10*Em(m,2))/10,  0.5, Sz0, Sz0 + 0.5)^2)/( Coefm*(Em(m,1)) - E0(1,1) );
            %Jmat(j,k) = Jmat(j,k) - (Lehmann_amplitudes(C0(:,n),Cm(:,m),B0,Bm,k,-1)*Lehmann_amplitudes(C0(:,n),Cm(:,m),B0,Bm,j,-1)/clebsch_gordan(0.5, S0, round(10*Em(m,2))/10, -0.5, Sz0, Sz0 + 0.5)^2)/( Coefm*(Em(m,1)) - E0(1,1) );

           end % k
        end % j
    end % m
end % n




end % end function calculate_J_matrix_Kondo_orbitals

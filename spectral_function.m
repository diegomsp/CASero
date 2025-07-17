function[DOS] = spectral_function(C,E,B,CM,EM,BM,Cm,Em,Bm,sites,eta,w,T,States_DOS_max,Charged_states_max)
% DOS via Lehmann representation of the Green function
kb = 0.000086173324;  % ev/K
beta = 1/(kb*T);
h = 0.6582; % ev*fs
DOS = 0*w;

if (T > 0)
% set up Boltzman factors
p = zeros(1,States_DOS_max);
Q = 0;
for n = 1:States_DOS_max % loop over states of the neutral system
    p(n) = exp(-beta*E(n,1));
    Q = Q + p(n);
end % n
p = p/Q;
for n = 1:States_DOS_max % loop over states of the neutral system
    for j = sites   % loop over sites of active space
        [Cc,Ca] = create_anihilate_e(C(:,n),abs(j),sign(j),B,BM,Bm); 
        
        for m = 1:Charged_states_max % loop over states of the +1 charged system
            
            DOS = DOS + p(n)*eta*(dot(CM(:,m),Cc).^2)./( (w - (EM(m,1)-E(n,1))).^2 + eta*eta );

            DOS = DOS + p(n)*eta*(dot(Cm(:,m),Ca).^2)./( (w - (E(n,1)-Em(m,1))).^2 + eta*eta );

        end % m
    end % j
end % n

else % if T == 0 

% set up Boltzman factors
indsGS = find(abs(E(:,1)-E(1,1))<0.000000001); % indices of Ground State
for n = 1:length(indsGS) % loop over states of the neutral system
    for j = sites   % loop over sites of active space
        [Cc,Ca] = create_anihilate_e(C(:,n),abs(j),sign(j),B,BM,Bm); 
        for m = 1:Charged_states_max % loop over states of the +1 charged system 
            DOS = DOS + eta*(dot(CM(:,m),Cc).^2)./( (w - (EM(m,1)-E(n,1))).^2 + eta*eta );
            DOS = DOS + eta*(dot(Cm(:,m),Ca).^2)./( (w - (E(n,1)-Em(m,1))).^2 + eta*eta );
        end % m
    end % j
end % n
DOS = DOS/length(indsGS);


end % if T > 0 

end % end function spectral_function
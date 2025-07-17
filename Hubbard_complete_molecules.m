function[E0,C0,H0,B0]=Hubbard_complete_molecules(G,t,U)
% MOLECULES PPP complete
%cycloocta_relaxed
%G=load('azulene_relaxed.xyz');
%G=load('pentalene_relaxed.xyz');
%G=load('cycloocta_relaxed.xyz');
%G = G(:,[2,3]);
H = zeros(length(G));
Hfile = [];
Hfile_mb = [];

for j = 1:length(G)

    for k = 1:length(G)

        if (j~=k)

           distance=norm(G(j,:)-G(k,:));


           if (distance<2)

                H(j,k) = t;

                if (j<k)

                    Hfile = [Hfile; t, j, 0, 0, k];

                end % if (j<k)

           end % if (distance<2)


        end % if j ~=k

    end % k

    Hfile_mb = [Hfile_mb; U, j, -j, -j, j];

end % j

Hfile_mb = clean_Hfile(Hfile_mb);
[C,E]=eig(H);
                                                                            %Nel      %Nst
[C0,n0u,n0d,E0,DOS0,w0,B0,H0] = Many_body_hamil_S2([Hfile;Hfile_mb],length(G),length(G),mod(length(G),2),2,0);
%[C0oe,n0uoe,n0doe,E0oe,DOS0oe,w0oe,B0oe,H0oe,Cnatoe,Enatoe] = Many_body_hamil([Hfile],1,length(G),1,2,0);
E0ppp=0;  C0ppp=0; E0ppp_2=0;  C0ppp_2=0; E0ppp_4=0;  C0ppp_4=0;
%[~,~,~,E0ppp,~,C0ppp,] = PPP_and_CAS(G,length(G));

%[~,~,~,E0ppp_2,~,C0ppp_2] = PPP_and_CAS(G,length(G)-2);
%[~,~,~,E0ppp_4,~,C0ppp_4] = PPP_and_CAS(G,length(G)-4);



[E0(1:4,1)-E0(1,1),E0(1:4,2)]


end %end function PPP compekte molecules

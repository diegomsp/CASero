function[E_st,C,E,AE,E0,E1,C0,C1,B0,B1,H0,H1,n0u,n0d,n1u,n1d,gap_oe,Cnat,Enat,Rho,Fm,FM,CM,BM,Cm,Bm,EM,Em,HM,Hm,nuM,ndM,num,ndm,nFM,nFm,L_spin_ex,L_spin_flip,DOS,w,Hfile,Jmat] = Huckel_Hubbard(G,t,U,M,varargin)
%INPUTS:
%G: Geometry file
%t : hopping for the tight-binding
%U: On-site Hubbard's U
%M: size of the active space

%C : eigenvectors
%E: energies
%Q: Chiral matrix

narginv = length(varargin);



n0d=0; n0u=n0d; Rho=n0u; Enat=Rho; Cnat=Enat; L_spin_flip=Cnat; L_spin_ex=L_spin_flip; nFm=L_spin_ex; nFM=nFm; FM=nFM; Fm=FM; H1=Fm; B1=H1; w1=B1; DOS1=w1; E1=DOS1; n1d=E1; n1u=n1d; C1=n1u; Hm=C1; Bm=Hm; wm=Bm; DOSm=wm; Em=DOSm; ndm=Em; num=ndm; Cm=num; HM=Cm; BM=HM; wM=BM; DOSM=wM; EM=DOSM; ndM=EM; nuM=ndM; CM=nuM;
Jmat = 0;

Nat = length(G); Nel = floor(Nat/2);  %This Nel is the number of electrons per spin in a closed shell calculation
%Nel is the number of occupied molecular orbitals
Nelu = Nel;  Neld = Nel;
if (narginv == 1) %specify number of electrons in the system

  Nel1 = varargin{1};
  Nel = ceil(Nel1/2);
  Nelu = Nel;
  if (mod(Nel1,2) ~= 0)
   Neld = Nel-1;
  else % if (mod(Nel1,2) ~= 0)
   Neld = Nel;
  end %if (mod(Nat,2) ~= 0)

else %if (narginv == 1)

  if (mod(Nat,2) ~= 0)
   Nelu = Nel + 1;
   Neld = Nel;
  end %if (mod(Nat,2) ~= 0)

end %if narginv == 1

H = zeros(Nat);
tole = t/1000; %tolerance for the detection of zero mode states
%first neighbours tight binding hamil
sublat = zeros(1,Nat);
sublat(1) = 1;   %we have sublattice 1 and -1.
%this array stores the sublattice index of each site.
%the neighbours of one site belong to the opposite sublattice.
for j = 1:Nat
   for k = j+1:Nat
      distance = norm(G(j,:)-G(k,:));
      if (distance < 1.7)
         if (sublat(j) ~= 0)
            sublat(k) = -sublat(j);
         elseif (sublat(k) ~= 0)
            sublat(j) = -sublat(k);
         end %if (sublat(j) ~= 0)
         H(j,k) = t;
         H(k,j) = t;
      end % if (distance < 2)
   end %k
end % j
%{
%finish calculation of sublattice map:
while (prod(sublat) == 0)
   for j = 1:Nat
   for k = j+1:Nat
      distance = norm(G(j,:)-G(k,:));
      if (distance <1.7)
         if (sublat(j) ~= 0)
            sublat(k) = -sublat(j);
         elseif (sublat(k) ~= 0)
            sublat(j) = -sublat(k);
         end %if (sublat(j) ~= 0)
      end % if (distance < 2)
   end %k
end % j
end % while sublat == 0
%}
[C,E] = eig(H);  E=diag(E);
display('Energies SUMOs...')
E(Nel)
E(Nel+1)
gap_oe = E(Nel+1)-E(Nel);

Umol=U*sum(C(:,Nel).^4);

SSB = abs(E(:)')./(U*sum(C.^4));
Q=SSB';
%E_st = J/2;  %??  Ojo, revisar!!



Nel; Nel_oe=Nel;
AE = zeros(Nat,M);
for ii = 1:M
  if (mod(M,2)==0)
     AE(:,ii) = C(:,Nel - 0.5*M + ii);
  else
     AE(:,ii) = C(:,Nel - 0.5*(M+1) + ii + 1);
  end % end if
end % ii


Sz = 0; %default
SzM = 1; Szm = 1;
nLUMO = (M/2)+1;

if (mod(Nelu+Neld,2) ~= 0)

Sz = 1;
Szm = 0; SzM = Szm;
nLUMO = (M+1)/2;

end % if (mod(Nelu+Neld,2) ~= 0)


if (size(M,2) > 1)
   AE = zeros(Nat,size(M));
   for a = 1:size(M,2)
      AE(:,a) = C(:,M(a));
   end %a
   M = size(M,2);
end % end if (size(M) > 1)

%for indora = 1: M

%end % indora

Hfile = [];

for j = 1:M
   for k = j:M
      W = 0;
      for a = 1:Nat
         for b = 1:Nat
            W = W + H(a,b)*(AE(a,j)*AE(b,k));
         end %b
      end % a
      %W=0.1;
      Hfile = [Hfile; W , j, 0, 0, k];
   end % k
end % j

%ADD SCISSOR OPERATOR

Scissor = -0*0.30;

Hfile = [Hfile; Scissor , nLUMO , 0, 0, nLUMO];

%}
%Next add the many-body terms
%%{

for j = 1:M
   for jp = 1:M
      for k = 1:M
         for kp = 1:M
            % many-body parameter for this interaction
            V = U*sum( AE(:,j).*AE(:,jp).*AE(:,k).*AE(:,kp) );
            Hfile = [Hfile;  V, j, -k, -kp, jp];
            %Hfile = [Hfile; V, -j, k, kp, -jp];
            %Hfile = [Hfile; V, j, k, kp, jp];
            %Hfile = [Hfile; V, -j, -k, -kp, -jp];
         end % kp
      end % k
   end % jp
end % j

%Define Ncore
if (mod(M,2)==0)
    Ncore = Nel - 0.5*M;
else
    Ncore = Nel - 0.5*(M+1) + 1;
end % end if
%ADD efective field of the "deep" orbitals
for j = 1:M
    for k = j:M
        W = 0;
        for n = 1:Ncore
                 W = W + U*sum( AE(:,j).*AE(:,k).*C(:,n).*C(:,n) );
        end % n
       Hfile = [Hfile; W, j, 0, 0, k];
    end % k
end % j

%Separate 1c, 2c, 3c & 4c

%we can try adding a many-body interaction between nearest neighbors
%{
J = 0.1*U;
Up = 0.6412*U;  %take Up (the many-body nearest neighbor interacttion) as a 10% of U
for j = 1:M
   for jp = 1:M
      for k = 1:M
         for kp = 1:M
            V = Vuu = Vdd = Vud = Vdu = 0;
            % many-body parameter for this interaction
            for mu = 1:Nat
               for nu = 1:Nat
                  if (mu ~= nu && abs(H(mu,nu)) > 1)
                     V = V + 0.5*Up*( AE(mu,j)*AE(nu,k)*AE(nu,kp)*AE(mu,jp) );
                     % We have yet to include spin-dependency!!
                     %Vuu = Vuu + 0.5*Up*( AE(mu,j)*AE(nu,k)*AE(nu,kp)*AE(mu,jp) );
                     %Vud = Vud + 0.5*Up*( AE(mu,j)*AE(nu,k)*AE(nu,kp)*AE(mu,jp) );
                     %Vdd = Vdd + 0.5*Up*( AE(mu,j)*AE(nu,k)*AE(nu,kp)*AE(mu,jp) );
                     %Vdu = Vdu + 0.5*Up*( AE(mu,j)*AE(nu,k)*AE(nu,kp)*AE(mu,jp) );
                  end % end if mu ~= nu && H(mu,nu) > 1
               end % nu
            end % mu
            Hfile = [Hfile; V, j, k, kp, jp];     %UU
            Hfile = [Hfile; V, -j, -k, -kp, -jp]; %DD
            Hfile = [Hfile; V, j, -k, -kp, jp];   %UD
            Hfile = [Hfile; V, -j, k, kp, -jp];   %DU
         end % kp
      end % k
   end % jp
end % j
%}



%[C0,n0u,n0d,E0,DOS0,w0,B0,H0,Cnat,Enat,Rho] = Many_body_hamil(Hfile,M,M,Sz,2,0);

% To calculate many-body DOS:

 %E_ini = E(Nel - M);  E_end = E(Nel + M); nwsteps = 10000; eta = 0.01;
 E_ini = 0;  E_end = 0; nwsteps = 10000; eta = 0.01;
 levelsR = [1:M]; spinR = [1];
 Nel=M; Nst=M; max_oc = 2;

 EXCTS = 0;

    Hfile = clean_Hfile(Hfile);
    [C0,n0u,n0d,E0,DOS,w,B0,H0,Cnat,Enat,Rho] = Many_body_hamil_S2(Hfile,Nel,Nst,Sz,max_oc,EXCTS,E_ini,E_end,nwsteps,eta,levelsR,spinR);
    E0=real(E0);
    E_st=E0(2)-E0(1);
    %E0(1:4,:)

%[a] = draw_orbital_save(G,C,Nel_oe,'HOMO');
%[a] = draw_orbital_save(G,C,Nel_oe+1,'LUMO');
%[a] = draw_orbital_save(G,C,Nel_oe-1,'HOMO1');
%[a] = draw_orbital_save(G,C,Nel_oe+2,'LUMO1');
%[a] = draw_orbital_save(G,C,Nel_oe-2,'HOMO2');
%[a] = draw_orbital_save(G,C,Nel_oe+3,'LUMO2');
%dlmwrite ("HOMO.dat", [C(:,Nel_oe)], "delimiter", "   ", "newline", "\n","precision",5);
%dlmwrite ("HOMO1.dat", [C(:,Nel_oe-1)], "delimiter", "   ", "newline", "\n","precision",5);
%dlmwrite ("HOMO2.dat", [C(:,Nel_oe-2)], "delimiter", "   ", "newline", "\n","precision",5);
%dlmwrite ("LUMO.dat", [C(:,Nel_oe+1)], "delimiter", "   ", "newline", "\n","precision",5);
%dlmwrite ("LUMO1.dat", [C(:,Nel_oe+2)], "delimiter", "   ", "newline", "\n","precision",5);
%dlmwrite ("LUMO2.dat", [C(:,Nel_oe+3)], "delimiter", "   ", "newline", "\n","precision",5);
%[E0(1:4,1)-E0(1,1),E0(1:4,2)]

 if (Cnat==0)
 else
        size(Cnat)
        size(AE)
        Cnat_ab = AE*Cnat;
        for nn = 1:size(Cnat_ab,2)
            %[a] = draw_orbital_save(G,Cnat_ab,nn,strcat('Natural_Orbital',num2str(nn)),Enat(nn));
        end % nn
 end

%[a] = draw_orbital_save(G,C,Nel_oe,'HOMO');
%[a] = draw_orbital_save(G,C,Nel_oe+1,'LUMO');
%[a] = draw_orbital_save(G,C,Nel_oe-1,'HOMO1');
%[a] = draw_orbital_save(G,C,Nel_oe+2,'LUMO1');
%f = figure('visible','off')
%Nel_oe=Nel_oe+1
%plot([(Nel_oe-3:Nel_oe+4)],E(Nel_oe-3:Nel_oe+4),'+','LineWidth',6); line([Nel_oe-3 Nel_oe+4],[0 0])
%set(gca, "linewidth", 4, "fontsize", 12)
%xlabel('State number','FontSize', 32)
%ylabel('Energy (eV)','FontSize', 32)
%saveas (f,'Huckel_spectrum_ZOOM','jpg')

%{
%[C1,n1u,n1d,E1,DOS1,w1,B1,H1] = Many_body_hamil_S2(Hfile,M,M,Sz+2,2,0);

%[E0(1:4,1)-E0(1,1),E0(1:4,2)]

[CM,nuM,ndM,EM,DOSM,wM,BM,HM] = Many_body_hamil_S2(Hfile,M+1,M,SzM,2,0); %[Cm,num,ndm,Em,DOSm,wm,Bm] = Many_body_hamil_S2(Hfile,M+1,M,SzM+2,2,0);
[Cm,num,ndm,Em,DOSm,wm,Bm,Hm] = Many_body_hamil_S2(Hfile,M-1,M,Szm,2,0);
 %% Kondo Orbitals
 Kondo_orb = 0; mu = -2.00;

 if (Kondo_orb==1)
   EM(:,1) = EM(:,1) + mu; Em(:,1) = Em(:,1) - mu;
   [Jmat]=calculate_J_matrix_Kondo_orbitals(C0,B0,E0,CM,BM,EM,Cm,Bm,Em,1,2);
   [Korb_mo,Keig] = eig(Jmat); Keig = diag(Keig); [Keig,iKs] = sort(Keig); Korb_mo = Korb_mo(:,iKs);
   for j = 1:size(Korb_mo,2)
     Korb(:,j) = AE*Korb_mo(:,j);
     plot_sts_map(G,Korb(:,j),0 , strcat('Kondo_orb_',num2str(j),'_dIdV'),'plot_bond',false,'h',0.9,'edge_space',5.0,'dx',0.1,'z_eff',3.25);
     [a] = draw_orbital_save(G,Korb(:,j),1, strcat('Kondo_orb_',num2str(j)),num2str(Keig(j)));
   end % j
 end
%%{

[C1,n1u,n1d,E1,DOS1,w1,B1,H1] = Many_body_hamil_S2(Hfile,M,M,Sz+2,2,0);
[a] = draw_orbital_save(G,C,Nel_oe,'HOMO');
[a] = draw_orbital_save(G,C,Nel_oe+1,'LUMO');
[a] = draw_orbital_save(G,C,Nel_oe-1,'HOMO1');
[a] = draw_orbital_save(G,C,Nel_oe+2,'LUMO1');
[a] = draw_orbital_save(G,C,Nel_oe-2,'HOMO2');
[a] = draw_orbital_save(G,C,Nel_oe+3,'LUMO2');

plot_sts_map(G, C(:,Nel_oe),0 , strcat('HOMO','_dIdV'),'plot_bond',false,'h',0.9,'edge_space',5.0,'dx',0.1,'z_eff',3.25)
plot_sts_map(G, C(:,Nel_oe-1),0 , strcat('HOMO1','_dIdV'),'plot_bond',false,'h',0.9,'edge_space',5.0,'dx',0.1,'z_eff',3.25)
plot_sts_map(G, C(:,Nel_oe-2),0 , strcat('HOMO2','_dIdV'),'plot_bond',false,'h',0.9,'edge_space',5.0,'dx',0.1,'z_eff',3.25)
plot_sts_map(G, C(:,Nel_oe+1),0 , strcat('LUMO','_dIdV'),'plot_bond',false,'h',0.9,'edge_space',5.0,'dx',0.1,'z_eff',3.25)
plot_sts_map(G, C(:,Nel_oe+2),0 , strcat('LUMO1','_dIdV'),'plot_bond',false,'h',0.9,'edge_space',5.0,'dx',0.1,'z_eff',3.25)
plot_sts_map(G, C(:,Nel_oe+3),0 , strcat('LUMO2','_dIdV'),'plot_bond',false,'h',0.9,'edge_space',5.0,'dx',0.1,'z_eff',3.25)

f = figure('visible','off')
%Nel_oe=Nel_oe+1
plot([(Nel_oe-2:Nel_oe+2)],E(Nel_oe-2:Nel_oe+2),'+','LineWidth',6); line([Nel_oe-2 Nel_oe+2],[0 0])
set(gca, "linewidth", 4, "fontsize", 12)
xlabel('State number','FontSize', 32)
ylabel('Energy (eV)','FontSize', 32)
saveas (f,'Huckel_spectrum_ZOOM','jpg')


%%{


   %  DYSON ORBITALS

display('calculate Dyson orbitals');

for n = 1:4

   [Fm,FM] = Dyson_Orbitals(C0(:,1),Cm(:,n),CM(:,n),B0,Bm,BM,AE);

   nFm(n) = norm(Fm)
   nFM(n) = norm(FM)

   FM1=FM(1:Nat); FM2=FM(Nat+1:2*Nat);
   if ( sum(abs(FM1)) > sum(abs(FM2)) )
      FM = FM1;
   else
      FM = FM2;
   end % if

   Fm1=Fm(1:Nat); Fm2=Fm(Nat+1:2*Nat);
   if ( sum(abs(Fm1)) > sum(abs(Fm2)) )
      Fm = Fm1;
   else
      Fm = Fm2;
   end % if

   titleorbM = strcat('DysonMa',num2str(n));
   titleorbm = strcat('Dysonmi',num2str(n));

   %if (nFM(n) > 0.0001)
      FM=FM/norm(FM);
      [a] = draw_orbital_save(G,FM,1,titleorbM,num2str(nFM(n)));
      infoheadff = [Nat, Nat, 0];
      ff = zeros(1,2*length(FM)); ff(1:2:end) = FM'; %ff = [0,ff];
      %writematrix(FM,strcat('DysonMa_',num2str(n),'.dat'),"delimiter"," ")
      %dlmwrite(strcat('DysonMa_',num2str(n),'.dat'),infoheadff,"delimiter"," ")
      %dlmwrite(strcat('DysonMa_',num2str(n),'.dat'),ff,"delimiter"," ","-append")
      %%writematrix(infoheadff,strcat('DysonMa_',num2str(n),'.dat'),"delimiter"," ")
      %%writematrix(ff,strcat('DysonMa_',num2str(n),'.dat'),"delimiter"," ","-append")
      %writecell({infoheadff; ff},strcat('DysonMa_',num2str(n),'.dat'),"delimiter"," ")
      plot_sts_map(G, FM,0 , strcat(titleorbM,'_dIdV'),'plot_bond',false,'h',0.9,'edge_space',5.0,'dx',0.1,'z_eff',3.25)
   %end % if
   %if (nFm(n) > 0.0001)
      Fm=Fm/norm(Fm);
      [a] = draw_orbital_save(G,Fm,1,titleorbm,num2str(nFm(n)));
      infoheadff = [Nat, Nat, 0];
      ff = zeros(1,2*length(Fm)); ff(1:2:end) = Fm'; %ff = [0,ff];

      %writematrix(Fm,strcat('Dysonmi_',num2str(n),'.dat'),"delimiter"," ")
      %dlmwrite(strcat('Dysonmi_',num2str(n),'.dat'),infoheadff,"delimiter"," ")
      %dlmwrite(strcat('Dysonmi_',num2str(n),'.dat'),ff,"delimiter"," ","-append")
      %%writematrix(infoheadff,strcat('Dysonmi_',num2str(n),'.dat'),"delimiter"," ")
      %%writematrix(ff,strcat('Dysonmi_',num2str(n),'.dat'),"delimiter"," ","-append")
      %writecell({infoheadff; ff},strcat('Dysonmi_',num2str(n),'.dat'),"delimiter"," ")
      plot_sts_map(G, Fm,0 , strcat(titleorbm,'_dIdV'),'plot_bond',false,'h',0.9,'edge_space',5.0,'dx',0.1,'z_eff',3.25)
   %end % if

end % n dyson orbitals

%NATURAL TRANSITION ORBITALS       %Draw these orbitals
display('Natural Transition Orbitals');
%[F,L_spin_ex] = Natural_Transition_Orbitals_general(C0,B0,C0,B0,AE,1,2,1,1); L_spin_ex = diag(L_spin_ex);
[F,L_spin_ex] = Natural_Transition_Orbitals_general(C0,B0,C1,B1,AE,1,1,1,-1); L_spin_ex = diag(L_spin_ex);
dlmwrite ("NTO_spin_ex.dat", [L_spin_ex';F], "delimiter", "   ", "newline", "\n","precision",5)
for j = 1:size(F,2)

   titlento = strcat('NTO_spin_ex_',num2str(j));
   if (norm(F(:,j)) > 0.0001 )
      [a] = draw_orbital_save(G,F,j,titlento,num2str(L_spin_ex(j)));
      plot_sts_map(G, F(:,j),0 , strcat(titlento,'_dIdV'),'plot_bond',false,'h',0.9,'edge_space',5.0,'dx',0.1,'z_eff',3.25)
    %writematrix(F(:,j),strcat('NTO_spin_ex_',num2str(j),'.dat'),"delimiter"," ")
    %%writematrix(M,'M_tab.txt','Delimiter','tab')

   end % if

end % j

L_spin_flip = 0;

if (Sz ~= 0)

   [C0o,n0uo,n0do,E0o,DOSo,wo,B0o] = Many_body_hamil_S2(Hfile,M,M,-Sz,2,0);

   [F,L_spin_flip] = Natural_Transition_Orbitals_general(C0,B0,C0o,B0o,AE,1,1,-1,1);  L_spin_flip = diag(L_spin_flip);
   dlmwrite ("NTO_spin_flip.dat", [L_spin_ex';F], "delimiter", "   ", "newline", "\n","precision",5)
   for j = 1:size(F,2)

      titlento = strcat('NTO_spin_flip_',num2str(j));
      if (norm(F(:,j)) > 0.0001 )
          titlento
         [a] = draw_orbital_save(G,F,j,titlento,num2str(L_spin_flip(j)));
         plot_sts_map(G, F(:,j),0 , strcat(titlento,'_dIdV'),'plot_bond',false,'h',0.9,'edge_space',5.0,'dx',0.1,'z_eff',3.25)
         %writematrix(F(:,j),strcat('NTO_spin_flip_',num2str(j),'.dat'),"delimiter"," ")
      end % if

   end % j


end % end if
%}
%}

end % end function tight_binding_model_and_parameters

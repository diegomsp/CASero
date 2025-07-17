function[C,E,AE,E00,C00,B00,H00,E02,C02,B02,H02,E0m2,C0m2,B0m2,H0m2,EM1,CM1,BM1,HM1,EMm1,CMm1,BMm1,HMm1,EM3,CM3,BM3,HM3,EMm3,CMm3,BMm3,HMm3,Em1,Cm1,Bm1,Hm1,Emm1,Cmm1,Bmm1,Hmm1,Em3,Cm3,Bm3,Hm3,Emm3,Cmm3,Bmm3,Hmm3] = Extract_Multiplets(G,t,U,M,varargin)
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


%Sz = 0; %default
%SzM = 1; Szm = 1;
nLUMO = (M/2)+1;

if (mod(Nelu+Neld,2) ~= 0)

%Sz = 1;
%Szm = 0; SzM = Szm;
nLUMO = (M+1)/2;

end % if (mod(Nelu+Neld,2) ~= 0)


Sz  = [-5:5];
SzM = [-5:5];
Szm = [-5:5];

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






% To calculate many-body DOS:

 %E_ini = E(Nel - M);  E_end = E(Nel + M); nwsteps = 10000; eta = 0.01;
 E_ini = 0;  E_end = 0; nwsteps = 10000; eta = 0.01;
 levelsR = [1:M]; spinR = [1];
 Nel=M; Nst=M; max_oc = 2;

 EXCTS = 0;

    Hfile = clean_Hfile(Hfile);


    [CM1,~,~,EM1,~,~,BM1,HM1] = Many_body_hamil_S2(Hfile,M+1,M,1,2,0);
    [CMm1,~,~,EMm1,~,~,BMm1,HMm1] = Many_body_hamil_S2(Hfile,M+1,M,-1,2,0);
    [CM3,~,~,EM3,~,~,BM3,HM3] = Many_body_hamil_S2(Hfile,M+1,M,3,2,0);
    [CMm3,~,~,EMm3,~,~,BMm3,HMm3] = Many_body_hamil_S2(Hfile,M+1,M,-3,2,0);

    [Cm1,~,~,Em1,~,~,Bm1,Hm1] = Many_body_hamil_S2(Hfile,M-1,M,1,2,0);
    [Cmm1,~,~,Emm1,~,~,Bmm1,Hmm1] = Many_body_hamil_S2(Hfile,M-1,M,-1,2,0);
    [Cm3,~,~,Em3,~,~,Bm3,Hm3] = Many_body_hamil_S2(Hfile,M-1,M,3,2,0);
    [Cmm3,~,~,Emm3,~,~,Bmm3,Hmm3] = Many_body_hamil_S2(Hfile,M-1,M,-3,2,0);

    [C00,~,~,E00,~,~,B00,H00] = Many_body_hamil_S2(Hfile,M,M,0,2,0);
    [C02,~,~,E02,~,~,B02,H02] = Many_body_hamil_S2(Hfile,M,M,2,2,0);
    [C0m2,~,~,E0m2,~,~,B0m2,H0m2] = Many_body_hamil_S2(Hfile,M,M,-2,2,0);





end % end function tight_binding_model_and_parameters

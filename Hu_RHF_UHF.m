function[Chu,Ehu,nhu,Chf,Ehf,nhf,Cu,Cd,Eu,Ed,nu,nd,Nel,SSB,Ehf_shift,Eu_shift] = Hu_RHF_UHF(G,t,U)
%DO HUCKEL, HARTREE FOCK RESTRICTED AND HARTREE FOCK UNRESTRICTED
%INPUTS:
%G: Geometry file
%t : hopping for the tight-binding
%U: On-site Hubbard's U


%C : eigenvectors
%E: energies
%Q: Chiral matrix

Nat = length(G); Nel = floor(Nat/2);  %This Nel is the number of electrons per spin in a closed shell calculation
%Nel is the number of occupied molecular orbitals
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
      if (distance < 1.8)
         if (sublat(j) ~= 0)
            sublat(k) = -sublat(j);
         elseif (sublat(k) ~= 0)
            sublat(j) = -sublat(k);
         end %if (sublat(j) ~= 0)

         H(j,k) = t;
         H(k,j) = t;
         % H(j,k) = -28.07749*exp(-1.65878*distance);
         % H(k,j) = -28.07749*exp(-1.65878*distance);
      end % if (distance < 2)
   end %k
end % j
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

[C,E] = eig(H);  E=diag(E);
display('Energies SUMOs...')
E(Nel)
E(Nel+1)
Chu=C; Ehu=E;
nhu = sum(Chu(:,1:Nel).^2,2);

SSB = (abs(E(:)')./(U*sum(C.^4)))';

% RESTRICTED HARTREE-FOCK CALCULATION

tole_scf = 0.000001;  %tolerance for scf calculation

beta = 0.70; %coefficient for the mixing algorithm
display('Hartree Fock (restricted)')

%Initial charges   !OJO   Initial guess important

n0hf = 0.5 + zeros(Nat,1);
n0hf(1) = n0hf(1) + 0.3; n0hf(Nat) = n0hf(Nat) - 0.3;

n0hf(2) = n0hf(2) + 0.5; n0hf(3) = n0hf(3) - 0.5;
 %(should be balanced, Sz=0)

error_scf = 1; %initial arbitrary large value for the error
while (error_scf > tole_scf)

   %Hhf = H + diag(U*n0hf.*(1-n0hf));
   Hhf = H + diag(U*n0hf) - 0.5*U*sum(n0hf.*n0hf)*eye(Nat);
   [Chf,Ehf] = eig(Hhf);  Ehf=diag(Ehf);  [Ehf,ishf]=sort(Ehf); Chf=Chf(:,ishf);

   n1hf = sum(Chf(:,1:Nel).^2,2); %n1hf(1)
   %n1hf(2)
   %n1hf(3)
   %n1hf(Nel)
   error_scf = norm(n1hf-n0hf);
   if (Nat ~= 2*Nel)
      error_scf = 0;
   end  %
   %Define next input charges (with linear mixing algorithm)
   n0hf = beta*n0hf + (1-beta)*n1hf;
   %error_scf=0;  %BORRAR
end % (error_scf > tole_scf)

nhf = n1hf;
Ehf_shift = Ehf - 0.5*(Ehf(Nel) + Ehf(Nel+1));


%   UNRESTRICTED HARTREE-FOCK CALCULATION

Nelu = Nel;  Neld = Nel;  %Number of up and down electrons.
                              %This can be changed
                              %Recall that Nel above was the number of electrons
                              %per spin in a closed shell calculation
  if (Nat ~= 2*Nel)
      Nelu = Nel + 1;
   end  %
tole_scf = 0.000001;  %tolerance for scf calculation

beta = 0.70; %coefficient for the mixing algorithm
display('Hartree Fock (Unrestricted)')


%Initial charges   !OJO   Initial guess important
%i1 = 2;  i2 = 48;   %rhomb4
i1 = 1;  i2 = 30;  %rhomb3
delta = 0.1;
m0 = 0.1;
n0u = 0.5 + zeros(Nat,1) + m0*((-1).^[1:Nat])';
n0d = 0.5 + zeros(Nat,1) - m0*((-1).^[1:Nat])';  %(should be balanced, Sz=0)
%n0u(i1) = 0.5 + delta; n0u(i2) = 0.5 - delta; n0d(i1) = 0.5 - delta; n0d(i2) = 0.5 + delta;
%How to (initially) break the ud symmetry preserving Sz=0 ?
error_scf = 1; %initial arbitrary large value for the error
while (error_scf > tole_scf)
   %create mean-field open shell hamiltonian
   %H2 = blkdiag(H + diag(U*n0d.*(1-n0u)),H + diag(U*n0u.*(1-n0d)));
   %[C2,E2] = eig(H2);  E2=diag(E2);  %How to sort them??
   %Sort eigenstates!!  First the ones of spin up, then the ones with spin down
   %We could diagonalize the two hamils separadtely
   %Hu = H + diag(U*n0d.*(1-n0u));
   Hu = H + diag(U*n0d) - 0.5*U*sum(n0d.*n0u)*eye(Nat);
   [Cu,Eu] = eig(Hu);  Eu=diag(Eu);  [Eu,isu]=sort(Eu); Cu=Cu(:,isu);
   %Hd = H + diag(U*n0u.*(1-n0d));
   Hd = H + diag(U*n0u) - 0.5*U*sum(n0d.*n0u)*eye(Nat);
   [Cd,Ed] = eig(Hd);  Ed=diag(Ed);  [Ed,isd]=sort(Ed); Cd=Cd(:,isd);
   %compute new charges
   %n1u = sum(C2(1:Nat,1:Nelu).^2,2); n1d = sum(C2(Nat+1:2*Nat,Neld+1:Neld+Neld).^2,2);
   n1u = sum(Cu(:,1:Nelu).^2,2); n1d = sum(Cd(:,1:Neld).^2,2);
   error_scf = norm(n1u-n0u) + norm(n1d-n0d);
   %Define next input charges (with linear mixing algorithm)
   n0u = beta*n0u + (1-beta)*n1u;
   n0d = beta*n0d + (1-beta)*n1d;
end % (error_scf > tole_scf)

nu = n1u;  nd = n1d;
Eu_shift = Eu - 0.5*(Eu(Nel) + Eu(Nel+1));
%%{
%PLOT ORBITALS

end  % end function Hu_RHF_UHF

function[F,L] = Natural_Transition_Orbitals_general(C0,B0,C1,B1,AE,n,m,s1,s2)
% n,m index of excited states we want to play with
%s1, s2 spin states of the creation/anihilation operators we want to play with:  1 means spin up, -1 means spin down.  0.5*(1-s)
s1 = 0.5*(1-s1); s2 = 0.5*(1-s2);
%[F] = Transition_Martin_Orbitals(C0p,B0p,C0m,B0m,AE,1,1);

C0 = C0(:,n);
C1 = C1(:,m);

Nst = numel(B0(1,:));
Nel0 = sum(abs(B0(1,:)));

Nelu0 = length(find(B0(1,:)==1)) + length(find(B0(1,:)==2));
Neld0 = length(find(B0(1,:)==-1)) + length(find(B0(1,:)==2));
Nelu1 = length(find(B1(1,:)==1)) + length(find(B1(1,:)==2));
Neld1 = length(find(B1(1,:)==-1)) + length(find(B1(1,:)==2));

%go over to 1's and 0's representation:
B0u = B0;  B0d = B0; B0u(B0u==2)=1; B0u(B0u==-1)=0;
B0d(B0d==1)=0; B0d(B0d==-1)=1; B0d(B0d==2)=1;
B0t = [B0u,B0d];

B1u = B1; B1d = B1; B1u(B1u==2)=1; B1u(B1u==-1)=0;
B1d(B1d==1)=0; B1d(B1d==-1)=1; B1d(B1d==2)=1;
B1t = [B1u,B1d];

%find(b==1)

tole = 0.000001;

C0_important = find(abs(C0)>tole);

%initiliaze variables

F2 = zeros(Nst);

F = zeros(Nst,1);

%recall that we use the convenion: first all spin up then all spin down
%  < C1  |  c'_As1 c_Bs2   |  C0 >

for A = 1:Nst

   for B = 1:Nst

      for J = 1:length(C0_important); %1:length(C0(:,1)) %Loop over Slater determinants
         J = C0_important(J);

         b = B0t(J,:);
         ib = find(b==1);

         for l = 1:Nel0

            if (ib(l) == B + s2*Nst )

               ibp = ib;
               ibp(l) = A + s1*Nst;
               bp = zeros(1,2*Nst);

               % sign of the permutation
               [seq1,iseq1] = sort(ibp);
               I = eye(length(iseq1));
               sg1 = det(I(:,iseq1));

               bp(ibp) = 1;


               [~,indbp] = ismember(int8(bp),B1t,'rows');

               if (indbp ~= 0)

                  F2(A,B) = F2(A,B) + sg1*C1(indbp)*C0(J);  %Corrected to include sign of the permutation (July 2023)

               end % end if (indbp != 0)

            end % end if ib(l) == B + Nst

         end % l

      end % end loop J %Loop over Slater determinants

   end % B

end % A

[U,L] = eig(F2);

F = AE*U;

end %end function Transition_Spin_Orbitals

function[F,L] = Transition_Martin_Orbitals_spinflip(C0,B0,C1,B1,AE,n,m)

%[F] = Transition_Martin_Orbitals(C0p,B0p,C0m,B0m,AE,1,1);

C0 = C0(:,n);
C1 = C1(:,m);

Nst = numel(B0(1,:));
Nel = sum(abs(B0(1,:)));

Nelu0 = length(find(B0(1,:)==1)) + length(find(B0(1,:)==2));
Neld0 = length(find(B0(1,:)==-1)) + length(find(B0(1,:)==2));
Nelu1 = length(find(B1(1,:)==1)) + length(find(B1(1,:)==2));
Neld1 = length(find(B1(1,:)==-1)) + length(find(B1(1,:)==2));

%go over to 1's and 0's representation:
B0u = B0; B0d = B0; B0u(B0u==2)=1; B0u(B0u==-1)=0; 
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

%recall that here we use the convenion: first all spin up then all spin down

for A = 1:Nst  %this with spin up
  
   for B = 1:Nst  %this with spin down
     
      for J = 1:length(C0_important); %1:length(C0(:,1)) %Loop over Slater determinants
         J = C0_important(J);
         %size(B0t)
         b = B0t(J,:);
         ib = find(b==1);
         %ib
         for l = 1:Nel
           
            if (ib(l) == B )
            %if (ib(l) == B)
              
               ibp = ib;
               ibp(l) = A + Nst;
               %ibp(l) = A + Nst;
               bp = zeros(1,2*Nst);
               bp(ibp) = 1;
               %Nst
               %size(B1t)
               %size(bp)
               [~,indbp] = ismember(int8(bp),B1t,'rows');
               if (indbp ~= 0)
                 
                  F2(A,B) = F2(A,B) + C1(indbp)*C0(J);
                 
               end % end if (indbp ~= 0)
                
            end % end if ib(l) == B + Nst
           
         end % l
         
      end % end loop J %Loop over Slater determinants
     
     %F(A) = F(A) + F2(A,B);
     
   end % B
  
end % A

[U,L] = eig(F2*F2');

F = AE*U;
          
end %end function Transition_Spin_Orbitals
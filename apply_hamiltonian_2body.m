function[R] = apply_hamiltonian_2body(Ri,Hfile,Nst) 
% Hfile is a N x 5 array (N is the number of terms in the hamiltonian)
% each row codifies a one-electron or a many-body hamiltonian terminal_size
% Each row of Hfile looks like:
% val i j k l
% val is the coefficient (hopping, on-site, Hubbard's U...)
%If k=l=0 this is a one-el term. If k and l are not 0, it is many body
% i j k l can be positive or negative. If negative, they correspond to spin down
%Ri is an input array. Each row is an element of the basis with a weight
%Ri(:,1) are the weights. 
%Ojo...!!  Check the signs!!
  %Nst = max(max(abs(Hfile(:,2:end))));  %change!!
  V = zeros(1,2*Nst);
  R=[];
  %HH = Hfile(:,2:end);
 %HH(HH<0) = abs(HH(HH<0)) + Nst;
 %Hfile = [Hfile(:,1),HH];
  for a = 1:size(Hfile,1)
      i = abs(Hfile(a,2)); j = abs(Hfile(a,3)); k = abs(Hfile(a,4)); l = abs(Hfile(a,5));
        %  MANY-BODY TERMS
        %%{
           si = sign(Hfile(a,2));  
           sj = sign(Hfile(a,3));
           sk = sign(Hfile(a,4));  
           sl = sign(Hfile(a,5));

           %transform b and b0 in 1's and 0's list. 
           if (si == -1)
             i = i + Nst;
           end
           if (sj == -1)
             j = j + Nst;
           end
           if (sk == -1)
             k = k + Nst;
           end
           if (sl == -1)
             l = l + Nst;
           end
           %}
     for el = 1:size(Ri,1)
        b0 = Ri(el,2:end); b = b0;    
        
           %}
           %create list of 1's and 0's. Order: first all spin up, from 1st level to
           % last level. Then all spin down, from 1st level to last level
           V1 = b; V1(V1==-1)=0; V1(V1==2)=1;
           V2 = b; V2(V2==1)=0; V2(V2==2)=1; V2=abs(V2);
           V = [V1,V2];
        
           [V,signo] = apply_operator_2body(i,j,k,l,V);
           if (~isempty(V))
          
       
              V1 = V(1:Nst);  V2 = 0.5*V(Nst+1:end); V = V1+V2;
              V(V==1.5)=2; V(V==0.5)=-1;
    
              R = [R; Ri(el,1)*Hfile(a,1)*signo,V];
           end % if length(V) > 0
           
           %ADD THE OTHER IMPLIED TERMS (If we have a REDUCED Hfile) (To do... 

     end % el = 1:size(Ri,1)
  end % a = 1:size(Hfile,1)
  
 
end %end function apply_hamiltonian

function[B,BM,Bm]=construct_basis2(Nst,Nel,max_oc,Sz,ifBM_Bm)
 
BM = zeros(2,2,2); Bm = zeros(2,2,2);
 
Nstps = Nst; 
Nst = 2*Nst;
  
 %Nel : total number of electrons in the system 
 %Nst : total number of  sites  in the system (counting spin)
 %Nstps : number of sites per spin
 %Distribute Net electrons in Nst sites in any way.
 %Basis element is encoded first in sequence of 1's and 0's [1 0...1 1...] 
 %of lenght Nst and Net 1's.  Then in a sequence of -1,0,1,2's of lenght Nst/2
 seed = zeros(1,Nst); %from this seed we define the many-body configurations of 
 %the basis by setting some of its components to 1.
 Nbas = 0;
 for in = 1:length(Nel)
    Nbas = Nbas + nchoosek(Nst,Nel(in)); %number of basis elements (and size of hamiltonian)
 end % j
 B = [];

for in = 1:length(Nel)
 basisel = nchoosek(1:Nst,Nel(in));  
 %basisel : collection of the different many-body configurations 
 if (size(max_oc)==1)
 for j = 1:nchoosek(Nst,Nel(in))
   vbasis = seed;
   vbasis(basisel(j,:))=1;  %an element of the basis; a many-body configuration
   vbasis = reshape(vbasis,2,Nstps);
   spinz = sum(vbasis(1,:))-sum(vbasis(2,:));
   vbasis = vbasis(1,:)-vbasis(2,:)+2*vbasis(1,:).*vbasis(2,:);
   %with this, we pass to representation 1 2 -1 0 ....
   %if ... add restrictions (typically to account for symmetries of the hamiltonian)
   if (ismember(spinz,[Sz; Sz]) && max(vbasis) < max_oc+1)
      B = [B; vbasis];
   end % end 
 end %j
 else % (size(max_oc)==1)  size(max_oc)==2
 %max_oc=[M,N]. This means that we want exactly N electrons in the levels 1,2,...,M
 Nlevel=max_oc(1); Nelec=max_oc(2);
 for j = 1:nchoosek(Nst,Nel(in))
   vbasis = seed;
   vbasis(basisel(j,:))=1;  %an element of the basis; a many-body configuration
   vbasis = reshape(vbasis,2,Nstps);
   spinz = sum(vbasis(1,:))-sum(vbasis(2,:));
   vbasis = vbasis(1,:)-vbasis(2,:)+2*vbasis(1,:).*vbasis(2,:);
   %with this, we pass to representation 1 2 -1 0 ....
   %if ... add restrictions (typically to account for symmetries of the hamiltonian)
   if (ismember(spinz,[Sz; Sz]) && sum(abs(vbasis(1:Nlevel)))==Nelec) 
       B = [B; vbasis];
   end % end 
 end %j
 end %if (size(max_oc)==1)
end % in





 
 
 if (ifBM_Bm == 1 && length(Nel)==1)  %Basis for : Create / anihilate electron with spin up
 %+++++++++++++++++++++++++++++++++++++++++
 %ONE MORE ELECTRON
 %+++++++++++++++++++++++++++++++++++++++++

 seed = zeros(1,Nst); %from this seed we define the many-body configurations of 
 %the basis by setting some of its components to 1.
 Nbas = nchoosek(Nst,Nel+1); %number of basis elements (and size of hamiltonian)
 basisel = nchoosek(1:Nst,Nel+1);   
 %basisel : collection of the different many-body configurations
 BM = zeros(Nbas,Nstps);
 ind_garb = [];
 for j = 1:Nbas
   vbasis = seed;
   vbasis(basisel(j,:))=1;  %an element of the basis; a many-body configuration
   vbasis = reshape(vbasis,2,Nstps);
   spinz = sum(vbasis(1,:))-sum(vbasis(2,:));
   vbasis = vbasis(1,:)-vbasis(2,:)+2*vbasis(1,:).*vbasis(2,:);
   %with this, we pass to representation 1 2 -1 0 ....
   %if ... add restrictions (typically to account for symmetries of the hamiltonian)
   if (ismember(spinz,[Sz+1,Sz-1]) )%&& sort(vbasis,'descend')(2) < max_oc+1) 
      BM(j,:) = vbasis;
   else
      ind_garb = [ind_garb,j];
   end % end 
 end %j
 BM(ind_garb,:) = [];
 
 %+++++++++++++++++++++++++++++++++++++++++
 %ONE LESS ELECTRON
 %+++++++++++++++++++++++++++++++++++++++++

 seed = zeros(1,Nst); %from this seed we define the many-body configurations of 
 %the basis by setting some of its components to 1.
 Nbas = nchoosek(Nst,Nel-1); %number of basis elements (and size of hamiltonian)
 basisel = nchoosek(1:Nst,Nel-1);   
 %basisel : collection of the different many-body configurations
 Bm = zeros(Nbas,Nstps);
 ind_garb = [];
 for j = 1:Nbas
   vbasis = seed;
   vbasis(basisel(j,:))=1;  %an element of the basis; a many-body configuration
   vbasis = reshape(vbasis,2,Nstps);
   spinz = sum(vbasis(1,:))-sum(vbasis(2,:));
   vbasis = vbasis(1,:)-vbasis(2,:)+2*vbasis(1,:).*vbasis(2,:);
   %with this, we pass to representation 1 2 -1 0 ....
   %if ... add restrictions (typically to account for symmetries of the hamiltonian)
   if (ismember(spinz,[Sz+1,Sz-1]) && max(vbasis) < max_oc+1) 
      Bm(j,:) = vbasis;
   else
      ind_garb = [ind_garb,j];
   end % end 
 end %j
 Bm(ind_garb,:) = [];

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



end % end function construct basis 
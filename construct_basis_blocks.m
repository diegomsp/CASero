function[B,BM,Bm]=construct_basis_blocks(Nst,Nel,max_oc,Sz,ifBM_Bm)
 
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
 Nbas = nchoosek(Nst,Nel); %number of basis elements (and size of hamiltonian)
 basisel = nchoosek(1:Nst,Nel);   
 %basisel : collection of the different many-body configurations
 B = zeros(Nbas,Nstps);
 ind_garb = [];
 
 nblocks_oc = size(max_oc,1);
 nblocks_Sz = size(Sz,1);
 
 for j = 1:Nbas
   ok_to_add = 1;
   vbasis = seed;
   vbasis(basisel(j,:))=1;  %an element of the basis; a many-body configuration
   vbasis = reshape(vbasis,2,Nstps);
   spinz = sum(vbasis(1,:))-sum(vbasis(2,:));
   vbasis = vbasis(1,:)-vbasis(2,:)+2*vbasis(1,:).*vbasis(2,:);
   %with this, we pass to representation 1 2 -1 0 ....
   %if ... add restrictions (typically to account for symmetries of the hamiltonian)
   
   for iblock_oc = 1:nblocks_oc
      occu = max_oc{iblock_oc,1}; block = max_oc{iblock_oc,2};
      occu_vbasis = sum(  abs(vbasis( block ))  );
      if (~ismember(occu_vbasis,occu))
         ok_to_add = 0;
         break
      end % occu
   end % iblock_oc
   if (ok_to_add) 
   for iblock_sz = 1:nblocks_Sz
      sz_allowed = Sz{iblock_sz,1}; block = Sz{iblock_sz,2};
      sz_vbasis = 0;   vbasis1=vbasis; vbasis1=vbasis1(block); vbasis1(vbasis1==2)=0; sz_vbasis = sum(vbasis1);   %%%%%%%%%% 
      if (~ismember(sz_vbasis,sz_allowed))
         ok_to_add = 0;
         break
      end % occu
   end % iblock_oc
   end % ok_to_add
   
   if (ok_to_add) 
      B(j,:) = vbasis;
   else
      ind_garb = [ind_garb,j];
   end % end 
   
 end %j
 B(ind_garb,:) = [];
 
 
 if (ifBM_Bm == 1)  %Basis for : Create / anihilate electron with spin up
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
   if (ismember(spinz,Sz+1) )%&& sort(vbasis,'descend')(2) < max_oc+1) 
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
   if (ismember(spinz,Sz-1) && max(vbasis) < max_oc+1) 
      Bm(j,:) = vbasis;
   else
      ind_garb = [ind_garb,j];
   end % end 
 end %j
 Bm(ind_garb,:) = [];

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



end % end function construct basis 
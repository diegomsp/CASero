function[H,HM,Hm] = create_hamiltonian_sparse_matrix_general_openshell(Hfile,B,BM,Bm,ifBM_Bm)

n = size(B,1);
nr = []; nc = []; nv = []; %structure arrays for the sparse hamiltonian matrix
for j = 1:n %loop over basis elements
  %j
  %tic
  Sc = B(j,:); 
  Ri = [1, Sc];
  [R] = apply_hamiltonian_general_openshell(Ri,Hfile);
  nn = size(R); nn = nn(1);
  Sr = R(:,2:end);
  %H(j,j) = R(1,1);
  for inn = 1:nn
    [~,row_indx] = ismember(int8(R(inn,2:end)),B,'rows');
    if (row_indx ~= 0)
       %H(j,m) = R(inn,1);
       nr = [nr; j]; nc = [nc; row_indx]; nv = [nv; R(inn,1)];
    end %if valm0
    %H(j,m) = R(inn,1);
  end %inn
  %toc
end %end for j = 1:n
H = sparse(nr,nc,nv);

%+++++++++++++++++ ONE LESS / ONE MORE ELECTRONS
if (ifBM_Bm == 1)
  
%Now hamiltonian for space with one less electron

n = size(Bm,1);
nr = []; nc = []; nv = []; %structure arrays for the sparse hamiltonian matrix
for j = 1:n %loop over basis elements
 % j
  %tic
  Sc = Bm(j,:); 
  Ri = [1, Sc];
  [R] =apply_hamiltonian_general_openshell(Ri,Hfile);
  nn = size(R); nn = nn(1);
  Sr = R(:,2:end);
  %H(j,j) = R(1,1);
  for inn = 1:nn
    [~,row_indx] = ismember(int8(R(inn,2:end)),Bm,'rows');
    if (row_indx ~= 0)
       %H(j,m) = R(inn,1);
       nr = [nr; j]; nc = [nc; row_indx]; nv = [nv; R(inn,1)];
    end %if valm0
    %H(j,m) = R(inn,1);
  end %inn
  %toc
end %end for j = 1:n
Hm = sparse(nr,nc,nv);

%Now hamiltonian for space with one more electron

n = size(BM,1);
nr = []; nc = []; nv = []; %structure arrays for the sparse hamiltonian matrix
for j = 1:n %loop over basis elements
 % j
  %tic
  Sc = BM(j,:); 
  Ri = [1, Sc];
  [R] = apply_hamiltonian_general_openshell(Ri,Hfile);
  nn = size(R); nn = nn(1);
  Sr = R(:,2:end);
  %H(j,j) = R(1,1);
  for inn = 1:nn
    [~,row_indx] = ismember(int8(R(inn,2:end)),BM,'rows');
    if (row_indx ~= 0)
       %H(j,m) = R(inn,1);
       nr = [nr; j]; nc = [nc; row_indx]; nv = [nv; R(inn,1)];
    end %if valm0
    %H(j,m) = R(inn,1);
  end %inn
  %toc
end %end for j = 1:n
HM = sparse(nr,nc,nv);

else
HM=zeros(2); Hm=zeros(2);
end % end if ifBM_Bm == 1

end  %end function create_hamiltonian_matrix
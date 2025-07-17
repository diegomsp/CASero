function[H,HM,Hm] = create_hamiltonian_sparse_matrix_general(Hfile1body,Hfile2body,B,BM,Bm,ifBM_Bm,Nst)

n = size(B,1);
nr = []; nc = []; nv = []; %structure arrays for the sparse hamiltonian matrix
%counte = 0;
%NNN=20*n; nr = zeros(NNN,1); nc=nr; nv=nr;
%n
%display('BOTTLENECK')
for j = 1:n %loop over basis elements
  %tic
  %j
  %tic
  Sc = B(j,:); 
  Ri = [1, Sc];
  %disp('apply hamil')
  %tic
  %disp('apply 1 hamil')
  %tic
  [R1] = apply_hamiltonian_1body(Ri,Hfile1body,Nst);
  %toc
  %disp('apply 2 hamil')
  %tic
  [R2] = apply_hamiltonian_2body(Ri,Hfile2body,Nst);
  %toc
  R = [R1;R2];
  %toc
  %R = [R1];
  %R = [R2];
  %disp('check for repeated rows and simplify')
 %tic
  %nrr = size(R,1);
  %inds = [];
  %for jj = 1:nrr
  %  for k = jj+1:nrr
  %    if R(jj,2:end)==R(k,2:end)
  %      inds =[inds,k];
  %      R(jj,1) = R(jj,1)+R(k,1);
  %    end % end if A(jj) == A(k)
  %  end % end for k = 1:nrr
  %end % end for jj = 1:nrr
  %R(inds,:) = [];
  %toc
  nn = size(R); nn = nn(1);
  %Sr = R(:,2:end);
  %H(j,j) = R(1,1);
  %disp('create sparse matrix structure')
  %tic
  %nn
  for inn = 1:nn
    %disp('ismember')
    %tic
    [~,row_indx] = ismember(int8(R(inn,2:end)),B,'rows');
    %toc
    %disp('create sparse matrix structure')
    %tic
    if (row_indx ~= 0)
       %H(j,m) = R(inn,1);
       %tic
       nr = [nr; j]; nc = [nc; row_indx]; nv = [nv; R(inn,1)];
       %toc
       %counte = counte + 1; 
       %nr(counte) = j; nc(counte) = row_indx; nv(counte) = R(inn,1);
    end %if valm0
    %toc
    %H(j,m) = R(inn,1);
  end %inn
  %toc
  %toc
  %toc
end %end for j = 1:n
%counte;
indf=find(nr==0); nr(indf)=[];  nc(indf)=[]; nv(indf)=[];  
H = sparse(nr,nc,nv,n,n);

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
  [R1] = apply_hamiltonian_1body(Ri,Hfile1body,Nst);
  [R2] = apply_hamiltonian_2body(Ri,Hfile2body,Nst);
  R = [R1;R2];
  %{
   %check for repeated rows and simplify
  nrr = size(R,1);
  inds = [];
  for jj = 1:nrr
    for k = jj+1:nrr
      if R(jj,2:end)==R(k,2:end)
        inds =[inds,k];
        R(jj,1) = R(jj,1)+R(k,1);
      end % end if A(jj) == A(k)
    end % end for k = 1:nrr
  end % end for jj = 1:nrr
  R(inds,:) = [];
  %}
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
Hm = sparse(nr,nc,nv,n,n);

%Now hamiltonian for space with one more electron

n = size(BM,1);
nr = []; nc = []; nv = []; %structure arrays for the sparse hamiltonian matrix
for j = 1:n %loop over basis elements
 % j
  %tic
  Sc = BM(j,:); 
  Ri = [1, Sc];
  [R1] = apply_hamiltonian_1body(Ri,Hfile1body,Nst);
  [R2] = apply_hamiltonian_2body(Ri,Hfile2body,Nst);
  R = [R1;R2];
  %R = [R1];
  %R = [R2];
   %{
   %check for repeated rows and simplify
  nrr = size(R,1);
  inds = [];
  for jj = 1:nrr
    for k = jj+1:nrr
      if R(jj,2:end)==R(k,2:end)
        inds =[inds,k];
        R(jj,1) = R(jj,1)+R(k,1);
      end % end if A(j) == A(k)
    end % end for k = 1:nrr
  end % end for jj = 1:nrr
  R(inds,:) = [];
   %}
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
HM = sparse(nr,nc,nv,n,n);

else
HM=zeros(2); Hm=zeros(2);
end % end if ifBM_Bm == 1

end  %end function create_hamiltonian_matrix
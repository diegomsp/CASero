function[Fm,FM,QPm,QPM] = Dyson_Orbitals(C,Cm,CM,B,Bm,BM,AE)

Nst = numel(B(1,:));
Nel = sum(abs(B(1,:)));

Fm = zeros(1,2*Nst)'; 

FM = zeros(1,2*Nst)';

%Fm

  for j = 1:2*Nst
            %for k = j:2*Nst
               vv = 1:2*Nst;  vv([j]) = []; 
               combss = nchoosek(vv,Nel-1);
               for se = 1:size(combss,1)
                  seq1 = [j,combss(se,:)]; seq2 = [combss(se,:)];
                  % find out sign of the permutations, sg1 and sg2
                  [seq1,iseq1] = sort(seq1); [seq2,iseq2] = sort(seq2);
                  I = eye(length(iseq1));
                  sg1 = det(I(:,iseq1));
                  I = eye(length(iseq2));
                  sg2 = det(I(:,iseq2));
                  % transform into representation [ 2 -1 0 1 2 ... ], b1 and b2
                  b1 = zeros(1,2*Nst);
                  b2 = zeros(1,2*Nst);
                  b1(seq1) = 1; b1u = b1(1:Nst); b1d = b1(Nst+1:2*Nst); 
                  b1 = b1u - 2*b1d; b1(b1==-1) = 2; b1(b1==-2) = -1;   
                  b2(seq2) = 1; b2u = b2(1:Nst); b2d = b2(Nst+1:2*Nst); 
                  b2 = b2u - 2*b2d; b2(b2==-1) = 2; b2(b2==-2) = -1;
                  % find corresponding indices, ib1 and ib2
                  [~,ib1] = ismember(int8(b1),B,'rows');  %OJO!! ismember SLOW
                  [~,ib2] = ismember(int8(b2),Bm,'rows');
                  %add element:  
                 
                  if (ib1 ~= 0 && ib2 ~= 0)
                     Fm(j) = Fm(j) + sg1*sg2*C(ib1)*Cm(ib2);
                  end % end if (ib1 ~= 0 && ib2 ~= 0)\
               end % se
            %end % k
  end % j
   prefac = 1;%(1/((factorial(Nel-1))));
   Fm = prefac*Fm;
         
%FM

  for j = 1:2*Nst
            %for k = j:2*Nst
               vv = 1:2*Nst;  vv([j]) = []; 
               combss = nchoosek(vv,Nel);   %!
               for se = 1:size(combss,1)
                  seq1 = [j,combss(se,:)]; seq2 = [combss(se,:)];
                  % find out sign of the permutations, sg1 and sg2
                  [seq1,iseq1] = sort(seq1); [seq2,iseq2] = sort(seq2);
                  I = eye(length(iseq1));
                  sg1 = det(I(:,iseq1));
                  I = eye(length(iseq2));
                  sg2 = det(I(:,iseq2));
                  % transform into representation [ 2 -1 0 1 2 ... ], b1 and b2
                  b1 = zeros(1,2*Nst);
                  b2 = zeros(1,2*Nst);
                  b1(seq1) = 1; b1u = b1(1:Nst); b1d = b1(Nst+1:2*Nst); 
                  b1 = b1u - 2*b1d; b1(b1==-1) = 2; b1(b1==-2) = -1;   
                  b2(seq2) = 1; b2u = b2(1:Nst); b2d = b2(Nst+1:2*Nst); 
                  b2 = b2u - 2*b2d; b2(b2==-1) = 2; b2(b2==-2) = -1;
                  % find corresponding indices, ib1 and ib2
                  [~,ib1] = ismember(int8(b1),BM,'rows');  %!
                  [~,ib2] = ismember(int8(b2),B,'rows');   %!
                  %add element:  
                  
                  if (ib1 ~= 0 && ib2 ~= 0)
                     FM(j) = FM(j) + sg1*sg2*CM(ib1)*C(ib2);  %!
                  end % end if (ib1 ~= 0 && ib2 ~= 0)\
               end % se
            %end % k
  end % j
 prefac = 1;%(1/((factorial(Nel))));
 FM = prefac*FM;
 
 %Now calculate the final Dyson Orbitals in the basis of atomic sites
 
 QPM=FM; QPm=Fm;

 FM = blkdiag(AE,AE)*FM;
 
 Fm = blkdiag(AE,AE)*Fm;
         
          
end %end function Dyson_Orbitals
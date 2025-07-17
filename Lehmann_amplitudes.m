function[F] = Lehmann_amplitudes(C1,C2,B1,B2,j,s)

% CHARGE = 1 OR -1
Nst = numel(B1(1,:));
Nel = sum(abs(B1(1,:)));
Nel2 = sum(abs(B2(1,:)));

if (Nel2 > Nel)
    CHARGE = 1;
else % Nel2 < Nel
    CHARGE = -1;
end

F = 0;
j = j + Nst*(s==-1);
  if (CHARGE==-1)
               vv = 1:2*Nst;  vv(j) = [];
               combss = nchoosek(vv,Nel-1);
               for se = 1:size(combss,1)
                  seq1 = [j,combss(se,:)]; seq2 = [combss(se,:)];
                  % find out sign of the permutation, sg1
                  [seq1,iseq1] = sort(seq1);
                  I = eye(length(iseq1));
                  sg1 = det(I(:,iseq1));
                  % transform into representation [ 2 -1 0 1 2 ... ], b1 and b2
                  b1 = zeros(1,2*Nst);
                  b2 = zeros(1,2*Nst);
                  b1(seq1) = 1; b1u = b1(1:Nst); b1d = b1(Nst+1:2*Nst);
                  b1 = b1u - 2*b1d; b1(b1==-1) = 2; b1(b1==-2) = -1;
                  b2(seq2) = 1; b2u = b2(1:Nst); b2d = b2(Nst+1:2*Nst);
                  b2 = b2u - 2*b2d; b2(b2==-1) = 2; b2(b2==-2) = -1;
                  % find corresponding indices, ib1 and ib2
                  [~,ib1] = ismember(int8(b1),B1,'rows');  %OJO!! ismember SLOW
                  [~,ib2] = ismember(int8(b2),B2,'rows');
                  %add element:

                  if (ib1 ~= 0 && ib2 ~= 0)
                     F = F + sg1*C1(ib1)*C2(ib2);
                  end % end if (ib1 ~= 0 && ib2 ~= 0)\
               end % se


  else  % CHARGE == 1

               vv = 1:2*Nst;  vv(j) = [];
               combss = nchoosek(vv,Nel);   %!
               for se = 1:size(combss,1)
                  seq1 = [j,combss(se,:)]; seq2 = [combss(se,:)];
                  % find out sign of the permutations, sg1 and sg2
                  [seq1,iseq1] = sort(seq1);
                  I = eye(length(iseq1));
                  sg1 = det(I(:,iseq1));
                  % transform into representation [ 2 -1 0 1 2 ... ], b1 and b2
                  b1 = zeros(1,2*Nst);
                  b2 = zeros(1,2*Nst);
                  b1(seq1) = 1; b1u = b1(1:Nst); b1d = b1(Nst+1:2*Nst);
                  b1 = b1u - 2*b1d; b1(b1==-1) = 2; b1(b1==-2) = -1;
                  b2(seq2) = 1; b2u = b2(1:Nst); b2d = b2(Nst+1:2*Nst);
                  b2 = b2u - 2*b2d; b2(b2==-1) = 2; b2(b2==-2) = -1;
                  % find corresponding indices, ib1 and ib2
                  [~,ib1] = ismember(int8(b1),B2,'rows');  %!
                  [~,ib2] = ismember(int8(b2),B1,'rows');   %!
                  %add element:

                  if (ib1 ~= 0 && ib2 ~= 0)
                     F = F + sg1*C2(ib1)*C1(ib2);  %!
                  end % end if (ib1 ~= 0 && ib2 ~= 0)\
               end % se

  end % end if CHARGE


end %end function Lehmann amplitudes

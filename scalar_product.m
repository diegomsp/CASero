function[R] = scalar_product(v1,v2)

if (isempty(v1) )
  R = 0;
elseif (isempty(v2) )
  R = 0;
else

  % v1 = [w1, S1]; v2 = [w2, S2]
  aux = v1(1,:); aux = aux(2:end); L = length(aux);
  % pre-cleaning
  v1(abs(v1(:,1))<0.00001,:) = [];
  v2(abs(v2(:,1))<0.00001,:) = [];
  n = size(v1); m = size(v2); n = n(1); m = m(1);
  R = 0;
  for j = 1:n
   for k = 1:m
     if v1(j,2:end)==v2(k,2:end)
       R = R + v1(j,1)*v2(k,1);
     end %end if v1(j,2:end)==v2(k,2:end)
   end % k
  end % j
end % end if
end % end function scalar_product

function[O3] = product_of_two_1e_Hfiles(O1,O2,varargin) 

%O1;  O2;  %PRODUCT OF 1-electron OPERATORS

n1 = size(O1,1); n2 = size(O2,1);

O3 = []; %here the product is stored

for l1 = 1:n1
 
   w1 = O1(l1,1);  c1 = O1(l1,2); d1 = O1(l1,5); 
   
   for l2 = 1:n2 
     
      w2 = O2(l2,1);  c2 = O2(l2,2); d2 = O2(l2,5);  
      
      if (d1 ~= c2 && d1 ~= d2)
         O3 = [O3; w1*w2, c1, c2, d2, d1];
      end
      
      if (d1 == c2)
         O3 = [O3; w1*w2, c1, 0, 0, d2];
         if (d1 ~= d2)
            O3 = [O3; -w1*w2, c1, c2, d1, d2];
         end 
      end % (d1 == c1)
     
   end % l2

end % l1

narg = length(varargin);
if (narg==1)
   coeff =   varargin{1};
   O3(:,1) = coeff*O3(:,1);
end


O3 = clean_Hfile(O3);

end % end function product_of_two_1e_Hfiles



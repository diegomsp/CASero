function[G] = Gfunction(x,beta)
  for j = 1:length(x)
   if (exp(-beta*x)==1)
     G(j) = 1/beta;
   else
     G =  x./(1 - exp(-beta*x));
   end
   if (isnan(G(j)) )
     G(j) = 1/beta;
   end % end if
   if (isinf(G(j)) )
     G(j) = 1/beta;
   end % end if
  end % j

end % end function Gfunction

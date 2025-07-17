function[On] = multiply_operator_by_scalar(O,x) 

   O(:,1) = x*O(:,1);
   
   On = O;

end % end function multiply_operator_by_scalar
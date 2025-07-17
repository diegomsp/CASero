function[V,signo] = apply_creation_operator(i,V)
%[i j k l]  %create destroy create destroy  if spin is down, just do i -> L+i
%V %list of ones and zeros
%L=length(V)/2;  
%Vu = V(1:L);   Vd = V(L+1:end);   %spin up and spin down

signo = 1;


% create  i
signo = signo*((-1)^( sum(V(1:i-1))  ));
V(i) = 1*(V(i)==0) + 2*(V(i)~=0);


if ( max(abs(V)) == 2)
   V = [];
   signo = 0;
end


end % function apply_creation_operator
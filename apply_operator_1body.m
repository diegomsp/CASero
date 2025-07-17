function[V,signo] = apply_operator_1body(i,l,V)
%[i j k l]  %create destroy create destroy  if spin is down, just do i -> L+i
%V %list of ones and zeros
%L=length(V)/2;  
%Vu = V(1:L);   Vd = V(L+1:end);   %spin up and spin down

signo = 1;


%destroy l
signo = signo*((-1)^( sum(V(1:l-1))  ));
V(l) = 0*(V(l)==1) + 2*(V(l)~=1);


% create  i
signo = signo*((-1)^( sum(V(1:i-1))  ));
V(i) = 1*(V(i)==0) + 2*(V(i)~=0);


if ( max(abs(V)) == 2)
   V = [];
end


end % function apply_operator_1body
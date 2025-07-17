function[value] = find_from_Hfile(Hfile,i,j,k,l)
%[Hfile]= clean_Hfile(Hfile);
value = 0;
 for a = 1:length(Hfile)
    ii = Hfile(a,2); jj = Hfile(a,3); kk = Hfile(a,4); ll = Hfile(a,5);
    if (ii==i)&&(jj==j)&&(kk==k)&&(ll==l)
        value = Hfile(a,1);
        break
    end
 end     
end
function[Hfile_clean]= clean_Hfile(Hfile)
indsgarb = find(abs(Hfile(:,1))<0.00000001); Hfile(indsgarb,:)=[];

a = 1;
while (a < size(Hfile,1))
  %a
     aaa=find(ismember(Hfile(:,2:5),Hfile(a,2:5),'rows'));
     for iaa = 2:length(aaa)
        aa = aaa(iaa);
        Hfile(aaa(1),1) = Hfile(aaa(1),1) + Hfile(aaa(iaa),1);
     end %iaa 
     Hfile(aaa(2:end),:) = [];
     a = a + 1;
end %a
indsgarb = find(abs(Hfile(:,1))<0.00000001); Hfile(indsgarb,:)=[];
Hfile_clean = Hfile;
end %end function clean_Hfile
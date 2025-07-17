function[Ope]=reduce_oe_terms_of_Hfile(Ope)

indxgr=[];
for j = 1:size(Ope,1)


    if (Ope(j,3)==0)
        if (Ope(j,2)<0 || Ope(j,2)>Ope(j,5))
            indxgr=[indxgr,j];
        end
    end % if 1e OE term


end % j
Ope(indxgr,:)=[];

end %end function reduce_oe_terms_of_Hfile
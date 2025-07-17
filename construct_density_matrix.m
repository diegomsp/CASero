 function[R]=construct_density_matrix(C1,C2,B)

    Nst = size(B,2);
    %spin up
    for j = 1:Nst
        [~,~,C2a,B2a]= apply_create_anihilate_onvector(C2,B,j);  % anihilate in j with spin up
        for k = 1:Nst
            [C2b,B2b]= apply_create_anihilate_onvector(C2a,B2a,k); % create in k with spin up
            R(j,k)=scalar_product_wfs_basis(C1,C2b,B,B2b);
        end % k
    end % j
     
    %spin down
    for j = Nst+1:2*Nst
        [~,~,C2a,B2a]= apply_create_anihilate_onvector(C2,B,j);  % anihilate in j with spin down
        for k = Nst+1:2*Nst
            [C2b,B2b]= apply_create_anihilate_onvector(C2a,B2a,k); % create in k with spin down
            R(j,k)=scalar_product_wfs_basis(C1,C2b,B,B2b);
        end % k
    end % j

end % end function
function[sign_traditional] = get_normal_order(V)

    if ( max(V)==2 || min(V)==-1 )
        V1 = V; V1(V1==-1)=0; V1(V1==2)=1;
        V2 = V; V2(V2==1)=0; V2(V2==2)=1; V2=abs(V2);
        V = [V1,V2];
    end % end if

    M = length(V)/2;

    iop_mine = find(V ~= 0);
    iop_traditional = iop_mine;

    for a = 1:length(iop_mine)
        if (iop_traditional(a) > M)
            iop_traditional_temp = iop_traditional;
            iop_traditional_temp(1:a-1) = mod(iop_traditional_temp(1:a-1),M);
            iop_traditional_temp(iop_traditional_temp==0) = M;
            level = iop_traditional(a) - M;
            isopl=find(iop_traditional_temp<=level, 1, 'last' ); if (isempty(isopl)); isopl=0; end
            iop_traditional = iop_traditional([1:isopl,a,isopl+1:a-1,a+1:M]);
        end % end if iop_mine(a) > M
    end % end for a

    %now we need to detect the sign of the permutation mapping iop_total
    %into iop_order
    [~,iseq1] = sort(iop_mine); [~,iseq2] = sort(iop_traditional);
    I = eye(length(iseq1));
    sg1 = det(I(:,iseq1));
    I = eye(length(iseq2));
    sg2 = det(I(:,iseq2));

    sign_traditional = sg1*sg2;


end % end function get_normal_order

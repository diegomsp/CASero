function[Cnew] = change_to_traditional_sign_convenion(C,B)
    Cnew = C;
    for j = 1:size(B,1)
        [sign_traditional] = get_normal_order(B(j,:));
        Cnew(j) = Cnew(j)*sign_traditional;
    end % j
end % end function change_to_traditional_sign_convenion



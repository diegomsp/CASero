function[R]=main_slater_det(C0,B0,state,threshold)


R = [B0(abs(C0(:,state))>threshold,:),  C0(abs(C0(:,state))>threshold,state)];

%[Cnew] = change_to_traditional_sign_convenion( C0(abs(C0(:,state))>threshold,state),B0(abs(C0(:,state))>threshold,:));

%R = [B0(abs(C0(:,state))>threshold,:),  Cnew];

[a,ir]=sort(abs(R(:,end)));
R=R(flipud(ir),:);

end %end function main_slater_det


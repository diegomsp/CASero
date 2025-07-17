N = 5;
t = -1;
U = 20;

%%%%%%
Occ = 2;
Sz = [-4:4];
Nel = N;

Hfile = [];
for j = 1:N-1

 Hfile = [Hfile;t, j, 0, 0, j+1];
 Hfile = [Hfile;U, j, -j, -j, j];

end % j
Hfile = [Hfile;t, 1, 0, 0, N];
Hfile = [Hfile;U, N, -N, -N, N];


M = N;

% odd

    [C01,~,~,E01,~,~,B01,H01] = Many_body_hamil_S2(Hfile,M,M,1,2,0);

    [CM0,~,~,EM0,~,~,BM0,HM0] = Many_body_hamil_S2(Hfile,M+1,M,0,2,0);






function[Hfile] = build_kanamori_hamiltonian
    Hfile_CAS = load('Hfile_porf_cobalto.txt');
    %%%%%%%%%%%%%%%%%% clean Hfile %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    indsgarb = find(abs(Hfile_CAS(:,1))<0.001); Hfile_CAS(indsgarb,:)=[];
    a = 1;
    while (a < size(Hfile_CAS,1))
    %a
        aaa=find(ismember(Hfile_CAS(:,2:5),Hfile_CAS(a,2:5),'rows'));
        for iaa = 2:length(aaa)
            aa = aaa(iaa);
            Hfile_CAS(aaa(1),1) = Hfile_CAS(aaa(1),1) + Hfile_CAS(aaa(iaa),1);
        end %iaa 
        Hfile_CAS(aaa(2:end),:) = [];
        a = a + 1;
    end %a
    indsgarb = find(abs(Hfile_CAS(:,1))<0.001); Hfile_CAS(indsgarb,:)=[];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %d = [3,5,11];   %orbitals d in Hfile_CAS
    %core = [1,2,4,6,7];   %orbitals that are not d but occupied
    d = [3,5,11];
    core = [1,2,4,6,7];
    Hfile = [];
    for id = 1:length(d)
       [U] = find_from_Hfile(Hfile_CAS,d(id),-d(id),-d(id),d(id));
       Hfile = [Hfile; U,  id, -id, -id,  id];
       [e] = find_from_Hfile(Hfile_CAS,d(id),0,0,d(id));
       %%%%%%%%%%%%%%%%%%% core effect %%%%%%%%%%%%%%%%%%%%%%%%
       sid = 0;
       sidx = 0;
       for j = core
           if (d(id)>j)
               sid = sid + find_from_Hfile(Hfile_CAS,d(id),j,j,d(id));
               sidx = sidx - find_from_Hfile(Hfile_CAS,j,d(id),j,d(id));
           else
               sid = sid + find_from_Hfile(Hfile_CAS,j,d(id),d(id),j);
               sidx = sidx - find_from_Hfile(Hfile_CAS,d(id),j,d(id),j);
           end
       end % j
       e = e + 2*sid +sidx;
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       Hfile = [Hfile; e,  id, 0, 0,  id];
       for id2 = id+1:length(d)
          [J] = find_from_Hfile(Hfile_CAS,d(id2),d(id),d(id),d(id2));
          Hfile = [Hfile; J,  id,  id2,  id2,  id];
          Hfile = [Hfile; J, -id, -id2, -id2, -id];
          Hfile = [Hfile; J,  id, -id2, -id2,  id];
          Hfile = [Hfile; J, -id,  id2,  id2, -id];
          % exchange terms
          [Jx] = find_from_Hfile(Hfile_CAS,d(id),d(id2),d(id),d(id2));
          Hfile = [Hfile; Jx,  id,  id2,  id,  id2];
          Hfile = [Hfile; Jx, -id, -id2, -id, -id2];
          % spin-flip
          Hfile = [Hfile; Jx,  id, -id2, -id,  id2];
          Hfile = [Hfile; Jx, -id,  id2,  id, -id2];
          % double-hopping
          Hfile = [Hfile; Jx,  id,  -id,   -id2,  id2];
          Hfile = [Hfile; Jx,  id2, -id2,  -id,   id];
          % single hopping
          [t] = find_from_Hfile(Hfile_CAS,d(id),0,0,d(id2));
          Hfile = [Hfile; t,  id,  0,   0,  id2];
       end % id2
    end % id
end % build_kanamori_hamiltonian
function[Cgc,Cga] = create_anihilate_e(Cg,atom,spin,B,BM,Bm)

      %create/anihilate electron:
      Cgc = zeros(length(BM),1);
      Cga = zeros(length(Bm),1);

      for j = 1:length(Cg)
        %j

        S = B(j,:);

        %Create
        if (S(atom) == -spin)
          S(atom) = 2;
          [~,row_indx] = ismember(S,BM,'rows');
          if (row_indx ~= 0)
             Cgc(row_indx) = Cg(j);
          end % if row_indx ~= 0
        elseif (S(atom) == 0)
          S(atom) = spin;
          [~,row_indx] = ismember(S,BM,'rows');
          if (row_indx ~= 0)
             Cgc(row_indx) = Cg(j);
          end % if row_indx ~= 0

        %Anihilate
        elseif (S(atom) == spin)
          S(atom) = 0;
          [~,row_indx] = ismember(S,Bm,'rows');
          if (row_indx ~= 0)
             Cga(row_indx) = Cg(j);
          end % if row_indx ~= 0
        elseif (S(atom) == 2)
          S(atom) = -spin;
          [~,row_indx] = ismember(S,Bm,'rows');
          if (row_indx ~= 0)
             Cga(row_indx) = Cg(j);
          end % if row_indx ~= 0
        else
          %S = [];
        end % end if
      end % j = 1:length(Cg) (create/anihilate electron)
end %end function create_anihilate_e

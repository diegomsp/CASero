function[SC] = calculate_spin_correl(C0,B0,nstate,StatesC)

%CALCULATE SPIN CORRELATION < Psi | Si Sj | Psi  >

%i j run over list StatesC

Px = 0.5*[ 0 1; 1 0];
Py = 0.5*[0 -1i; 1i 0];
Pz = 0.5*[1 0; 0 -1];

Nst = numel(B0(1,:)); %number of sites
F = C0(:,nstate); %state of which we want to calculate the spin
B = B0;  %many-body basis set
Bu = B; Bd = B;
Bu(Bu==2) = 1; Bu(Bu==-1) = 0; 
Bd(Bd==1) = 0; Bd(Bd==2) = 1;  Bd(Bd==-1) = 1;
B = [Bu,Bd];

S = 0;

%We want to loop only over non-zero coefficients of the wave-function
tol = 0.001;
indsF = find((F.^2)>tol);

%%{


     for j = StatesC
       
       Sx(j) =  0;  Sy(j) = 0;  Sz(j) = 0;
       
      for s = 0:1
         for sp = 0:1
            
                  
                  %for b = 1:size(B,1)
                  for b = 1:numel(indsF)
                     b;
                     %for bp = 1:size(B,1)
                     b = indsF(b);
                     for bp = 1:numel(indsF)
                        bp = indsF(bp);
                        v1 = B(b,:); 
                        v2 = B(bp,:);
                        % <v1 | O | v2 >
                        %for j = 1:Nst
                           %for k = 1:Nst   
                              R = expected_value_string_of_operators([fliplr(-find(v1==1)),(j + sp*Nst),-(j + s*Nst),find(v2==1)]);
                              Sx(j) = Sx(j) + Px(s+1,sp+1)*R*F(bp)*F(b); Sy(j) = Sy(j) + Py(s+1,sp+1)*R*F(bp)*F(b); Sz(j) = Sz(j) + Pz(s+1,sp+1)*R*F(bp)*F(b);
                              %R = expected_value_string_of_operators([fliplr(-find(v1==1)),(k + sp*Nst),-(k + s*Nst),find(v2==1)]);
                              %A4 = A4 + Px(s,sp)*R*F(bp)*F(b); A5 = A5 + Py(s,sp)*R*F(bp)*F(b); A6 = A6 + Pz(s,sp)*R*F(bp)*F(b);
                           %end %k
                         %end % j
                     end %bp
                  end %b
             
         end %sp
      end % s  
    end %j


     %display('k loop')
     % 
for j = StatesC
for k = StatesC
  
  SC(j,k) = 0;
     
      for s = 0:1
         for sp = 0:1
            for t = 0:1
               for tp = 0:1
                  P = Px(s+1,sp+1)*Px(t+1,tp+1) + Py(s+1,sp+1)*Py(t+1,tp+1) + Pz(s+1,sp+1)*Pz(t+1,tp+1);
                  %for b = 1:size(B,1)
                  for b = 1:numel(indsF)
                     b;
                     %for bp = 1:size(B,1)
                     b = indsF(b);
                     for bp = 1:numel(indsF)
                        bp = indsF(bp);
                        v1 = B(b,:); 
                        v2 = B(bp,:);
                        % <v1 | O | v2 >
                        %for j = 1:Nst
                           %for k = 1:Nst   
                              R = expected_value_string_of_operators([fliplr(-find(v1==1)),(j + sp*Nst),-(j + s*Nst),k + t*Nst,-(k + tp*Nst),find(v2==1)]);
                              %S = S + P*R*F(bp)*F(b); 
                              SC(j,k) = SC(j,k) + P*R*F(bp)*F(b);
                           %end %k
                         %end % j
                     end %bp
                  end %b
               end %tp
            end %t
         end %sp
      end % s   
      
      SC(j,k) = SC(j,k) - Sx(j)*Sx(k) - Sy(j)*Sy(k) - Sz(j)*Sz(k);
      
  end %k
 end %j
      
      
       
      %}
           
      %Next calculate <Si> and <Sj>
      

  
  
      
      
end %end function calculate_total_spin
   
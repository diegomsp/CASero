function varargout = Many_body_hamil_S2(Hfile,varargin)
%function[Ct, charges_up, charges_down, E,   DOS,   w,  B,  H,  Cnat,  Enat,  Rho] = Many_body_hamil(Hfile,varargin)  
%         v1    v2          v3          v4    v5   v6   v7  v8   v9    v10    v11 
% cases varargout:
   % special cases with reduced number of outputs:
   
   %function[E] = Many_body_hamil_S2(Hfile,varargin) 
   
   %function[E,Ct,B] = Many_body_hamil_S2(Hfile,varargin) 
   
   %function[B,H] = Many_body_hamil_S2(Hfile,varargin) 

%[Ct,charges_up,charges_down,E,DOS] = Many_body_hamil_S2(Hfile,Nel,Nst,Sz,max_oc,ifBM_Bm);
%If we want DOS:
%[Ct,charges_up,charges_down,E,DOS,w] = Many_body_hamil_S2(Hfile,Nel,Nel,0,2,1);

%the optional inputs go as Nel, Nst, Sz, max_oc, ifBM_Bm, H, B, HM, BM, Hm, Bm
%Example:
%[Ct,charges_up,charges_down,E,DOS,w] = Many_body_hamil_S2(Hfile,Nel,Nst,Sz,max_oc,excts,H,B); excts = 2

% To calculate many-body DOS:
%[Ct,charges_up,charges_down,E,DOS,w,B,H,Cnat,Enat,Rho] = Many_body_hamil_S2(Hfile,Nel,Nst,Sz,max_oc,1,E_ini,E_end,nwsteps,eta,levelsR,spinR)
%++++++++++++++++++++++++++++++++++++++++++++++++++++
%            DESCRIPTION OF THE VARIABLES:
%Nel : number of electrons; 
%Nst: number of sites; 
%Sz: total z spin: can be a vector (and then all spins whose value is an element
%   of Sz will be considered)
%max_oc : maximum occupation on each levels (usually 2, 1 en models such as 
%   heisenberg chains etc
%ifBM_Bm: set to 1 if we want to calculate Density of States

%eliminate superfluos terms from the hamiltonian
%indsgarb = find(Hfile(:,1)==0); Hfile(indsgarb,:)=[];

%%%%%%%%%%%% Is Hfile a FCIDUMP file?? %%%%%%%%%%
if (Hfile(end,2)==0) % Then this is a FCIDUMP file!
    Hfile(end,:) = [];
    [Hfile]=convert_FCIDUMP_Hfile(Hfile);
end % end if This is a FCIDUMP file!

indsgarb = find(abs(Hfile(:,1))<0.00000001); Hfile(indsgarb,:)=[];
[Hfile_clean]= clean_Hfile(Hfile);

narg = length(varargin);
%nargout = length(varargout);
noutputs = max(nargout,1);
%Default values of variables:
Nst = max(max(abs(Hfile(:,2:end))));
Nel = Nst;
Sz = 0;
max_oc = 2;
ifBM_Bm = 0; 

Lanczos = 1;

% Variables for temporal propagation (if dt=0, we do not do it)
tf = 1;
dt = 0;
Cin = 0;
%......................

if (narg >=1)
   Nel = varargin{1}; 
end 
if (narg >=2)
   Nst = varargin{2}; 
end 
if (narg >=3)
   Sz = varargin{3}; 
end
if (narg >=4)
   max_oc = varargin{4}; 
end  
if (narg >=5)
   ifBM_Bm = varargin{5}; 
end 
if (narg >=6)
   H = varargin{6}; 
end 
if (narg >=7)
   B = varargin{7}; 
end 
if (narg >=8)
   HM = varargin{8}; 
end 
if (narg >=9)
   BM = varargin{9}; 
end 
if (narg >=10)
   Hm = varargin{10}; 
end
if (narg >=11)
   Bm = varargin{11}; 
end  

%In case this is a one-electron (Tight-binding calculation...):
%if (Nel == 1)  %We remove all two-body terms
%   indsgarb = find(Hfile(:,3)~=0); Hfile(indsgarb,:)=[];
%else % separate 1 body and 2 body terms in the Hfile
   inds1e = find(Hfile(:,3)==0);
   Hfile1body = Hfile(inds1e,:);
   Hfile2body = Hfile;
   Hfile2body(inds1e,:) = [];
%end % end if Nel == 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Generate Hfile of S^2 operator
SxT = []; SyT = []; SzT = [];
for site = 1:Nst
    [Sx1,Sy1,Sz1] = Pauli_matrices_Hfile(site);
    SxT = [SxT; Sx1]; SyT = [SyT; Sy1]; SzT = [SzT; Sz1];
end % site
S2file = [product_of_two_1e_Hfiles(SxT,SxT);product_of_two_1e_Hfiles(SyT,SyT);product_of_two_1e_Hfiles(SzT,SzT)];
S2file = clean_Hfile(S2file);
%S2file = [0 1 0 0 1];
inds1e = find(S2file(:,3)==0);
S2file1body = S2file(inds1e,:);
S2file2body = S2file;
S2file2body(inds1e,:) = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%The Hfile : Does it come from an open shell or closed shell computation?
os = 0; %default value
Hsites = Hfile(:,2:5);
indoe = find(Hsites(:,2)==0);
Hsitesoe = Hsites(indoe,1);
%os = 1;   %change this by hand or let this code decide by analyzing the Hfile
if (min(Hsitesoe) < 0)
   os = 1;
else 
   os = 0;
end %if (min(Hsitesoe) < 0) 


if (narg ~= 7 || ifBM_Bm == 1) 
  
   if (iscell(max_oc))
     [B,BM,Bm]=construct_basis_blocks(Nst,Nel,max_oc,Sz,ifBM_Bm);
   else 
     [B,BM,Bm]=construct_basis2(Nst,Nel,max_oc,Sz,ifBM_Bm);
   end
   
   if (os == 0)
     display('Now we build the hamiltonian matrix')
     tic
      [H,HM,Hm] = create_hamiltonian_sparse_matrix_general(Hfile1body,Hfile2body,B,BM,Bm,ifBM_Bm,Nst);
      [S2] = create_hamiltonian_sparse_matrix_general_nonreduced(S2file1body,S2file2body,B,BM,Bm,0,Nst);
      %[H,HM,Hm] = create_hamiltonian_sparse_matrix_NEW(Hfile,B,BM,Bm,ifBM_Bm,Nst);
      %[H,HM,Hm] = create_hamiltonian_sparse_matrix_general_nonreduced(Hfile,B,BM,Bm,ifBM_Bm,Nst);
      toc
      display('Done with the hamiltonian matrix')
   else 
       display('Now we build the hamiltonian matrix')
       tic
      [H,HM,Hm] = create_hamiltonian_sparse_matrix_general_nonreduced(Hfile1body,Hfile2body,B,BM,Bm,ifBM_Bm,Nst);
      [S2] = create_hamiltonian_sparse_matrix_general_nonreduced(S2file1body,S2file2body,B,BM,Bm,0,Nst);
      %[H,HM,Hm] = create_hamiltonian_sparse_matrix_NEW(Hfile,B,BM,Bm,ifBM_Bm,Nst);
      %[H,HM,Hm] = create_hamiltonian_sparse_matrix_general_nonreduced(Hfile,B,BM,Bm,ifBM_Bm,Nst);
      toc
      display('Done with the hamiltonian matrix')  
   end %end if (os == 0)
else
    if (os == 0)
      display('Now we build the hamiltonian matrix')
      tic
      BM=0; Bm=0;
      [Hplus,HM,Hm] = create_hamiltonian_sparse_matrix_general(Hfile1body,Hfile2body,B,BM,Bm,ifBM_Bm,Nst);
      [S2] = create_hamiltonian_sparse_matrix_general_nonreduced(S2file1body,S2file2body,B,BM,Bm,0,Nst);
      H = H + Hplus;
      toc
      display('Done with the hamiltonian matrix')
    else 
        display('Now we build the hamiltonian matrix')
      tic
      BM=0; Bm=0;
      [Hplus,HM,Hm] = create_hamiltonian_sparse_matrix_general_nonreduced(Hfile1body,Hfile2body,B,BM,Bm,ifBM_Bm,Nst);
      [S2] = create_hamiltonian_sparse_matri_general_nonreduced(S2file1body,S2file2body,B,BM,Bm,0,Nst);
      H = H + Hplus;
      toc
      display('Done with the hamiltonian matrix')
    end %end if (os == 0)
end %end if nargin < 6

if (narg == 7 && ifBM_Bm == 2)
  %add something new to the hamiltonian (recycle the previous one)
 
    BM = [1]; Bm = [1];
    if (os == 0)
        %B
        %size(H)
        %size(B)
        %Hfile
       [Hplus] = create_hamiltonian_sparse_matrix_general(Hfile1body,Hfile2body,B,BM,Bm,ifBM_Bm,Nst);
       %size(Hplus)
    else 
       [Hplus] = create_hamiltonian_sparse_matrix_general_openshell(Hfile,B,BM,Bm,ifBM_Bm);
    end % end if (os == 0)
    if length(Hplus) > 0
       H = H + Hplus;
    end % end if length(Hplus) > 0
end % end if narg == 7

if max(max(abs(H - H'))) > 1e-7
    max(max(abs(H - H')))
    error('Input matrix is not symmetric to working precision!');
else
    H = 0.5*(H + H');
end

if (dt == 0)  
  
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  % D I A G O N A L I Z E     H A M I L T O N I A N
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if (Lanczos == 1) 
    %H = sparse(H);
    nit_Lanczos = 300; %other possible input
    nit = min(nit_Lanczos,length(H));  %number of Lanczos iterations
    display('Now Lanczos:');
    [T, Q] = Lanczos_repo(H,nit);
    display('Lanczos done');
    [CT,E] = eig(T);
    E = diag(E);  C = Q*CT;
    Cg = C(:,1); Eg = E(1);
 else
    [C,E] = eig(H);
    E = diag(E);
    Cg = C(:,1); Eg = E(1);
 end
   %Calculate occupations:
   charges_up = zeros(1,Nst);
   charges_down = zeros(1,Nst);
  
   %Calculate spin states and add these values next to the energies in the
   %same variable (E)
   %E = [E,zeros(length(E),1)];
   eigenS = diag(C\S2*C);
   sval = 0.5*(-1+sqrt(1+4*eigenS));
   E = [E,sval];

   if (noutputs > 3)
      for k = 1:length(Cg)
         bv = B(k,:);  
         for jatom = 1:Nst
            charges_up(jatom)   = charges_up(jatom)   + ((bv(jatom)==1) +(bv(jatom)==2))*(abs(Cg(k))^2);
            charges_down(jatom) = charges_down(jatom) + ((bv(jatom)==-1)+(bv(jatom)==2))*(abs(Cg(k))^2);
         end %jatom
      end %end k =1:length(Cin)
      varargout{2} = charges_up; 
      varargout{3} = charges_down;
      varargout{1} = C;
      varargout{4} = E;
   end % end if (length(varargout) > 3)
 
   if (noutputs == 1)
      varargout{1} = E;
   end %end if (length(varargout) == 1)
   if (noutputs == 3)
      varargout{1} = E;
      varargout{2} = C;
      varargout{3} = B;
   end %end if (length(varargout) == 3)
   if (noutputs == 2)
      varargout{1} = B;
      varargout{2} = H;
   end %end if (length(varargout) == 3)
 
 Ct = C;
 %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  % E N D     O F     D I A G O N A L I Z E     H A M I L T O N I A N
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
 
 else % else dt > 0
 
 %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     % T E M P O R A L     P R O P A G A T I O N     W I T H     L A N C Z O S
 %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   
    
     Nt = floor(tf/dt)+1;
     Ct(:,1) = Cin/norm(Cin);
        
        
    charges_up = zeros(Nst,Nt);
    charges_down = zeros(Nst,Nt);
    
    %initial charges
    jt = 1;
 
  for jt = 2:Nt   %loop over time
    nit_Lanczos = 300;
    nit = min(nit_Lanczos,length(H));  %number of Lanczos iterations
    [Tl,Q] = Lanczos_repo(H,nit,Cin);
    [CT,E] = eig(Tl);
     E = diag(E);  C = Q*CT;
    %Cg = C(:,1); Eg = E(1);
 
    % The key Lanczos evolution!:
    Ct(:,jt) = Q*(expm(-1i*dt*Tl)*(Q'*Ct(:,jt-1)));  %this is the evolved state
    
      % Change Hamiltonian:  If the hamiltonian evolves with time, we update it here:
  
      % End of change hamiltonian  
end % jt 

 for k = 1:length(Cg)
   bv = B(k,:);  
   for jatom = 1:Nst
        charges_up(jatom,:) = charges_up(jatom,:) + ((bv(jatom)==1)+(bv(jatom)==2))*(abs(Ct(k,:))^2);
        charges_down(jatom,:) = charges_down(jatom,:) + ((bv(jatom)==-1)+(bv(jatom)==2))*(abs(Ct(k,:))^2);
   end %jatom
 end %end k =1:length(Cin)
   
end % end if dt > 0  
      %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      % E N D    O F   T E M P O R A L     P R O P A G A T I O N 
      %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
  
  %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      %      L A N C Z O S   S P E C T R A L    F U N C T I O N S
      %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
      if (ifBM_Bm == 1)
        display('Now we calculate the DOS')
          Cg = Ct(:,1);
          % +++++++++++++++++++  Calculate LDOS over list of atoms provided
          % CAREFUL! For now, only spin-up is possible... (if the system is not
          % magnetic this is irrelevant)
          % Options for DOS 
          Eini= -4; Efin = -2; nwsteps = 100000;
          
          atomR=[1 2 3 4]; spinR=[1];  
          eta = 0.01;
          
          %  READ OPTIONS FROM INPUT LIST 
          
          if (narg >=6)
             Eini = varargin{6}; 
          end 
          if (narg >=7)
             Efin = varargin{7}; 
          end 
          if (narg >=8)
             nwsteps = varargin{8}; 
          end 
          if (narg >=9)
             eta = varargin{9}; 
          end 
          if (narg >=10)
             atomR = varargin{10}; 
          end 
          if (narg >=11)
             spinR = varargin{11}; 
          end
          w = linspace(Eini,Efin,nwsteps);  % OJO!! depends strongly on the system...
         
          %%%% done reading options
          G = 0*w;
          
          for iatom = 1:length(atomR)       
             for ispin = 1:length(spinR)
            
                atom=atomR(iatom);  spin = spinR(ispin);
                %++++++++ CHANGE THIS FOR THE GENERAL CASE!!!!!!!!
                [Cgc,Cga] = create_anihilate_e(Cg,atom,spin,B,BM,Bm);
                % Generalize for the implicit lanczos case!!!
                
                 mLanczos_M = min(80,length(HM)); %length(HM);
                 mLanczos_m = min(80,length(Hm));
         
                 [TM,QM] = Lanczos_repo(HM,mLanczos_M,Cgc);
                 [Tm,Qm] = Lanczos_repo(Hm,mLanczos_m,Cga);
         
                  aM = diag(TM); bM = diag(TM,1); bM = [bM;0];
                  am = diag(Tm); bm = diag(Tm,1); bm = [bm;0];
      
         fraccM = 0;
         fraccm = 0;
          
         
         for k = 1:length(TM)-1
           
            fraccM = (bM(length(TM)-k)^2)./(w+Eg+1i*eta-aM(length(TM)-k+1)-fraccM);
            
         end % k = 1:length(TM)
       
         for k = 1:length(Tm)-1
           
           fraccm = (bm(length(Tm)-k)^2)./(w-Eg+1i*eta+am(length(Tm)-k+1)-fraccm);
           
         end % k = 1:length(Tm)
         %This is atually just the trace of the Green function. We will need 
         %more work to actually compute the Green function
         nu = charges_up(atom);
         G = G + (nu)./(w-Eg+1i*eta+am(1)-fraccm) + (1-nu)./(w+Eg+1i*eta-aM(1)-fraccM);
       end %ispin
     end %iatom
         DOS = (-1/pi)*imag(G);
         trapz(w,DOS)
         %DOS
      else % else if exts
         DOS=0; w=0;
      end % end if excts == 1
      if (noutputs >= 5)
         varargout{5} = DOS;
         varargout{6} = w;
      end %
      
      %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      %     E N D   O F   L A N C Z O S   S P E C T R A L    F U N C T I O N S
      %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
      
      
      %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      %   N A T U R A L   O R B I T A L S    
      %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      %E(2)-E(1)
      Cg = C(:,1);
        Cnat = 0; Enat = 0; Rho = 0;
        Calculate_Rho_NOrb = (noutputs >=9);
        %if (Sz ~= 0)
           %Calculate_Rho_NOrb = 0;
        %end
        if (Calculate_Rho_NOrb == 1)
        %%{
        display('Calculating 1P Density Matrix');
         Rho = zeros(2*Nst);
         
         for j = 1:2*Nst
            for k = j:2*Nst
               vv = 1:2*Nst;  vv([j,k]) = []; 
               combss = nchoosek(vv,Nel(1)-1);
               for se = 1:size(combss,1)
                  seq1 = [j,combss(se,:)]; seq2 = [k,combss(se,:)];
                  % find out sign of the permutations, sg1 and sg2
                  [seq1,iseq1] = sort(seq1); [seq2,iseq2] = sort(seq2);
                  I = eye(length(iseq1));
                  sg1 = det(I(:,iseq1));
                  I = eye(length(iseq2));
                  sg2 = det(I(:,iseq2));
                  % transform into representation [ 2 -1 0 1 2 ... ], b1 and b2
                  b1 = zeros(1,2*Nst); b2 = zeros(1,2*Nst);
                  b1(seq1) = 1; b1u = b1(1:Nst); b1d = b1(Nst+1:2*Nst); 
                  b1 = b1u - 2*b1d; b1(b1==-1) = 2; b1(b1==-2) = -1;   
                  b2(seq2) = 1; b2u = b2(1:Nst); b2d = b2(Nst+1:2*Nst); 
                  b2 = b2u - 2*b2d; b2(b2==-1) = 2; b2(b2==-2) = -1;
                  % find corresponding indices, ib1 and ib2
                  [~,ib1] = ismember(int8(b1),B,'rows');  %OJO!! ismember SLOW
                  [~,ib2] = ismember(int8(b2),B,'rows');
                  %add element:  
                  prefac = 1;%(factorial(Nel-1));%*(1/((factorial(Nel))^2));
                  if (ib1 ~= 0 && ib2 ~= 0)
                     Rho(j,k) = Rho(j,k) + sg1*sg2*Cg(ib1)*Cg(ib2);
                     Rho(k,j) = Rho(j,k);
                  end % end if (ib1 ~= 0 && ib2 ~= 0)\
               end % se
            end % k
         end % j
         Rho = prefac*Rho;
         %Construct RHO. Second Option
         %Rho = zeros(2*Nst);
         
         %{
         for ib = 1 : length(Cg)
            %indxs = zeros(1,Nel);
            b = B(:,ib); 
            bu=b; bd=b; bu(bu==2)=1; 
            bu(bu==-1)=0; 
            bd(bd==1)=0; bd(bd==2)=1; bd(bd==-1)==1; 
            b = [bu,bd];
            indxs = find(b==1);
         end % end for ib = 1 : length(Cg)
         %}
         
         
          display('Construct Natural Orbitals');
          %[Cnat,Enat] = eig(Rho);  Enat = diag(Enat);
          [Cnat,Enat] = eig(Rho(1:Nst,1:Nst)+Rho(Nst+1:2*Nst,Nst+1:2*Nst));  Enat = diag(Enat);
          end % end if (Calculate_Rho_NOrb == 1)
         %}
          %B,  H,  Cnat,  Enat,  Rho
       if (noutputs >= 7)
         varargout{7} = B;
       end %
       if (noutputs >= 8)
         varargout{8} = H;
       end %
       if (noutputs >= 9)
         varargout{9} = Cnat;
       end %
       if (noutputs >= 10)
         varargout{10} = Enat;
       end %
       if (noutputs >= 11)
         varargout{11} = Rho;
       end %
      %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      %   E N D   O F   N A T U R A L   O R B I T A L S    
      %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
end % end function Many_body_hamil
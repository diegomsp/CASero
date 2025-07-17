
t = -3;  U = 4;
% how many occupied and unoccupied states we plot in the spectra:
statesplot_unoccupied = 3;
statesplot_occupied   = 3;

[Chu,Ehu,nhu,Chf,Ehf,nhf,Cu,Cd,Eu,Ed,nu,nd,Nel,SSB,Ehf_shift,Eu_shift] = Hu_RHF_UHF(G,t,U);


% Huckel orbitals
[a] = draw_orbital_save(G,Chu,Nel,'HOMOhu');
[a] = draw_orbital_save(G,Chu,Nel+1,'LUMOhu');
[a] = draw_orbital_save(G,Chu,Nel-1,'HOMO1hu');
[a] = draw_orbital_save(G,Chu,Nel+2,'LUMO1hu');
[a] = draw_orbital_save(G,Chu,Nel-2,'HOMO2hu');
[a] = draw_orbital_save(G,Chu,Nel+3,'LUMO2hu');

% restricted Hartree Fock orbitals
[a] = draw_orbital_save(G,Chf,Nel,'HOMOhf');
[a] = draw_orbital_save(G,Chf,Nel+1,'LUMOhf');
[a] = draw_orbital_save(G,Chf,Nel-1,'HOMO1hf');
[a] = draw_orbital_save(G,Chf,Nel+2,'LUMO1hf');
[a] = draw_orbital_save(G,Chf,Nel-2,'HOMO2hf');
[a] = draw_orbital_save(G,Chf,Nel+3,'LUMO2hf');

% unrestricted Hartree Fock orbitals
 % spin up channel
[a] = draw_orbital_save(G,Cu,Nel,'HOMOu');
[a] = draw_orbital_save(G,Cu,Nel+1,'LUMOu');
[a] = draw_orbital_save(G,Cu,Nel-1,'HOMO1u');
[a] = draw_orbital_save(G,Cu,Nel+2,'LUMO1u');
[a] = draw_orbital_save(G,Cu,Nel-2,'HOMO2u');
[a] = draw_orbital_save(G,Cu,Nel+3,'LUMO2u');
 % spin down channel
[a] = draw_orbital_save(G,Cd,Nel,'HOMOd');
[a] = draw_orbital_save(G,Cd,Nel+1,'LUMOd');
[a] = draw_orbital_save(G,Cd,Nel-1,'HOMO1d');
[a] = draw_orbital_save(G,Cd,Nel+2,'LUMO1d');
[a] = draw_orbital_save(G,Cd,Nel-2,'HOMO2d');
[a] = draw_orbital_save(G,Cd,Nel+3,'LUMO2d');

 % spin density
 SD = zeros(length(G),1);
 for j = 1:Nel
   SD = SD + abs(Cu(:,j) - Cd(:,j)).^2;
 end % j
 [a] = draw_orbital_lessradius_save(G,SD,1,'SpinDensity');

% Spectra:

f = figure('visible','off')
%Nel_oe=Nel_oe+1
plot([(Nel-statesplot_occupied:Nel+statesplot_unoccupied)],Ehu(Nel-statesplot_occupied:Nel+statesplot_unoccupied),'+','LineWidth',6); line([Nel-statesplot_occupied Nel+statesplot_unoccupied],[0 0])
set(gca, "linewidth", 4, "fontsize", 12)
xlabel('State number','FontSize', 32)
ylabel('Energy (eV)','FontSize', 32)
saveas (f,'Huckel_spectrum_ZOOM','jpg')

f = figure('visible','off')
%Nel_oe=Nel_oe+1
plot([(Nel-statesplot_occupied:Nel+statesplot_unoccupied)],Ehf(Nel-statesplot_occupied:Nel+statesplot_unoccupied),'+','LineWidth',6); line([Nel-statesplot_occupied Nel+statesplot_unoccupied],[0 0])
set(gca, "linewidth", 4, "fontsize", 12)
xlabel('State number','FontSize', 32)
ylabel('Energy (eV)','FontSize', 32)
saveas (f,'RHF_spectrum_ZOOM','jpg')

f = figure('visible','off')
%Nel_oe=Nel_oe+1
plot([(Nel-statesplot_occupied:Nel+statesplot_unoccupied)],Eu(Nel-statesplot_occupied:Nel+statesplot_unoccupied),'+','LineWidth',6); line([Nel-statesplot_occupied Nel+statesplot_unoccupied],[0 0])
hold on; plot([(Nel-statesplot_occupied:Nel+statesplot_unoccupied)],Ed(Nel-statesplot_occupied:Nel+statesplot_unoccupied),'+','LineWidth',6);
set(gca, "linewidth", 4, "fontsize", 12)
xlabel('State number','FontSize', 32)
ylabel('Energy (eV)','FontSize', 32)
saveas (f,'UHF_spectrum_ZOOM','jpg')

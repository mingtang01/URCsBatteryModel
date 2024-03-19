% URCs_TrialScript.m provides two application examples of the URCs
%      analytical model described in:
%      "A Fast Analytical Model for Predicting the Rate Performance of
%      Battery Cells Under Mixed Kinetic Control", Hongxuan Wang, Fan Wang,
%      Harsh Agarwal, Ryan Stephens, and Ming Tang
%
% Example 1: Full Cell 4um Rate Performance
%      Generates a simplified version of Figure 3d in the above paper,
%      which predicts the rate performance of NMC/Gr full cells (r_cat=4um)
%      with varying cathode thickness underoing galvanostatic discharge at
%      different C rates.
%
% Example 2: Half Cell 4um Qw/Ew      
%      Generates simplified versions of Figures 5a and 5c in the above
%      paper, which plots the specific discharge capacity Qw and discharge
%      energy Ew as a function of cathode thickness and porosity for
%      NMC/Gr half cells (r_cat=4um) undergoing 1C discharge.
%
% Created by: Hongxuan Wang, Fan Wang, and Ming Tang
%             Department of Materials Science & NanoEngineering 
%             Rice University, Houston, TX, USA
% Contact: mingtang@rice.edu
% Date: 03/01/2024

%% Full Cell 4um Rate Performance
clear
param=DefineParam_NMC_Graphite;
anode_type='Gr';
Lc_list = [70, 100, 150, 200, 250, 300]*1e-6;
logC_list = linspace(-1, 1, 50);
C_list = 10.^logC_list;
RatePerf = zeros(size(Lc_list, 1)*size(C_list, 1), 3);
epsc = 0.25;
epsa = 1-1.1*(1-epsc)*0.55*49761/31507/1.15;
counter = 0;

for Lc = Lc_list
    La = 1.15*Lc;
    for Crate = C_list
        Itot = param.F*(param.csmaxc-param.cs0c)*(1-epsc)*Lc/(3600/Crate); 
        [DoD, exitflag, Uout, ~, ~, ~, ~] = URCs_Model_Discharge_ConstEps(anode_type, param, Itot, Lc, epsc, La, epsa);
        counter = counter + 1;
        RatePerf(counter, :) = [Lc, Crate, DoD];
    end
end

% plot
figure
for Lc = Lc_list
   RatePerf_Lc_curr = RatePerf(RatePerf(:, 1)==Lc, :);
   semilogx(RatePerf_Lc_curr(:, 2), RatePerf_Lc_curr(:, 3), LineWidth=3);
   hold on
end
xlabel('{\it C}')
ylabel('DoDf')
title('Rate Performance: NMC/Gr Full Cell {\it r_{cat}}=4\mum')

%% Half Cell 4um Qw/Ew
clear
param=DefineParam_NMC_Graphite;
anode_type='Li';
n_Lc = 100;
n_epsc= 100;
Lc_list = linspace(50e-6, 400e-6, n_Lc);
epsc_list = linspace(0.1, 0.6, n_epsc);
[Lc_mesh, epsc_mesh] = meshgrid(Lc_list, epsc_list);
Qw_mesh = zeros(size(Lc_mesh));
Ew_mesh = zeros(size(Lc_mesh));

Crate = 1;

% constant cell params for specific Q/E
rho_NMC = 4.77;
rho_sep_s = 0.946;
rho_e = 1.3;
rho_Cu = 8.96;
rho_Al = 2.7;
rho_Li = 0.534;
rho_Gra = 2.27;
L_cc = 15e-6;
L_sep = param.Ls;
epss = param.epss;
rho_sep = epss*rho_e+(1-epss)*rho_sep_s;

for i = 1:n_Lc
    Lc = Lc_list(i);
    La = 1.15*Lc;
    for j = 1:n_epsc
        epsc = epsc_list(j);
        epsa = 1-1.1*(1-epsc)*0.55*49761/31507/1.15;
       
        % variable cell params for specific Q/E
        L_Li = 0.25*(1-epsc)*Lc*0.3557;   
        rho_cat = epsc*rho_e+(1-epsc)*rho_NMC;
        rho_an = epsa*rho_e+(1-epsa)*rho_Gra;
        
        Itot = param.F*(param.csmaxc-param.cs0c)*(1-epsc)*Lc/(3600/Crate); 
        [DoD, exitflag, Uout, ~, ~, ~, ~] = URCs_Model_Discharge_ConstEps(anode_type, param, Itot, Lc, epsc, La, epsa);
        Ecell = integral(@(x) Uout(x)*(3600/Crate)*Itot, 0, DoD)/3600; % Wh/m^2
        
        if strcmp(anode_type, 'Li')==1
            Qw = DoD*153.78*(rho_NMC*Lc*(1-epsc))/(rho_cat*Lc+rho_Li*L_Li+...
                rho_sep*L_sep+(rho_Al+rho_Cu)*L_cc/2); % mAh/g
            Ew = Ecell/(1000*(rho_cat*Lc+rho_Li*L_Li+...
                rho_sep*L_sep+(rho_Al+rho_Cu)*L_cc/2)); % Wh/kg
        elseif strcmp(anode_type, 'Gr')==1
            Qw = DoD*153.78*(rho_NMC*Lc*(1-epsc))/(rho_cat*Lc+rho_an*La+...
                rho_sep*L_sep+(rho_Al+rho_Cu)*L_cc/2); % mAh/g
            Ew = Ecell/(1000*(rho_cat*Lc+rho_an*La+...
                rho_sep*L_sep+(rho_Al+rho_Cu)*L_cc/2)); % Wh/kg
        end
        Qw_mesh(j, i) = Qw;
        Ew_mesh(j, i) = Ew;
    end
end

% plot
figure
contourf(Lc_mesh, epsc_mesh, Qw_mesh, ShowText="on")
xlabel('{\it L_{cat}} [m]')
ylabel('{\it \epsilon_{cat}}')
title('{\it Q_{w}} [mAh/g]: NMC/Li Half Cell {\it r_{cat}}=4\mum')

figure
contourf(Lc_mesh, epsc_mesh, Ew_mesh, ShowText="on")
xlabel('{\it L_{cat}} [m]')
ylabel('{\it \epsilon_{cat}}')
title('{\it E_{w}} [Wh/kg]: NMC/Li Half Cell {\it r_{cat}}=4\mum')
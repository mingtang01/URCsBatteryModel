function param = DefineParam_NMC_Graphite
%
% DefineParam_NMC_Graphite returns a MATLAB structure array that contains
%      the parameters of NMC/Li half and NMC/Graphite full cells
%      used in the discharge performance predictions based on
%      the URCs analytical model described in:
%      "A Fast Analytical Model for Predicting the Rate Performance of
%      Battery Cells Under Mixed Kinetic Control", Hongxuan Wang, Fan Wang,
%      Harsh Agarwal, Ryan Stephens, and Ming Tang
%
% Created by: Hongxuan Wang, Fan Wang, and Ming Tang
%             Department of Materials Science & NanoEngineering 
%             Rice University, Houston, TX, USA
% Contact: mingtang@rice.edu
% Date: 03/01/2024


% Constants
param.F = 96485; % Faraday constant [C/mol]
param.R = 8.314; % Gas constant [J/(mol*K)]
param.T = 298;   % Temperature (K)

% Mesh size used in numerical integration (cathode, separator, anode)
param.Nc = 51; param.Ns = 11; param.Na = 51;

%%% NMC cathode properties
param.csmaxc = 49761; % Max. Li concentration in cathode acitve material [mol/m^3]
param.cs0c = 0.45 * param.csmaxc; % Initial Li concentration in cathode active material in fully charged state [mol/m^3]
param.ftauc = @(x) (1.3*x.^(-0.8)); % cathode tortuosity function
param.Rc = 4e-6; % NMC particle size [m]
param.Dsc = 1e-14; % NMC Li diffusivity [m^2/s]
param.k0c = 3e-11; % NMC exchange current rate constant [mol·m^-2·s^-1·(mol·m-3)^-1.5]
param.f_cat = 1; % volume fraction of active material among cathode solid phase
% Cathode OCV
% interpolated from table
data = table2array( readtable('NMC111_OCV.csv','HeaderLines',1) );
param.Ueqc = griddedInterpolant(data(:,1), data(:,2), 'linear');
% Inverse function of Ueq(c): c(Ueq)
param.fInvUeqc = griddedInterpolant(flipud(data(:,2)), flipud(data(:,1)), 'linear');

%%% Separator properties
param.Ls = 25e-6; % Separator thickness [m]
param.epss = 0.55; % Separator porosity
param.ftaus = @(x) (x.^(-0.5)); % separator tortuosity function

%%% Graphite anode properties
param.csmaxa = 31507; % Max. Li concentration in anode active material [mol/m^3]
param.cs0a = 28986; % Initial Li concentration in anode active material in fully charged state [mol/m^3]
param.ftaua = @(x) (2.8*x.^(-1)); % anode tortuosity function
param.Ra = 1e-6; % graphite particle radius [m]
param.Dsa = 9e-14; % Li diffusivity in graphite [m^2/s]
param.k0a = 3e-11; % Graphite exchange current rate constant [mol·m^-2·s^-1·(mol·m-3)^-1.5]
param.f_an = 1; % volume fraction of active material among anode solid phase
% Anode OCV
% analytical expression
param.Ueqa = @(x) 0.124+1.5*exp(-70*x)-0.0351*tanh((x-0.286)/0.083)-0.0045*tanh((x-0.9)/0.119)-0.035*tanh((x-0.99)/0.05)...
    -0.0147*tanh((x-0.5)/0.034)-0.102*tanh((x-0.194)/0.142)-0.022*tanh((x-0.98)/0.0164)-0.011*tanh((x-0.124)/0.0226)...
    +0.0155*tanh((x-0.105)/0.029);
% Inverse function of Ueq(c): c(Ueq)
clist=linspace(0.9999,0.0001,100);
param.fInvUeqa = griddedInterpolant(param.Ueqa(clist), clist, 'linear');

%%% Lithium metal anode properties
param.i0_Li = 20; % exchange current density [A/m^2]

%%% Electrolyte properties
param.c0 = 1000;  % lithium salt concentration [mol/m^3]
T = param.T; % temperature 
% Ambipolar diffusivity
% concentration unit: mol/m^3
param.fDamb = @(c) ((c>=0)&(c<=4000)).*(1e-4*10.^(-4.43-54./(T-(229+5e-3.*c.*((c>=0)&(c<=4000))))-0.22e-3.*c.*((c>=0)&(c<=4000))))+...
   3.8722e-11.*(c>4000)+...
   6.1290e-10.*(c<0); % [m^2/s]
% Li transference number
param.ftLi = @(c) 0.38*ones(size(c));
% Ionic conductivity
param.fkappa = @(c) (1/10*(c/1000).*(-10.5+0.668*(c/1000)+...
     0.494*(c/1000).^2+0.0740*T-0.0178*(c/1000)*T-8.86e-4*(c/1000).^2*...
     T-6.96e-5*T^2+2.8e-5*(c/1000)*T^2).^2).*((c>=0)&(c<=4000))+...
     0.0814.*(c>4000) + 0.*(c<0); % [S/m]
% Thermodynamic factor
param.fTDF = @(c) 1*ones(size(c));
% Define a function G(c) = integrate( Damb(x)/(1-tLi(x)), x=[0, c])                
param.fIDamb_tLi = @(c)integral(@(x) param.fDamb(x)./(1-param.ftLi(x)),0,c);
IDamb_list = (0:50:6e3)'; % list of electrolyte concentrations;
IDamb_list(:,2) = arrayfun(param.fIDamb_tLi,IDamb_list); % list of G(c) 
%define the inverse function of G(c) through interpolation 
param.FinvIDamb_tLi = griddedInterpolant(IDamb_list(:,2), IDamb_list(:,1),'linear');


%%% Parameters used in cell energy density calculation
param.rho_AMc = 4.77; % Cathode active material density [g/cm^3]
param.rho_AMa = 2.27; % Anode active material density [g/cm^3]
param.rho_Li = 0.534; % Lithium anode density [g/cm^3]
param.rho_CCc = 2.7; % Cathode current collector density [g/cm^3]
param.rho_CCa = 8.96; % Anode active material density [g/cm^3]
param.L_CCc = 7.5e-6; % Cathode current collector thickness [m]
param.L_CCa = 7.5e-6; % Anode current collector thickness [m]
param.rho_electrolyte = 1.3; % Electrolyte density [g/cm^3]
param.rho_sep = 0.946; % Separator density [g/cm^3]


%%% calculate the roots of tan(lm) = lm
lms = arrayfun( @(lm0) fzero( @(x) tan(x)-x, ([0.4, 0.5-1e-9]+lm0)*pi), 1:20); 
param.lmssq = lms.^2;
function [DoD, exitflag, Uout, Lpz, cesol, Philsol, cssol] = URCs_Model_Discharge_ConstEps(anode_type, param, I, Lc, epsc, La, epsa)
%
% URCs_Model_Discharge_ConstEps(anode_type, param, I, Lc, epsc, La, epsa)
%      implements the URCs analytical model for a galvanostatic
%      discharge process, assuming electrodes with constant porosity 
%      The URCs model is described in:
%      "A Fast Analytical Model for Predicting the Rate Performance of
%      Battery Cells Under Mixed Kinetic Control", Hongxuan Wang, Fan Wang,
%      Harsh Agarwal, Ryan Stephens, and Ming Tang
%      
% Output arguments:
%   DoD      - Normalized galvanostatic discharge capacity
%   exitflag - electrolyte salt availability under pseudo-steady-state
%              0: no salt depletiion (Lpz = Lc)
%              1: partial salt depletion (Lc > Lpz > 0)
%              2: complete salt depletion (Lpz = 0)
%   Uout     - Cell output voltage as a function of DoD (unit: V)
%   Lpz      - Size of salt penetration zone (unit: m)
%   cesol    - salt concentration along electrode thickness (unit: mol/m^3)
%   Philsol  - electrolyte potential along electrode thickness (unit: V)
%   cssol    - Li concentration on active material particle surface
%              along electrode thickness (unit: mol/m^3)
%              
%                
% Input arguments:
%   anode_type - Type of anode
%                'Li' (Li metal half cells)
%                'Gr' (Graphite full cells)
%   param    - MATLAB structure array containing cell parameter values.
%              See 'DefineParam_NMC_Graphite.m' for an example.
%   I        - Applied discharge current density (unit: A/m^2)
%   Lc       - Cathode thickness (unit: m)
%   epsc     - Cathode porosity
%   La       - Anode thickness (unit: m)
%   epsa     - Anode porosity
%
% Created by: Hongxuan Wang, Fan Wang, and Ming Tang
%             Department of Materials Science & NanoEngineering 
%             Rice University, Houston, TX, USA
% Contact: mingtang@rice.edu
% Date: 03/01/2024


%%% Load Parameters from input
F = param.F;    
T = param.T;
R = param.R;
% Electrolyte parameters
c0 = param.c0;
fG = param.fIDamb_tLi; % function G(c)
FinvG = param.FinvIDamb_tLi; % Inverse function of G(c) 
fkappa = param.fkappa; % ionic conductivity
fomega = @(c) 1+2*R*T/F^2*param.fkappa(c).*(1-param.ftLi(c)).^2.*...
        param.fTDF(c)./(param.fDamb(c).*c); % a term in electrolyte potential

% Electrode and separator parameters
% separator
Ls = param.Ls; 
epss = param.epss;
taus = param.ftaus(epss);
% cathode 
tauc = param.ftauc(epsc); 
f_cat = param.f_cat; % volume fraction of cathode active material in solid
cs0c = param.cs0c; % Li concentration in cathode active material
csmaxc = param.csmaxc;
fInvUeqc = param.fInvUeqc; % inverse of cathode OCV
% anode 
if strcmp(anode_type,'Gr')==1
    taua = param.ftaua(epsa);
end
f_an = param.f_an; % volume fraction of anode active material in solid
cs0a = param.cs0a; % Li concentration in anode active material
csmaxa = param.csmaxa;
fInvUeqa = param.fInvUeqa; % inverse of anode OCV

%%% Discretize domain
Xc = linspace(0, Lc, param.Nc); 
Xs = linspace(Lc, Lc+Ls, param.Ns); 
if strcmp(anode_type,'Gr')==1
    Xa = linspace(Lc+Ls, Lc+Ls+La, param.Na);
end

%%% Electrolyte Transport Module
% Integration of eps*c(x) in cathode, separator and anode regions
fintc = @(Lpz) epsc*trapz(Xc, FinvG( (tauc./epsc.*I./(2*F*Lpz).*(Xc-Lc+Lpz).^2).*(Xc>Lc-Lpz) ), 2);
fints = @(Lpz) epss*trapz(Xs, FinvG( (tauc/epsc*I/(2*F)).*Lpz + (taus/epss*I/F).*(Xs-Lc) ), 2);
if strcmp(anode_type,'Li')==1
    % objective function for finding Lpz using salt conservation 
    fIc = @(Lpz) (fintc(Lpz) + fints(Lpz) - (epsc*Lc+epss*Ls)*c0 );
elseif strcmp(anode_type,'Gr')==1
    finta = @(Lpz) epsa*trapz(Xa, FinvG( (tauc/epsc*I/(2*F)).*Lpz + taus/epss*I/F*Ls + (taua/epsa*I/F).*(Xa-Lc-Ls)-(taua/epsa*I/(2*F)/La).*(Xa-Lc-Ls).^2 ), 2);
    fIc = @(Lpz) (fintc(Lpz) + fints(Lpz) + finta(Lpz) - (epsc*Lc+epss*Ls+epsa*La)*c0);
else
    error('Unrecognized anode_type');
end

% Calculate electrolyte penetration depth Lpz
Ic0 = fIc(Lc); % If Ic0<0, then Lpz > Lc, i.e. c @ current collector>0
minPZ = 0.01; % set a minimum Lpz to avoid numerical issue.
Ic1 = fIc(minPZ*Lc); % Ic1>0 if Lpz<minPZ*Lc 

if Ic0 < 0
    Lpz = Lc;      
    exitflag = 0;
    % Integration of eps*c(x) in cathode and separator regions in the
    % case of full electrolyte penetration Lpz = Lc
    fintcFP = @(cc) epsc*trapz(Xc, FinvG( fG(cc) + tauc./epsc.*I./(2*F*Lc).*Xc.^2 ), 2); 
    fintsFP = @(cc) epss*trapz(Xs, FinvG( fG(cc) + tauc/epsc*I/(2*F).*Lc+(taus/epss*I/F).*(Xs-Lc) ), 2);
    % Find the salt concentration at current collector, cc
    if strcmp(anode_type,'Li')==1
        cc = fzero(@(cc) (fintcFP(cc) + fintsFP(cc)- (epsc*Lc+epss*Ls)*c0), [0,c0]);
    elseif strcmp(anode_type,'Gr')==1
        fintaFP = @(cc) epsa*trapz(Xa, FinvG( fG(cc) + tauc./epsc*I/(2*F).*Lc+(taus/epss*I/F).*Ls ...
            + (taua/epsa*I/F).*(Xa-Lc-Ls)-(taua/epsa*I/(2*F)/La).*(Xa-Lc-Ls).^2) ); 
        cc = fzero(@(cc) (fintcFP(cc) + fintsFP(cc) + fintaFP(cc) - (epsc*Lc+epss*Ls+epsa*La)*c0), [0,c0]);
    end
    % calculate salt concentration profile
    cesolc = FinvG( fG(cc) + tauc./epsc.*I./(2*F*Lc).*Xc.^2 );
    cesols = FinvG( fG(cc) + (tauc/epsc*I/(2*F)).*Lc+(taus/epss*I/F).*(Xs-Lc) );
    if strcmp(anode_type,'Li')==1
        cesol = griddedInterpolant([Xc(1:end-1), Xs], [cesolc(1:end-1), cesols], 'linear');
    elseif strcmp(anode_type,'Gr')==1
        cesola = FinvG( fG(cc) + (tauc/epsc*I/(2*F)).*Lc+taus/epss*I/F*Ls ...
            + (taua/epsa*I/F).*(Xa-Lc-Ls)-(taua/epsa*I/(2*F)/La).*(Xa-Lc-Ls).^2);
        cesol = griddedInterpolant([Xc(1:end-1), Xs(1:end-1), Xa], [cesolc(1:end-1), cesols(1:end-1), cesola], 'linear');
    end
elseif Ic1 > 0
    Lpz = minPZ*Lc; % lowest bound of DoD
    exitflag = 2;
    cesol=@(x) c0*ones(size(x));
elseif Ic0 >= 0 && Ic1 < 0
    options = optimset('TolX', 1e-8);
    Lpz = fzero(fIc, [Lc*minPZ Lc], options);
    exitflag = 1;
    % output electrolyte concentration solution
    cesolc = FinvG( (tauc./epsc.*I./(2*F*Lpz).*(Xc-Lc+Lpz).^2).*(Xc>Lc-Lpz) );
    cesols = FinvG( (tauc/epsc*I/(2*F)).*Lpz+(taus/epss*I/F).*(Xs-Lc) );
    if strcmp(anode_type,'Li')==1
        cesol = griddedInterpolant([Xc(1:end-1), Xs], [cesolc(1:end-1), cesols], 'linear');
    elseif strcmp(anode_type,'Gr')==1
        cesola = FinvG( (tauc/epsc*I/(2*F)).*Lpz+(taus/epss*I/F)*Ls ...
            + (taua/epsa*I/F).*(Xa-Lc-Ls)-(taua/epsa*I/(2*F)/La).*(Xa-Lc-Ls).^2);
        cesol = griddedInterpolant([Xc(1:end-1), Xs(1:end-1), Xa], [cesolc(1:end-1), cesols(1:end-1), cesola], 'linear');
    end
end


%%% Solid phase module
Rc = param.Rc; % cathode particle size
Dsc = param.Dsc;  % cathode Li diffusivity
lmssq = param.lmssq(:); % square of the roots of tan(lambda)=lambda; lmssq should be a column vector 
a_p = 3*(1-epsc)/Rc;  % volumetric particle surface area
j_in = I./(F*a_p.*Lpz); % Li flux at cathode particle surface 
tmaxPZ = (csmaxc-cs0c)*Rc./(3*j_in); % max time to fully discharge particles in PZ at current density I
tmax0 = tmaxPZ*Lc/Lpz; % nominal time to fully discharge cell assuming all particles are reacting
% Li surface concentration at cathode particle surface, t should be a scalar or row vector
csurfc = @(t, j_in) cs0c + j_in.*(3.*t/Rc + Rc./(5*Dsc) - 2*Rc/Dsc.*sum(exp(-lmssq*Dsc*t./Rc^2)./lmssq, 1)); 

% Calculate average Li concentration in cathode particle as a function of Li surface concentration
tlistcat=linspace(0, tmaxPZ, 100);
csurfcatlist = csurfc(tlistcat, j_in); 
t_vs_csurfcat = griddedInterpolant(csurfcatlist, tlistcat, 'linear');
csmeancat = @(csurf) cs0c + (j_in*3/Rc)*t_vs_csurfcat(csurf); 
if strcmp(anode_type,'Gr')==1
    Ra = param.Ra; % cathode particle size
    Dsa = param.Dsa;  % anode Li diffusivity
    a_pa = 3*(1-epsa)/Ra;  % volumetric particle surface area
    j_out = I./(F*a_pa.*La); % Li flux at anode particle surface   
    % Li surface concentration at anode particle surface
    csurfa = @(t, j_out) cs0a - j_out.*(3.*t/Ra + Ra./(5*Dsa) - 2*Ra/Dsa.*sum(exp(-lmssq*Dsa*t./Ra^2)./lmssq, 1));
    % Calculate average Li concentration in anode particle as a function of Li surface concentration
    tlistan=linspace(tmax0, 0, 100);
    csurfanlist = csurfa(tlistan, j_out); 
    t_vs_csurfan = griddedInterpolant(csurfanlist, tlistan, 'linear');
    csmeanan = @(csurf) cs0a - (j_out*3/Ra)*t_vs_csurfan(csurf); 
end

% Calculate cathode surface overpotential (negative value)
min_cl = eps; % salt concentration floor in DZ against numerical error
j0c = param.k0c/2*sqrt(max(cesol(Xc),min_cl)) * sqrt(param.csmaxc.^2 - param.cs0c.^2);
eta_cat = -2*R*T/F*asinh(j_in./(2*j0c)); 
% Calculate Li or graphtie anode surface overpotential 
if strcmp(anode_type,'Li')==1
    % calculate Li metal anode surface overpotential (positive value)
    eta_Li = 2*R*T/F*asinh(I/(2*param.i0_Li));
elseif strcmp(anode_type,'Gr')==1
    eta_Li = 0;
    % calculate graphite anode surface overpotential (positive value) 
    j0a = param.k0a/2*sqrt(max(cesol(Xa),min_cl))*param.cs0a;
    eta_an = 2*R*T/F*asinh(j_out./(2*j0a)); 
end

% Calculate liquid potential in separator and cathode
% Zero liquid potential reference location: Half cell - Li metal, full cell - separator/anode interface
% Note that eta_Li is set to 0 when anode_type=='Gr'
if exitflag==1 || exitflag==0  % Lpz < Lcat  
    % separator
    fPhils = @(x) -eta_Li-taus/epss*I*trapz(Xs, (Xs>x)./fkappa(max(min_cl, cesol(Xs))).*fomega(max(min_cl, cesol(Xs))), 2);
    Philslist = arrayfun(fPhils, Xs);
    dPhils = taus/epss*I*trapz(Xs, 1./fkappa(max(min_cl, cesol(Xs))).*fomega(max(min_cl, cesol(Xs))), 2); % liquid potential drop across separator
    % cathode
    fPhilc = @(x) -eta_Li - dPhils - tauc./epsc.*I./Lpz.*trapz(Xc, (Xc>x).*(Xc>(Lc-Lpz)).*(Xc-Lc+Lpz)./fkappa(max(min_cl, cesol(Xc))).*fomega(max(min_cl, cesol(Xc))), 2);
    Philclist = arrayfun(fPhilc, Xc);
else % Zero PZ
    Philslist = zeros(size(Xs));
    Philclist = zeros(size(Xc));
end

% Output liquid potential profile
if strcmp(anode_type,'Li')==1
    Philsol = griddedInterpolant([Xc(1:end-1), Xs], [Philclist(1:end-1), Philslist], 'linear');
elseif strcmp(anode_type,'Gr')==1
    if exitflag==1 || exitflag==0
        % calculate electrolyte potential drop across anode
        fPhila = @(x) taua./epsa.*I./La.*trapz(Xa, (Xa<x).*(Lc+Ls+La-Xa)./fkappa(max(min_cl, cesol(Xa))).*fomega(max(min_cl, cesol(Xa))), 2);
        Philalist = arrayfun(fPhila,Xa);
    else
        Philalist = zeros(size(Xa));
    end
    Philsol = griddedInterpolant([Xc(1:end-1), Xs(1:end-1), Xa], [Philclist(1:end-1), Philslist(1:end-1), Philalist], 'linear');
end

% Calculate solid Li concentration and DoD.
NU = 100; % discretize voltage
if strcmp(anode_type,'Li')==1
    Ucutoff = 3;
elseif strcmp(anode_type,'Gr')==1
    Ucutoff = 2.8;
end
Ucatlist = linspace(4.2, Ucutoff, NU);

% calculate the amount of Li intercalated into cathode as a function of Ucat
fintcscat = @(Ucat) (1-epsc)*f_cat*trapz(Xc, (Xc>Lc-Lpz).*(csmeancat(fInvUeqc(Ucat-Philclist-eta_cat)*csmaxc)-cs0c), 2);
% calculate cathode DoD corresponding to the applied Ucat
dodcatlist = arrayfun(@(U) fintcscat(U)/((1-epsc)*f_cat*(csmaxc-cs0c).*Lc), Ucatlist);
% Cathode potential as a function of DoD
Ucat= griddedInterpolant(dodcatlist(:), Ucatlist, 'linear');
if strcmp(anode_type,'Li')==1
    Uout = Ucat;
    DoD = min(1, max(0, dodcatlist(end)));
    % calculate solid concentration profile from Ueq(c(x)) = Uout - Phil(x) - eta_cat
    cssol = griddedInterpolant(Xc, (Xc>Lc-Lpz).*max(fInvUeqc(Ucutoff-Philclist-eta_cat)*csmaxc, cs0c) + (Xc<=Lc-Lpz).*cs0c, 'linear');    
elseif strcmp(anode_type,'Gr')==1
    Uanlist = linspace(0,1.5,NU);
    % claculate the amount of Li deintercalated from anode
    fintcsan = @(Uan) (1-epsa)*f_an*trapz(Xa, cs0a-csmeanan(fInvUeqa(Uan-Philalist-eta_an)*csmaxa), 2);
    % calculate cell DoD based on Li removed from anode
    dodanlist = arrayfun(@(U) fintcsan(U)/((1-epsc)*f_cat*(csmaxc-cs0c).*Lc), Uanlist);
    Uan = griddedInterpolant(dodanlist(:), Uanlist, 'linear');
    Uout = @(dod) Ucat(dod)-Uan(dod);
    DoD = min(1, max(0, fzero(@(x) Uout(x)-Ucutoff, [-0.5, 1]))); % lower bound below 0 against numerical error
    % calculate solid concentration profile in cathode and anode at the end of discharge from Ueq(c(x)) = Uout - Phil(x) - eta_cat
    cssol{1} = griddedInterpolant(Xc, (Xc>Lc-Lpz).*max(fInvUeqc(Ucat(DoD)-Philclist-eta_cat)*csmaxc, cs0c) + (Xc<=Lc-Lpz).*cs0c, 'linear');
    cssol{2} = griddedInterpolant(Xa, fInvUeqa(Uan(DoD)-Philalist-eta_an)*csmaxa, 'linear');
end
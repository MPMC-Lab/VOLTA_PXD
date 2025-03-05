%==================================================================================================
%> @file       parameter_functions.m
%> @brief      Defines various parameter functions used in the electrochemical-thermal model.
%>
%> @details    This module provides function handles to compute:
%>             - Stoichiometric coefficient,
%>             - Solid and electrolyte diffusivities,
%>             - Solid and electrolyte conductivities,
%>             - Reaction rate for intercalation,
%>             - Transference number,
%>             - Open circuit potentials (reference and adjusted),
%>             - Their derivatives (e.g., OCP derivative and entropy change),
%>             - Specific surface area.
%>
%>             These functions take in the valid index structure (vind), base parameters (param0),
%>             and sometimes additional variables or function handles (pfunc) to calculate model-specific
%>             quantities.
%>
%> @authors    Sanghyun Kim (shnkim@yonsei.ac.kr), Jung-Il Choi (jic@yonsei.ac.kr)
%> @date       March 2025
%> @version    1.0
%> @license    MIT License
%==================================================================================================

function pfunc = parameter_functions(vind, param0)
    % Define function handles for various parameter calculations
    pfunc.theta = @(cse) Stoichiometric_coefficient(vind, param0, cse);
    pfunc.Ds = @(cse, T) SolidDiffusivity(vind, param0, pfunc, cse, T);
    pfunc.De = @(epse, ce, T) ElectrolyteDiffusivity(vind, param0, epse, T, ce);

    pfunc.kappa = @(epse, ce, T) ElectrolyteConductivity(vind, param0, epse, T, ce);
    pfunc.kint = @(cse, T) ReactionRate_intercalation(vind, param0, pfunc, cse, T);
    
    pfunc.tplus = @(ce, T) Transferece_number(vind, T, ce);
    pfunc.Uref = @(cse) OpenCircuitPotential_ref(vind, pfunc, cse, param0);
    pfunc.dUdT = @(cse) EntropyChange(vind, pfunc, cse, param0);
    pfunc.Uint = @(cse, T) OpenCircuitPotential(vind, pfunc, param0, cse, T);
    pfunc.eta = @(cse, phis, phie, T) OverPotential_int(pfunc, cse, phis, phie, T);

    pfunc.as = @(epss, Rp) specific_surface(vind, epss, Rp, param0);
    pfunc.sigma = @(epss, epse, T) SolidConductivity(vind, param0, epss, epse, T);

    pfunc.dUdcse = @(cse, T) OCP_derivative(vind, pfunc, cse, param0, T);
end

%% Stoichiometric coefficient based on surface concentration
function theta = Stoichiometric_coefficient(vind, param0, cse)
    ind = vind.cs;
    csmax = param0.csmax(ind);
    cse = cse(ind);
    
    theta = zeros(vind.nmax, 1);
    theta(ind) = cse ./ csmax;
end

%% Solid phase diffusivity computation
function Ds = SolidDiffusivity(vind,param0,pfunc,cse,T)
    theta = pfunc.theta(cse);
    Tref = param0.Tref;
    R = param0.R;
    
    indne = vind.ne;
    indpe = vind.pe;
    Ds = zeros(vind.nmax, 1);

    % Set diffusivity for negative and positive electrodes
    Ds(indne) = 2.30e-14;
    Ds(indpe) = 4.80e-14;
    
    % Adjust diffusivity for the active material
    ind = vind.cs;
    T = T(ind);
    X = theta(ind);
    EaD = param0.EaD(ind);
    Ds(ind) = Ds(ind) .* exp(-EaD/R .* (1./T - 1/Tref)) .* exp(-6*(X-0.1).^5);
end

%% Electrolyte diffusivity computation
function De = ElectrolyteDiffusivity(vind,param0,epse,T,ce)
    ind = vind.ce;
    Tref = param0.Tref;
    brug = param0.brug(ind);
    epse = epse(ind);
    T = T(ind);
    ce = ce(ind);

    % Data for electrolyte diffusivity interpolation
    data = [200   3.9e-10/(1 - 200*59e-6)
            500   4.12e-10/(1 - 500*59e-6)
            800   4e-10/(1 - 800*59e-6)
            1000  3.8e-10/(1 - 1000*59e-6)
            1200  3.50e-10/(1 - 1200*59e-6)
            1600  2.68e-10/(1 - 1600*59e-6)
            2000  1.9e-10/(1 - 2000*59e-6)];
    Deint = interp1(data(:,1), data(:,2), ce, "pchip");
    
    De = zeros(vind.nmax, 1);
    De(ind) = (epse.^brug) .* Deint .* exp(-16500/8.314 .* (1./T - 1/Tref));
end

%% Solid phase conductivity computation
function sigma = SolidConductivity(vind,param0,epss,epse,T0)
    ind = vind.phis;
    sigma0 = param0.sigma;
    epss = epss(ind);
    brug = param0.brug(ind);
    
    sigma = zeros(vind.nmax, 1);
    sigma(ind) = sigma0(ind) .* (epss.^brug);
end

%% Electrolyte conductivity computation
function kappa = ElectrolyteConductivity(vind,param0,epse,T,ce)
    ind = vind.phie;
    Tref = param0.Tref;
    brug = param0.brug(ind);
    epse = epse(ind);
    T = T(ind);
    ce = ce(ind);

    % Data for electrolyte conductivity interpolation
    data = [0    1e-6
            200  0.455
            500  0.783
            800  0.935
            1000 0.95
            1200 0.927
            1600 0.78
            2000 0.60
            2200 0.515];
    kappaint = interp1(data(:,1), data(:,2), ce, "pchip");
    
    kappa = zeros(vind.nmax, 1);
    kappa(ind) = (epse.^brug) .* kappaint .* exp(-4000/8.314 .* (1./T - 1/Tref));

end

%% Reaction rate for intercalation
function kint = ReactionRate_intercalation(vind,param0,pfunc,cse,T)
    indne = vind.ne;
    indpe = vind.pe;
    Tref = param0.Tref;
    R = param0.R;
    
    kint = zeros(vind.nmax, 1);
    kint(indne) = 1.70e-11;
    kint(indpe) = 1.00e-11;

    % Adjust reaction rate for intercalation using activation energy
    ind = vind.jint;
    Eak = param0.Eak(ind);
    T = T(ind);
    kint(ind) = kint(ind) .* exp(-Eak/R .* (1./T - 1/Tref));

end

%% Transference number for electrolyte
function tplus = Transferece_number(vind,T,ce)
    ind = vind.phie;
    T = T(ind);
    ce = ce(ind);
    
    % Data for transference number interpolation
    data = [200  0.37
            500  0.322
            800  0.27
            1000 0.251
            1200 0.248
            1600 0.236
            2000 0.11];
    tplus = zeros(vind.nmax, 1);
    tplus(ind) = interp1(data(:,1), data(:,2), ce, "pchip");
end

%% Open circuit potential reference (OCP) computation
function Uref = OpenCircuitPotential_ref(vind,pfunc,cse,param0)
    theta = pfunc.theta(cse);
    indne = vind.ne;
    indpe = vind.pe;
    
    % Data for negative electrode OCP interpolation
    data = [0	1.5
            0.006649624	0.626
            0.013299247	0.508
            0.019948872	0.433
            0.026598496	0.376
            0.03324812	0.331
            0.039897743	0.294
            0.046547367	0.262
            0.053196992	0.235
            0.059846616	0.213
            0.066496239	0.202
            0.076249021	0.198
            0.088661653	0.193
            0.098414434	0.189
            0.1072806	0.184
            0.11481684	0.18
            0.122353081	0.175
            0.129002704	0.171
            0.136095637	0.166
            0.14274526	0.161
            0.150281501	0.155
            0.156931124	0.151
            0.165353982	0.146
            0.174220146	0.141
            0.183086312	0.137
            0.192395786	0.132
            0.201261951	0.128
            0.211014733	0.123
            0.223870672	0.119
            0.238499845	0.115
            0.274851122	0.11
            0.388781346	0.106
            0.431782248	0.101
            0.468133525	0.097
            0.486752471	0.092
            0.503598186	0.088
            0.516897433	0.083
            0.541279388	0.079
            0.668065551	0.074
            0.712396377	0.071
            0.799284796	0.067
            0.839182539	0.062
            0.870214118	0.058
            0.889719682	0.053
            0.903462238	0.049
            0.914544945	0.044
            0.92341111	0.04
            0.931833968	0.035
            0.938483591	0.031
            0.945133215	0.026
            0.951782839	0.02
            0.958432462	0.014
            0.965082086	0.008
            0.967741935	0.005 ];
    X = theta(indne);
    Uref_n = interp1(data(:,1), data(:,2), X, "pchip");
    
    % Data for positive electrode OCP interpolation
    data = [0.3	4.3
            0.30234274	4.29621
            0.310411138	4.27183
            0.32173126	4.2425
            0.33332314	4.21508
            0.349544275	4.17816
            0.366899296	4.14007
            0.401197015	4.07053
            0.422272308	4.03122
            0.452437434	3.97922
            0.470504648	3.94957
            0.492498295	3.91719
            0.51709707	3.88519
            0.544825745	3.857
            0.564092443	3.84099
            0.588691217	3.82042
            0.608229672	3.80867
            0.63880712	3.79233
            0.682616366	3.76812
            0.727250257	3.74332
            0.758661721	3.72577
            0.778284514	3.71004
            0.79130078	3.70218
            0.804120256	3.6925
            0.815065539	3.67859
            0.825607871	3.66468
            0.836553154	3.64835
            0.843993698	3.63082
            0.850394065	3.61147
            0.855969787	3.5891
            0.859690059	3.55887
            0.862585686	3.52683
            0.864647297	3.47847
            0.867121231	3.40834
            0.869791955	3.29106
            0.872059728	3.16169
            0.873287324	3.04563
            0.873949851	3       ];
    X = theta(indpe);
    Uref_p = interp1(data(:,1), data(:,2), X, "pchip");

    Uref = zeros(vind.nmax, 1);
    Uref(indne) = Uref_n;
    Uref(indpe) = Uref_p;
end

%% Entropy change computation for OCP temperature dependence
function dUdT = EntropyChange(vind,pfunc,cse,param0)
    theta = pfunc.theta(cse);
    indne = vind.ne;
    indpe = vind.pe;
    
    % Data for negative electrode entropy change interpolation
    data = [0	3.00E-04
            0.17	0
            0.24	-6.00E-05
            0.28	-1.60E-04
            0.5	-1.60E-04
            0.54	-9.00E-05
            0.71	-9.00E-05
            0.85	-1.00E-04
            1	-1.20E-04       ];
    X = theta(indne);
    dUdT_n = interp1(data(:,1), data(:,2), X, "pchip");
    
    % Set constant entropy change for positive electrode
    X = theta(indpe);
    dUdT_p = -10/96485 * ones(size(X));
    
    dUdT = zeros(vind.nmax, 1);
    dUdT(indne) = dUdT_n;
    dUdT(indpe) = dUdT_p;
end

%% OCP derivative with respect to the surface concentration
function dUdcse = OCP_derivative(vind,pfunc,cse,param0,T)
    ind = vind.cs;
    
    perturb = 1e-4;
    csep = cse * (1 + perturb);
    csem = cse * (1 - perturb);
    dcse = csep - csem;
    dU = pfunc.Uint(csep, T) - pfunc.Uint(csem, T);
    dUdcse = zeros(vind.nmax, 1);
    dUdcse(ind) = dU(ind) ./ dcse(ind);
end

%% Open Circuit Potential (OCP) calculation
function U = OpenCircuitPotential(vind,pfunc,param0,cse,T)
    Uref = pfunc.Uref(cse);
    dUdT = pfunc.dUdT(cse);
    
    ind = vind.cs;
    Uref = Uref(ind);
    dUdT = dUdT(ind);
    T = T(ind);
    Tref = param0.Tref;
    
    U = zeros(vind.nmax, 1);
    U(ind) = Uref + (T - Tref) .* dUdT;
end

%% Interfacial overpotential for intercalation
function eta = OverPotential_int(pfunc, cse, phis, phie, T)
    Uint = pfunc.Uint(cse, T);
    eta = phis - phie - Uint;
end

%% Specific surface area calculation
function as = specific_surface(vind, epss, Rp, param0)
    ind = vind.cs;
    epss = epss(ind);
    Rp = Rp(ind);
    
    as = zeros(vind.nmax, 1);
    as(ind) = 3 * epss ./ Rp;
end
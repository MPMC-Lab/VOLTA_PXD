%==================================================================================================
%> @file       init_Parameters.m
%> @brief      Initialization of parameters and initial conditions for the electrochemical-thermal model.
%>
%> @details    This module defines constants and model parameters for the electrochemical-thermal model.
%>             It includes definitions of physical constants (e.g., Faraday and gas constants), reference 
%>             and ambient temperatures, electrochemical-thermal model parameters, and initial conditions 
%>             for temperature, concentration, and potentials. Additional functions are provided for 
%>             initializing concentration and potential profiles, as well as for naming variable domains.
%>
%> @authors    Sanghyun Kim (shnkim@yonsei.ac.kr), Jung-Il Choi (jic@yonsei.ac.kr)
%> @date       March 2025
%> @version    1.0
%> @license    MIT License
%==================================================================================================

function init = init_Parameters
    % Define physical constants
    param.F = 96485;          % Faraday constant [C/mol]
    param.R = 8.314;          % Gas constant [J/(mol K)]
    
    % Temperature settings
    param.Tref = 298.15;      % Reference temperature [K]
    param.Tamb = 298.65;      % Ambient temperature [K]
    param.Tinit = 298.65;     % Initial temperature [K]
    
    % Electrochemical-Thermal model parameters
    param.h = 0.28;                                     % Heat exchange coefficient [W/(m^2 K)]
    param.tplus = [0, 0.364, 0.364, 0.364, 0]';           % Transference number
    param.alpha_int_a = 0.5;                            % Charge transfer coefficient (anodic)
    param.alpha_int_c = 0.5;                            % Charge transfer coefficient (cathodic)
    param.brug   = [0 1.5 1.5 1.5 0]';                  % Bruggeman coefficient
    param.Lx = [15 65 22 60 15]' * 1e-6;                % Thickness [m]
    param.Ly = 227 * 1e-3;                              % Height [m]
    param.Lz = [22 79 45 79 22]' * 1e-3;                % Width [m] / single-sided tab
    % param.Lz = [20 60 20 0 0]'*1e-3;                  % Width [m] / double-sided tab
    param.csmax  = [0, 31360, 0, 48867, 0]';            % Maximum concentration in the solid phase [mol/m^3]
    param.Rp     = [0, 5, 0, 7.5, 0]' * 1e-6;           % Solid particle radius [m]
    param.rho    = [8963, 1347.33, 1170.98, 2873, 2780]'; % Density [kg/m^3]
    param.Cp     = [390, 1457, 1978.16, 810, 910]';       % Specific heat capacity [J/(kg K)]
    param.Cp = (param.Cp + param.Cp / 50 * (param.Tinit - 298.15));  % Adjusted specific heat capacity
    param.lambdax  = 0.89724 * ones(1,5)';               % Thermal conductivity in x-direction [W/(m K)]
    param.lambdayz = 29.557 * ones(1,5)';                % Thermal conductivity in y and z directions [W/(m K)]
    param.EaD    = [0, 35000/3, 0, 31556/3, 0]';         % Activation energy for solid phase diffusion [J/mol]
    param.Eak    = [0, 10000, 0, 10000, 0]';             % Activation energy for reaction constant [J/mol]
    param.Ds     = [0, 2.30, 0, 4.80, 0]' * 1e-14;       % Solid diffusion coefficient [m^2/s]
    param.De     = [0, 7.5, 7.5, 7.5, 0]' * 1e-10;       % Electrolyte diffusion coefficient [m^2/s]
    param.sigma  = [5.998e7, 100, 0, 100, 3.774e7]';      % Solid phase conductivity [S/m]
    param.kappa  = [0, 0, 0, 0, 0]';                      % Electrolyte phase conductivity [S/m]
    param.kint   = [0, 1.70, 0, 1.00, 0]' * 1e-11;       % Reaction rate constant [m^(2.5)/(mol^0.5 s)]
    param.epse   = [0, 0.3558, 0.3700, 0.5767, 0]';        % Porosity
    param.epss   = [1, 0.4890, 0.0000, 0.4233, 1]';        % Volume fraction of active maaterial
    param.stoi_range.ne = [0.0066 0.6745];                % Stoichiometric range for negative electrode
    param.stoi_range.pe = [0.3088 0.8455];                % Stoichiometric range for positive electrode
    param.SOC0 = 1.00;                                  % Initial state of charge
    
    % Aging related parameters
    param.RSEI = 0;                                     % SEI resistance [Ohm m^2]

    % Initial conditions
    variable.T    = param.Tinit * ones(1,5)';            % Temperature [K]
    variable.cs   = init_concentration(param.stoi_range, param.SOC0, param.csmax); % Solid phase concentration [mol/m^3]
    variable.ce   = [0, 1150, 1150, 1150, 0]';           % Electrolyte concentration [mol/m^3]
    variable.phis = init_potential(variable.cs, param.csmax, param.Tinit);         % Solid phase potential [V]
    variable.phie = [0, 1e-6, 1e-6, 1e-6, 0]';           % Electrolyte potential [V]
    variable.jint = [0, 1e-16, 0, 1e-16, 0]';            % Pore-wall flux [mol/(m^2 s)]
    variable.jSEI = [0, 1e-12, 0, 0, 0]';                % SEI flux [mol/(m^2 s)]
    variable.jLP  = [0, 1e-12, 0, 0, 0]';                % Lithium plating flux [mol/(m^2 s)]
    variable.deltaSEI = [0, 1e-12, 0, 0, 0]';            % SEI thickness change [m/s]
    variable.deltaLP = [0, 1e-12, 0, 0, 0]';             % Lithium plating thickness change [m/s]
    variable.q = [0, 1e-12, 0, 1e-12, 0]';               % Mean ion flux in particles [mol/(m^2 s)]
    variable = name_domain(variable);

    % Package parameters and variables into the output structure
    init.param = param;
    init.variable = variable;
end

function variable = name_domain(variable)
    % Assign domain names for each variable field based on nonzero entries
    vnames = fieldnames(variable);
    vdomains = cell(size(vnames));
    for iv = 1:length(vnames)
        vn = vnames{iv};
        vdomains{iv} = find(variable.(vn));
    end
    variable.names = vnames;
    variable.domains = vdomains;
end

function cs0 = init_concentration(stoi_range,SOC0,csmax)
    % Initialize solid phase concentration based on stoichiometric ranges and SOC
    stoi0_ne = diff(stoi_range.ne) * SOC0 + stoi_range.ne(1);
    stoi0_pe = diff(stoi_range.pe) * (1 - SOC0) + stoi_range.pe(1);
    cs0  = [0, csmax(2) * stoi0_ne, 0, csmax(4) * stoi0_pe, 0]';
end

function phis0 = init_potential(cse0,csmax,T0)
    % Initialize solid phase potential based on concentration and temperature
    x = cse0(2) / csmax(2);
    y = cse0(4) / csmax(4);

    % Negative electrode potential data interpolation
    data = [0, 1.5
            0.006649624, 0.626
            0.013299247, 0.508
            0.019948872, 0.433
            0.026598496, 0.376
            0.03324812,  0.331
            0.039897743, 0.294
            0.046547367, 0.262
            0.053196992, 0.235
            0.059846616, 0.213
            0.066496239, 0.202
            0.076249021, 0.198
            0.088661653, 0.193
            0.098414434, 0.189
            0.1072806,   0.184
            0.11481684,  0.18
            0.122353081, 0.175
            0.129002704, 0.171
            0.136095637, 0.166
            0.14274526,  0.161
            0.150281501, 0.155
            0.156931124, 0.151
            0.165353982, 0.146
            0.174220146, 0.141
            0.183086312, 0.137
            0.192395786, 0.132
            0.201261951, 0.128
            0.211014733, 0.123
            0.223870672, 0.119
            0.238499845, 0.115
            0.274851122, 0.11
            0.388781346, 0.106
            0.431782248, 0.101
            0.468133525, 0.097
            0.486752471, 0.092
            0.503598186, 0.088
            0.516897433, 0.083
            0.541279388, 0.079
            0.668065551, 0.074
            0.712396377, 0.071
            0.799284796, 0.067
            0.839182539, 0.062
            0.870214118, 0.058
            0.889719682, 0.053
            0.903462238, 0.049
            0.914544945, 0.044
            0.92341111,  0.04
            0.931833968, 0.035
            0.938483591, 0.031
            0.945133215, 0.026
            0.951782839, 0.02
            0.958432462, 0.014
            0.965082086, 0.008
            0.967741935, 0.005];
    Uref_n = interp1(data(:,1), data(:,2), x, "pchip");

    % Negative electrode entropy change data interpolation
    data = [0,       3.00E-04
            0.17,    0
            0.24,   -6.00E-05
            0.28,   -1.60E-04
            0.5,    -1.60E-04
            0.54,   -9.00E-05
            0.71,   -9.00E-05
            0.85,   -1.00E-04
            1,      -1.20E-04];
    dUdT_n = interp1(data(:,1),data(:,2),x,"pchip");
    
    % Positive electrode potential data interpolation
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
    Uref_p = interp1(data(:,1),data(:,2),y,"pchip");

    % Positive electrode entropy change
    dUdT_p = -10/96485;

    % Compute potentials adjusted for temperature
    U_ne = Uref_n+(T0-298.15)*dUdT_n;
    U_pe = Uref_p+(T0-298.15)*dUdT_p;
    
    % Set initial potential profile for the cell
    phis0 = [U_ne,U_ne,0,U_pe,U_pe]';
end

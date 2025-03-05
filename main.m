%==================================================================================================
%> @file       main.m
%> @brief      Main script for running electrochemical-thermal model simulations.
%>
%> @details    This script runs simulations using both a 2D (P2D) and a 4D (P4D) model for
%>             different C-rates. The simulation results (time, cell voltage, and temperature) are
%>             plotted for comparison. The script also defines functions for time marching,
%>             adaptive time stepping, and termination conditions.
%>
%> @authors    Sanghyun Kim (shnkim@yonsei.ac.kr), Jung-Il Choi (jic@yonsei.ac.kr)
%> @date       March 2025
%> @version    1.0
%> @license    MIT License
%==================================================================================================

clear; warning off;

% Set up figure for simulation results
figure(1); clf; 
nexttile(1); hold on; nexttile(2); hold on;
crates = [0.5, 1, 2, 5];
results = cell(2, length(crates));

for i = 1:length(crates)
    crate = crates(i);

    % Run 2D model simulation (P2D)
    [mysim, time, Vcell, T, elapsed_time] = run_P2D(crate);
    figure(1);
    nexttile(1); plot(time, Vcell, ":", "LineWidth", 3, "DisplayName", "VOLTA_P2D")
    nexttile(2); plot(time, mean(T, 1) - 273.15, ":", "LineWidth", 3, "DisplayName", "VOLTA_P2D")
    results{1, i} = struct('mysim', mysim, 'time', time, 'Vcell', Vcell, 'T', T, 'elapsed_time', elapsed_time);

    % Run 4D model simulation (P4D)
    [mysim, time, Vcell, T, elapsed_time] = run_P4D(crate);    
    figure(1);
    nexttile(1); hold on; plot(time, Vcell, "-", "LineWidth", 2, "DisplayName", "VOLTA_P4D")
    nexttile(2); hold on; plot(time, T(1, :) - 273.15, "-", "LineWidth", 2, "DisplayName", "VOLTA_P4D")
    results{2, i} = struct('mysim', mysim, 'time', time, 'Vcell', Vcell, 'T', T, 'elapsed_time', elapsed_time);

end
%% Run P2D Simulation
function [mysim, time, Vcell, T, elapsed_time] = run_P2D(Crate)
    tic;
    init = init_Parameters;
    opt = opt_Simulation;
    opt.nx = [3, 13, 4, 12, 3]';
    opt.ny = 1;
    opt.nz0 = [1, 0, 0, 0, 0]';
    opt.nz = sum(opt.nz0);
    opt.Crate = Crate;
    opt.Iapp = opt.Iapp * Crate;

    mysim = init_Simulation(init, opt);
    mysim.isinit = 1;
    
    % Perform time marching until reaching a minimum of 30 steps
    while mysim.nt < 30
        mysim = time_marching(mysim);
    end
    Vcell = mysim.Vcell(mysim.nt);
    mysim.time(1) = 0;
    mysim.Vcell(1) = Vcell;
    isEnd = terminate_condition(mysim);
    
    % Switch to normal simulation settings
    mysim.isinit = 0;
    mysim.nt = 1;
    mysim.dt = mysim.opt.dt_normal;
    mysim.dt0 = mysim.dt;
    mysim.result.sys1{1} = mysim.sys1.x;  % Save system1 initial state
    mysim.result.sys1T{1} = mysim.sys2.T;   % Save system2 initial state
    
    for i = 1:2
        mysim = time_marching(mysim);
        isEnd = terminate_condition(mysim);
        mysim.result.sys1T{mysim.nt} = mysim.sys2.T;
    end
    while ~isEnd
        mysim = time_marching(mysim);
        mysim = time_stepping(mysim);
        isEnd = terminate_condition(mysim);
        mysim.result.sys1T{mysim.nt} = mysim.sys2.T;
    end
    
    elapsed_time = toc;
    
    time = mysim.time(1:mysim.nt);
    Vcell = mysim.Vcell(1:mysim.nt);
    T = cell2mat(mysim.result.sys1T);
    T = T(:, 1:mysim.nt);
end

%% Run P4D Simulation
function [mysim, time, Vcell, T, elapsed_time] = run_P4D(Crate)
    tic;
    init = init_Parameters;
    opt = opt_Simulation;
    opt.nx = [3, 13, 4, 12, 3]';
    opt.ny = 12;
    opt.nz0 = [1, 4, 2, 4, 1]';
    opt.nz = sum(opt.nz0);
    opt.Crate = Crate;
    opt.Iapp = opt.Iapp * Crate;

    mysim = init_Simulation(init, opt);
    mysim.isinit = 1;
    
    % Perform time marching until reaching a minimum of 30 steps
    while mysim.nt < 30
        mysim = time_marching(mysim);
    end
    Vcell = mysim.Vcell(mysim.nt);
    mysim.time(1) = 0;
    mysim.Vcell(1) = Vcell;
    isEnd = terminate_condition(mysim);
    
    % Switch to normal simulation settings
    mysim.isinit = 0;
    mysim.nt = 1;
    mysim.dt = mysim.opt.dt_normal;
    mysim.dt0 = mysim.dt;
    mysim.result.sys1{1} = mysim.sys1.x;  % Save system1 initial state
    mysim.result.sys1T{1} = mysim.sys2.T;   % Save system2 initial state
    
    for i = 1:2
        mysim = time_marching(mysim);
        isEnd = terminate_condition(mysim);
        mysim.result.sys1T{mysim.nt} = mysim.sys2.T;
    end
    while ~isEnd
        mysim = time_marching(mysim);
        mysim = time_stepping(mysim);
        isEnd = terminate_condition(mysim);
        mysim.result.sys1T{mysim.nt} = mysim.sys2.T;
    end
    
    elapsed_time = toc;
    
    time = mysim.time(1:mysim.nt);
    Vcell = mysim.Vcell(1:mysim.nt);
    T = cell2mat(mysim.result.sys1T);
    T = T(:, 1:mysim.nt);
end

%% Time Marching Function
function mysim = time_marching(mysim)
    mysim = sys2_update_parameter(mysim);
    mysim = sys2_update_variable(mysim);
    mysim = sys1_update_parameter(mysim);
    mysim = sys1_update_variable(mysim);
end

%% Termination Condition Function
function isEnd = terminate_condition(mysim)
    Vcell = mysim.Vcell(mysim.nt);
    cond1 = Vcell < mysim.opt.CutoffV;
    cond2 = ~isreal(Vcell);
    cond3 = isnan(Vcell);
    cond4 = Vcell > mysim.opt.CutoverV;
    isEnd = cond1 + cond2 + cond3 + cond4;
end

%% Adaptive Time Stepping Function
function mysim = time_stepping(mysim)
    tol = mysim.controller.TOL;
    b1 = mysim.controller.b1;
    b2 = mysim.controller.b2;
    b3 = mysim.controller.b3;
    eps0 = mysim.controller.eps;

    nt = mysim.nt;
    y = mysim.Vcell(nt-2:nt);
    t = mysim.time(nt-2:nt);
    y = (y - 2.0) / (3.6 - 2.0);

    dydt = diff(y) ./ diff(t);
    e1 = norm((dydt(2) - dydt(1)) * (t(3) - t(1)) / 2);
    e2 = norm(y(3) - y(2));
    eps = max(e1, e2);

    dt0 = mysim.dt;
    if eps <= tol
        % Adaptive time stepping: increase or decrease dt based on the error
        dt = (tol / eps)^0.5 * dt0;  % Coefficients chosen as b1=0.25, b2=0.25, b3=-0.25
    else
        dt = dt0 / 2;
    end
    dt = max(1, dt);
    dt = min(60, dt);

    mysim.dt0 = mysim.dt;
    mysim.dt = dt;
    mysim.controller.eps = eps;
end
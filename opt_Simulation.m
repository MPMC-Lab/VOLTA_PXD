%==================================================================================================
%> @file       opt_Simulation.m
%> @brief      Sets simulation options and initial conditions.
%>
%> @details    This function defines simulation options for the electrochemical-thermal model.
%>             It includes settings for time stepping, grid discretization, applied current, 
%>             solid phase diffusion model, voltage cutoff conditions, and adaptive time stepping.
%>
%> @authors    Sanghyun Kim (shnkim@yonsei.ac.kr), Jung-Il Choi (jic@yonsei.ac.kr)
%> @date       March 2025
%> @version    1.0
%> @license    MIT License
%==================================================================================================

function opt = opt_Simulation
    % Initial time step sizes [s]
    opt.dt_init = 1e-8;
    opt.dt_normal = 1;
    
    % Number of control volumes for each region: [ncc, ne, sep, pe, pcc]
    % For P4D settings (uncomment if needed)
    % opt.nx = [3, 13, 4, 12, 3]';
    % opt.ny = 12;
    % opt.nz0 = [1, 4, 2, 4, 1]';   % single-sided tab
    % opt.nz0 = [2, 8, 2, 0, 0]';   % double-sided tab
    % For P2D settings
    opt.nx = [3, 13, 4, 12, 3]';
    opt.ny = 1;
    opt.nz0 = [1, 0, 0, 0, 0]';

    opt.nz = sum(opt.nz0);
    
    % Applied current settings
    opt.Crate = 1;
    opt.Q = 52.4 / 26 / 2;
    opt.Iapp = -opt.Crate * opt.Q;
    opt.Iapp = opt.Iapp;
    
    % Solid phase diffusion model selection:
    % 1 - Full order (to be updated)
    % 2 - Two parameters (quadratic distribution)
    % 3 - Three parameters (biquadratic distribution)
    opt.SolidDiffusion = 3;

    % Number of grid points for full-order solid phase diffusion
    opt.nr = 10;
    if opt.SolidDiffusion > 1
        opt.nr = 1;
    end

    % Voltage cutoff conditions
    opt.CutoffV = 3.0;
    opt.CutoverV = 4.2;
    
    % Additional model options
    % 0 - Disabled
    % 1 - Enabled
    opt.modelInclude.Thermal = 1;
    opt.modelInclude.Aging = 0; % To be updated
    
    % System variables for simulation
    opt.sys1 = {'cs','ce','phis','phie','jint'}';
    opt.sys2 = {'T'}';
    if opt.SolidDiffusion == 3
        opt.sys2{end+1} = 'q';    % Mean flux
    end
    
    % Adaptive time stepping options
    opt.adaptive_time.on = 1;
    if opt.adaptive_time.on == true
        controller.TOL = 0.01;
        controller.eps = controller.TOL;
        controller.b1 = 0.25;
        controller.b2 = 0.25;
        controller.b3 = -0.25;
        opt.adaptive_time.controller = controller;
    end
end
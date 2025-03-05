%==================================================================================================
%> @file       sys2_update_parameter.m
%> @brief      Updates system 2 parameters based on the current state variables.
%>
%> @details    This function updates the parameters used in system 2 of the simulation. 
%>             It interpolates the temperature between the current and previous states, and then 
%>             recalculates several key parameters (e.g., diffusivity, open-circuit potential, 
%>             overpotential, transference number, specific surface area, and conductivity) using 
%>             the provided parameter functions.
%>
%> @authors    Sanghyun Kim (shnkim@yonsei.ac.kr), Jung-Il Choi (jic@yonsei.ac.kr)
%> @date       March 2025
%> @version    1.0
%> @license    MIT License
%==================================================================================================

function mysim = sys2_update_parameter(mysim)
    % Retrieve function handles and current state variables
    pfunc = mysim.pfunc;
    sys1 = mysim.sys1;
    sys2 = mysim.sys2;
    param = mysim.param;
    
    % Get the current state of system 1 and interpolate temperature in system 2
    x = sys1.x;
    w = (-mysim.dt / 2) / mysim.dt0;
    T = (1 - w) * sys2.T + w * sys2.To;
    
    % Retrieve state variables from system 1
    cs = sys1.cs(x);
    ce = sys1.ce(x);
    phis = sys1.phis(x);
    phie = sys1.phie(x);
    cse = cs;  % Assuming surface concentration equals solid concentration
    
    % Retrieve additional parameters for the electrolyte and solid phases
    epse = param.epse;
    epss = param.epss;
    Rp = param.Rp;
    
    % Update parameters using the corresponding parameter functions
    param.Ds = pfunc.Ds(cse, T);
    param.theta = pfunc.theta(cse);
    param.kappa = pfunc.kappa(epse, ce, T);
    param.Uref = pfunc.Uref(cse);
    param.dUdT = pfunc.dUdT(cse);
    param.Uint = pfunc.Uint(cse, T);
    param.eta = pfunc.eta(cse, phis, phie, T);
    param.tplus = pfunc.tplus(ce, T);
    param.as = pfunc.as(epss, Rp);
    param.sigma = pfunc.sigma(epss, epse, T);
    
    % Save updated parameters back into the simulation structure
    mysim.param = param;
end

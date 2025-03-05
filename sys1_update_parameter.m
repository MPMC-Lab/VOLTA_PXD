%==================================================================================================
%> @file       sys1_update_parameter.m
%> @brief      Updates system 1 parameters based on the current state variables.
%>
%> @details    This function updates the parameters used in system 1 of the simulation.
%>             It interpolates between the current and previous states, then calculates updated
%>             parameters (e.g., diffusivities, conductivities, reaction rates, open-circuit
%>             potential, and overpotential) using the provided parameter functions.
%>
%> @authors    Sanghyun Kim (shnkim@yonsei.ac.kr), Jung-Il Choi (jic@yonsei.ac.kr)
%> @date       March 2025
%> @version    1.0
%> @license    MIT License
%==================================================================================================

function mysim = sys1_update_parameter(mysim)
    % Retrieve parameter function handles and current system state
    pfunc = mysim.pfunc;
    sys1 = mysim.sys1;
    sys2 = mysim.sys2;
    param = mysim.param;
    
    % Compute weighting factor for interpolation between current and initial state
    w = (-mysim.dt / 2) / mysim.dt0;
    x = (1 - w) * sys1.x + w * sys1.xo;
    
    % Update state variables for system 1 using the interpolated value
    cs = sys1.cs(x);
    ce = sys1.ce(x);
    phis = sys1.phis(x);
    phie = sys1.phie(x);
    T = sys2.T;
    cse = cs;  % Assume surface concentration equals solid concentration

    % Retrieve additional parameters for electrolyte and solid phases
    epse = param.epse;
    epss = param.epss;
    Rp = param.Rp;
    
    % Update parameters using the corresponding parameter functions
    param.theta = pfunc.theta(cse);
    param.Ds = pfunc.Ds(cse, T);
    param.De = pfunc.De(epse, ce, T);
    param.kappa = pfunc.kappa(epse, ce, T);
    param.kint = pfunc.kint(cse, T);
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

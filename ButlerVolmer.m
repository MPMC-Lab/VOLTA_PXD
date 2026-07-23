%==================================================================================================
%> @file       ButlerVolmer.m
%> @brief      Calculates pore-wall flux and its derivatives.
%>
%> @details    This function computes the pore-wall flux using a formulation similar to the 
%>             Butler-Volmer equation. It takes into account the overpotential, concentration 
%>             terms, temperature, and internal resistance effects. Additionally, the function 
%>             computes the partial derivatives of the flux with respect to surface potentials and 
%>             concentrations, which are used in sensitivity analysis and parameter updates.
%>
%> @authors    Sanghyun Kim (shnkim@yonsei.ac.kr), Jung-Il Choi (jic@yonsei.ac.kr)
%> @date       March 2025
%> @version    1.0
%> @license    MIT License
%==================================================================================================

function BV = ButlerVolmer(vind, param0, pfunc, kint, cse, ce, phis, phie, T, j)
    % Extract constants from parameters
    R = param0.R;
    F = param0.F;
    RSEI = param0.RSEI;
    csmax = param0.csmax;
    
    % Calculate derivative of the surface potential with respect to cse
    dUdcse = pfunc.dUdcse(cse, T);
    
    % Calculate overpotential using the provided parameter function
    eta = pfunc.eta(cse, phis, phie, T);
    
    % Adjust overpotential for non-electrolyte regions using internal resistance
    indne = vind.ne;
    eta(indne) = eta(indne) - j(indne) .* RSEI;
    
    % Compute concentration difference and its square root term for the reaction
    cdiff = csmax - cse;
    csqrt = sqrt(ce .* cdiff .* cse);
    
    % Calculate the baseline pore-wall flux (analogous to the exchange current density)
    j0 = kint .* csqrt;
    
    % Compute exponential term based on the overpotential
    jexp = 2 * sinh((F ./ (2 * R * T)) .* eta);
    
    % Compute pore-wall flux
    BV.jint = j0 .* jexp;
    
    
    % Derivative of the exponential term with respect to overpotential
    djexpdeta = (F ./ (2 * R * T)) * 2 .* cosh((F ./ (2 * R * T)) .* eta);

    % Calculate the derivatives of the pore-wall flux with respect to surface potentials
    % and concentrations. The concentration sensitivities of the exchange-current
    % prefactor are not carried in the single-pass Jacobian; the linearization uses
    % the overpotential sensitivities and the open-circuit-potential slope dU/dcse.
    BV.djdphis = j0 .* djexpdeta;
    BV.djdphie = -BV.djdphis;
    BV.djdce = zeros(size(j0));
    BV.djdcse = -j0 .* djexpdeta .* dUdcse;
end

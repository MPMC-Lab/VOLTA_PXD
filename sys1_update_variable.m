%==================================================================================================
%> @file       sys1_update_variable.m
%> @brief      Updates system 1 state variables (solid and electrolyte properties, potentials, and fluxes)
%>             based on the current simulation state.
%>
%> @details    This function updates the state variables for system 1 by computing time derivatives,
%>             diffusion terms, and charge conservation equations for the solid and electrolyte phases.
%>             It then assembles the overall system Jacobian and nonlinear residual to solve for the 
%>             updated state vector, updates the cell voltage, and stores the results.
%>
%>             The function calls several nested functions to compute:
%>             - Time derivative of solid concentration (cs_time_derivative)
%>             - Diffusion in the solid phase (cs_diffusion)
%>             - Time derivative of electrolyte concentration (ce_time_derivative)
%>             - Diffusion in the electrolyte (ce_diffusion)
%>             - Charge conservation for solid phase potential (phis_charge_conservation)
%>             - Charge conservation for electrolyte potential (phie_charge_conservation)
%>             - Nonlinear terms for electrolyte potential (phie_nonlinear)
%>             - Self contribution of intercalation flux (jint_self)
%>             - Intercalation flux contribution via Butler-Volmer kinetics (jint_intercalation)
%>
%> @authors    Sanghyun Kim (shnkim@yonsei.ac.kr), Jung-Il Choi (jic@yonsei.ac.kr)
%> @date       March 2025
%> @version    1.0
%> @license    MIT License
%==================================================================================================

function mysim = sys1_update_variable(mysim)
    % Retrieve simulation options, finite difference operators, parameter functions, and grid info
    SD = mysim.opt.SolidDiffusion;
    FDM = mysim.FDM;
    pfunc = mysim.pfunc;
    vind = mysim.variable.vind;
    vdx = mysim.grid.vdx;
    dy = mysim.grid.dy;
    param = mysim.param;
    sys1 = mysim.sys1;
    sys2 = mysim.sys2;
    
    % Update state history: store previous state for time interpolation
    sys1.xoo = sys1.xo;
    sys1.xo = sys1.x;

    % Extract previous state variables at time t0 (x0)
    x0 = sys1.xo;
    cs = sys1.cs(x0);
    ce = sys1.ce(x0);
    phis = sys1.phis(x0);
    phie = sys1.phie(x0);
    j = sys1.jint(x0);
    T = sys2.T;
    if SD==1

    else
        cse = cs;
    end
    
    Iapp = mysim.Iapp;
    dt = mysim.dt;
    
    % Solid phase concentration time derivative and diffusion
    Mdx_1_1 = cs_time_derivative(dt, sys1);
    [Mdx_1_2, Mx_1, C_1] = cs_diffusion(vind, dt, sys1, sys2, param, SD);
    Mdx_1 = Mdx_1_1 + Mdx_1_2;
    
    % Electrolyte concentration time derivative and diffusion
    Mdx_2 = ce_time_derivative(dt, sys1);
    [Mx_2, bc_2] = ce_diffusion(vind, FDM, sys1, param);
    
    % Solid phase potential (charge conservation)
    [Mx1_3, bc_3] = phis_charge_conservation(vind, vdx, FDM, sys1, param, Iapp, dy);
    
    % Electrolyte potential (charge conservation)
    [Mx1_4, bc_4] = phie_charge_conservation(vind, vdx, FDM, sys1, param);
    [Jnl_4, nl_4] = phie_nonlinear(vind, FDM, sys1, param, T, ce);
    
    % Intercalation flux: Butler-Volmer kinetics
    Mx1_5 = jint_self(sys1);
    [Jnl_5, nl_5] = jint_intercalation(vind, sys1, param, pfunc, cse, ce, phis, phie, T, j);
    
    % Assemble contributions to the system Jacobian and nonlinear residual
    Mdx_Jacobian = Jnl_4 + Jnl_5;
    Mdx = Mdx_1 + Mdx_2 + Mdx_Jacobian;
    Mx  = Mx_1 + Mx_2;
    Mx1 = Mx1_3 + Mx1_4 + Mx1_5;
    nl = nl_4 + nl_5;
    bc = bc_2 + bc_3 + bc_4;
    C = C_1;

    % Form the total system: L1*x1 + L0*x0 + nl + bc + C = 0, solve for updated state x1
    L1 = Mdx + Mx1 + Mx / 2;
    L0 = -Mdx + Mx / 2;
    x1 = -L1 \ (L0 * x0 + nl + bc + C);
    
    % Update system state and simulation history
    sys1.x = x1;
    mysim.sys1 = sys1;
    
    % Update cell voltage and simulation time for the next time step
    nt = mysim.nt + 1;
    mysim.nt = nt;
    mysim.time(nt) = mysim.time(nt-1) + dt;
    mysim.Vcell(nt) = update_Vcell(vind, sys1.phis(x1), T, Iapp);

    % Store updated variables and diagnostic outputs
    mysim.result.sys1{nt} = x1;
    mysim.result.sys1T{nt} = T;
end

% Time derivative for solid concentration (cs)
function Mdx = cs_time_derivative(dt, sys1)
    row = 'cs';
    indr = sys1.xind.(row);
    col = 'cs';
    indc = sys1.xind.(col);

    N = length(sys1.x);
    A0 = sys1.speye(row, col);
    Mdx = sparse(N, N);
    Mdx(indr, indc) = A0 / dt;
end

% Diffusion in the solid phase (cs) and intercalation flux coupling (jint)
function [Mdx, Mx, C] = cs_diffusion(vind, dt, sys1, sys2, param, SolidDiffusion)
    row = 'cs';
    indr = sys1.xind.(row);
    col = 'jint';
    indc = sys1.xind.(col);

    ind = vind.(row);
    Rp = param.Rp(ind);
    Ds = param.Ds(ind);

    N = length(sys1.x);
    A0 = sys1.speye(row, col);

    Mdx = sparse(N, N);
    Mx = sparse(N, N);
    C = zeros(N, 1);

    if SolidDiffusion == 2
        Mdx0 = A0 .* Rp / 5 / dt ./ Ds;
        Mx0 = A0 * 3 ./ Rp;
    elseif SolidDiffusion == 3
        q = sys2.q(ind);
        Mdx0 = A0 .* Rp / 35 / dt ./ Ds;
        Mx0 = A0 * 57 / 7 ./ Rp;
        C(indr) = 48 / 7 * Ds ./ Rp .* q;
    end

    Mdx(indr, indc) = Mdx0;
    Mx(indr, indc) = Mx0;
end

% Time derivative for electrolyte concentration (ce)
function Mdx = ce_time_derivative(dt, sys1)
    row = 'ce';
    indr = sys1.xind.(row);
    col = 'ce';
    indc = sys1.xind.(col);
    
    N = length(sys1.x);
    A0 = sys1.speye(row, col);
    Mdx = sparse(N, N);
    Mdx(indr, indc) = A0 / dt;
end

% Diffusion for electrolyte concentration (ce)
function [Mx,bc] = ce_diffusion(vind,FDM,sys1,param)
    row = 'ce';
    indr = sys1.xind.(row);
    L = @(C) FDM.Diffusion(C, row).all;
    ind = vind.(row);
    
    De = param.De;
    tplus = param.tplus(ind);
    as = param.as(ind);
    epse = param.epse(ind);
    
    N = length(sys1.x);
    Mx = sparse(N, N);
    bc = zeros(N, 1);
    
    col = 'ce';
    indc = sys1.xind.(col);
    Mx(indr, indc) = -L(De) ./ epse;
    
    col = 'jint';
    indc = sys1.xind.(col);
    A0 = sys1.speye(row, col);
    Mx(indr, indc) = A0 .* (-as) .* (1 - tplus) ./ epse;
end

% Charge conservation for solid phase potential (phis)
function [Mx1, bc] = phis_charge_conservation(vind, vdx, FDM, sys1, param, Iapp, dy)
    row = 'phis';
    indr = sys1.xind.(row);
    L = @(C) FDM.Diffusion(C, row).all;
    ind = vind.(row);
    
    sigma = param.sigma;
    F = param.F;
    as = param.as(ind);
    N = length(sys1.x);
    Mx1 = sparse(N, N);
    bc = zeros(N, 1);
    
    col = 'phis';
    indc = sys1.xind.(col);
    Mx1(indr, indc) = L(sigma);
    
    indLB = indr(ismember(ind, vind.itabn));
    indRB = indr(ismember(ind, vind.itabp));
    Itabn = -Iapp.I / Iapp.Atabn;
    Itabp = Iapp.I / Iapp.Atabp;
    bc(indLB) = Itabn / Iapp.nneg;
    bc(indRB) = Itabp / Iapp.npos;
    
    col = 'jint';
    indc = sys1.xind.(col);
    A0 = sys1.speye(row, col);
    Mx1(indr, indc) = A0 .* (-as) * F;
end

% Charge conservation for electrolyte potential (phie)
function [Mx1, bc] = phie_charge_conservation(vind, vdx, FDM, sys1, param)
    row = 'phie';
    indr = sys1.xind.(row);
    L = @(C) FDM.Diffusion(C, row).all;
    ind = vind.(row);
    
    kappa = param.kappa;
    F = param.F;
    as = param.as(ind);
    N = length(sys1.x);
    Mx1 = sparse(N, N);
    bc = zeros(N, 1);
    
    col = 'phie';
    indc = sys1.xind.(col);
    Mx1(indr, indc) = L(kappa);
    
    dx = vdx.(row);
    nx = length(dx);
    indLB = (indr(1):nx:indr(end))';
    kappa = kappa(ind);
    indLB0 = (1:nx:length(ind))';
    Mx1(indLB, indLB) = Mx1(indLB, indLB) - 2 * diag(kappa(indLB0)) / dx(1)^2;
    
    col = 'jint';
    indc = sys1.xind.(col);
    A0 = sys1.speye(row, col);
    Mx1(indr, indc) = A0 .* as * F;
end

% Nonlinear contribution for electrolyte potential (phie)
function [Jnl, nl] = phie_nonlinear(vind, FDM, sys1, param, T, ce)
    row = 'phie';
    indr = sys1.xind.(row);
    col = 'ce';
    L = @(C) FDM.Diffusion(C, col).all;
    ind = vind.(row);
    
    F = param.F;
    R = param.R;
    tplus = param.tplus;
    nu = 1 - tplus;             % nu=(1-tplus)*(1+dlnfdlnce)
    kappa = param.kappa;
    kappaD = 2 * R / F .* nu .* T .* kappa;
    L0 = L(kappaD);
    
    N = length(sys1.x);
    Jnl = sparse(N, N);
    nl = zeros(N, 1);
    
    col = 'ce';
    ceinv = 1 ./ ce(ind);
    indc = sys1.xind.(col);
    A0 = sys1.speye(row, col);
    Jnl(indr, indc) = -A0 .* (L0 * ceinv);
    
    lnce = log(ce(ind));
    lnce(isinf(lnce)) = 0;
    nl(indr) = -L0 * lnce;
end

% Self-contribution of intercalation flux (jint)
function Mx1 = jint_self(sys1)
    row = 'jint';
    indr = sys1.xind.(row);
    col = 'jint';
    indc = sys1.xind.(col);
    
    N = length(sys1.x);
    A0 = sys1.speye(row, col);
    Mx1 = sparse(N, N);
    Mx1(indr, indc) = A0;
end

% Intercalation flux contribution via Butler-Volmer kineticsfunction [Jnl, nl] = jint_intercalation(vind, sys1, param, pfunc, cse, ce, phis, phie, T, j)
function [Jnl, nl] = jint_intercalation(vind, sys1, param, pfunc, cse, ce, phis, phie, T, j)
    kint = param.kint;
    BV = ButlerVolmer(vind, param, pfunc, kint, cse, ce, phis, phie, T, j);
    row = 'jint';
    indr = sys1.xind.(row);
    
    N = length(sys1.x);
    Jnl = sparse(N, N);
    nl = zeros(N, 1);
    
    ind = vind.(row);
    nl0 = BV.jint;
    nl(indr) = -nl0(ind);
    
    col = 'cs';
    indc = sys1.xind.(col);
    A0 = sys1.speye(row, col);
    Jnl(indr, indc) = -A0 .* BV.djdcse(ind);
    
    col = 'ce';
    indc = sys1.xind.(col);
    A0 = sys1.speye(row, col);
    Jnl(indr, indc) = -A0 .* BV.djdce(ind);
    
    col = 'phis';
    indc = sys1.xind.(col);
    A0 = sys1.speye(row, col);
    Jnl(indr, indc) = -A0 .* BV.djdphis(ind);
    
    col = 'phie';
    indc = sys1.xind.(col);
    A0 = sys1.speye(row, col);
    Jnl(indr, indc) = -A0 .* BV.djdphie(ind);
end

% Update cell voltage
function [Vcell] = update_Vcell(vind, phis, T, Iapp)
    itabn = vind.itabn;
    itabp = vind.itabp;
    
    Ttabn = mean(T(itabn));
    Ttabp = mean(T(itabp));
    sigman = (-0.04889 * Ttabn.^3 + 54.65 * Ttabn.^2 - 218 * Ttabn + 3.52e6) * 100;
    sigmap = (-0.0325 * Ttabp.^3 + 37.07 * Ttabp.^2 - 15000 * Ttabp + 2.408e6) * 100;
    Lytabn = 30e-3;
    Lytabp = 30e-3;
    Itabn = -Iapp.I / Iapp.Atabn;
    Itabp = Iapp.I / Iapp.Atabp;
    Vtabn = mean(phis(itabn)) + Itabn / sigman * Lytabn;
    Vtabp = mean(phis(itabp)) + Itabp / sigmap * Lytabp;
    Vcell = Vtabp - Vtabn;
end
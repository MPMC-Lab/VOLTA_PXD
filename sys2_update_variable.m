%==================================================================================================
%> @file       sys2_update_variable.m
%> @brief      Updates system 2 state variables (temperature and ion flux) based on the current simulation state.
%>
%> @details    This function updates the temperature field in system 2 using a thermal model and,
%>             if enabled, updates the mean ion flux in the solid phase. It calls nested functions to
%>             compute the temperature evolution via conduction, ohmic, entropic, and reaction heat.
%>
%> @authors    Sanghyun Kim (shnkim@yonsei.ac.kr), Jung-Il Choi (jic@yonsei.ac.kr)
%> @date       March 2025
%> @version    1.0
%> @license    MIT License
%==================================================================================================

function mysim = sys2_update_variable(mysim)
    % Update temperature if the Thermal model is enabled
    if mysim.opt.modelInclude.Thermal == 1
        mysim.sys2.Too = mysim.sys2.To;  % Save previous temperature field
        mysim.sys2.To = mysim.sys2.T;     % Update old temperature field
        mysim.sys2.T = thermal_model(mysim);
    end

    % Update mean ion flux if SolidDiffusion model 3 is enabled
    if mysim.opt.SolidDiffusion == 3
        mysim.sys2.q = mean_ion_flux(mysim);
    end
end

%% Thermal model for temperature update
function T = thermal_model(mysim)
    FDM = mysim.FDM;
    vind = mysim.variable.vind;
    vdx = mysim.grid.vdx;
    dy = mysim.grid.dy;
    dz = mysim.grid.dz;
    param = mysim.param;
    sys1 = mysim.sys1;
    sys2 = mysim.sys2;
    
    % Retrieve current state variables from system 1
    x = sys1.x;
    ce = sys1.ce(x);
    phis = sys1.phis(x);
    phie = sys1.phie(x);
    jint = sys1.jint(x);
    T0 = sys2.To;
    
    Iapp = mysim.Iapp;
    dt = mysim.dt;
    
    % Compute time derivative matrix based on thermal capacity
    Mdx = T_time_derivative(vind, dt, param);
    
    % Compute conduction matrix and boundary contributions
    [Mx, bc] = T_conduction(vind, vdx, dy, dz, param, FDM);
    
    % Compute ohmic heating contribution
    [Q1, QT1] = ohmic_heat(vind, vdx, param, FDM, ce, phis, phie, Iapp);
    
    % Compute entropic heating contribution
    QT2 = entropic_heat(vind, param, jint);
    
    % Compute reaction heat contribution
    Q3 = reaction_heat(vind, param, jint);
    
    % Total heat contributions
    Q = Q1 + Q3;
    QTv = QT1 + QT2;
    
    QT = spdiags(QTv, 0, Mdx);
    L1 =  Mdx - (Mx + QT) / 2;
    L0 = -Mdx - (Mx + QT) / 2;
    
    % Solve for updated temperature
    T = -L1 \ (L0 * T0 - bc - Q);
    
    % Nested function: Time derivative matrix for temperature
    function Mdx = T_time_derivative(vind, dt, param)
        ind = vind.T;
        rho = param.rho(ind);
        Cp = param.Cp(ind);
        Mdx = speye(length(ind)) / dt .* rho .* Cp;
    end

    % Nested function: Temperature conduction matrix and boundary conditions
    function [Mx, bc] = T_conduction(vind, vdx, dy, dz, param, FDM)
        ind = vind.T;
        dx = vdx.T;
        nx = length(dx);
        ny = length(dy);
        nz = length(dz);

        indices = reshape(1:nx*ny*nz, nx, ny, nz);
        indB1 = reshape(indices(1, :, :), [], 1, 1);
        indB2 = reshape(indices(nx, :, :), [], 1, 1);
        indB3 = reshape(indices(:, 1, :), [], 1, 1);
        indB4 = reshape(indices(:, ny, :), [], 1, 1);
        indB5 = reshape(indices(:, :, 1), [], 1, 1);
        indB6 = reshape(indices(:, :, nz), [], 1, 1);
    
        Diffusion = @(C) FDM.Diffusion(C, 'T');
        Tamb = param.Tamb;
        h = param.h;
        lambdax = param.lambdax;
        lambdayz = param.lambdayz;
        Dlambdax = Diffusion(lambdax);
        Dlambdayz = Diffusion(lambdayz);
        Mx = Dlambdax.x + Dlambdayz.y + Dlambdayz.z;

        % Compute boundary convection terms for each face
        h1 = h ./ (1 + h/2 * dx(1) ./ lambdax(indB1));
        h2 = h ./ (1 + h/2 * dx(nx) ./ lambdax(indB2));
        h3 = h ./ (1 + h/2 * dy(1) ./ lambdayz(indB3));
        h4 = h ./ (1 + h/2 * dy(ny) ./ lambdayz(indB4));
        h5 = h ./ (1 + h/2 * dz(1) ./ lambdayz(indB5));
        h6 = h ./ (1 + h/2 * dz(nz) ./ lambdayz(indB6));
        
        bc = zeros(length(ind), 1);
        bc(indB1) = bc(indB1) + h1 / dx(1) * Tamb;
        bc(indB2) = bc(indB2) + h2 / dx(nx) * Tamb;
        bc(indB3) = bc(indB3) + h3 / dy(1) * Tamb;
        bc(indB4) = bc(indB4) + h4 / dy(ny) * Tamb;
        bc(indB5) = bc(indB5) + h5 / dz(1) * Tamb;
        bc(indB6) = bc(indB6) + h6 / dz(nz) * Tamb;
        
        % Adjust conduction matrix at boundaries
        Mx(indB1, indB1) = Mx(indB1, indB1) - diag(h1 / dx(1) .* ones(length(indB1), 1));
        Mx(indB2, indB2) = Mx(indB2, indB2) - diag(h2 / dx(nx) .* ones(length(indB2), 1));
        Mx(indB3, indB3) = Mx(indB3, indB3) - diag(h3 / dy(1) .* ones(length(indB3), 1));
        Mx(indB4, indB4) = Mx(indB4, indB4) - diag(h4 / dy(ny) .* ones(length(indB4), 1));
        Mx(indB5, indB5) = Mx(indB5, indB5) - diag(h5 / dz(1) .* ones(length(indB5), 1));
        Mx(indB6, indB6) = Mx(indB6, indB6) - diag(h6 / dz(nz) .* ones(length(indB6), 1));
    end
    
    % Nested function: Ohmic heat calculation
    function [Q,QT] = ohmic_heat(vind,vdx,param,FDM,ce,phis,phie,Iapp)
        hmean = @(C) FDM.harmmean(C, 'T');
        G = FDM.Grad.T;
    
        ind = vind.T;
        indLB = ind(ismember(ind, vind.itabn));
        indRB = ind(ismember(ind, vind.itabp));
    
        phis = phis(ind);
        phie = phie(ind);
        ce = ce(ind);
        lnce = log(ce);
        lnce(isinf(lnce)) = 0;
    
        R = param.R;
        F = param.F;
        tplus = param.tplus(ind);
        nu = 1 - tplus;             % nu=(1-tplus)*(1+dlnfdlnce)
        sigma = param.sigma(ind);
        kappa = param.kappa(ind);
        gamma = -2 * R .* nu / F;
        
        fsigma = hmean(sigma);
        fkappa = hmean(kappa);
        
        Q01 = zeros(vind.nmax, 1);
        Q02 = zeros(vind.nmax, 1);
        Q03 = zeros(vind.nmax, 1);
        for i = 1:6
            Q01(ind) = Q01(ind) + (fsigma{i} .* (G{i} * phis).^2) / 2;
            Q02(ind) = Q02(ind) + (fkappa{i} .* (G{i} * phie).^2) / 2;
            Q03(ind) = Q03(ind) + (fkappa{i} .* (G{i} * phie)) .* (G{i} * lnce) / 2;
        end
        Itabn = -Iapp.I / Iapp.Atabn;
        Itabp = Iapp.I / Iapp.Atabp;
        Q01(indLB) = Q01(indLB) + fsigma{1}(indLB) .* ((-Itabn ./ fsigma{1}(indLB)).^2) / 2;
        Q01(indRB) = Q01(indRB) + fsigma{2}(indRB) .* ((-Itabp ./ fsigma{2}(indRB)).^2) / 2;
        Q03(ind) = Q03(ind) .* gamma;
    
        Q = Q01 + Q02;
        QT = Q03;
    end
    
    % Nested function: Entropic heat calculation
    function QT = entropic_heat(vind, param, jint)
        ind = vind.T;
        jint = jint(ind);
        F = param.F;
        as = param.as(ind);
        dUdT = param.dUdT(ind);
        QT = F * as .* jint .* dUdT;
    end
    
    % Nested function: Reaction heat calculation
    function Q = reaction_heat(vind, param, jint)
        ind = vind.T;
        jint = jint(ind);
        F = param.F;
        as = param.as(ind);
        eta = param.eta(ind);
        Q = F * as .* jint .* eta;
    end
end

%% Compute mean ion flux using a finite difference update for the ion flux in particles
function q = mean_ion_flux(mysim)
    vind = mysim.variable.vind;
    param = mysim.param;
    sys1 = mysim.sys1;
    
    x = sys1.x;
    jint = sys1.jint(x);
    
    q0 = mysim.sys2.q;
    dt = mysim.dt;
    
    ind = vind.q;
    Ds = param.Ds(ind);
    Rpsq = (param.Rp(ind)).^2;
    
    c1 = 1 / dt;
    c2 = 30 * Ds ./ Rpsq;
    c3 = 45 / 2 ./ Rpsq;
    
    n = length(ind);
    L1 = spdiags(c1 + c2 / 2, 0, n, n);
    L0 = spdiags(-c1 + c2 / 2, 0, n, n);
    C = c3 .* jint(ind);
    
    q = zeros(vind.nmax, 1);
    q(ind) = -L1 \ (L0 * q0(ind) + C);
end
%==================================================================================================
%> @file       init_Simulation.m
%> @brief      Initializes simulation settings, grids, and operators for the electrochemical-thermal model.
%>
%> @details    This module initializes the simulation by setting up numerical operators, grid domains,
%>             variable indices, and initial parameters/variables. It also builds system components for
%>             simulation (system1 and system2) and computes the applied current density.
%>
%> @authors    Sanghyun Kim (shnkim@yonsei.ac.kr), Jung-Il Choi (jic@yonsei.ac.kr)
%> @date       March 2025
%> @version    1.0
%> @license    MIT License
%==================================================================================================

function mysim = init_Simulation(init, opt)
    % Create numerical operators instance
    FDM0 = Numerical_operators;

    mysim.init = init;
    mysim.opt = opt;
    
    % Copy variable names and domains from initialization
    mysim.variable.names = init.variable.names;
    mysim.variable.domains = variable_domain(mysim);
    
    % Create the computational grid based on the domain settings
    mysim.grid = domain_grid(mysim);
    
    % Determine valid indices for each variable based on the grid
    mysim.variable.vind = valid_index(mysim);
    
    % Set up finite difference operators for numerical computations
    mysim.FDM = Operators(mysim, FDM0);
    
    % Initialize parameters and variables
    mysim.param0 = initialize_params(mysim);
    mysim.param = mysim.param0;
    mysim.pfunc = parameter_functions(mysim.variable.vind, mysim.param0);
    mysim.var0 = initialize_variables(mysim);
    
    % Build system components for the simulation
    mysim.sys1 = build_system1(mysim);
    mysim.sys2 = build_system2(mysim);
    
    % Set applied current density and time step settings
    mysim.Iapp = mysim.opt.Iapp;
    mysim.Iapp = current_density(mysim);
    mysim.dt = mysim.opt.dt_init;
    mysim.dt0 = mysim.dt;
    mysim.controller = mysim.opt.adaptive_time.controller;
    
    % Preallocate arrays for time and cell voltage recording
    mysim.time = zeros(100000, 1);
    mysim.Vcell = zeros(100000, 1);
    mysim.nt = 1;
end

% Compute applied current density based on grid dimensions and electrode areas.
function iapp = current_density(mysim)
    I = mysim.opt.Iapp;
    
    dx = mysim.grid.dx;
    dy = mysim.grid.dy;
    dz = mysim.grid.dz;

    if length(dz) > 1
        Atabn = mysim.init.param.Lz(2) * mysim.init.param.Lx(1);
        Atabp = mysim.init.param.Lz(4) * mysim.init.param.Lx(5); % single-sided tab
        % Atabp = mysim.init.param.Lz(2)*mysim.init.param.Lx(5); % double-sided tab
        normaln = dy(end);
        normalp = dy(end);
    else
        if length(dy) > 1
            Atabn = sum(mysim.init.param.Lz) * mysim.init.param.Lx(1);
            Atabp = sum(mysim.init.param.Lz) * mysim.init.param.Lx(5);
            normaln = dy(end);
            normalp = dy(end);
        else
            Atabn = sum(mysim.init.param.Lz) * mysim.init.param.Ly;
            Atabp = sum(mysim.init.param.Lz) * mysim.init.param.Ly;
            normaln = dx(1);
            normalp = dx(end);
        end
    end

    iapp.I = I;
    iapp.Atabn = Atabn;
    iapp.Atabp = Atabp;
    iapp.nneg = normaln;
    iapp.npos = normalp;
end

% Set variable domains based on nonzero entries for each variable field.
function vdomain = variable_domain(mysim)
    vnames = mysim.init.variable.names;
    vdomains = mysim.init.variable.domains;
    for iv = 1:numel(vnames)
        d0 = vdomains{iv};
        ddiff = diff(d0);
        idend = [find(ddiff > 1); numel(d0)];
        idstart = [1; idend(1:end-1) + 1];
        d1 = cell(numel(idend), 1);
        for i = 1:numel(idend)
            d1{i} = d0(idstart(i):idend(i));
        end
        vdomain.vdomain0.(vnames{iv}) = d0;
        vdomain.(vnames{iv}) = d1;
    end
end

% Create computational grid and determine index positions.
function grid = domain_grid(mysim)
    Lx = mysim.init.param.Lx;
    Ly = mysim.init.param.Ly;
    Lz = mysim.init.param.Lz;
    nx = mysim.opt.nx;
    ny = mysim.opt.ny;
    nz0 = mysim.opt.nz0;
    nz = mysim.opt.nz;
    Lx1 = cumsum(Lx);
    Lx0 = [0; Lx1(1:4)];

    % Set grid ticks based on cumulative dimensions
    x_tick = [0; cumsum(Lx)];
    z_tick = [0; cumsum(Lz)];
    x = zeros(sum(nx) + 1, 1);
    z = zeros(sum(nz0) + 1, 1);
    ii = 0; kk = 0;
    for i = 1:5
        x(ii+1:ii+nx(i)+1) = linspace(x_tick(i), x_tick(i+1), nx(i)+1);
        z(kk+1:kk+nz0(i)+1) = linspace(z_tick(i), z_tick(i+1), nz0(i)+1);
        ii = ii + nx(i);
        kk = kk + nz0(i);
    end
    y = linspace(0, Ly, ny+1)';

    dx = diff(x);
    dy = diff(y);
    dz = diff(z);

    % Specify indices for different regions
    xend = cumsum(nx);
    xstart = [1; xend(1:end-1) + 1];
    indx = cell(5, 1);
    indxyz = cell(5, 1);
    ind0 = reshape(1:sum(nx) * ny * nz, [sum(nx), ny, nz]);
    for i = 1:5
        indx{i} = (xstart(i):xend(i))';
        indxyz{i} = reshape(ind0(xstart(i):xend(i), :, :), [], 1, 1);
    end

    if nz > 1
        nz0csum = cumsum(nz0);
        iztabn = (nz0csum(1) + 1):nz0csum(2);
        iztabp = (nz0csum(3) + 1):nz0csum(4);   % single-sided tab
        % iztabp = (nz0csum(1)+1):nz0csum(2);   % double-sided tab
        itabn = reshape(ind0(xstart(1):xend(1), ny, iztabn), [], 1, 1);
        itabp = reshape(ind0(xstart(5):xend(5), ny, iztabp), [], 1, 1);     % single-sided tab
        % itabp = reshape(ind0(xstart(5):xend(5), ny, iztabp), [], 1, 1);   % double-sided tab
    else
        if ny > 1
            iztabn = 1;
            iztabp = 1;
            itabn = reshape(ind0(xstart(1):xend(1), ny, iztabn), [], 1, 1);
            itabp = reshape(ind0(xstart(5):xend(5), ny, iztabp), [], 1, 1); % single-sided tab
            % itabp =  reshape(ind0(xstart(5):xend(5),1,iztabp),[],1,1);    % double-sided tab
        else
            itabn = xstart(1);
            itabp = xend(5);
        end
    end

    grid.x = x;
    grid.y = y;
    grid.z = z;
    grid.dx = dx;
    grid.dy = dy;
    grid.dz = dz;
    grid.indx = indx;
    grid.indxyz = indxyz;
    grid.itabn = itabn;
    grid.itabp = itabp;

    % Compute grid spacing for each variable domain
    vnames = mysim.variable.names;
    vdomain = mysim.variable.domains;
    for i = 1:length(vnames)
        vn = vnames{i};
        grid.vdx.(vn) = dx(sort(cell2mat(indx(cell2mat(vdomain.(vn))))));
    end
end

% Determine valid index ranges for each variable.
function vind = valid_index(mysim)
    grid = mysim.grid;
    indx = grid.indx;
    indxyz = grid.indxyz;

    vnames = mysim.variable.names;
    vdomain = mysim.variable.domains;
    for i = 1:length(vnames)
        vn = vnames{i};
        vind.vindx.split.(vn) = cell(length(vdomain.(vn)), 1);
        vind.split.(vn) = cell(length(vdomain.(vn)), 1);
        for di = 1:length(vdomain.(vn))
            vind.vindx.split.(vn){di} = sort(cell2mat(indx(vdomain.(vn){di})));
            vind.split.(vn){di} = sort(cell2mat(indxyz(vdomain.(vn){di})));
        end
        vind.vindx.(vn) = sort(cell2mat(indx(cell2mat(vdomain.(vn)))));
        vind.(vn) = sort(cell2mat(indxyz(cell2mat(vdomain.(vn)))));
    end

    % Assign indices for different regions
    vind.ncc = indxyz{1};
    vind.ne = indxyz{2};
    vind.sep = indxyz{3};
    vind.pe = indxyz{4};
    vind.pcc = indxyz{5};
    vind.nmax = length(cell2mat(indxyz));
    vind.itabn = grid.itabn;
    vind.itabp = grid.itabp;
end

% Set up numerical operators for the simulation based on the finite difference method.
function FDM = Operators(mysim, FDM0)
    vnames = mysim.variable.names;
    vdomain = mysim.variable.domains;
    grid = mysim.grid;
    vind = mysim.variable.vind;

    FDM.Grad = gradient(vnames, vdomain, grid, vind, FDM0);
    FDM.harmmean = @(C, vname) harmmean(C, vname, vdomain, grid, vind, FDM0);
    FDM.Diffusion = @(C, vname) Diffusion_Neumann(C, vname, grid, FDM.Grad, FDM.harmmean);

    % Compute gradient at faces
    function Grad = gradient(vnames, vdomain, grid, vind, FDM0)
        dx0 = grid.dx;
        dy = grid.dy;
        dz = grid.dz;
        ny = length(dy);
        nz = length(dz);

        for iv = 1:length(vnames)
            vn = vnames{iv};
            ind0 = vind.split.(vn);
            indx0 = vind.vindx.split.(vn);

            Grad.(vn) = cell(6, 1);
            Grad0 = cell(6, length(vdomain.(vn)));
            ind = cell(length(vdomain.(vn)), 1);
            nx = zeros(length(vdomain.(vn)), 1);
            Bout = cell(6, length(vdomain.(vn)));
            for di = 1:length(vdomain.(vn))
                ind{di} = reshape(ind0{di}, [], ny, nz);
                dx = dx0(indx0{di});
                nx(di) = length(dx);
                Grad0(:, di) = FDM0.Gradient(dx, dy, dz);

                stepsize = reshape(repmat([1, nx(di), nx(di)*ny], 2, 1), [], 1);
                for i = 1:6
                    id = [-1; 0; 1] * stepsize(i);
                    Bout{i, di} = spdiags(Grad0{i, di}, id);
                end
            end
            
            ind0 = cat(1, ind{:});
            stepsize = reshape(repmat([1, sum(nx), sum(nx)*ny], 2, 1), [], 1);
            for i = 1:6
                id = [-1; 0; 1] * stepsize(i);
                Bin = zeros(sum(nx) * ny * nz, 3);
                for di = 1:length(vdomain.(vn))
                    Bin(ismember(ind0, ind{di}), :) = Bout{i, di};
                end
                Grad.(vn){i} = spdiags(Bin, id, sum(nx)*ny*nz, sum(nx)*ny*nz);
            end
        end
    end

    % Compute harmonic mean at faces
    function Cface = harmmean(C0, vname, vdomain, grid, vind, FDM0)
        vn = vname;
        dx0 = grid.dx;
        dy = grid.dy;
        dz = grid.dz;
        ind0 = vind.split.(vn);
        indx0 = vind.vindx.split.(vn);
        ny = length(dy);
        nz = length(dz);
    
        Cface = cell(6, 1);
        Cface0 = cell(6, length(vdomain.(vn)));
        for di = 1:length(vdomain.(vn))
            dx = dx0(indx0{di});
            C = C0(ind0{di});
            Cface0(:, di) = FDM0.harmmean(C, dx, dy, dz);
            for i = 1:6
                Cface0{i, di} = reshape(Cface0{i, di}, [], ny, nz);
            end
        end
        for i = 1:6
            Cface{i} = reshape(cat(1, Cface0{i, :}), [], 1, 1);
        end
    end

    % Diffusion operator with Neumann boundary conditions
    function L = Diffusion_Neumann(C0, vname, grid, grad, hmean)
        vn = vname;
        dx = grid.vdx.(vn);
        dy = grid.dy;
        dz = grid.dz;
        nx = length(dx);
        ny = length(dy);
        nz = length(dz);
        dx = reshape(repmat(reshape(dx, [], 1, 1), [1, ny, nz]), [], 1);
        dy = reshape(repmat(reshape(dy, 1, [], 1), [nx, 1, nz]), [], 1);
        dz = reshape(repmat(reshape(dz, 1, 1, []), [nx, ny, 1]), [], 1);
        
        Gface = grad.(vn);
        Cface = hmean(C0, vn);
        CG = cell(6, 1);
        for i = 1:6
            CG{i} = Cface{i} .* Gface{i};
        end
        L.x = (CG{2} - CG{1}) ./ dx;
        L.y = (CG{4} - CG{3}) ./ dy;
        L.z = (CG{6} - CG{5}) ./ dz;
        L.all = L.x + L.y + L.z;
    end
end

% Initialize parameters for simulation using the input parameters.
function param = initialize_params(mysim)
    nx = mysim.opt.nx;
    ny = mysim.opt.ny;
    nz = mysim.opt.nz;
    
    param0 = mysim.init.param;
    pnames = fieldnames(param0);
    pones = mat2cell(ones(sum(nx), ny, nz), nx, ny, nz);
    for ip = 1:length(pnames)
        pname = pnames{ip};
        pvals = param0.(pname);
        if length(pvals) > 1
            pmat = cell(size(pones));
            for i = 1:length(nx)
                pmat{i} = pvals(i) * pones{i};
            end
            param.(pname) = reshape(cell2mat(pmat), [], 1, 1);
        else
            param.(pname) = pvals;
        end
    end
end

% Initialize variables for simulation using the input variable structure.
function var = initialize_variables(mysim)
    nx = mysim.opt.nx;
    ny = mysim.opt.ny;
    nz = mysim.opt.nz;

    var0 = mysim.init.variable;
    vnames = var0.names;
    vones = mat2cell(ones(sum(nx), ny, nz), nx, ny, nz);
    for iv = 1:length(vnames)
        vname = vnames{iv};
        vvals = var0.(vname);
        vmat = cell(size(vones));
        for i = 1:length(nx)
            vmat{i} = vvals(i) * vones{i};
        end
        var.(vname) = reshape(cell2mat(vmat), [], 1, 1);
    end
end

% Build system1 (n+1 step) for the simulation.
function sys1 = build_system1(mysim)
    vdomain = mysim.variable.domains;
    vind = mysim.variable.vind;
    nx = mysim.opt.nx;
    ny = mysim.opt.ny;
    nz = mysim.opt.nz;

    vnames = mysim.opt.sys1;
    for iv = 1:length(vnames)
        vn = vnames{iv};
        sys1.(vn) = mysim.var0.(vn);
    end

    iend = 0;
    for iv = 1:length(vnames)
        vn = vnames{iv};
        vval = sys1.(vn);
        ind0 = vind.(vn);
        istart = iend + 1;
        iend = iend + length(ind0);
        x(istart:iend) = vval(ind0);
        xind.(vn) = (istart:iend)';
    end
    
    for iv = 1:length(vnames)
        vn = vnames{iv};
        sys1.(vn) = @(x) x2v(x, xind, vn, vind);
    end

    sys1.names = vnames;
    sys1.x = reshape(x, [], 1);
    sys1.xo = sys1.x;
    sys1.xoo = sys1.xo;
    sys1.xind = xind;
    sys1.speye = @(vr, vc) sys1_speye(vdomain, vr, vc, nx, ny, nz);

    function v = x2v(x, xind, vn, vind)
        v = zeros(vind.nmax, 1);
        ind = vind.(vn);
        v(ind) = x(xind.(vn));
    end

    function A = sys1_speye(vdomain, vr, vc, nx, ny, nz)
        dr = vdomain.vdomain0.(vr);
        dc = vdomain.vdomain0.(vc);

        nr = zeros(5, 1);
        nc = zeros(5, 1);
        A0 = cell(5, 1);
        for i = 1:length(dr)
            nr(dr(i)) = nx(dr(i));
        end
        for i = 1:length(dc)
            nc(dc(i)) = nx(dc(i));
        end
        
        isparse = find(~(nr .* nc));
        ispeye = find(nr .* nc);
    
        for i = 1:length(isparse)
            A0{isparse(i)} = sparse(nr(isparse(i)), nc(isparse(i)));
        end
        for i = 1:length(ispeye)
            A0{ispeye(i)} = speye(nr(ispeye(i)), nc(ispeye(i)));
        end
    
        A = kron(speye(ny * nz), blkdiag(A0{:}));
    end
end

% Build system2 (n+1/2 step) for the simulation.
function sys2 = build_system2(mysim)
    vnames = mysim.opt.sys2;
    for iv = 1:length(vnames)
        vn = vnames{iv};
        sys2.(vn) = mysim.var0.(vn);
    end

    sys2.names = vnames;
    sys2.To = sys2.T;
    sys2.Too = sys2.To;
end
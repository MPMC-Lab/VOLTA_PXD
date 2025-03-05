%==================================================================================================
%> @file       Numerical_operators.m
%> @brief      Provides numerical operators for finite difference methods.
%>
%> @details    This module defines functions to compute the harmonic mean at cell faces,
%>             the gradient matrices, and the diffusion operator with Neumann boundary
%>             conditions. These operators are used in the discretization of partial
%>             differential equations within the simulation.
%>
%> @authors    Sanghyun Kim (shnkim@yonsei.ac.kr), Jung-Il Choi (jic@yonsei.ac.kr)
%> @date       March 2025
%> @version    1.0
%> @license    MIT License
%==================================================================================================

function FDM = Numerical_operators
    % Assign function handles for numerical operations
    FDM.harmmean = @(C, dx, dy, dz) face_harmmean(C, dx, dy, dz);
    FDM.Gradient = @(dx, dy, dz) face_Gradient(dx, dy, dz);
    FDM.Diffusion = @(C, dx, dy, dz) Diffusion_Neumann(C, dx, dy, dz);
end

% Compute harmonic mean at cell faces for a given property C.
function Cface = face_harmmean(C, dx, dy, dz)
    nx = length(dx);
    ny = length(dy);
    nz = length(dz);
    dx = reshape(dx, [nx, 1, 1]);
    dy = reshape(dy, [1, ny, 1]);
    dz = reshape(dz, [1, 1, nz]);

    % Replicate dx, dy, dz to match cell dimensions
    ddx = repmat(dx, [1, ny, nz]);
    ddy = repmat(dy, [nx, 1, nz]);
    ddz = repmat(dz, [nx, ny, 1]);

    % Initialize Cface as a cell array for the six faces
    Cface = cell(6, 1);
    for i = 1:length(Cface)
        Cface{i} = C;
    end

    % Define indices for interior faces
    ind0 = reshape(1:nx*ny*nz, nx, ny, nz);
    ind1 = reshape(ind0(1:nx-1, :, :), [], 1);
    ind2 = reshape(ind0(2:nx, :, :), [], 1);
    ind3 = reshape(ind0(:, 1:ny-1, :), [], 1);
    ind4 = reshape(ind0(:, 2:ny, :), [], 1);
    ind5 = reshape(ind0(:, :, 1:nz-1), [], 1);
    ind6 = reshape(ind0(:, :, 2:nz), [], 1);

    % Compute harmonic mean for x-direction faces
    C12 = (ddx(ind1) + ddx(ind2)) ./ (ddx(ind2)./C(ind1) + ddx(ind1)./C(ind2));
    % Compute harmonic mean for y-direction faces
    C34 = (ddy(ind3) + ddy(ind4)) ./ (ddy(ind4)./C(ind3) + ddy(ind3)./C(ind4));
    % Compute harmonic mean for z-direction faces
    C56 = (ddz(ind5) + ddz(ind6)) ./ (ddz(ind6)./C(ind5) + ddz(ind5)./C(ind6));

    % Assign computed values to corresponding faces
    Cface{1}(ind2) = C12;
    Cface{2}(ind1) = C12;
    Cface{3}(ind4) = C34;
    Cface{4}(ind3) = C34;
    Cface{5}(ind6) = C56;
    Cface{6}(ind5) = C56;

    % Set boundary values to the original values of C
    if nx > 1
        indB1 = reshape(ind0(1, :, :), [], 1);
        indB2 = reshape(ind0(nx, :, :), [], 1);
        Cface{1}(indB1) = C(indB1);
        Cface{2}(indB2) = C(indB2);
    end
    if ny > 1
        indB3 = reshape(ind0(:, 1, :), [], 1);
        indB4 = reshape(ind0(:, ny, :), [], 1);
        Cface{3}(indB3) = C(indB3);
        Cface{4}(indB4) = C(indB4);
    end
    if nz > 1
        indB5 = reshape(ind0(:, :, 1), [], 1);
        indB6 = reshape(ind0(:, :, nz), [], 1);
        Cface{5}(indB5) = C(indB5);
        Cface{6}(indB6) = C(indB6);
    end
    % Replace any NaN values with zeros
    for i = 1:6
        indNaN = isnan(Cface{i});
        Cface{i}(indNaN) = 0;
    end
end

% Compute gradient matrices at cell faces based on grid spacings.
function Gface = face_Gradient(dx, dy, dz)
    nx = length(dx);
    ny = length(dy);
    nz = length(dz);
    dx = reshape(dx, [], 1);
    dy = reshape(dy, [], 1);
    dz = reshape(dz, [], 1);
    G0 = cell(6, 1);
    Gface = cell(6, 1);
    
    % Compute first-order gradient matrices in each spatial direction
    [G0{1}, G0{2}] = first_order_gradient(dx, nx);  % x-direction gradients
    [G0{3}, G0{4}] = first_order_gradient(dy, ny);  % y-direction gradients
    [G0{5}, G0{6}] = first_order_gradient(dz, nz);  % z-direction gradients
    
    % Define repetition factors for constructing the full gradient matrix
    nrepI = reshape(repmat([1, nx, nx*ny], 2, 1), [], 1);   % Inner repetition
    nrepO = reshape(repmat([ny*nz, nz, 1], 2, 1), [], 1);     % Outer repetition
    
    % Assemble gradient matrices for each face
    for i = 1:length(G0)
        [Bout, id] = spdiags(G0{i});
        id = id * nrepI(i);
        Bin = zeros(nx*ny*nz, length(id));
        for ii = 1:length(id)
            Bin(:, ii) = repmat(reshape(repmat(Bout(:, ii)', nrepI(i), 1), [], 1), nrepO(i), 1);
        end
        Gface{i} = spdiags(Bin, id, nx*ny*nz, nx*ny*nz);
    end

    % Nested function: first-order gradient using finite differences.
    function [Gl, Gr] = first_order_gradient(dx, nx)
        if nx > 1
            h = (dx(1:nx-1) + dx(2:nx)) / 2;
            m1 = -1 ./ h;
            d0 = -[0; m1];
            Gl = sparse(diag(m1, -1)) + sparse(diag(d0, 0));
            p1 = 1 ./ h;
            d0 = -[p1; 0];
            Gr = sparse(diag(p1, 1)) + sparse(diag(d0, 0));
        else
            Gl = 0;
            Gr = 0;
        end
    end

    % Nested function (not used): second-order gradient (provided for reference)
    function [Gl, Gr] = second_order_gradient(dx, nx)
        if nx > 1
            h1 = -dx(1:nx-2)/2; h2 = dx(2:nx-1)/2; h3 = dx(2:nx-1) + dx(3:nx)/2;
            denom = (h3-h2) .* (h2-h1) .* (h1-h3);
            p1 = [0; (h2.^2-h1.^2) ./ denom];
            m1 = [(h3.^2-h2.^2) ./ denom; -2/(dx(nx-1) + dx(nx))];
            d0 = -([p1; 0] + [0; m1]);
            Gl = sparse(diag(m1, -1)) + sparse(diag(p1, 1)) + sparse(diag(d0, 0));
            
            h1 = -dx(1:nx-2)/2 - dx(2:nx-1); h2 = -dx(2:nx-1)/2; h3 = dx(3:nx)/2;
            denom = (h3-h2) .* (h2-h1) .* (h1-h3);
            p1 = [2/(dx(1)+dx(2)); (h2.^2-h1.^2) ./ denom];
            m1 = [(h3.^2-h2.^2) ./ denom; 0];
            d0 = -([p1; 0] + [0; m1]);
            Gr = sparse(diag(m1, -1)) + sparse(diag(p1, 1)) + sparse(diag(d0, 0));
        else
            Gl = 0;
            Gr = 0;
        end
    end
end

% Compute diffusion operator with Neumann boundary conditions.
function L = Diffusion_Neumann(C, dx, dy, dz)
    % Compute harmonic means at faces and gradient matrices
    Cface = face_harmmean(C, dx, dy, dz);
    Gface = face_Gradient(dx, dy, dz);

    nx = length(dx);
    ny = length(dy);
    nz = length(dz);
    % Replicate grid spacings to match the number of cells
    dx = repmat(reshape(dx, [nx, 1, 1]), [1, ny, nz]);
    dy = repmat(reshape(dy, [1, ny, 1]), [nx, 1, nz]);
    dz = repmat(reshape(dz, [1, 1, nz]), [nx, ny, 1]);
    dx = dx(:);
    dy = dy(:);
    dz = dz(:);

    % Compute product of harmonic mean and gradient for each face
    CG = cell(6, 1);
    for i = 1:6
        CG{i} = Cface{i} .* Gface{i};
    end
    % Assemble diffusion contributions in each spatial direction
    L.x = (CG{2} - CG{1}) ./ dx;
    L.y = (CG{4} - CG{3}) ./ dy;
    L.z = (CG{6} - CG{5}) ./ dz;
    L.all = L.x + L.y + L.z;
end
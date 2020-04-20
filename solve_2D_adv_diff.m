% 2-D advection diffusion solution

% parameters required:
% N = number of CVs
% L = physical length of geometry (m)
% gamma = diffusion coefficient (kg/(m*s))
% rho = density of fluid (kg/m^3)
% u = velocity (m/s)
% phi_i = left-end boundary condition of quantity being solved (SI unit)
% phi_f = right-end boundary condition of quantity being solved (SI unit)
% solver_type = type of solver to be used for the calculations:
%   1) 'central_diff'
%   2) 'exponential'
%   3) 'power_law'

% return output:
% phi = distribution of the quantity over CVs

function [phi] = solve_2D_adv_diff(Nx, Ny, L, d, h, k, Cp, rho, u, v, phi_i, phi_f, solver_type, omega)
    % define phi for all nodes
    phi = zeros(Nx + 2, Ny + 2);
    
    % assign boundary conditions
    phi(1, :) = phi_i;
    phi(end, :) = phi_f;
    
    % calculate diffusion length array
    dx = L / Nx;
    dy = d / Ny;
    gamma = k / Cp;
    
    Dw = (dy * gamma / dx) * ones(Nx, Ny+2);
    Dw(1, :) = (dy * gamma /(dx / 2));
    De = (dy * gamma / dx) * ones(Nx, Ny+2);
    De(end, :) = (dy * gamma / (dx / 2));
    Ds = (dx * gamma / dy) * ones(Nx, Ny+2);
    Ds(:, 1) = 0;
    Ds(:, 2) = (dx * gamma / (dy/2));
    Ds(:, end) = (dx * gamma / (dy/2));
    Dn = (dx * gamma / dy) * ones(Nx, Ny+2);
    Dn(:, 1) = (dx * gamma / (dy/2));
    Dn(:, end-1) = (dx * gamma / (dy/2));
    Dn(:, end) = 0;
    
    % calculate flow strength
    Fw = rho * u * dy * ones(Nx, Ny+2);
    Fw(:, 1) = 0;
    Fw(:, end) = 0;
    Fe = rho * u * dy * ones(Nx, Ny+2);
    Fe(:, 1) = 0;
    Fe(:, end) = 0;
    Fs = rho * v * dx * ones(Nx, Ny+2);
    Fn = rho * v * dx * ones(Nx, Ny+2);
    
    % calculate Peclet number
    Pw = Fw.* Dw.^-1;
    Pe = Fe.* De.^-1;
    Ps = Fs.* Ds.^-1;
    Pn = Fn.* Dn.^-1;
    
    % calculate convective term
    q_conv = (h / Cp) * (2 * dx * dy / d);
    
    % compute coefficient arrays
    Aw = Dw.*peclet_function(abs(Pw), solver_type) + max(Fw, 0);
    Ae = De.*peclet_function(abs(Pe), solver_type) + max(-Fe, 0);
    As = Ds.*peclet_function(abs(Ps), solver_type) + max(Fs, 0);
    An = Dn.*peclet_function(abs(Pn), solver_type) + max(-Fn, 0);
    Ap = Aw + Ae + As + An + q_conv;
    
    % compute the value of quantity in steady state for each CV
    error = 1;
    while error > 10^-6
        phi_old = phi;

        % x-sweep
        for i = 2:Nx+1
            A = zeros(Ny+2, 3);
            A(2:end, 1) = -As(i-1, 2:end);
            A(:, 2) = Ap(i-1, :) / omega;
            A(1:end-1, 3) = -An(i-1, 1:end-1);
            b = Aw(i-1, :).* phi(i-1, :);
            b = b + Ae(i-1, :).* phi(i+1, :);
            b = b + Ap(i-1, :).* phi(i, :) * ((1-omega)/omega);
            b(1) = b(1) + As(i-1, 1) * phi(i, 1);
            b(end) = b(end) + An(i-1, end) * phi(i, end);

            phi(i, :) = TDMA_solver(Ny+2, A, b);
        end
        
        % y-sweep
        % updating bottom boundary layer
        A = zeros(Nx, 3);
        A(2:end, 1) = -Aw(2:end, 1);
        A(:, 2) = Ap(:, 1) / omega;
        A(1:end-1, 3) = -Ae(1:end-1, 1);
        b = An(:, 1).* phi(2:end-1, 2);
        b = b + Ap(:, 1).* phi(2:end-1, 1) * ((1-omega)/omega);
        b(1) = b(1) + Aw(1, 1) * phi(1, 1);
        b(end) = b(end) + Ae(end, 1) * phi(end, 1);
        
        phi(2:end-1, 1) = TDMA_solver(Nx, A, b);
        
        % sweeping through interior cells
        for j = 2:Ny+1
            A = zeros(Nx, 3);
            A(2:end, 1) = -Aw(2:end, j);
            A(:, 2) = Ap(:, j) / omega;
            A(1:end-1, 3) = -Ae(1:end-1, j);
            b = As(:, j).* phi(2:end-1, j-1);
            b = b + An(:, j).* phi(2:end-1, j+1);
            b = b + Ap(:, j).* phi(2:end-1, j) * ((1-omega)/omega);
            b(1) = b(1) + Aw(1, j) * phi(1, j);
            b(end) = b(end) + Ae(end, j) * phi(end, j);
            
            phi(2:end-1, j) = TDMA_solver(Nx, A, b);
        end
        
        % updating top boundary layer
        A = zeros(Nx, 3);
        A(2:end, 1) = -Aw(2:end, Ny+2);
        A(:, 2) = Ap(:, Ny+2) / omega;
        A(1:end-1, 3) = -Ae(1:end-1, Ny+2);
        b = As(:, Ny+2).* phi(2:end-1, Ny+1);
        b = b + Ap(:, Ny+2).* phi(2:end-1, Ny+2) * ((1-omega)/omega);
        b(1) = b(1) + Aw(1, Ny+2) * phi(1, Ny+2);
        b(end) = b(end) + Ae(end, Ny+2) * phi(end, Ny+2);
        
        phi(2:end-1, Ny+2) = TDMA_solver(Nx, A, b);
        
        % calculate error
        error = sum(sum(abs(phi - phi_old)));
    end
end
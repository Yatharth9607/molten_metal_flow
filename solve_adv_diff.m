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

function [phi] = solve_adv_diff(Nx, Ny, L, d, h, k, Cp, rho, u, v, phi_i, phi_f, solver_type, omega)
    % define phi for all nodes
    phi = zeros(Nx + 2, Ny + 2);
    
    % assign boundary conditions
    phi(1, :) = phi_i;
    phi(end, :) = phi_f;
    
    % calculate diffusion length array
    dx = L / Nx;
    dy = d / Ny;
    gamma = k / Cp;
    
    Dw = (dy * gamma / dx) * ones(Nx, Ny);
    Dw(1, :) = (dy * gamma /(dx / 2));
    De = (dy * gamma / dx) * ones(Nx, Ny);
    De(end, :) = (dy * gamma / (dx / 2));
    Ds = (dx * gamma / dy) * ones(Nx, Ny);
    Ds(:, 1) = (dx * gamma / (dy/2));
    Dn = (dx * gamma / dy) * ones(Nx, Ny);
    Dn(:, end) = (dx * gamma / (dy/2));
    
    % calculate flow strength
    Fw = rho * u * dy;
    Fe = rho * u * dy;
    Fs = rho * v * dx;
    Fn = rho * v * dx;
    
    % calculate Peclet number
    Pw = Fw * Dw.^-1;
    Pe = Fe * De.^-1;
    Ps = Fs * Ds.^-1;
    Pn = Fn * Dn.^-1;
    
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
        A = zeros(Ny+2, 3);
        b = zeros(Ny+2, 1);
        
        % equations for CVs adjacent to left inlet (i = 2)
        % equation for bottom boundary node
        A(1, 2) = ((k/dx) + (2*k/dy) + h) / omega;
        A(1, 3) = -(2*k/dy);
        b(1) = (k/dx) * phi(3, 1);
        b(1) = b(1) + ((k/dx) + (2*k/dy) + h) * phi(2, 1) * (1-omega) / omega;

        % equation for interior cells
        A(2:end-1, 1) = -As(1, :);
        A(2:end-1, 2) = Ap(1, :) / omega;
        A(2:end-1, 3) = -An(1, :);
        b(2:end-1) = Aw(1, :).* phi(1, 2:end-1) + Ae(1, :).* phi(3, 2:end-1) + Ap(1, :).* phi(2, 2:end-1) * ((1-omega)/omega);

        % equation for top boundary node
        A(end, 1) = -(2*k/dy);
        A(end, 2) = ((k/dx) + (2*k/dy) + h) / omega;
        b(end) = (k/dx) * phi(3, end);
        b(end) = b(end) + ((k/dx) + (2*k/dy) + h) * phi(2, end) * (1-omega) / omega;

        phi(2, :) = TDMA_solver(Ny+2, A, b);

        % equation for interior CVs (i = 3 to ITMAX-2)
        for i = 3:Nx
            % equation for bottom boundary node
            A(1, 2) = ((2*k/dx) + (2*k/dy) + h) / omega;
            A(1, 3) = -(2*k/dy);
            b(1) = (k/dx) * phi(i-1, 1) + (k/dx) * phi(i+1, 1);
            b(1) = b(1) + ((2*k/dx) + (2*k/dy) + h) * phi(i, 1) * (1-omega) / omega;
            
            % equation for interior cells
            A(2:end-1, 1) = -As(i-1, :);
            A(2:end-1, 2) = Ap(i-1, :) / omega;
            A(2:end-1, 3) = -An(i-1, :);
            b(2:end-1) = Aw(i-1, :).* phi(i-1, 2:end-1) + Ae(i-1, :).* phi(i+1, 2:end-1) + Ap(i-1, :).* phi(i, 2:end-1) * ((1-omega)/omega);
            
            % equation for top boundary node
            A(end, 1) = -(2*k/dy);
            A(end, 2) = ((2*k/dx) + (2*k/dy) + h) / omega;
            b(end) = (k/dx) * phi(i-1, end) + (k/dx) * phi(i+1, end);
            b(end) = b(end) + ((2*k/dx) + (2*k/dy) + h) * phi(i, end) * (1-omega) / omega;
            
            phi(i, :) = TDMA_solver(Ny+2, A, b);
        end

        % equations for CVs adjacent to left inlet (i = ITMAX - 1)
        % equation for bottom boundary node
        A(1, 2) = ((k/dx) + (2*k/dy) + h) / omega;
        A(1, 3) = -(2*k/dy);
        b(1) = (k/dx) * phi(end-2, 1);
        b(1) = b(1) + ((k/dx) + (2*k/dy) + h) * phi(end-1, 1) * (1-omega) / omega;

        % equation for interior cells
        A(2:end-1, 1) = -As(1, :);
        A(2:end-1, 2) = Ap(1, :) / omega;
        A(2:end-1, 3) = -An(1, :);
        b(2:end-1) = Aw(1, :).* phi(end-2, 2:end-1) + Ae(1, :).* phi(end, 2:end-1) + Ap(1, :).* phi(end-1, 2:end-1) * ((1-omega)/omega);

        % equation for top boundary node
        A(end, 1) = -(2*k/dy);
        A(end, 2) = ((k/dx) + (2*k/dy) + h) / omega;
        b(end) = (k/dx) * phi(end-2, end);
        b(end) = b(end) + ((k/dx) + (2*k/dy) + h) * phi(end-1, end) * (1-omega) / omega;

        phi(end-1, :) = TDMA_solver(Ny+2, A, b);

        % y-sweep
        A = zeros(Nx, 3);
        b = zeros(Nx, 1);

        % updating bottom boundary layer
        % equation for left most boundary node
        A(1, 2) = ((k/dx) + (2*k/dy) + h) / omega;
        A(1, 3) = -(k/dx);
        b(1) = (2*k/dy) * phi(2, 2);
        b(1) = b(1) + ((k/dx) + (2*k/dy) + h) * phi(2, 1) * (1-omega) / omega;
        
        % equation for interior boundary nodes
        A(2:end-1, 1) = -(k/dx);
        A(2:end-1, 2) = ((2*k/dx) + (2*k/dy) + h) / omega;
        A(2:end-1, 3) = -(k/dx);
        b(2:end-1) = (2*k/dy) * phi(3:end-2, 2);
        b(2:end-1) = b(2:end-1) + ((2*k/dx) + (2*k/dy) + h) * phi(3:end-2, 1) * (1-omega) / omega;
        
        % equation for right most boundary node
        A(end, 1) = -(k/dx);
        A(end, 2) = ((k/dx) + (2*k/dy) + h) / omega;
        b(end) = (2*k/dy) * phi(end-1, 2);
        b(end) = b(1) + ((k/dx) + (2*k/dy) + h) * phi(end-1, 1) * (1-omega) / omega;

        phi(2:end-1, 1) = TDMA_solver(Nx, A, b);
        
        % sweeping through interior cells
        for j = 2:Ny+1
            A(2:end, 1) = -Aw(2:end, j-1);
            A(:, 2) = Ap(:, j-1) / omega;
            A(1:end-1, 3) = -Ae(1:end-1, j-1);
            b = As(:, j-1).* phi(2:end-1, j-1);
            b = b + An(:, j-1).* phi(2:end-1, j+1);
            b = b + Ap(:, j-1).* phi(2:end-1, j) * ((1-omega)/omega);
            b(1) = b(1) + Aw(1, j-1) * phi(1, j);
            b(end) = b(end) + Ae(end, j-1) * phi(end, j);
            
            phi(2:end-1, j) = TDMA_solver(Nx, A, b);
        end
        
        % updating top boundary layer
        % equation for left most boundary node
        A(1, 2) = ((k/dx) + (2*k/dy) + h) / omega;
        A(1, 3) = -(k/dx);
        b(1) = (2*k/dy) * phi(2, end-1);
        b(1) = b(1) + ((k/dx) + (2*k/dy) + h) * phi(2, end) * (1-omega) / omega;
        
        % equation for interior boundary nodes
        A(2:end-1, 1) = -(k/dx);
        A(2:end-1, 2) = ((2*k/dx) + (2*k/dy) + h) / omega;
        A(2:end-1, 3) = -(k/dx);
        b(2:end-1) = (2*k/dy) * phi(3:end-2, end-1);
        b(2:end-1) = b(2:end-1) + ((2*k/dx) + (2*k/dy) + h) * phi(3:end-2, end) * (1-omega) / omega;
        
        % equation for right most boundary node
        A(end, 1) = -(k/dx);
        A(end, 2) = ((k/dx) + (2*k/dy) + h) / omega;
        b(end) = (2*k/dy) * phi(end-1, end-1);
        b(end) = b(1) + ((k/dx) + (2*k/dy) + h) * phi(end-1, end) * (1-omega) / omega;

        phi(2:end-1, 1) = TDMA_solver(Nx, A, b);
        
        % calculate error
        error = sum(sum(abs(phi - phi_old)));
    end
end
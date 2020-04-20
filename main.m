% homework 5

% given data
Nx = 100;       % number of control volumes in stream-wise direction
Ny = 41;        % number of control volumes in cross-stream direction
L = 2;          % length of the channel (m)
d = 0.05;       % height of the channel (m)
rho = 10000;    % density of liquid lead (kg/m^3)
Cp = 140;       % specific heat capacity of liquid lead (J/kg C)
k = 21;         % thermal conductivity of liquid lead (W/m C)
h = 2000;       % convective heat transfer coefficient (W/m^2 C)
T_i = 1000;     % temperature at entrance (Celcius)
T_f = 700;      % temperature at exit (Celcius)
T_amb = 27;     % ambient temperature (Celcius)
u = 0.2;        % velocity in positive x direction (m/s)
v = 0;          % velocity in positive y direction (m/s)

omega = 0.95;   % over-relaxation factor

% solution:
dx = L / Nx;
dy = d / Ny;
x_axis = linspace(-L/(2*Nx), L + (L/(2*Nx)), Nx+2);
x_axis(1) = 0;
x_axis(end) = L;
y_axis = linspace(-d/(2*Ny), d + (d/(2*Ny)), Ny+2);
y_axis(1) = 0;
y_axis(end) = d;

% Power law scheme
theta = solve_adv_diff(Nx, Ny, L, d, h, k, Cp, rho, u, v, T_i-T_amb, T_f-T_amb, "power_law", omega);
Temperature = theta + T_amb;
num = sum(rho*Cp*u.*Temperature(2:end-1, 2:end-1), 2);
denom = rho*Cp*u*Ny;
avg_temp = num / denom;

% plot temperature distribution
figure (1)
plot(x_axis, Temperature(:, (Ny+3)/2)')
hold on
plot(x_axis, Temperature(:, 1)')
hold on
plot(x_axis(2:end-1), avg_temp, '--')

lgd1 = legend('center line', 'wall', 'average along channel');
title(lgd1, 'Temperature at -')
xlabel('location [m]')
ylabel('Temperature [degree C]')
title('Temperature distribution with heat convection')
hold off

% plot Nusselt number
numerator = (Temperature(2:end-1, end) - Temperature(2:end-1, end-1)) / (dy/2);
denominator = (Temperature(2:end-1, end) - avg_temp) / (L);
nu = numerator./ denominator;
figure (2)
plot(x_axis(2:end-1), nu)
title('Nusselt number across the length of the channel')
xlabel('location [m]')
ylabel('Nusselt number')

% plot contour graph
figure (3)
contourf(x_axis, y_axis, Temperature', 50, 'edgecolor', 'none')
title('Temperature contour plot')
xlabel('Length of channel [m]')
ylabel('Height of channel [m]')
c = colorbar;
c.Label.String = 'Temperature [degree C]';

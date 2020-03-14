clc
clear all

N = 100;

A = 2;
B = 4;
CC = 12;

e = A / 100;
p = B / 10; 
t = CC / 100;
phi = linspace(0, pi, N/2+1);
count = 1;

alpha = 8 * (pi/180); % angle of attack in radians

x = 0.5*(1-cos(phi));

[T, ybar, dTdx, dybardx] = naca4(e, p, t, x);

y_up = ybar+T/2;
y_low = ybar-T/2;

x_plot = [flip(x), x(2:end)];
y_plot = [flip(y_low), y_up(2:end)];

% PLOTTING THE AIRFOIL
figure(1)
plot(x_plot, y_plot, 'b')
axis('equal')

% SOLVING FOR THE FOURIER SERIES COEFFICIENTS
A_0 = alpha - (1/pi) * trapz(phi, dybardx);
for R = 1:10
    A_n(R) = (2/pi) * trapz(phi, dybardx.*cos(R*phi));
end

% Calculation for C_ell 
m = 2*pi;
alpha_max = 20;
alpha_0 = (1/pi) * trapz(phi, dybardx.*(1-cos(phi)));
X = [0:alpha_max] .* (pi/180); % list of angle of attacks
Y_C_ell = m.*(X-alpha_0);
X_C_ell = [0:alpha_max];

% Calculation for C_mac
Y_C_mac = -pi/4*(A_n(1) - A_n(2));
X_C_mac = [0:alpha_max];
Y_C_mac = Y_C_mac * ones(1,length(X_C_mac));

% PLOTTING C_ell and C_mac as functions of angle of attack
figure(2)
plot(X_C_ell, Y_C_ell, 'b')
hold on
plot(X_C_mac, Y_C_mac, 'r')
xlabel('Angle of Attack (in degrees)')
ylabel('Coefficients')
legend('Coefficient of Lift', 'Pitching Moment Coefficient',...
    'Location', 'northwest')
grid on

% Calculation for C_p
V_inf = 10; % m/s (freestream velocity)
% compute fourier series coefficients
A_0 = alpha - (1/pi) * ...
    trapz(phi, dybardx); % angle of attack is 8 degrees (C_ell is about 1)

sum = 0;
for R = 1:length(A_n)
    sum = sum + A_n(R)*sin(R.*phi);
end

gamma = 2 * V_inf * ...
    (A_0.*(1+cos(phi))./sin(phi) + sum);


% compute u_c, u_t, v_c, and v_t
dTdx(1) = 0; 
q = V_inf .* dTdx;
for R = 1:size(x,2)
    integrand = q./(x(R) - x);
    integrand(R) = 0.0;
    u_t(R) = 1.0/(2*pi).*trapz(x, integrand);
end

for R = 1:size(x, 2)
    integrand = gamma./(x(R)-x);
    integrand(R) = 0.0;
    v_c(R) = -1.0/(2*pi).*trapz(x, integrand);
end
v_c = V_inf.*(dybardx - alpha);

v_t_up = q./ 2.0;
v_t_low = -q./ 2.0;
u_c_up = gamma./ 2.0;
u_c_low = -gamma./ 2.0;

% compute V and V_mod
V_up = sqrt((V_inf*cos(alpha) + u_t + u_c_up).^2 + ...
    (V_inf*sin(alpha) + v_t_up + v_c).^2);

V_low = sqrt((V_inf*cos(alpha) + u_t + u_c_low).^2 + ...
    (V_inf*sin(alpha) + v_t_low + v_c).^2);

riegel_up = 1./sqrt(1+(dybardx+0.5.*dTdx).^2);
riegel_low = 1./sqrt(1+(dybardx-0.5.*dTdx).^2);

V_mod_up = V_up.*riegel_up;
V_mod_low = V_low.*riegel_low;

C_p_up = 1 - (V_mod_up ./ V_inf).^2;
C_p_low = 1 - (V_mod_low ./ V_inf).^2;

figure(3)
plot(x, -C_p_up, 'b');
hold on
plot(x, -C_p_low, 'b');
grid on
ylim([-1, 1])
xlabel('x/c')
ylabel('-Cp')
% axis('equal')

%---------------FUNCTIONS---------------
function [T, ybar, dTdx, dybardx] = naca4(e, p, t, x)

    T = 10*t*(0.2969*sqrt(x) - 0.126*x - 0.3536*x.^2 + ...
        0.2843*x.^3 - 0.1015*x.^4);
    dTdx = 10*t*(0.2969*0.5./sqrt(x) - 0.126 - 0.3537*2*x + ...
        0.2843*3*x.^2 - 0.1015*4*x.^3);

    n = length(x);
    ybar = zeros(1, n);
    dybardx = zeros(1, n);
    
    for i = 1:n
        if x(i) <= p
            ybar(i) = e/p^2 * (2*p*x(i) - x(i)^2);
            dybardx(i) = e/p^2 * (2*p - 2*x(i));
        else
            ybar(i) = e/(1-p)^2 * (1 - 2*p + 2*p*x(i) - x(i)^2);
            dybardx(i) = e/(1-p)^2 * (2*p - 2*x(i));
        end
    end

end
clc
clear all

[x_cfd, y_cfd] = csvimport('pressure coefficient.csv', 'columns', {'X', 'Y'});
% N = 100;
% 
% A = 2;
% B = 4;
% CC = 12;
% 
% e = A / 100;
% p = B / 10; 
% t = CC / 100;
% phi = linspace(0, pi, N/2+1);
% count = 1;
% 
% alpha = 8 * (pi/180); % angle of attack in radians
% 
% x = 0.5*(1-cos(phi));
% 
% [T, ybar, dTdx, dybardx] = naca4(e, p, t, x);
% 
% y_up = ybar+T/2;
% y_low = ybar-T/2;
% 
% x_plot = [flip(x), x(2:end)];
% y_plot = [flip(y_low), y_up(2:end)];
% 
% % SOLVING FOR THE FOURIER SERIES COEFFICIENTS
% A_0 = alpha - (1/pi) * trapz(phi, dybardx);
% for R = 1:10
%     A_n(R) = (2/pi) * trapz(phi, dybardx.*cos(R*phi));
% end
% 
% % Calculation for C_ell 
% m = 2*pi;
% alpha_max = 20;
% alpha_0 = (1/pi) * trapz(phi, dybardx.*(1-cos(phi)));
% X = [0:alpha_max] .* (pi/180); % list of angle of attacks
% Y_C_ell = m.*(X-alpha_0);
% X_C_ell = [0:alpha_max];
% 
% % Calculation for C_mac
% Y_C_mac = -pi/4*(A_n(1) - A_n(2));
% X_C_mac = [0:alpha_max];
% Y_C_mac = Y_C_mac * ones(1,length(X_C_mac));
% 
% % compare taf and panel
% x_plot = [5, 8];
% y_panel_Cl = [0.8566, 1.2149];
% y_panel_Cm = [-0.0538, -0.0540];
% figure(1)
% plot(X_C_ell, Y_C_ell, 'b')
% hold on
% plot(X_C_mac, Y_C_mac, 'r')
% hold on
% plot(x_plot, y_panel_Cl, 'ob', 'MarkerSize', 5)
% hold on
% plot(x_plot, y_panel_Cm, 'or', 'MarkerSize', 5)
% xlabel('Angle of Attack (in degrees)')
% ylabel('Coefficients')
% legend('Cl (Thin Airfoil)', 'Cm (Thin Airfoil)',...
%     'Cl (Panel Method)', 'Cm (Panel Method)',...
%     'Location', 'northwest')
% grid on
% 
% % compare all three
% x_plot = [5, 8];
% y_cfd_Cl = [0.8195918, 1.16060531];
% y_cfd_Cm = [-0.0599837, -0.0636735];
% y_panel_Cl = [0.8566, 1.2149];
% y_panel_Cm = [-0.0538, -0.0540];
% figure(2)
% plot(X_C_ell, Y_C_ell, 'b')
% hold on
% plot(X_C_mac, Y_C_mac, 'r')
% hold on
% plot(x_plot, y_cfd_Cl, 'xb', 'MarkerSize', 10)
% hold on
% plot(x_plot, y_cfd_Cm, 'xr', 'MarkerSize', 10)
% hold on
% plot(x_plot, y_panel_Cl, 'ob', 'MarkerSize', 5)
% hold on
% plot(x_plot, y_panel_Cm, 'or', 'MarkerSize', 5)
% xlabel('Angle of Attack (in degrees)')
% ylabel('Coefficients')
% legend('Cl (Thin Airfoil)', 'Cm (Thin Airfoil)',...
%     'Cl (CFD)', 'Cm (CFD)', ...
%     'Cl (Panel Method)', 'Cm (Panel Method)', ...
%     'Location', 'northwest')
% grid on
% 
% 
% function [T, ybar, dTdx, dybardx] = naca4(e, p, t, x)
% 
%     T = 10*t*(0.2969*sqrt(x) - 0.126*x - 0.3536*x.^2 + ...
%         0.2843*x.^3 - 0.1015*x.^4);
%     dTdx = 10*t*(0.2969*0.5./sqrt(x) - 0.126 - 0.3537*2*x + ...
%         0.2843*3*x.^2 - 0.1015*4*x.^3);
% 
%     n = length(x);
%     ybar = zeros(1, n);
%     dybardx = zeros(1, n);
%     
%     for i = 1:n
%         if x(i) <= p
%             ybar(i) = e/p^2 * (2*p*x(i) - x(i)^2);
%             dybardx(i) = e/p^2 * (2*p - 2*x(i));
%         else
%             ybar(i) = e/(1-p)^2 * (1 - 2*p + 2*p*x(i) - x(i)^2);
%             dybardx(i) = e/(1-p)^2 * (2*p - 2*x(i));
%         end
%     end
% 
% end



N = 400;

A = 2;
B = 4;
CC = 12;
chord = 1;

e = A / 100;
p = B / 10; 
t = CC / 100;
phi = linspace(0, pi, N/2+1);
count = 1;
rho = 1.225; % kg/m^3

alpha = 8 * (pi/180); % angle of attack in radians

x = 0.5*(1-cos(phi));

[T, ybar, dTdx, dybardx] = naca4(e, p, t, x);

y_up = ybar+T/2;
y_low = ybar-T/2;

x = [flip(x), x(2:end)];
y = [flip(y_low), y_up(2:end)];


% control points
for R = 1:size(x,2)-1
    x_bar(R) = (x(R) + x(R+1))/2;
    y_bar(R) = (y(R) + y(R+1))/2;
end

% Flow tangency boundary condition
V_inf = 10;
for I = 1:N
    theta_i(I) = atan2((y(I+1) - y(I)) , ...
        (x(I+1) - x(I)));
    
    for J = 1:N
        r_1(I,J) = sqrt((x_bar(I)-x(J+1))^2 + ...
            (y_bar(I)-y(J+1))^2);
        
        theta_j(J) = atan2((y(J+1)-y(J)) , ...
            (x(J+1) - x(J)));
        
        r(I,J) = sqrt((x_bar(I)-x(J))^2 + ...
            (y_bar(I)-y(J))^2);        
        
        
        if I ~= J
            const_1 = x(J)-x_bar(I);
            const_2 = y(J+1)-y_bar(I);
            const_3 = y(J)-y_bar(I);
            const_4 = x(J+1)-x_bar(I);
            
            beta(I,J) = atan2(((const_1*const_2) - ...
                            (const_3*const_4)) , ...
                            ((const_1*const_4) + ...
                            (const_3*const_2)));             
        else
            beta(I,J) = pi;
        end       
    end
end

for I = 1:N
    A(I, N+1) = 0.0;
    for J = 1:N
        A(I,J) = log(r_1(I,J)/r(I,J))*sin(theta_i(I)-theta_j(J)) + ...
            beta(I,J)*cos(theta_i(I)-theta_j(J));
        A(I,N+1) = A(I,N+1) + ...
            log(r_1(I,J)/r(I,J))*cos(theta_i(I)-theta_j(J)) - ...
            beta(I,J)*sin(theta_i(I)-theta_j(J));
    end
    b(I) = 2*pi*V_inf*sin(theta_i(I) - alpha);        
end

k = 1; 
A_k_sum = 0.0;
A_N_sum = 0.0;
    
% Kutta condition
for J = 1:N        
    A_k = beta(k,J)*sin(theta_i(k)-theta_j(J)) - ...
        log(r_1(k,J)/r(k,J))*cos(theta_i(k)-theta_j(J));
    A_N = beta(N,J)*sin(theta_i(N)-theta_j(J)) - ...
        log(r_1(N,J)/r(N,J))*cos(theta_i(N)-theta_j(J));
    
    A_k_sum = A_k_sum + beta(k,J)*cos(theta_i(k)-theta_j(J)) + ...
        log(r_1(k,J)/r(k,J))*sin(theta_i(k)-theta_j(J));
    A_N_sum = A_N_sum + beta(N,J)*cos(theta_i(N)-theta_j(J)) + ...
        log(r_1(N,J)/r(N,J))*sin(theta_i(N)-theta_j(J));
    
    A(N+1,J) = A_k + A_N; 
end
A(N+1, N+1) = A_k_sum + A_N_sum;
b(N+1) = -2*pi*V_inf * ...
    (cos(theta_i(k)-alpha)+cos(theta_i(N)-alpha));

q_gamma = inv(A)*b'; % solve for source and vortex strengths

% Solve for V_ti
for I = 1:N
    vel_src = 0.0;
    vel_vtx = 0.0;
    vel_free = V_inf*cos(theta_i(I)-alpha);
    
    for J = 1:N      
        vel_src = vel_src + ...
            q_gamma(J)*(beta(I,J)*sin(theta_i(I)-theta_j(J)) - ...
            log(r_1(I,J)/r(I,J))*cos(theta_i(I)-theta_j(J)));
        vel_vtx = vel_vtx + ...
            q_gamma(N+1)*(beta(I,J)*cos(theta_i(I)-theta_j(J)) + ...
            log(r_1(I,J)/r(I,J))*sin(theta_i(I)-theta_j(J)));
    end
    
    V_ti(I) = vel_free + ...
        (1/(2*pi))*vel_src + ...
        (1/(2*pi))*vel_vtx;
end

% Calculate Cp
Cp = 1 - (V_ti/V_inf).^2;
figure(1)
plot(x_bar, -Cp, 'r')
hold on
xlabel('x/c')
ylabel('-Cp')
grid on

% Calculate Cm
P = 0.5*rho*V_inf^2.*Cp;
for R = 1:N
    F(R) = P(R)*sqrt((x(R+1)-x(R))^2+(y(R+1)-y(R))^2);
    x_dist(R) = 0.25 - x(R); % quarter chord
    y_dist(R) = 0.0 - y(R);
end
F_x = F.*sin(theta_i);
F_y = F.*cos(theta_i);
F = [F_x;F_y;zeros(1,length(F_x))]';
distance = [x_dist;y_dist;zeros(1,length(F_x))]';

moment = cross(F,distance,2);
moment = moment(:,3);
moment_sum = 0.0;
for R = 1:N
    moment_sum = moment_sum + moment(R);
end

Cm = moment_sum / (0.5*rho*V_inf^2*chord);

% Calculate Cl
gamma_sum = 0;
for R = 1:N
    dist = sqrt((x(R+1)-x(R))^2 + (y(R+1)-y(R))^2);
    gamma_sum = gamma_sum + (q_gamma(end)*dist);
end

lift = rho*V_inf*gamma_sum;
Cl = lift / (0.5*rho*V_inf^2*chord);

% Calculate Cd
Cd = 0.0; % because drag is 0



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

C_p = [-C_p_up,-C_p_low];
x = [x,x];
plot(x, C_p, 'b');
hold on
plot(x_cfd, -y_cfd, '.k');
legend ('Panel Methods', 'Thin Airfoil Theory', 'CFD')
xlim([0,1])




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


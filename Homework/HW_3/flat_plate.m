clc
clear all

% Parameters
V_e = 10; % m/s
mu = 0.00001802;  % dynamic viscosity
rho = 1.225;
nu = mu / rho;  % kinematic viscosity
R_e = 1 * 10^6;
L = R_e*nu/V_e;  % length of plate
% N = 20;  % number of iterations
% x = L;
dVedx = 0;  % for flat plate

% Thwaite's method
x_0 = 10^-6;
R_e_x_0 = rho*V_e*x_0 / mu;
theta_0 = 0.664*x_0 / sqrt(R_e_x_0);

[t_thw, theta_thw] = ode45(@(x, theta_thw) dthetadx(x, theta_thw, nu, V_e, dVedx), [x_0 L], theta_0);
gamma = theta_thw.^2 * dVedx ./ nu;

for R = 1:length(t_thw)
    if gamma(R) > 0.1
        gamma(R) = 0.1;
    end

    if gamma(R) >= 0 && gamma(R) <= 0.1
        ell_thw(R) = 0.22 + 1.57*gamma(R) - 1.8*gamma(R)^2;
        H_thw(R) = 2.61 - 3.75*gamma(R) - 5.24*gamma(R)^2;
        
    elseif gamma(R) > -0.1 && gamma(R) < 0
        ell_thw(R) = 0.22 + 1.402*gamma(R) + 0.018*gamma(R) ./ ...
            (0.107 + gamma(R));
        H_thw(R) = 0.2088 + 0.0731 / (0.14 + gamma(R));
    end
    
    delta_star_thw(R) = H_thw(R) * theta_thw(R);
    c_f_thw(R) = 2*nu*ell_thw(R) / (theta_thw(R)*V_e);
end

C_f_thw = 1/L * trapz(t_thw,c_f_thw);

% Blasius Solution for Flat Plate (laminar and incompressible)
R_e_bla = rho*V_e*t_thw ./ mu;

delta_star_bla = 1.72 .* t_thw ./ sqrt(R_e_bla);
theta_bla = 0.664 .* t_thw ./ sqrt(R_e_bla);
c_f_bla = 0.664 ./ sqrt(R_e_bla);
H_bla = delta_star_bla ./ theta_bla;

C_f_bla = 1/L * trapz(t_thw, c_f_bla);


% Head's Method
% initial conditions
R_e_x_0 = rho*V_e*x_0 ./ mu;
theta_head_0 = 0.036 * x_0 / R_e_x_0^0.2;
H1_0 = 10.6;

y0 = [theta_head_0, H1_0];

[t_head, y_head] = ode45(@(x, y_head) dydt(x, y_head, V_e, dVedx, nu), [x_0 L], y0);
theta_head = y_head(:, 1);
H1 = y_head(:, 2);
for R = 1:length(t_head)
    H_head(R) = H_solver(H1(R));
    c_f_head(R) = cf_solver(theta_head(R), H_head(R), V_e, nu);
    delta_star_head(R) = H_head(R) * theta_head(R);
end
C_f_head = 1/L * trapz(t_head,c_f_head);


% Schlichting Solutions
R_e_x_sch = rho*V_e.*t_head ./ mu;

delta_star_sch = 0.046 .* t_head ./ R_e_x_sch.^0.2;
theta_sch = 0.036 .* t_head ./ R_e_x_sch.^0.2;
c_f_sch = 0.0592 ./ R_e_x_sch.^0.2;
H_sch = delta_star_sch ./ theta_sch;

C_f_sch = 0.074 / (R_e_x_sch(length(t_head))^0.2);

%Plotting delta star
figure(1)
plot(t_thw,delta_star_bla,'b')
hold on
plot(t_thw,delta_star_thw,'r')
hold on
plot(t_head,delta_star_sch,'k')
hold on
plot(t_head,delta_star_head,'g')
legend('Blasius','Thwaite','Schlichting','Head', 'Location', 'northwest');
xlabel('x distance on flat plate'); 
ylabel('Displacement Thickness');
% title('Delta Star')
 

%---------------------------------------------------------------
%Plotting theta
figure(2)
plot(t_thw,theta_bla,'b')
hold on
plot(t_thw,theta_thw,'r')
hold on
plot(t_head,theta_sch,'k')
hold on
plot(t_head,theta_head,'g')
legend('Blasius','Thwaite','Schlichting','Head', 'Location', 'northwest');
xlabel('x distance on flat plate');
ylabel('Momentum Thickness');
% title('Theta')
 

%---------------------------------------------------------------
%Plotting H
figure(3)
plot(t_thw,H_bla,'b')
hold on
plot(t_thw,H_thw,'r')
hold on
plot(t_head,H_sch,'k')
hold on
plot(t_head,H_head,'g')
legend('Blasius','Thwaite','Schlichting','Head');
xlabel('x distance on flat plate'); 
ylabel('Shape Factor');
ylim([1 3.5])
% title('H')
 

%---------------------------------------------------------------
%Plotting cf
figure(4)
plot(t_thw,c_f_bla,'b')
hold on
plot(t_thw,c_f_thw,'r')
hold on
plot(t_head,c_f_sch,'k')
hold on
plot(t_head,c_f_head,'g')
legend('Blasius','Thwaite','Shlichting','Head');
xlabel('x distance on flat plate'); 
ylabel('Local Coefficient of Skin Friction Drag');
% title('cf')
ylim([0 .015])


function y = dthetadx(x, theta, nu, V_e, dVedx) % for Thwaite's Method
    % y is new theta
    y = 0.225*nu / (V_e * theta) - 3*theta*dVedx / V_e;  
end

function y = dydt(x, y_head, V_e, dVedx, nu) % for Head's Method
    % y_head(1) = theta, y_head(2) = H1
    y = zeros(2,1);
    H = H_solver(y_head(2));
    cf = cf_solver(y_head(1), H, V_e, nu);
    
    y(1) = cf/2 - dVedx*y_head(1)/V_e * (H + 2); % dthetadx
    temp = y_head(2) - 3;    
    y(2) = (0.0306*(temp^-0.6169))/y_head(1) - ...
        dVedx*y_head(2)/V_e - (y(1)*y_head(2))/y_head(1); % dH1dx
end


function y = H_solver(H1)
    if H1 < 3.3
        y = 3.0;        
    end
    temp = H1 - 3.3;
    if H1 >= 5.3
        y = 0.86*(temp^-0.777) + 1.1;
    elseif H1 < 5.3
        y = 1.1538*(temp^-0.326) + 0.6778;
    end
end

function y = cf_solver(theta, H, V_e, nu)
    R_e = V_e * theta / nu;
    y = 0.246 * 10^(-0.678*H) * R_e^-0.268;
end

function y = H1_solver(H)
    temp1 = H - 1.1;
    temp2 = H - 0.6778;
    if H <= 1.6
        y = 0.8234*(temp1^-1.287) + 3.3;
    elseif H > 1.6
        y = 1.5501*(temp2^-3.064) + 3.3;
    end
end
clc 
clear all

N = 200;

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
mu = 0.00001802;  % dynamic viscosity
nu = mu / rho;
R_e = 1 * 10^6;
V_inf = R_e*nu/chord;  % freestream velocity

alpha = 5 * (pi/180); % angle of attack in radians

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
% figure(2)
% plot(x_bar, -Cp, 'b')
% xlabel('x/c')
% ylabel('-Cp')

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
% inviscid setup
[xu, Veu, xl, Vel] = parseairfoil(V_ti, x_bar, y_bar);
xu(1) = 10^-6;
xl(1) = 10^-6;
% initial conditions for Thwaite's
[Ve_upper_0, dVeudx_0] = velocity(xu, Veu, xu(1));
[Ve_lower_0, dVeldx_0] = velocity(xl, Vel, xl(1));

thetau_0_thw = sqrt((0.075*nu) / dVeudx_0);
thetal_0_thw = sqrt((0.075*nu) / dVeldx_0);


[t_u_thw, thetau_thw] =  ode45(@(x, thetau_thw) dthetadx(x, thetau_thw, nu, xu, Veu), ...
                                xu, thetau_0_thw);
[t_l_thw, thetal_thw] =  ode45(@(x, thetal_thw) dthetadx(x, thetal_thw, nu, xl, Vel), ...
                                xl, thetal_0_thw);
                                                        
for R = 1:length(xu)
    [Ve_upper(R), dVeudx(R)] = velocity(xu, Veu, xu(R));
end

[ell_u_thw, H_u_thw, ds_u_thw, cf_u_thw] = calculate_thw(thetau_thw, Ve_upper, ...
                                            dVeudx, t_u_thw, nu);

for R = 1:length(xl)
    [Ve_lower(R), dVeldx(R)] = velocity(xl, Vel, xl(R));
end

[ell_l_thw, H_l_thw, ds_l_thw, cf_l_thw] = calculate_thw(thetal_thw, Ve_lower, ...
                                            dVeldx, t_l_thw, nu);

% Check for transition
for R = 1:length(H_u_thw)
    check_Re = -40.4557 + 64.8066*H_u_thw(R) - 26.7538*H_u_thw(R)^2 + ...
        3.3819*H_u_thw(R)^3;
    Re_x = Ve_upper(R) * xu(R) / nu;
    if H_u_thw(R) > 2.1 && log10(Re_x) > check_Re % && H_u_thw(R) < 2.8
        transu_idx = R;
        break;
    end
end

for R = 1:length(H_l_thw)
    check_Re = -40.4557 + 64.8066*H_l_thw(R) - 26.7538*H_l_thw(R)^2 + ...
        3.3819*H_l_thw(R)^3;
    Re_x = Ve_lower(R) * xl(R) / nu;
    if H_l_thw(R) > 2.1 && log10(Re_x) > check_Re && H_l_thw(R) < 2.8 
        transl_idx = R;
        break;
    end
end

% initial conditions for Head's
thetau_0_head = thetau_thw(transu_idx);
thetal_0_head = thetal_thw(transl_idx);
H_one_0 = 10.6; 

y0_upper = [thetau_0_head, H_one_0];
y0_lower = [thetal_0_head, H_one_0];


[t_u_head, y_u_head] = ode45(@(x, y_u_head) dydt(x, y_u_head, xu(transu_idx:end), ...
                        Veu(transu_idx:end), nu), xu(transu_idx:end), y0_upper);
                    
[t_l_head, y_l_head] = ode45(@(x, y_l_head) dydt(x, y_l_head, xl(transl_idx:end), ...
                        Vel(transl_idx:end), nu), xl(transl_idx:end), y0_lower);                    
                    
theta_u_head = y_u_head(:,1);
H1_u_head = y_u_head(:,2);
theta_l_head = y_l_head(:,1);
H1_l_head = y_l_head(:,2);

for R = 1:length(t_u_head)
    H_u_head(R) = H_solver(H1_u_head(R));
    cf_u_head(R) = cf_solver(theta_u_head(R), H_u_head(R), Veu(R), nu);
    ds_u_head(R) = H_u_head(R) * theta_u_head(R);
end

for R = 1:length(t_l_head)
    H_l_head(R) = H_solver(H1_l_head(R));
    cf_l_head(R) = cf_solver(theta_l_head(R), H_l_head(R), Vel(R), nu);
    ds_l_head(R) = H_l_head(R) * theta_l_head(R);
end

% separation
for R = 1:length(t_u_head)
    if H_u_head(R) >= 3
        sepu = R-1;
        break;
    else
        sepu = R;
    end
end
for R = 1:length(t_l_head)
    if H_l_head(R) >= 3
        sepl = R - 1;
        break;
    else
        sepl = R;
    end
end

% concatenate thwaite and head
% upper
H_u = [H_u_thw(1:transu_idx), H_u_head(1:sepu)];
theta_u = [thetau_thw(1:transu_idx)', theta_u_head(1:sepu)'];
delta_star_u = [ds_u_thw(1:transu_idx), ds_u_head(1:sepu)];
cf_u = [cf_u_thw(1:transu_idx), cf_u_head(1:sepu)];
t_u = [t_u_thw(1:transu_idx)', t_u_head(1:sepu)'];
% Cf_u = 1/chord * trapz(t_u, cf_u);

% lower
H_l = [H_l_thw(1:transl_idx), H_l_head(1:sepl-1)];
theta_l = [thetal_thw(1:transl_idx)', theta_l_head(1:sepl-1)'];
delta_star_l = [ds_l_thw(1:transl_idx), ds_l_head(1:sepl-1)];
cf_l = [cf_l_thw(1:transl_idx), cf_l_head(1:sepl-1)];
t_l = [t_l_thw(1:transl_idx)', t_l_head(1:sepl-1)'];
% Cf_l = 1/chord * trapz(t_l, cf_l);


Cd_u = (2*theta_u(end)/chord)*(Veu(end)/V_inf)^((H_u(end)+5)/2);
Cd_l = (2*theta_l(end)/chord)*(Vel(end)/V_inf)^((H_l(end)+5)/2);

Cd = Cd_u + Cd_l;

% plot boundary layer properties
figure(3)
plot(t_u, H_u, 'b', t_l, H_l, 'r', t_u(transu_idx), H_u(transu_idx), 'bx', ...
    t_l(transl_idx), H_l(transl_idx), 'rx', 'MarkerSize', 10)
xlabel('Distance from Stagnation Point')
ylabel('Shape factor (H)')
ylim([1.2 3.5])
legend('Upper Surface', 'Lower Surface', 'Upper Surface Transition', 'Lower Surface Transition')
figure(4)
plot(t_u, theta_u, 'b', t_l, theta_l, 'r', t_u(transu_idx), theta_u(transu_idx), 'bx', ...
    t_l(transl_idx), theta_l(transl_idx), 'rx', 'MarkerSize', 10)
xlabel('Distance from Stagnation Point')
ylabel('Momentum Thickness')
legend('Upper Surface', 'Lower Surface', 'Upper Surface Transition', 'Lower Surface Transition', 'Location', 'northwest')
figure(5)
plot(t_u, delta_star_u, 'b', t_l, delta_star_l, 'r', t_u(transu_idx), delta_star_u(transu_idx), 'bx', ...
    t_l(transl_idx), delta_star_l(transl_idx), 'rx', 'MarkerSize', 10)
xlabel('Distance from Stagnation Point')
ylabel('Displacement thickness')
legend('Upper Surface', 'Lower Surface', 'Upper Surface Transition', 'Lower Surface Transition', 'Location', 'northwest')
figure(6)
plot(t_u, cf_u, 'b', t_l, cf_l, 'r', t_u(transu_idx), cf_u(transu_idx), 'bx', ...
    t_l(transl_idx), cf_l(transl_idx), 'rx', 'MarkerSize', 10)
xlabel('Distance from Stagnation Point')
ylabel('Local Coefficient of Friction')
ylim([0 0.02])
legend('Upper Surface', 'Lower Surface', 'Upper Surface Transition', 'Lower Surface Transition')

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

%--------------Boundary Conditions-------------
function y = dthetadx(x, theta, nu, xv, Vev) % for Thwaite's Method
    % y is new theta
    [V_e, dVedx] = velocity(xv, Vev, x);
    y = 0.225*nu / (V_e * theta) - 3*theta*dVedx / V_e;
end

function y = dydt(x, y_head, xv, Vev, nu) % for Head's Method
    % y_head(1) = theta, y_head(2) = H1
    [V_e, dVedx] = velocity(xv, Vev, x);
    y = zeros(2,1);
    H = H_solver(y_head(2));
    cf = cf_solver(y_head(1), H, V_e, nu);
    
    y(1) = cf/2 - dVedx*y_head(1)/V_e * (H + 2); % dthetadx
    if y_head(2) < 3
        y(2) = 0;
    else
        temp = y_head(2) - 3;    
        y(2) = (0.0306*temp^-0.6169)/y_head(1) - ...
            (dVedx*y_head(2))/V_e - (y(1)*y_head(2))/y_head(1); % dH1dx
    end
end


function y = H_solver(H1)
    temp = H1 - 3.3;
    
    if H1 < 3.3
        y = 3.0;
    elseif H1 >= 5.3
        y = 0.86*temp^-0.777 + 1.1;
    elseif H1 < 5.3
        y = 1.1538*temp^-0.326 + 0.6778;
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
        y = 0.8234*(temp1)^-1.287 + 3.3;
    elseif H > 1.6
        y = 1.5501*(temp2)^-3.064 + 3.3;
    end
end

function [ell, H, delta_star, cf] = calculate_thw(theta, Ve, dVedx, t, nu)
    for R = 1:length(t)
        gamma = theta(R)^2 * dVedx(R) / nu;
        if gamma >= 0 && gamma <= 0.1
            ell(R) = 0.22 + 1.57*gamma - 1.8*gamma^2;
            H(R) = 2.61 - 3.75*gamma - 5.24*gamma^2;

        elseif gamma > -0.1 && gamma < 0
            ell(R) = 0.22 + 1.402*gamma + 0.018*gamma / ...
                (0.107 + gamma);
            H(R) = 2.088 + (0.0731 / (0.14 + gamma));
        elseif gamma > 0.1
            gamma = 0.1;
        else
            break; % laminar separation
        end
        delta_star(R) = H(R) * theta(R);
        cf(R) = 2*nu*ell(R) / (theta(R)*Ve(R));
    end
end

%-----------Ve and dVEdx---------------
function [xu, Veu, xl, Vel] = parseairfoil(Vt, x, y)
    % find stagnation point
    idx = find(Vt > 0, 1, 'first');
    
    % separate into upper and lower surfaces
    xu = x(idx:end);
    yu = y(idx:end);
    Veu = Vt(idx:end);

    xl = x(idx-1:-1:1);
    yl = y(idx-1:-1:1);
    Vel = -Vt(idx-1:-1:1);  %negative sign because of definition of Vt

    % compute distance along curved path around airfoil
    nu = length(xu);
    su = zeros(nu,1);
    for i = 1:nu-1
        su(i+1) = su(i) + sqrt((xu(i+1) - xu(i))^2 + (yu(i+1) - yu(i))^2);
    end

    nl = length(xl);
    sl = zeros(nl,1);
    for i = 1:nl-1
        sl(i+1) = sl(i) + sqrt((xl(i+1) - xl(i))^2 + (yl(i+1) - yl(i))^2);
    end

    % rename to x (boundary layer definition for x)
    xu = su';
    xl = sl';
end

function [Ve, dVedx] = velocity(xv, Vev, x)
    if x < xv(1)
        disp("provided x value is below xv range")
    elseif x > xv(end)
        disp("provided x value is above xv range")
    end

    % find index for linear interpolation
    idx = find(x >= xv, 1, 'last');
    if idx == length(xv)
        idx = idx - 1;
    end

    frac = (x - xv(idx)) / (xv(idx+1) - xv(idx));

    Ve = Vev(idx) +  frac*(Vev(idx+1) - Vev(idx));

    dVedx = (Vev(idx+1) - Vev(idx)) / (xv(idx+1) - xv(idx));
end

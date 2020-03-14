clc 
clear all

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


% PLOTTING THE AIRFOIL
figure(1)
plot(x, y, 'b')
hold on
plot(x_bar, y_bar, 'xk')
axis('equal')

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
figure(2)
plot(x_bar, -Cp, 'b')
xlabel('x/c')
ylabel('-Cp')

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
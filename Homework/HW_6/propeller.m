clc
clear all

%6.1-----------------------------------------------
xfoil = importdata("naca4412_FINAL.txt");
exp_data = importdata("exp_data.txt");
J_exp = exp_data(:, 1);
CT_exp = exp_data(:, 2);
CP_exp = exp_data(:, 3);
eta_exp = exp_data(:, 4);
alpha = xfoil(:, 1);
cl = xfoil(:, 2);
cd = xfoil(:, 3);
D_prop = 0.254; % meters
c_R = 0.128; % from UIUC data
AR = 1 / c_R;
R_prop = D_prop/2;
B = 2;  %number of blades
c_R_tot = [0.149;0.173;0.189;0.197;0.201;0.200;0.194;0.186;....
    0.174;0.160;0.145;0.128;0.112;0.096;0.081;0.061];
c = c_R_tot*R_prop; % chord length at any location of the propeller

% Viterna Method-----------------------------------------------
% Upper angles of attack
cd_max = 1.11 + 0.18 * AR;
B_1 = cd_max;
A_1 = B_1 / 2;
cl_stall = cl(end);
cd_stall = cd(end);
alpha_stall = alpha(end) * (pi/180);
A_2 = (cl_stall - cd_max * sin(alpha_stall) * cos(alpha_stall)) * ...
    sin(alpha_stall)/(cos(alpha_stall))^2;
B_2 = cd_stall - (cd_max*(sin(alpha_stall))^2 / cos(alpha_stall));
alpha_up = (16.8:0.1:180) .* (pi/180);
C_D_up = B_1 .* (sin(alpha_up)).^2 + B_2 .* cos(alpha_up);
C_L_up = A_1 .* sin(2 .* alpha_up) + A_2 .* ((cos(alpha_up)).^2 ./ sin(alpha_up));
 
% Lower angles of attack
cl_stall = cl(1);
cd_stall = cd(1);
alpha_stall = alpha(1) * (pi/180);
A_2 = (cl_stall - cd_max * sin(alpha_stall) * cos(alpha_stall)) * ...
    sin(alpha_stall)/(cos(alpha_stall))^2;
B_2 = cd_stall - (cd_max*(sin(alpha_stall))^2 / cos(alpha_stall));
alpha_low = (-180:0.1:-14.4) .* (pi/180);
C_D_low = B_1 .* (sin(alpha_low)).^2 + B_2 .* cos(alpha_low);
C_L_low = A_1 .* sin(2 .* alpha_low) + A_2 .* ((cos(alpha_low)).^2 ./ sin(alpha_low));
 
%Plot Cd and Cl
plotting1(alpha,cl,C_L_up,C_L_low,cd,C_D_up,C_D_low,alpha_up,alpha_low);
 
%6.2 BEM-----------------------------------------------

%setup
%Concantenate cl, cd, and alpha values (from xfoil and vaterni)----
cl_0 = horzcat(C_L_low, cl', C_L_up);
cd_0 = horzcat(C_D_low, cd', C_D_up);
alpha = horzcat(alpha_low, alpha' .* (pi/180), alpha_up);
N = 100; % number of points
J = linspace(0.08, 0.695, N);
omega = 1;
rho = 1;
r_R = [0.20;0.25;0.30;0.35;0.40;0.45;0.50;0.55;0.60;...
    0.65;0.70;0.75;0.80;0.85;0.90;0.95]; % position on propeller
r = r_R .* R_prop;
 
for R = 1:N             %advance ratios
    for K = 1:length(r)     %across blade
        V_inf = J(R) * omega * D_prop / (2 * pi); 
        sigma_prime = B*c(K)/(2*pi*r(K));
 
        %solving for twist (theta)
        [theta] = twist(K);
 
        %solving residual equation
        % x = phi variable
        fun = @(x) residuals(x,omega,sigma_prime,theta,r(K),V_inf,alpha,R_prop,cl_0,cd_0,c(K),B);
        x0 = [(1*10^-6) pi/2];    
        
        phi = fzero(fun,x0);
        
        %solve for values using known phi
        alpha_fin = theta - phi;
        [cn_fin, ct_fin] = corrections(cl_0, cd_0, alpha_fin,omega,R_prop,alpha,r(K),V_inf,phi,c(K));
        
        %Tip modificaitons
        R_hub = 0.15*R_prop;
        ftip = B/2*((R_prop-r(K))/(r(K)*abs(sin(phi))));
        Ftip = (2/pi)*(acos(exp(-ftip)));
        fhub = (B/2)*((r(K)-R_hub)/(R_hub*abs(sin(phi))));
        Fhub = (2/pi)*(acos(exp(-fhub)));
        F = Ftip*Fhub;
 
        %induction factors
        a_prime = sigma_prime*ct_fin / (4*F*sin(phi)*cos(phi)+sigma_prime*ct_fin);
        a = sigma_prime * cn_fin/(4*F*(sin(phi))^2-sigma_prime*cn_fin);
        
        W_sq = (V_inf*(1+a))^2 + (omega*r(K)*(1-a_prime))^2;
        
        T_prime(K) = B*cn_fin*0.5*rho*W_sq*c(K);
        Q_prime(K) = B*r(K)*ct_fin*0.5*rho*W_sq*c(K);
    end
   
    T = trapz(r,T_prime);
    Q = trapz(r,Q_prime); 
    P = Q*omega;
        
    n = omega/(2*pi);
    CT(R) = T/(rho*n^2*D_prop^4);
    CQ(R) = Q/(rho*n^2*D_prop^5);
    CP(R) = P/(rho*n^3*D_prop^5);
        
    eff(R) = J(R)*(CT(R)/CP(R));
end 

%Plot CT,CP,Eff vs J
plotting2(J,CT,CP,eff,J_exp,CT_exp,CP_exp,eta_exp)
 
%*******************Functions************************
function [cn, ct] = corrections(cl_0, cd_0, alpha_curr,omega,R_prop,alpha,r,V_inf,x,c)
    a = 343.0; % speed of sound in m/s
    % call interpolation
    alpha_max = 30 * pi /180;
    alpha_min = -30 * pi / 180;
    
    if alpha_curr > alpha_max
        alpha_curr = alpha_max;
    elseif alpha_curr < alpha_min
        alpha_curr = alpha_min;
    end
    
    cl_0 = interp1(alpha, cl_0, alpha_curr);
    cd_0 = interp1(alpha, cd_0, alpha_curr);
    
    %V_l or W
    V_l = sqrt((omega*r)^2+(V_inf)^2);

    % Prandtl-Glauert
    M = V_l / a; % local mach number
    cl = cl_0 ./ (1 - M^2);
    cd = cd_0 ./ (1 - M^2);
    
    cn_0 = cl * cos(x) + cd * sin(x);
    ct = cl * sin(x) - cd * cos(x);
    
    %Rotational correction (only for c_n)
    cn = cn_0 + 1.5*(c/R_prop)^2*(2*pi*alpha_curr-cl_0)*(omega*R_prop/V_l)^2;
end
 
function y = residuals(x,omega,sigma_prime,theta,r,V_inf,alpha,R_prop,cl_0,cd_0,c,B)
    %alpha
    alpha_curr = theta - x;

    %cn and ct using alpha
    [cn, ct] = corrections(cl_0, cd_0, alpha_curr,omega,R_prop,alpha,r,V_inf,x,c);

    %Tip modificaitons
    R_hub = 0.15*R_prop;
    ftip = (B/2)*((R_prop-r)/(r*abs(sin(x))));
    Ftip = (2/pi)*(acos(exp(-ftip)));
    fhub = (B/2)*((r-R_hub)/(R_hub*abs(sin(x))));
    Fhub = (2/pi)*acos(exp(-fhub));
    F = Ftip*Fhub;

    %induction factors
    a_prime = sigma_prime*ct / (4*F*sin(x)*cos(x)+sigma_prime*ct);
    a = sigma_prime * cn/(4*F*(sin(x))^2-sigma_prime*cn);

    %residual equation
    y = sin(x)/(1+a) - V_inf*cos(x)/(omega*r*(1-a_prime));
end
 
function [theta] = twist(J)
    %Twist Values
    twist_value = [37.19;33.54;29.25;25.64;22.54;20.27;...
        18.46;17.05;15.97;14.87;14.09;13.39;12.84;12.25;11.37;...
        10.19];
    
    theta = twist_value(J) * (pi/180); % radians
end
 
function plotting1(alpha,cl,C_L_up,C_L_low,cd,C_D_up,C_D_low,alpha_up,alpha_low)
    %Plotting
    figure(1)
    plot(alpha, cl, 'r--')
    hold on
    plot(alpha_up .* (180/pi), C_L_up, 'b--')
    hold on 
    plot(alpha_low .* (180/pi), C_L_low, 'b--')
    xlabel('Angle of Attack')
    ylabel('Coefficient of Lift')
    xlim([-30 30])
    ylim([-1.5 2.0])
    legend("XFOIL", "Viterna", "Location", "northwest")

    figure(2)
    plot(alpha, cd, 'r--')
    hold on
    plot(alpha_up .* (180/pi), C_D_up, 'b--')
    hold on 
    plot(alpha_low .* (180/pi), C_D_low, 'b--')
    xlabel('Angle of Attack')
    ylabel('Coefficient of Drag')
    xlim([-30 30])
    legend("XFOIL", "Viterna")
    set(gca, 'YScale', 'log')
    end

    function plotting2(J,CT,CP,eff,J_exp,CT_exp,CP_exp,eta_exp)
    figure(3)
    plot(J,CT,'r')
    hold on
    plot(J_exp,CT_exp,'b.','MarkerSize', 10)
    xlabel('Advance Ratio (J)')
    ylabel('Thrust Coefficient')
    xlim([0.0 0.8])
    ylim([0.0 0.1])
    legend("BEM", "Experimental Data")

    figure(4)
    plot(J,CP,'r')
    hold on
    plot(J_exp,CP_exp,'b.','MarkerSize', 10)
    xlabel('Advance Ratio (J)')
    ylabel('Power Coefficient')
    xlim([0.0 0.8])
    ylim([0.0 0.1])
    legend("BEM", "Experimental Data")

    figure(5)
    plot(J,eff,'r')
    hold on
    plot(J_exp,eta_exp,'b.','MarkerSize', 10)
    xlabel('Advance Ratio (J)')
    ylabel('Efficency')
    xlim([0.0 0.8])
    ylim([0.0 1.0])
    legend("BEM", "Experimental Data", "Location", "northwest")
end
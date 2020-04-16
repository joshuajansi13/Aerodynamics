clc
clear all

N = 50;
t_c = 0.05;
alpha = linspace(-2, 2, N);
alpha = alpha .* (pi / 180); % radians
M_1 = 2.0;
theta_up = atan(t_c) - alpha; % radians
theta_low = atan(t_c) + alpha;
gamma = 1.4;

% supersonic thin airfoil theory
c_n_taf = 4 * alpha ./ (sqrt(M_1^2 - 1));
c_a_taf = (4 / (sqrt(M_1^2 - 1))) * (t_c)^2;
c_l_taf = c_n_taf;
c_d_taf = c_n_taf .* alpha + c_a_taf;
c_mle_taf = -2 .* alpha ./ (sqrt(M_1^2 - 1));

% shock-expansion theory
for R = 1:length(alpha)
    % First shockwave
    if theta_up(R) < 0 % Prandtl-Meyer
        fun = @(x) f_fan(x, gamma, M_1, -theta_up(R));
        x0 = [1 5];
        M_2_up = fzero(fun, x0);
        P2_P1_up(R) = ((1+((gamma-1)/2)*M_1^2) / ...
            (1+((gamma-1)/2)*M_2_up^2))^(gamma/(gamma-1));
        
    else % Shock-Exp
        fun = @(x) f_shock(x, gamma, M_1, theta_up(R));
        x0 = pi/4; % initial guess
        beta_up = fzero(fun, x0); % radians

        % upper surface
        Mn_1_up = M_1 * sin(beta_up);
        Mn_2_up = sqrt((1+((gamma-1)*Mn_1_up^2)/2) ./ ...
            (gamma*Mn_1_up^2 - (gamma-1)/2));
        M_2_up = Mn_2_up / sin(beta_up-theta_up(R));
        P2_P1_up(R) = 1 + 2*gamma*(Mn_1_up^2 - 1) / (gamma+1);
    end
    
    fun = @(x) f_shock(x, gamma, M_1, theta_low(R));
    x0 = pi/4; % initial guess
    beta_low = fzero(fun, x0); % radians
    
    % lower surface
    Mn_1_low = M_1 * sin(beta_low);
    Mn_2_low = sqrt((1+((gamma-1)*Mn_1_low^2)/2) ./ ...
        (gamma*Mn_1_low.^2 - (gamma-1)/2));
    M_2_low = Mn_2_low / sin(beta_low-theta_low(R));
    P2_P1_low(R) = 1 + 2*gamma*(Mn_1_low^2 - 1) ./ (gamma+1);
    
    % Second wave (Prandtl-Meyer: fan expansion)
    theta = 2 * atan(t_c);
    
    % upper surface
    fun = @(x) f_fan(x, gamma, M_2_up, theta);
    x0 = [1 5];
    M_3_up = fzero(fun, x0);
    P3_P2_up(R) = ((1+((gamma-1)/2)*M_2_up^2) / ...
        (1+((gamma-1)/2)*M_3_up^2))^(gamma/(gamma-1));
    
    % lower surface
    fun = @(x) f_fan(x, gamma, M_2_low, theta);
    x0 = [1 5];
    M_3_low = fzero(fun, x0);
    P3_P2_low(R) = ((1+((gamma-1)/2)*M_2_low^2) / ...
        (1+((gamma-1)/2)*M_3_low^2))^(gamma/(gamma-1));
end

P3_P1_up = P3_P2_up .* P2_P1_up;
P3_P1_low = P3_P2_low .* P2_P1_low;

c_mle_shock = (1/(4*gamma*M_1^2)) .* (P2_P1_up - P2_P1_low + ...
    3 .* P3_P1_up - 3 .* P3_P1_low) + (t_c^2/(4*gamma*M_1^2)) .* ...
    (P2_P1_up - P2_P1_low - P3_P1_up + P3_P1_low);

c_n_shock = (1/(gamma*M_1^2)) .* (P2_P1_low + P3_P1_low - ...
    P2_P1_up - P3_P1_up);

c_a_shock = (t_c/(gamma*M_1^2)) .* (P2_P1_up - P3_P1_up + ...
    P2_P1_low - P3_P1_low);

c_l_shock = c_n_shock .* cos(alpha) - c_a_shock .* sin(alpha);
c_d_shock = c_n_shock .* sin(alpha) + c_a_shock .* cos(alpha);

% figure(1)
% plot(alpha .* (180/pi), c_mle_taf, alpha .* (180/pi), c_mle_shock)
% legend("Thin Airfoil Theory", "Shock-Expansion Theory")
% xlabel("Angle of Attack")
% ylabel("Moment Coefficient")
% figure(2)
% plot(alpha .* (180/pi), c_l_taf, alpha .* (180/pi), c_l_shock)
% legend("Thin Airfoil Theory", "Shock-Expansion Theory", "Location", "northwest")
% xlabel("Angle of Attack")
% ylabel("Lift Coefficient")
% figure(3)
% plot(alpha .* (180/pi), c_d_taf, alpha .* (180/pi), c_d_shock)
% legend("Thin Airfoil Theory", "Shock-Expansion Theory", "Location", "northwest")
% xlabel("Angle of Attack")
% ylabel("Drag Coefficient")

% line of best fit for supersonic TAF
p_lift = polyfit(alpha, c_l_shock, 1); % p(1) is the predicted slope of the curve
p_moment = polyfit(alpha, c_mle_shock, 1);
function y = f_fan(x, gamma, M, theta)
    nu_M = sqrt((gamma+1)/(gamma-1)) * ...
        atan(sqrt((gamma-1)*(M^2-1)/(gamma+1))) - ...
        atan(sqrt(M^2-1));
    
    y = sqrt((gamma+1)/(gamma-1)) * ...
        atan(sqrt((gamma-1)*(x^2-1)/(gamma+1))) - ...
        atan(sqrt(x^2-1)) - nu_M - theta;
end

function y = f_shock(x, gamma, M, theta)
    y = (2*cot(x)*(M^2*(sin(x))^2 - 1)) / ...
        (M^2 * (gamma+cos(2*x)) + 2) - ...
        tan(theta); % function
end


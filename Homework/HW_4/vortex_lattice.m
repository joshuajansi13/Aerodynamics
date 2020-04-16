clc
clear all

points = 20;
alpha = linspace(-3,15,points);

for X = 1:points
    N = 100; % number of panels 
    b = 180; % 1 % inches; wing span 
    alpha_curr = alpha(X); % make into an array of varying angles
    beta = 0;
    phi = zeros(1,N); % dihedral angle 
    theta = linspace(2,0,N); % [0, 0, 0]; %angle of twist
    sweep = 30; % 45
    max_span = 0.5 * b;
    lat_spacing = max_span / N; 
    center_spacing = lat_spacing / 2;
    rho = 1.225;
    c_cp = linspace(28.57, 11.43, N); % linspace((0.2*b),(0.2*b),N);
    c_bv = linspace(28.57, 11.43, N+1); % linspace((0.2*b), (0.2*b), N+1));
    c_avg = (28.57+11.43)/2; % (0.2 * b);
    S = b*c_avg;
    x_hat = [1; 0; 0];

    V_inf = [cosd(alpha_curr)*cosd(beta);
        -sind(beta);
        sind(alpha_curr)*cosd(beta)];

    n_hat = [sind(theta); ...
        -cosd(theta).*sind(phi); ...
        cosd(theta).*cosd(phi)]; % normal vector

    % control points and vortex points
    % s = starboard; p = port
    y_s(1) = 0;
    x(1) = 0;
    for R = 1:N
        if R == 1
            ybar_s(R) = center_spacing;
        else
            ybar_s(R) = ybar_s(R-1) + lat_spacing;
        end
        y_s(R+1) = y_s(R) + lat_spacing;
    end

    xbar = (tand(sweep)*ybar_s)+((3/4)*c_cp);
    x = (tand(sweep)*y_s)+((1/4)*c_bv);
    y_p = -y_s;

    r_s = [x; y_s; zeros(1,N+1)];
    r_p = [x; y_p; zeros(1,N+1)];
    rcp = [xbar; ybar_s; zeros(1,N)];

    for I = 1:N
        for J = 1:N
            r1 = rcp(:,I) - r_s(:,J);
            r2 = rcp(:,I) - r_s(:,J+1);
            const = 1 / (4*pi);
            temp1 = cross(r1, r2) / (norm(r1)*norm(r2) + dot(r1, r2));
            temp2 = (1/norm(r1)) + (1/norm(r2));
            temp3 = cross(r1,x_hat)/((norm(r1)-r1(1)))*(1/norm(r1));
            temp4 = cross(r2,x_hat)/((norm(r2)-r2(1)))*(1/norm(r2));
            V_hat = const*(temp1*temp2 + temp3 - temp4);
            AIC_s = dot(V_hat, n_hat(:,I));

            r1 = rcp(:,I) - r_p(:,J);
            r2 = rcp(:,I) - r_p(:,J+1);
            temp1 = cross(r1, r2) / (norm(r1)*norm(r2) + dot(r1, r2));
            temp2 = (1/norm(r1)) + (1/norm(r2));
            temp3 = cross(r1,x_hat)/((norm(r1)-r1(1)))*(1/norm(r1));
            temp4 = cross(r2,x_hat)/((norm(r2)-r2(1)))*(1/norm(r2));
            V_hat = const*(temp1*temp2 + temp3 - temp4);
            AIC_p = dot(V_hat, n_hat(:,I));

            AIC(I,J) = AIC_s - AIC_p;
        end
        B(I) = dot(-V_inf,n_hat(:,I));
    end

    gamma = AIC \ B';

    % Kutta-Joukowski Theorem
    L = 0;
    for R = 1:N
        L = L + gamma(R)  * (y_s(R+1) - y_s(R));
    end

    L = 2 * L * 1 * rho;
    cL(X) = L / (0.5 * rho * 1^2 * S);

    % % make elliptical lift distribution
    % L_ellip = 0;
    % for R = 1:N
    %     gamma_ellip(R) = sqrt(1 - ((2*ybar_s(R))/b)^2);
    %     L_ellip = L_ellip + gamma_ellip(R) * (y_s(R+1) - y_s(R));
    % end
    % 
    % L_ellip = 2 * L_ellip * rho * 1;

    % solve for kij
    D = 0;
    temp = 0;
    for I = 1:N
        for J = 1:N+1
            kij_s = ((y_s(J)-ybar_s(I)) * (y_s(I+1)-y_s(I))) / ...
                (y_s(J) - ybar_s(I))^2;
            kij_p = ((y_p(J)-ybar_s(I)) * (y_p(I+1)-y_p(I))) / ...
                (y_p(J) - ybar_s(I))^2;
            if J == 1
                small_gamma = -gamma(1);
            elseif J == N+1
                small_gamma = gamma(N);
            else
                small_gamma = gamma(J-1) - gamma(J);
            end
            D = D + gamma(I) * small_gamma * (kij_s + kij_p);
        end
    end
    D = (rho / (2*pi)) * D;
    cD(X) = D / (0.5 * rho * 1^2 * S);
    q = 0.5 * rho * 1;
    e_inv = L^2 / (q * pi * b^2 * D);
end

CD_cL = [0.007482262815025453, -0.18027073690930484;
    0.0061416363187905, -0.08790679453883343;
    0.005630921463081973, 0.004422326364203588;
    0.004709603470539586, 0.09257287117507929;
    0.005015742205069423, 0.1974551310883157;
    0.007414070774632547, 0.28546639002945273;
    0.01065101635158075, 0.36505085385139924;
    0.013879256561670245, 0.4530272913251019;
    0.01752680527545228, 0.5367903312392088;
    0.021580604442638873, 0.6289279340713548;
    0.02688362375404435, 0.7168173178764708;
    0.032588540835425034, 0.8172772514255037;
    0.039555736111312606, 0.9009010054698721];
CD0_cL = [0.006661022468251384, -0.16507109519157703
    0.006134809508303485, -0.08118962335829805;
    0.0050198632367307065, 0.019971779007923463;
    0.004326929338977531, 0.10402691848474976;
    0.003799413871703028, 0.19207641376316054;
    0.004604363399544121, 0.28290459133832635;
    0.006078801693259525, 0.36470205144903933;
    0.006551611852816673, 0.4517095408661673;
    0.00735786388798437, 0.5383696949962009;
    0.008082057961576034, 0.6209486595028763;
    0.008636926082709214, 0.7120373385433625;
    0.009862585477043311, 0.7899272766742644;
    0.011331813741452302, 0.8883968305655051];

alpha_cL = [-3.0811361981063428, -0.1573780043699926;
    -2.10050983248361, -0.06438455935906795;
    -0.9541150764748778, 0.032920611798980426;
    0.11070648215585521, 0.12177713037144944;
    1.0071376547705668, 0.21890750182083019;
    2.156154406409314, 0.30362709395484333;
    3.2209759650400542, 0.39248361252731256;
    4.284923525127454, 0.4855353241077933;
    5.264675892206846, 0.5827239621267297;
    6.248798252002906, 0.6589366351056082;
    7.395193008011653, 0.7562418062636563;
    8.457392571012377, 0.8576839038601602;
    9.522214129643118, 0.9465404224326294;
    10.50371449380917, 1.0353386744355428;
    12.22636562272396, 1.1666132556445739;
    13.128040786598689, 1.238572469045885;
    13.626219956300076, 1.2473124544792427;
    14.147997086671523, 1.1427822286962854];

cD_exp = CD_cL(:,1) - CD0_cL(:,1);
figure(1)
plot(cD_exp, CD_cL(:,2), 'ko')
hold on
plot(cD, cL, 'b--')
hold on
plot(cD(end), cL(end), 'rx', 'MarkerSize', 10)
xlabel('Coefficient of Drag')
ylabel('Coefficient of Lift')
set(gca,'XMinorTick','on','YMinorTick','on')
legend('Experimental Data', 'Calculated Data', ...
    'Max Cl', 'Location', 'northwest')
xlim([-0.01, 0.06])

figure(2)
plot(alpha_cL(:,1), alpha_cL(:,2), 'ko')
hold on
plot(alpha, cL, 'b--')
hold on
plot(alpha(end), cL(end), 'rx', 'MarkerSize', 10)
xlabel('Angle of Attack')
ylabel('Coefficient of Lift')
set(gca,'XMinorTick','on','YMinorTick','on')
legend('Experimental Data', 'Calculated Data', ...
    'Max Cl', 'Location', 'northwest')

% e_inv_plot = [1.0884, 1.0431, 1.0084, 1.0042,...
%     1.0021, 1.0008, 1.0004]; % 5,10,50,100,200,500,1000
% figure(1)
% plot([5,10,50,100,200,500,1000], e_inv_plot, 'b.', 'MarkerSize', 20)
% ylim([0.95,1.10])
% xlabel('Number of Panels')
% ylabel('Inviscid Span Efficiency')
% set(gca,'XMinorTick','on','YMinorTick','on')

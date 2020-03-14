clc
clear all

% Vortex 1-4 (starting at bottom left moving clockwise)
d = 1.0;
delta_t = 0.01;
t_init = 0.0;
t_final = 60.0;

gamma1 = [0, 0, -1]; % same strength for each vortex
gamma2 = [0, 0, 1];
gamma3 = [0, 0, 1];
gamma4 = [0, 0, -1];

pos_1_old = [0, -0.5, 0]; % initialize positions of vortices
pos_2_old = [0, 0.5, 0];
pos_3_old = [1.0, 0.5, 0];
pos_4_old = [1.0, -0.5, 0];

vel_1_old = [0, 0, 0]; % initialize velocities
vel_2_old = [0, 0, 0];
vel_3_old = [0, 0, 0];
vel_4_old = [0, 0, 0];

vortex_1 = zeros(4000,2); % creates empty matrices
vortex_2 = zeros(4000,2);
vortex_3 = zeros(4000,2);
vortex_4 = zeros(4000,2);

count = 1;

for R = t_init:delta_t:t_final
    % computes the x and y distances between each vortex
    dist_12_x = pos_2_old(1)-pos_1_old(1);
    dist_13_x = pos_3_old(1)-pos_1_old(1);
    dist_14_x = pos_4_old(1)-pos_1_old(1);
    dist_21_x = pos_1_old(1)-pos_2_old(1);
    dist_23_x = pos_3_old(1)-pos_2_old(1);
    dist_24_x = pos_4_old(1)-pos_2_old(1);
    dist_31_x = pos_1_old(1)-pos_3_old(1);
    dist_32_x = pos_2_old(1)-pos_3_old(1);
    dist_34_x = pos_4_old(1)-pos_3_old(1);
    dist_41_x = pos_1_old(1)-pos_4_old(1);
    dist_42_x = pos_2_old(1)-pos_4_old(1);
    dist_43_x = pos_3_old(1)-pos_4_old(1);
    dist_12_y = pos_2_old(2)-pos_1_old(2);
    dist_13_y = pos_3_old(2)-pos_1_old(2);
    dist_14_y = pos_4_old(2)-pos_1_old(2);
    dist_21_y = pos_1_old(2)-pos_2_old(2);
    dist_23_y = pos_3_old(2)-pos_2_old(2);
    dist_24_y = pos_4_old(2)-pos_2_old(2);
    dist_31_y = pos_1_old(2)-pos_3_old(2);
    dist_32_y = pos_2_old(2)-pos_3_old(2);
    dist_34_y = pos_4_old(2)-pos_3_old(2);
    dist_41_y = pos_1_old(2)-pos_4_old(2);
    dist_42_y = pos_2_old(2)-pos_4_old(2);
    dist_43_y = pos_3_old(2)-pos_4_old(2);
    
    % computes the induced velocity of each vortex
    vel_12 = cross(gamma1, [dist_12_x,dist_12_y,0])/(2*pi*(dist_12_x^2+dist_12_y^2));
    vel_13 = cross(gamma1, [dist_13_x,dist_13_y,0])/(2*pi*(dist_13_x^2+dist_13_y^2));
    vel_14 = cross(gamma1, [dist_14_x,dist_14_y,0])/(2*pi*(dist_14_x^2+dist_14_y^2));
    vel_21 = cross(gamma2, [dist_21_x,dist_21_y,0])/(2*pi*(dist_21_x^2+dist_21_y^2));
    vel_23 = cross(gamma2, [dist_23_x,dist_23_y,0])/(2*pi*(dist_23_x^2+dist_23_y^2));
    vel_24 = cross(gamma2, [dist_24_x,dist_24_y,0])/(2*pi*(dist_24_x^2+dist_24_y^2));
    vel_31 = cross(gamma3, [dist_31_x,dist_31_y,0])/(2*pi*(dist_31_x^2+dist_31_y^2));
    vel_32 = cross(gamma3, [dist_32_x,dist_32_y,0])/(2*pi*(dist_32_x^2+dist_32_y^2));
    vel_34 = cross(gamma3, [dist_34_x,dist_34_y,0])/(2*pi*(dist_34_x^2+dist_34_y^2));
    vel_41 = cross(gamma4, [dist_41_x,dist_41_y,0])/(2*pi*(dist_41_x^2+dist_41_y^2));
    vel_42 = cross(gamma4, [dist_42_x,dist_42_y,0])/(2*pi*(dist_42_x^2+dist_42_y^2));
    vel_43 = cross(gamma4, [dist_43_x,dist_43_y,0])/(2*pi*(dist_43_x^2+dist_43_y^2));
    
    if R == 0
        vel_1_new = [0,0,0];
        vel_2_new = [0,0,0];
        vel_3_new = [0,0,0];
        vel_4_new = [0,0,0];
    else
        % updates the velocities of each vortex
        vel_1_new = vel_21 + vel_31 + vel_41;
        vel_2_new = vel_12 + vel_32 + vel_42;
        vel_3_new = vel_13 + vel_23 + vel_43;
        vel_4_new = vel_14 + vel_24 + vel_34;
    end

    vel_1_old = vel_1_new;
    vel_2_old = vel_2_new;
    vel_3_old = vel_3_new;
    vel_4_old = vel_4_new;
    
    % updates the position of each vortex
    pos_1_new = pos_1_old + (vel_1_old * 0.01);
    pos_2_new = pos_2_old + (vel_2_old * 0.01);
    pos_3_new = pos_3_old + (vel_3_old * 0.01);
    pos_4_new = pos_4_old + (vel_4_old * 0.01);
    
    pos_1_old = pos_1_new;
    pos_2_old = pos_2_new;
    pos_3_old = pos_3_new;
    pos_4_old = pos_4_new; 
    
    % saves the position in its respective matrix (for plotting purposes)
    vortex_1(count,:) = [pos_1_old(1), pos_1_old(2)];
    vortex_2(count,:) = [pos_2_old(1), pos_2_old(2)];
    vortex_3(count,:) = [pos_3_old(1), pos_3_old(2)];
    vortex_4(count,:) = [pos_4_old(1), pos_4_old(2)];
    
    count = count + 1;
end

% creates plot for the simulation
figure(1)
plot(vortex_1(:,1), vortex_1(:,2), 'g')
hold on
plot(vortex_2(:,1), vortex_2(:,2), 'r')
plot(vortex_3(:,1), vortex_3(:,2), 'b')
plot(vortex_4(:,1), vortex_4(:,2), 'k')
legend('Vortex 1', 'Vortex 2', 'Vortex 3', 'Vortex 4')
xlabel('X-direction')
ylabel('Y-direction')
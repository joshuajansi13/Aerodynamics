clc
clear all

N = 4; % number of panels 
b = 1; % wing span
c = 0.2 * b;
max_span = 0.5 * b; 
lat_spacing = max_span / N;
center_spacing = lat_spacing / 2;

% control points and vortex points
% s = starboard; p = port
for R = 1:N
    if R == 1
        ym_s(R) = center_spacing;
        y1_n_s(R) = 0;
        y2_n_s(R) = lat_spacing;
    else
        ym_s(R) = ym_s(R-1) + lat_spacing;
        y1_n_s(R) = y1_n_s(R-1) + lat_spacing;
        y2_n_s(R) = y2_n_s(R-1) + lat_spacing;
    end
    xm(R)= ym_s(R) + 3 / 4 * c;
    x1_n(R) = y1_n_s(R) + 1 / 4 * c;
    x2_n(R) = y2_n_s(R) + 1 / 4 * c;
end
ym_p = -ym_s;
y1_n_p = -y1_n_s;
y2_n_p = -y2_n_s;

% calculate downwash coefficients for each panel (starboard and port)
for R = 1:N
    for I = 1:N
        temp_1 = 1 / ((xm(R)-x1_n(I))*(ym_s(R)-y2_n_s(I)) - (xm(R)-x2_n(I))*(ym_s(R)-y1_n_s(I)));
        temp_2 = ((x2_n(I)-x1_n(I))*(xm(R)-x1_n(I)) + (y2_n_s(I)-y1_n_s(I))*(ym_s(R)-y1_n_s(I))) / ...
            sqrt((xm(R)-x1_n(I))^2 + (ym_s(R)-y1_n_s(I))^2);
        temp_3 = ((x2_n(I)-x1_n(I))*(xm(R)-x2_n(I)) + (y2_n_s(I)-y1_n_s(I))*(ym_s(R)-y2_n_s(I))) / ...
            sqrt((xm(R)-x2_n(I))^2 + (ym_s(R)-y2_n_s(I))^2);
        temp_4 = (1 / (y1_n_s(I)-ym_s(R))) * (1 + (xm(R)-x1_n(I)) / ...
            sqrt((xm(R)-x1_n(I))^2 + (ym_s(R)-y1_n_s(I))^2));
        temp_5 = (1 / (y2_n_s(I)-ym_s(R))) * (1 + (xm(R)-x2_n(I)) / ...
            sqrt((xm(R)-x2_n(I))^2 + (ym_s(R)-y2_n_s(I))^2));
        wmn_s(R, I) = temp_1 * (temp_2 - temp_3) + temp_4 - temp_5; % starboard side
    end
end

for R = 1:N
    for I = 1:N
        temp_1 = 1 / ((xm(R)-x1_n(I))*(ym_s(R)-y2_n_p(I)) - (xm(R)-x2_n(I))*(ym_s(R)-y1_n_p(I)));
        temp_2 = ((x2_n(I)-x1_n(I))*(xm(R)-x1_n(I)) + (y2_n_p(I)-y1_n_p(I))*(ym_s(R)-y1_n_p(I))) / ...
            sqrt((xm(R)-x1_n(I))^2 + (ym_s(R)-y1_n_p(I))^2);
        temp_3 = ((x2_n(I)-x1_n(I))*(xm(R)-x2_n(I)) + (y2_n_p(I)-y1_n_p(I))*(ym_s(R)-y2_n_p(I))) / ...
            sqrt((xm(R)-x2_n(I))^2 + (ym_s(R)-y2_n_p(I))^2);
        temp_4 = (1 / (y1_n_p(I)-ym_s(R))) * (1 + (xm(R)-x1_n(I)) / ...
            sqrt((xm(R)-x1_n(I))^2 + (ym_s(R)-y1_n_p(I))^2));
        temp_5 = (1 / (y2_n_p(I)-ym_s(R))) * (1 + (xm(R)-x2_n(I)) / ...
            sqrt((xm(R)-x2_n(I))^2 + (ym_s(R)-y2_n_p(I))^2));
        wmn_p(R, I) = temp_1 * (temp_2 - temp_3) + temp_4 - temp_5; % port side
    end
end


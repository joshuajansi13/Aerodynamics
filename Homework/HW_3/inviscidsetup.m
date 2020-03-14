% inputs:
% Vt: vector of tangential velocities around airfoil
% x: corresponding x locations
% y: corresponding y locations
% 
% assumes x, y is ordered starting from trailing edge going clockwise
% around the airfoil back to the trailing edge 
% (the convention we used for Hess-Smith)
% 
% also assumes that the tangential velocity from the panel method is
% clockwise: positive on lower surface is upstream
% (also Hess-Smith convention)


function out = parseairfoil(Vt, x, y)
    % find stagnation point

    for R = 1:length(Vt)
        if Vt(R) > 0
            idx = R;
            break
        end
    end
    
    % separate into upper and lower surfaces
    xu = x(idx:end);
    yu = y(idx:end);
    Veu = Vt(idx:end);

    xl = x(idx-1:-1:1);
    yl = y(idx-1:-1:1);
    Vel = -Vt(idx-1:-1:1);  %negative sign because of definition of Vt

    % compute distance along curved path around airfoil
    nu = length(xu);
    su = zeros(nu);
    for i = 1:nu-1
        su(i+1) = su(i) + sqrt((xu(i+1) - xu(i))^2 + (yu(i+1) - yu(i))^2);
    end

    nl = length(xl);
    sl = zeros(nl);
    for i = 1:nl-1
        sl(i+1) = sl(i) + sqrt((xl(i+1) - xl(i))^2 + (yl(i+1) - yl(i))^2);
    end

    % rename to x (boundary layer definition for x)
    xu = su;
    xl = sl;

    out = [xu, Veu, xl, Vel];
end


% inputs
% xv: vector of x locations along boundary layer
% Vev: corresponding vector of edge velocities
% x: point to evaluate Ve and dVedx at
% 
% uses linear interpolation
% 
% outputs
% Ve: evaluated at x
% dVedx: evaluated at x


function out = velocity(xv, Vev, x)

    if x < xv(1)
        disp("provided x value is below xv range")
    elseif x > xv(end)
        disp("provided x value is above xv range")
    end

    % find index for linear interpolation
    for R = 1:length(x)
        if x >= xv
            idx = R;
        end
    end


    frac = (x - xv(idx)) / (xv(idx+1) - xv(idx));

    Ve = Vev(idx) +  frac*(Vev(idx+1) - Vev(idx));

    dVedx = (Vev(idx+1) - Vev(idx)) / (xv(idx+1) - xv(idx));

    out = [Ve, dVedx];
end
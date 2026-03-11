function [PRVs] = DCMtoPRV(C)
    % DCMtoPRV  Convert a DCM to the corresponding PRV
    %   [PRVs] = DCMtoPRV(C) calculates the rotation axes (unit vector)
    %   and the rotation angles (radians) from a 3x3 DCM. The PRV is a 4
    %   element vector containing the 3 axes and the angle, in that order.
    %
    % 
    
    arguments
        C (3,3) {mustBeNumeric}
    end

    if abs(C'*C - eye(3)) > 1e-3
        error("Please input an orthogonal matrix");
    end
    cosPhi = (trace(C) - 1) / 2;
    phi = acos(max(min(cosPhi, 1), -1));

    if phi < 1e-9
        axis = [1; 0; 0];
        phi = 0;
        return;
    end

    u_x_sin = [C(3,2) - C(2,3); 
               C(1,3) - C(3,1); 
               C(2,1) - C(1,2)];
    
    axis = u_x_sin / (2 * sin(phi));
    
    axis = axis / norm(axis);

    PRVs = {[axis' phi], [-axis' -phi], [axis' phi+2*pi], [-axis' 2*pi-phi]};
end
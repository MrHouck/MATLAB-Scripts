function C = PRVtoDCM(axis, phi)
    % PRVtoDCM  Convert a PRV to the corresponding DCM
    %   C = PRVtoDCM(axis, angle) calculates the direction cosine
    %   matrix, C, for a given axis and angle of rotation. The axis must be
    %   a vector with 3 elements. The rotation angle is input in radians.
    
    arguments
        axis (3,1) {mustBeNumeric}
        phi (1,1) {mustBeNumeric}
    end
    if abs(norm(axis)-1) > 1e-4
        warning("Axis is not a unit vector, correcting now")
        axis = axis/norm(axis);
    end
    c = cos(phi);
    s = sin(phi);
    u = axis(1);
    v = axis(2);
    w = axis(3);
    
    C = [(u^2*(1-c)+c),     (u*v*(1-c)-w*s),   (u*w*(1-c)+v*s);
         (v*u*(1-c)+w*s),   (v^2*(1-c)+c),     (v*w*(1-c)-u*s);
         (w*u*(1-c)-v*s),   (w*v*(1-c)+u*s),   (w^2*(1-c)+c)];
end
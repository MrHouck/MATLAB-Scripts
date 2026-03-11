function [yaw, pitch, roll] = DCMtoEA323(C)
    arguments
        C (3,3) {mustBeNumeric}
    end
    if any(abs(C'*C-eye(3)) > 1e-6, 'all')
        error("Input an orthonormal matrix")
    end
    pitch = acosd(C(3,3));
    yaw = atan2d(C(3, 2), C(3, 1));
    roll = atan2d(C(2,3),-C(1,3));
    
end
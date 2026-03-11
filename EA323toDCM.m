function C = EA323toDCM(yaw, pitch, roll)
    % EA323toDCM  Convert a 3-2-3 sequence of Euler Angles to a DCM
    % [C] = EA323toDCM(yaw, pitch, roll) converts the yaw, pitch, and roll
    % assuming a 3-2-3 rotation sequence. Yaw, pitch, and roll in radians.
    arguments
        yaw (1,1) {mustBeNumeric}
        pitch (1,1) {mustBeNumeric}
        roll (1,1) {mustBeNumeric}
    end

    c1 = cos(roll);
    c2 = cos(pitch);
    c3 = cos(yaw);
    s1 = sin(roll);
    s2 = sin(pitch);
    s3 = sin(yaw);

    C = [
        c3*c2*c1-s3*s1 c3*c2*s1+s3*c1 -c3*s2;
        -s3*c2*c1-c3*s1 -s3*c2*s1+c3*c1 s3*s2;
        s2*c1 s2*s1 c2
    ];

    

end
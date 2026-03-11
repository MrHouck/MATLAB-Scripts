function C = EPtoDCM(q)
    %assumes the convention where the quaternion vector is the first 3
    %components of the input vector, i.e. [q1:3, q4] rather than [q0, q1:3]
    %or other notation
    q1=q(1);
    q2=q(2);
    q3=q(3);
    q4=q(4);
    C = [
        1-(2*q2^2)-(2*q3^2) 2*(q1*q2+q3*q4) 2*(q1*q3 - q2*q4);
        2*(q1*q2-q3*q4) 1-(2*q1^2)-(2*q3^2) 2*(q2*q3 + q1*q4);
        2*(q1*q3+q2*q4) 2*(q2*q3 - q1*q4) 1-(2*q1^2)-(2*q2^2);
    ];
end
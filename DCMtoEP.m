function q = DCMtoEP(C)
    %shepperds method :sad:

    if abs(C'*C - eye(size(C))) > 1e-4
        error("C must be an orthonormal matrix.")
    end

    q12 = 0.25*(1+2*C(1,1) - trace(C));
    q22 = 0.25*(1+2*C(2,2) - trace(C));
    q32 = 0.25*(1+2*C(3,3) - trace(C));
    q42 = 0.25*(1+trace(C));

    biggestQ = max([q12 q22 q32 q42]);

    if biggestQ == q12
        q1 = sqrt(q12);
        q2 = ((C(1,2)+C(2,1))/4)/q1;
        q3 = ((C(3,1)+C(1,3))/4)/q1;
        q4 = ((C(2,3)-C(3,2))/4)/q1;
    elseif biggestQ == q22
        q2 = sqrt(q22);
        q1 = ((C(1,2)+C(2,1))/4)/q2;
        q3 = ((C(2,3)+C(3,2))/4)/q2;
        q4 = ((C(3,1)-C(1,3))/4)/q2;
    elseif biggestQ == q32
        q3 = sqrt(q32);
        q1 = ((C(3,1)+C(1,3))/4)/q3;
        q2 = ((C(2,3)+C(3,2))/4)/q3;
        q4 = ((C(1,2)-C(2,1))/4)/q3;
    elseif biggestQ == q42
        q4 = sqrt(q42);
        q1 = ((C(2,3)-C(3,2))/4)/q4;
        q2 = ((C(3,1)-C(1,3))/4)/q4;
        q3 = ((C(1,2)-C(2,1))/4)/q4;
    end
    q = [q1 q2 q3 q4];
end
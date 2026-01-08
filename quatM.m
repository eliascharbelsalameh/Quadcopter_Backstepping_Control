function M = quatM(q)
    q0 = q(1);
    qv = q(2:4).';
    sign_q0 = sign(q0);
    if sign_q0 == 0; sign_q0 = 1; end
    M = [sign_q0 * qv; 
        q0 * eye(3) + skew_symmetric(qv)];
end

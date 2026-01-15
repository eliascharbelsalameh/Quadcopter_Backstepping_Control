%This function allow to compiute the matrix M for rotational control

function M = quatM(q)
    q = q(:);
    q0 = q(1);
    qv = q(2:4); 

    qx = qv(1); qy = qv(2); qz = qv(3);

    S = [  0   -qz   qy;
          qz    0   -qx;
         -qy   qx    0  ];

    sign_q0 = sign(q0);

    if sign_q0 == 0; sign_q0 = 1; end
    M = [sign_q0 * qv.'; 
        q0 * eye(3) + S]; 

    
end

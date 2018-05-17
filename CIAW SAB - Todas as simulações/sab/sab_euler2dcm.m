% Rotation matrix from Euler angles [phi theta psi] <=> [roll pitch yaw]
% rotation sequence is 'ZYX'

function dcm = sab_euler2dcm( euler )
    c = cos( euler );
    s = sin( euler );
    
    %     [          cy*cz,          cy*sz,            -sy]
    %     [ sy*sx*cz-sz*cx, sy*sx*sz+cz*cx,          cy*sx]
    %     [ sy*cx*cz+sz*sx, sy*cx*sz-cz*sx,          cy*cx]

    dcm = zeros(3,3);
    dcm(1,1) = c(2)*c(3);
    dcm(1,2) = c(2)*s(3);
    dcm(1,3) = -s(2);
    dcm(2,1) = s(1)*s(2)*c(3) - c(1)*s(3);
    dcm(2,2) = s(1)*s(2)*s(3) + c(1)*c(3);
    dcm(2,3) = s(1)*c(2);
    dcm(3,1) = c(1)*s(2)*c(3) + s(1)*s(3);
    dcm(3,2) = c(1)*s(2)*s(3) - s(1)*c(3);
    dcm(3,3) = c(1)*c(2);
end
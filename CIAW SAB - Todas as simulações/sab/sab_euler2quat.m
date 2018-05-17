% Quaternion from Euler angles [phi theta psi] <=> [roll pitch yaw]
% rotation sequence is 'ZYX'

function q = sab_euler2quat(euler)
    
    c = cos( euler/2 );
    s = sin( euler/2 );
    
    q = [ c(1)*c(2)*c(3) + s(1)*s(2)*s(3), ...
          c(1)*c(2)*s(3) - s(1)*s(2)*c(3), ...
          c(1)*s(2)*c(3) + s(1)*c(2)*s(3), ...
          s(1)*c(2)*c(3) - c(1)*s(2)*s(3)];
end
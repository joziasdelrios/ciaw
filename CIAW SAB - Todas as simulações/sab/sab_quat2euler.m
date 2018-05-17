function euler = sab_quat2euler( q )
 
    q = sab_quatnormalize( q );
    
    % The default rotation sequence is 'ZYX' where the order of rotation
    % angles for the default rotation are:
    %   e(1) = phi   - x axis rotation, 
    %   e(2) = theta - y axis rotation, 
    %   e(3) = psi   - z axis rotation, 

    euler = [0 0 0]';
    euler(3) = atan2( 2*( q(2)*q(3) + q(1)*q(4)),  q(1)^2 + q(2)^2 - q(3)^2 - q(4)^2 );
    euler(2) = asin( -2*( q(2)*q(4) - q(1)*q(3)) );
    euler(1) = atan2( 2*( q(3)*q(4) + q(1)*q(2)),  q(1)^2 - q(2)^2 - q(3)^2 + q(4)^2 );
end
    
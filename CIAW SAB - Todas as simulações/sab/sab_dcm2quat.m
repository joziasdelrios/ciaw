%  dcm2quat 
%
%   Convert direction cosine matrix (DCM 3x3 matrix) to Quaternion.
%   Quaternion is a row vector, scalar part is the first element.

function q = sab_dcm2quat( dcm )
    q =  [0 0 0 0]; 
    
    tr = trace(dcm);

    if (tr > 0)
        sqtrp1 = sqrt( tr + 1.0 );
        
        q(1) = 0.5*sqtrp1; 
        q(2) = (dcm(2,3) - dcm(3,2)) / (2.0*sqtrp1);
        q(3) = (dcm(3,1) - dcm(1,3)) / (2.0*sqtrp1); 
        q(4) = (dcm(1,2) - dcm(2,1)) / (2.0*sqtrp1); 
    else
        d = diag(dcm);
        
        if (d(2) > d(1)) && (d(2) > d(3))
            % max value at dcm(2,2)
            sqdip1 = sqrt(d(2) - d(1) - d(3) + 1.0 );
            
            q(3) = 0.5*sqdip1; 
            
            if ( sqdip1 ~= 0 )
                sqdip1 = 0.5/sqdip1;
            end
            
            q(1) = (dcm(3,1) - dcm(1,3))*sqdip1; 
            q(2) = (dcm(1,2) + dcm(2,1))*sqdip1; 
            q(4) = (dcm(2,3) + dcm(3,2))*sqdip1; 
        
        elseif d(3) > d(1)
            % max value at dcm(3,3,i)
            sqdip1 = sqrt(d(3) - d(1) - d(2) + 1.0 );
            
            q(4) = 0.5*sqdip1; 
            
            if ( sqdip1 ~= 0 )
                sqdip1 = 0.5/sqdip1;
            end
            
            q(1) = (dcm(1,2) - dcm(2,1))*sqdip1;
            q(2) = (dcm(3,1) + dcm(1,3))*sqdip1; 
            q(3) = (dcm(2,3) + dcm(3,2))*sqdip1; 
        else
            % max value at dcm(1,1)
            sqdip1 = sqrt(d(1) - d(2) - d(3) + 1.0 );
            
            q(2) = 0.5*sqdip1; 
            
            if ( sqdip1 ~= 0 )
                sqdip1 = 0.5/sqdip1;
            end
            
            q(1) = (dcm(2,3) - dcm(3,2))*sqdip1; 
            q(3) = (dcm(1,2) + dcm(2,1))*sqdip1; 
            q(4) = (dcm(3,1) + dcm(1,3))*sqdip1; 
        end
    end
end




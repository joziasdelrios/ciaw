% MATLAB: fun��o de rejei��o de u em v

function r = vecrej(u, v)
    vn = vecdir(v);
    r = cross(vn, cross(u, vn) );
end

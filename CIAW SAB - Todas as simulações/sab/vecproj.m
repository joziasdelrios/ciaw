% MATLAB: função de projeção de u em v

function p = vecproj(u, v)
    vn = vecdir(v);
    p = vn * dot(u, vn);
end

function qout = sab_quatnormalize(q)
    qout = q ./ norm(q);
end
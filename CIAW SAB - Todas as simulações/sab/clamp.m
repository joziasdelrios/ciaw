% clamp: satura um valor num intervalo [-limite,+limite]

function y = clamp(x, lim)
    y = min( +lim, max( -lim, x ) );
end
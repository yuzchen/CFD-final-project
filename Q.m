function result = Q(x,delta)
    result = abs(x) .* (abs(x)>=delta) + (x.^2 + delta^2) / 2 / delta .* (abs(x)<delta);
end
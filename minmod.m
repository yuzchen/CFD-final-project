function phi = minmod(vector)
    if  sum(sign(vector)) == size(vector,2)
        phi = min(vector);
    elseif sum(sign(vector)) == - size(vector,2)
        phi = max(vector);
    else 
        phi = 0;
    end
end

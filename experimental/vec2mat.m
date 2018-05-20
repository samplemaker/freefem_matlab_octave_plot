function [M] = vec2mat (V, c)
    r = length (V) / c;
    M = reshape (V, c, r);
end

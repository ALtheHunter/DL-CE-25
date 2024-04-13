function [out]=damping(old, new, m)
    
    out = m .* old + ( 1 - m ) .* new;
end
function r = mobius_on_point(T, z)
% MOBIUS_ON_POINT - calculates the image of a Möbius map applied to a
% complex number
%
% T - a 2x2 matrix representing the Möbius map. T = [a, b; c, d] means
%     T(z) = (a*z + b) / (c*z + d)
%
% z - the complex number to apply the map to
%
% Conley May 2013, adapted from page 75 of Indra's Pearls by Mumford et al.

if (z == Inf)                % z is infinity
    if (T(2,1) ~= 0)
        r = T(1,1) / T(2,1); % c ~= 0, return a/c
    else
        r = Inf;             % c == 0, return inifinity
    end
else    
    n = T(1,1)*z + T(1,2);
    d = T(2,1)*z + T(2,2);
    if (d == 0)
        r = Inf;             % if c*z + d == 0, return infinity
    else
        r = n / d;
    end
end

end
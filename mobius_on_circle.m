function [new_center, new_radius] = mobius_on_circle(T, center, radius)
% MOBIUS_ON_CIRCLE - calculates the image of a Möbius map applied to a
% circle in the complex plane
%
% T - a 2x2 matrix representing the Möbius map. T = [a, b; c, d] means
%     T(z) = (a*z + b) / (c*z + d)
%
% center - the center of the circle
%
% radius - the radius of the circle
%
% Conley May 2013, adapted from page 91 of Indra's Pearls by Mumford et al.

z = center - (radius^2)/conj((T(2,2)/T(2,1)) + center);
new_center = mobius_on_point(T,z);
new_radius = abs(new_center - mobius_on_point(T, center + radius));

end
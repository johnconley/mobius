function dfs(max_level)
% DFS - implementation of depth first search for fractal plotting
%
% max_level - maximum depth to plot (50 looks nice but it takes a while)
%
% consulted
% http://www.mcs.uvawise.edu/msh3e/resources/MobiusDFS/MobiusDFS.pde for
% use of repet array
%
% Conley May 2013. Adapted from page 148 of Indra's Pearls by Mumford et al


% % From the website
% ta = 1.95 + 0.2*1i;
% tb = 2.01 + 0.0*1i;

% Apollonian gasket
% axis([-1.2 1.2 -1.2 1.2]);
% axis([-.2 .2 0 .4]); for close-up
ta = 2;
tb = 2;

% % Jordan curve from page 238
% % axis([-1.2 1.2 -1.2, 1.2]);
% ta = 1.91 + .05*1i;
% tb = 3;

% % Jordan curve from page 226 (modified)
% % axis([-1.5 1.5 -2.5, 2.5]);
% ta = 1.87 + .5*1i;
% tb = 1.87 - .5*1i;

% given ta and tb, calculate tab and z0
tab = (ta*tb - sqrt((ta*tb)^2 - 4*((ta)^2 + (tb)^2)))/2;
z0 = ((tab - 2)*tb)/(tb*tab - 2*ta + 2*1i*tab);

% from ta, tb, tab, and z0, calculate the generator matrices
gens{1} = [(ta/2), ((ta*tab-2*tb+4*1i)/((2*tab+4)*z0)); (((ta*tab-2*tb-4*1i)*z0)/(2*tab-4)), (ta/2)];
gens{2} = [((tb-2*1i)/2), (tb/2); (tb/2), ((tb+2*1i)/2)];
gens{3} = matrix_inverse(gens{1});
gens{4} = matrix_inverse(gens{2});

% generate the repet array
repet(1,1) = find_fix(gens{2}*gens{3}*gens{4}*gens{1});
repet(1,2) = find_fix(gens{1});
repet(1,3) = find_fix(gens{4}*gens{3}*gens{2}*gens{1});
repet(2,1) = find_fix(gens{3}*gens{4}*gens{1}*gens{2});
repet(2,2) = find_fix(gens{2});
repet(2,3) = find_fix(gens{1}*gens{4}*gens{3}*gens{2});
repet(3,1) = find_fix(gens{4}*gens{1}*gens{2}*gens{3});
repet(3,2) = find_fix(gens{3});
repet(3,3) = find_fix(gens{2}*gens{1}*gens{4}*gens{3});
repet(4,1) = find_fix(gens{1}*gens{2}*gens{3}*gens{4});
repet(4,2) = find_fix(gens{2});
repet(4,3) = find_fix(gens{3}*gens{2}*gens{1}*gens{4});

tags = zeros(1, 10000);
word = cell(1, 10000);

level = 1;
epsilon = .001;
tags(1) = 1;
tags(2) = 1;
word{1} = gens{1};

% set the first oldpoint for plotting
p = gens{1}*gens{2}*gens{3}*gens{4};
fp = find_fix(p);
oldpoint = fp;

% formatting for the plot. see commented code near traces for what to set
% the axis values to
axis([-1.2 1.2 -1.2 1.2]);
axis equal;
set(gca, 'visible', 'off');

% main loop of DFS
while (true)
    [done, oldpoint] = end_of_branch(level, epsilon, tags, word, oldpoint, max_level, repet);
    while (~done)
        [level, tags, word] = go_forward(level, tags, word, gens);
        
        [done, oldpoint] = end_of_branch(level, epsilon, tags, word, oldpoint, max_level, repet);
    end

    % go backward until we can go forward or we are at level 0
    level = go_backward(level);
    while ((level ~= 0) && (~turn_available(level, tags)))
        level = go_backward(level);
    end
    
    % if we are at level 0 and there are no more turns to make, break
    if ((level <= 0) && (tags(1) == 2))
        break;
    end
    
    % there are still turns to make, so turn_and_go_forward
    [level, tags, word] = turn_and_go_forward(level, tags, word, gens);
    
end

fprintf('done\n\n');

end

function [lev, tags, word] = go_forward(level, tags, word, gens)
% GO_FORWARD - step forward one level and update tags and word

lev = level + 1;
tags(lev) = mod(tags(lev - 1), 4) + 1;
word{lev} = word{lev - 1} * gens{tags(lev)};

end

function lev = go_backward(level)
% GO_BACKWARD - step back one level

lev = level - 1;

end

function b = turn_available(level, tags)
% TURN_AVAILABLE - return true if there is an unexplored turn at the
% current level, else return false

if ( mod(tags(level + 1) - 1, 4) == mod(tags(level) + 2, 4))
    b = false;
else
    b = true;
end

end

function [lev, tags, word] = turn_and_go_forward(level, tags, word, gens)
% TURN_AND_GO_FORWARD - move ot the left on the current level and go
% forward one level from there, updating tags and word

tags(level + 1) = mod(tags(level + 1) - 1, 4);
if (tags(level + 1) == 0)
    tags(level + 1) = 4;
end

if (level == 0)
    word{1} = gens{tags(1)};
else
    word{level + 1} = word{level} * gens{tags(level + 1)};
end
lev = level + 1;

end

function [b, new_old] = end_of_branch(level, epsilon, tags, word, oldpoint, max_level, repet)
% END_OF_BRANCH - if the the points we are plotting are close enough or we
% have reached the maximum depth, plot the points, update the oldpoint, and
% return true. Else return false

% for use with repet
newpoint = mobius_on_point(word{level}, oldpoint);
dist = abs(newpoint - oldpoint);
z = zeros(1,3);
for i = 1:3
    z(i) = mobius_on_point(word{level}, repet(tags(level),i));
end
if (((dist < epsilon) && (abs(z(1) - z(2)) < epsilon) && (abs(z(2) - z(3)) < epsilon)) || (level >= max_level))
    line([real(oldpoint) real(newpoint)], [real(imag(oldpoint)) real(imag(newpoint))], 'Color', [1 0 0]);
    hold on;
    b = true;
    new_old = newpoint;
else
    b = false;
    new_old = oldpoint;
end

end

function z = find_fix(p)
% FIND_FIX - find the fixed point of a Mobius transformation

z = ((p(1,1)-p(2,2)) - sqrt((p(1,1)-p(2,2))^2 + 4*p(1,2)*p(2,1)))/(2*p(2,1));

end

function A = matrix_inverse(a)
% MATRIX_INVERSE - calculate the inverse of a 2x2, determinant-1 matrix

A(1,1) = a(2,2);
A(2,2) = a(1,1);
A(1,2) = -a(1,2);
A(2,1) = -a(2,1);

end


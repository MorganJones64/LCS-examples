function [dy] = doublegyreVEC(t,yin,A,eps,om)
x = yin(1,:,:);
y = yin(2,:,:);

u = zeros(size(x)); v = u;

a = eps * sin(om * t);
b = 1 - 2 * a;

f = a * x.^2 + b * x;
df = 2 * a * x + b;

u = -pi * A * sin(pi * f) .* cos(pi * y);
v =  pi * A * cos(pi * f) .* sin(pi * y) .* df;

dy = [u;v];
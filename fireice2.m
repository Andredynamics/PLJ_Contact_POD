function cmap = bbw(m)

% LightCyan-Cyan-Blue-Black-Red-Yellow-LightYellow
clrs = [0 0 .5; 0 0 1;0 1 1; ...
    0 0 0; 1 1 0; 1 0 0; 0.5 0 0];

% clrs = [0 0 0; 0 0 1; 0 1 1;...
%         0 1 1; 1 1 0; 1 0 0; 1 1 1];


y = -3:3;
if mod(m,2)
    delta = min(1,6/(m-1));
    half = (m-1)/2;
    yi = delta*(-half:half)';
else
    delta = min(1,6/m);
    half = m/2;
    yi = delta*nonzeros(-half:half);
end
cmap = interp2(1:3,y,clrs,1:3,yi);

end
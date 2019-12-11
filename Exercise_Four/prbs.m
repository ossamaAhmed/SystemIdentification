%Initialization
clc;
close all;

% Screen size used to place plots
screensize = get(groot, 'ScreenSize');
screenwidth = screensize(3);
screenheight = screensize(4);

% Generate PRBS using idinput
order = [5, 6, 7, 8];
size_order = size(order, 2);
N = (2 .^ order) - 1;
Band = [0,1];
Range = [-5,5];
u = cell(size_order,1);
u_P = cell(size_order,1);
Wn = cell(size_order,1);
figure(1);
fig = gcf;
fig.Position = [0, 0, screenwidth, 12/13 * screenheight];
for i = 1 : size_order
    u{i} = idinput([N(i), 1, 1], 'prbs', Band, Range)';
    u_P{i} = periodogram(u{i});
    s = size(u_P{i},2);
    omega = (2*pi/s)*[0:s-1];
    idx = find(omega > 0 & omega < pi);
    Wn{i} = omega(idx);
    figure(1);
    Pplot = u_P{i};
    semilogx(Wn{i}, Pplot(idx));
    hold on;
end
function plotAntennaResponse(fname,varargin)

figure;
angles = h5read(fname,'/angles');
power = h5read(fname,'/antennaPower');

h = polar(angles,power/max(power));
set(h,'LineWidth',2)

if nargin == 2
    tstring = varargin{1};
    title(tstring);
end

set(gca,'FontWeight','bold');
set(gca,'LineWidth',2);
set(gca,'FontSize',12);

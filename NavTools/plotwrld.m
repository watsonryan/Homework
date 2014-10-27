%
%   PLOTWRLD.M
%
%   This routine loads the world geographic data and
%   sets up the color map for imaging the globe
%
%   The world topographic data (GTOPO30) was made available
%   courtesy of the U.S. Geological Survey, EROS
%   Data Center, Sioux Falls, South Dakota
%
%	M. & S. Braasch 6-98
%	Copyright (c) 1998 by GPSoft LLC
%	All Rights Reserved.
%
clf
close
load wrldmatr

wrldcmap = ...
   [0 0 1;
    0 .86 0;
    0 .87 0;
   .5 .88 0;
   .6 .89 0;
   .70 .90 0;
   .72 .91 0;
   .74 .92 0;
   .76 .6 0;
   .78 .65 0;
   .8 .7 0;
   .82 .75 0;
   .84 .8 0;
   .85 .85 0;
   .9 .9 0;
   .96 .96 0;
   .97 .97 0;
   .98 .98 0;
   .99 .99 0;
   1 1 0];

[x,y,z]=sphere(24);
h=surface(x,y,z,'FaceColor','texture','CData',wrldmatr);
view(0,25)
colormap(wrldcmap)
set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1], ...
	'NextPlot','add', ...
   'Visible','off');


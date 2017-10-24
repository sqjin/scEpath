function generate_three_dimension_ball(x0,y0,z0,color,x1range,x2range)
%% plot an ellipsoid with center (x0,y0,z0) and semi-axis lengths (xr,yr,zr)
% you should change xr, yr, zr based on your data, otherwise the ball will look like an oval
% xr = 1*(max(x1range)-min(x1range))*0.045/10; yr = 1*(max(x2range)-min(x2range))*0.045/10; zr = 0.075*1;
xr = 1*(max(x1range)-min(x1range))*0.045/3; yr = 1*(max(x2range)-min(x2range))*0.045/3; zr = 0.075*0.2;
[x,y,z]  = ellipsoid(x0,y0,z0,xr,yr,zr);
surf(x,y,z,'FaceColor',color,'EdgeColor','none')
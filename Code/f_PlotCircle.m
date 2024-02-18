function f_PlotCircle(Radius,Center,Normal)
%https://it.mathworks.com/matlabcentral/answers/111144-2d-circle-3d-plot

phi = atan2(Normal(2),Normal(1)); %azimuth angle, in [-pi, pi]
theta = atan2(sqrt(Normal(1)^2 + Normal(2)^2) ,Normal(3));% zenith angle, in [0,pi]    
t = 0:pi/32:2*pi;
x = Center(1)- Radius*( cos(t)*sin(phi) + sin(t)*cos(theta)*cos(phi) );
y = Center(2)+ Radius*( cos(t)*cos(phi) - sin(t)*cos(theta)*sin(phi) );
z = Center(3)+ Radius*sin(t)*sin(theta);
plot3(x,y,z,'LineWidth', 2,'Color','b')

end
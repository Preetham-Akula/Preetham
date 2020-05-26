clc; clear all;
filename = uigetfile('*.txt');
TLEdata = readmatrix(filename);
i = TLEdata(2,3); %inclination of orbit in degrees
%disp(i);
o = TLEdata(2,4); %RAAN of satellite in degrees
%disp(o);
e = TLEdata(2,5)*10^(-7); % Eccentricity of orbit
%disp(e);
w = TLEdata(2,6); % Argument of perigee in degrees
%disp(w);
M0 = TLEdata(2,7); % Mean Anamoly in degrees
%disp(M0);
n = TLEdata(2,8); % revoluitons per day

%constants used in plot
G = 398600.4418; % Geocentric gravitational constant
Re = 6378; % Radius of earth at equator

% Unit conversions
o = o*pi/180; % degress to radians
i = i*pi/180; % degress to radians
M0 = M0*pi/180; % degress to radians
w = w*pi/180; % degress to radians

% Derived calculations
N = n*2*pi/86164; % Revolutions per minute to radians per second
a = G^(1/3)/N^(2/3); % Semi-major axis
disp(a);
P = a*(1-e^2); % Semi latus rectum
rp = a*(1-e); % Radius at perigee
ra = a*(1+e); % radius at apogee
vp = sqrt(G*(2/rp-1/a)); %Velocity at perigee in km/s
va = sqrt(G*(2/ra-1/a)); % Velocity at apogee in km/s
T = 2*pi/N; % Period of satellite
hours = floor(T/3600); %Hours
minutes = floor((T-hours*3600)/60); %Minutes
seconds = floor(T-hours*3600-minutes*60); %Seconds

%true;
% fprintf('\n Period is %3d h: %3d m: %3d s', hours,minutes, seconds);
% fprintf('\n Radius at Perigee %10.3f km    Radius at Perigee %10.3f km', rp, rp-Re);
% fprintf('\n Altitude at Apogee %10.3f km    Altitude at Perigee %10.3f km', ra, ra-Re);
% fprintf('\n Velocity at Perigee %6.4f km/s    Velocity at Apogee %6.4f km/s', vp, va);
% disp(i);
% disp(e);

% Plotting the orbit
n0=1; %number of orbits to plot
t0=0; %intial time
tn = n0*T; %Final time
s = 30; %steps in secods
t=t0:s:tn+s; %time vector

%Ecentricity anomoly
M=M0+n*(t-t0);
E = zeros(size(t,2),1);
for j=1:size(t,2)
    E(j) = anom_ecc(M(j),e);
end

%Calculation of radius
sinv = (sqrt(1-e.^2).*sin(E))./(1-e.*cos(E)); % Sine of true anomoly
cosv = (cos(E)-e)./(1-e.*cos(E)); % cosine of true anomoly
v = atan2(sinv, cosv); %true anomoly
u = v + w; % Argument of latitude
r = (a*(1-e^2))/(1+e*cos(v)); % Radius

% ECI Earth centered Inertial
xp = r.*cos(u);
yp = r.*sin(u);
xs = xp.*cos(o)-yp.*cos(i).*sin(o);  %ECI x
ys = xp.*sin(o)+yp.*cos(i).*cos(o); %ECI y
zs = yp.*sin(i); %ECI z

hold on;
[x,y,z] = ellipsoid(0,0,0,Re, Re, Re, 40);
hsurface = surf(x,y,z);
set(hsurface, 'Facecolor', [0 0 1], 'Facealpha', 0.9, 'Facelighting', 'gouraud', 'Edgecolor', 'k')
axis equal
plot3(xs, ys, zs, 'r*')
angle_eq = linspace(0,2*pi, 361);
xeq = (Re*1.0001).*cos(angle_eq);
yeq = (Re*1.0001).*sin(angle_eq);
zeq = zeros(1,size(angle_eq,2));
plot3([0,2*ra],[0,0],[0,0],'-.k','Linewidth',1); %X
plot3([0,0],[0,2*ra],[0,0],'-.k','Linewidth',1); %Y
plot3([0,0],[0,0],[0,2*ra],'-.k','Linewidth',1); %Z
text(2*ra+120,10,0,texlabel('X'),'Color','k','Fontsize',16);
text(0,2*ra+120,0,texlabel('Y'),'Color','k','Fontsize',16);
text(0,0,2*ra+120,texlabel('Z'),'Color','k','Fontsize',16);
plot3(xeq, yeq, zeq, '--y','Linewidth',1); %Equator
title('Output: Interactive Orbit View')
grid on;

view(135, 45)








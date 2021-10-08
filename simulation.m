%% Approach 1
clear all; clc;
% Code for Rectangular Integration based Kinematic model

T = 2; % Sampling time
ti = 0; % Initial time
tf = 100; % Final time
td = ti:T:tf; % Discrete time defined
tf = td(end);
% Due to discrete steps, t_d(end) is not necessarily equal to tf.
% That is why we need to recompute tf by assigning t_d(end) to it.
sl = length(td); % Number of samples in discrete time frame
wR = 0.3*cos(0.002*td); % Right wheel angular speed
wL = sin(0.002*td); % Left wheel angular speed
r = 5; % Radius of wheels [cm]
L = 10; % Width of robot platform [cm]
xd = zeros(1,sl); % Vector to store global x_n
yd = zeros(1,sl); % Vector to store global y_n
xd1 = zeros(1,sl); % Vector to store global x_n for trapezoidal integration
yd1 = zeros(1,sl); % Vector to store global y_n for trapezoidal integration
xd2 = zeros(1,sl); % Vector to store global x_n for exact integration
yd2 = zeros(1,sl); % Vector to store global y_n for exact integration
xd3 = zeros(1,sl); % Vector to store global x_n for geometry based 1 integration
yd3 = zeros(1,sl); % Vector to store global y_n for geometry based 1 integration
xd4 = zeros(1,sl); % Vector to store global x_n for taylor integration
yd4 = zeros(1,sl); % Vector to store global y_n for taylor integration


phd = zeros(1,sl); % Vector to store global \phi_n
phd_exact = zeros(1,sl); % Vector to store global \phi_n for exact method

% continuous model

dt = 0.01; % Time increment for continuous scenario
tc = ti:dt:tf; % Continuous time (CT) frame
wR_c = 0.3*cos(0.002*tc); % Right wheel angular speed
wL_c = sin(0.002*tc); % Left wheel angular speed
sl_c = length(tc); % Number of steps in CT simulation
xc = zeros(1,sl_c); % Vector to store CT x(t)
yc = zeros(1,sl_c); % Vector to store CT y(t)
phc = zeros(1,sl_c); % % Vector to store CT \phi(t)

for ii=1:sl_c-1 % <= one less simulation for the reason mentioned above
vR_c=r*wR_c(ii);
vL_c=r*wL_c(ii);
v_c=(vR_c+vL_c)/2;
w_c=(vR_c-vL_c)/L;
dph_c=dt*w_c;
xc(ii+1)=xc(ii)+v_c*dt*cos(phc(ii));
yc(ii+1)=yc(ii)+v_c*dt*sin(phc(ii));
phc(ii+1)=phc(ii)+dph_c;
end

for ii=2:sl-1
% note that we need to run this loop for one less iteration because
% the pose quantities have time index (ii+1)
vR=r*wR(ii);
vL=r*wL(ii);
v=(vR+vL)/2;
w=(vR-vL)/L;
dph=T*w;
xd(ii+1)=xd(ii)+v*T*cos(phd(ii));
yd(ii+1)=yd(ii)+v*T*sin(phd(ii));

xd1(ii+1)=xd1(ii)+(v*T./2)*(cos(phd(ii)) + cos(phd(ii) + T*w));
yd1(ii+1)=yd1(ii)+(v*T./2)*(sin(phd(ii)) + sin(phd(ii) + T*w));

xd2(ii)=xd2(ii-1)+(v./w)*(sin(phd_exact(ii-1)+dph) - sin(phd_exact(ii-1)));
yd2(ii)=yd2(ii-1)-(v./w)*(cos(phd_exact(ii-1)+dph) - cos(phd_exact(ii-1)));
phd_exact(ii) = phd_exact(ii-1)+dph;

xd3(ii+1)=xd3(ii)+v*T*(cos(phd(ii)+dph/2) );
yd3(ii+1)=yd3(ii)+v*T*(sin(phd(ii)+dph/2) );

phd(ii+1)=phd(ii)+dph;
end
hold on;
plot(xc,yc,'k','linewidth',1)
plot(xd,yd,'ks:','linewidth',1);
plot(xd1,yd1,'-o','linewidth',1);
plot(xd2,yd2,'v','linewidth',1);
plot(xd3,yd3,'x','linewidth',1);
legend('Çontinuous','Rectangular','Trapezoidal','Éxact','Geometry-Based 1')
sgtitle('Discrete-time Kinematic models')

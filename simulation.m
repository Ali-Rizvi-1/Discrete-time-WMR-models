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

phd = zeros(1,sl); % Vector to store global \phi_n

for ii=1:sl-1
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

phd(ii+1)=phd(ii)+dph;
end
hold on;
plot(xd,yd,'ks:','linewidth',1);
plot(xd1,yd1,'-o','linewidth',1);
legend('Rectangular','Trapezoidal')
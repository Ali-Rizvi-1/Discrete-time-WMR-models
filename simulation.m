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
xd5 = zeros(1,sl); % Vector to store global x_n for RK4 integration
yd5 = zeros(1,sl); % Vector to store global y_n for RK4 integration

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
xd(ii+1)=xd(ii)+v*T*cos(phd(ii)); % rectangular model
yd(ii+1)=yd(ii)+v*T*sin(phd(ii));

xd1(ii+1)=xd1(ii)+(v*T./2)*(cos(phd(ii)) + cos(phd(ii) + T*w)); % Trapezoidal
yd1(ii+1)=yd1(ii)+(v*T./2)*(sin(phd(ii)) + sin(phd(ii) + T*w));

xd2(ii)=xd2(ii-1)+(v./w)*(sin(phd_exact(ii-1)+dph) - sin(phd_exact(ii-1))); % Exact
yd2(ii)=yd2(ii-1)-(v./w)*(cos(phd_exact(ii-1)+dph) - cos(phd_exact(ii-1)));
phd_exact(ii) = phd_exact(ii-1)+dph;

xd3(ii+1)=xd3(ii)+v*T*(cos(phd(ii)+dph/2) ); % Geometry Based I
yd3(ii+1)=yd3(ii)+v*T*(sin(phd(ii)+dph/2) );

xd4(ii)=xd4(ii-1)+((v*dph)./w)*(cos(phd_exact(ii-1)) - 0.5*sin(phd_exact(ii-1))*dph); % Taylor
yd4(ii)=yd4(ii-1)+((v*dph)./w)*(sin(phd_exact(ii-1)) + 0.5*cos(phd_exact(ii-1))*dph);
phd_exact(ii) = phd_exact(ii-1)+dph;

xd5(ii+1) = xd5(ii) + ((T*v)/6)*( cos(phd(ii)) + 2*cos(phd(ii)+dph/2) ... % RK4
    + 2*cos(phd(ii)+dph) + 2*cos(phd(ii)+(3*dph)/2) ...
    + cos(phd(ii)+4*dph) );
yd5(ii+1) = yd5(ii) + ((T*v)/6)*( sin(phd(ii)) + 2*sin(phd(ii)+dph/2) ...
    + 2*sin(phd(ii)+dph) + 2*sin(phd(ii)+(3*dph)/2) ...
    + sin(phd(ii)+4*dph) );

phd(ii+1)=phd(ii)+dph;
end
hold on;
plot(xc,yc,'k','linewidth',2);
plot(xd,yd,'ks:','linewidth',1.5);
plot(xd1,yd1,'-o','linewidth',1.5);
plot(xd2,yd2,'v','linewidth',1.5);
plot(xd3,yd3,'x','linewidth',1.5);
plot(xd4,yd4,'^','linewidth',1.5);
plot(xd5,yd5,'--','linewidth',1.5);
legend('Çontinuous','Rectangular','Trapezoidal','Éxact','Geometry-Based 1',...
    'Taylor','RK4','fontsize',16,'location','best')
xlabel('x','fontsize',16); ylabel('y','fontsize',16);
sgtitle('Discrete-time Kinematic models','fontsize',18)

%% Experiment 2
clear all; clc;

%Compute trajectory for continuous-time for dt = 0:01
ti = 0; % Initial time
tf = 100; % Final time
dt = 0.01; % Time increment for continuous scenario
tc = ti:dt:tf; % Continuous time (CT) frame
r = 5; % Radius of wheels [cm]
L = 10; % Width of robot platform [cm]
wR_c = 0.3*cos(0.002*tc); % Right wheel angular speed
wL_c = sin(0.002*tc); % Left wheel angular speed
sl_c = length(tc); % Number of steps in CT simulation
xc = zeros(1,sl_c); % Vector to store CT x(t)
yc = zeros(1,sl_c); % Vector to store CT y(t)
phc = zeros(1,sl_c); % % Vector to store CT \phi(t)

for ii=2:sl_c-1 % <= one less simulation for the reason mentioned above
vR_c=r*wR_c(ii);
vL_c=r*wL_c(ii);
v_c=(vR_c+vL_c)/2;
w_c=(vR_c-vL_c)/L;
dph_c=dt*w_c;
xc(ii+1)=xc(ii)+v_c*dt*cos(phc(ii));
yc(ii+1)=yc(ii)+v_c*dt*sin(phc(ii));
phc(ii+1)=phc(ii)+dph_c;
end

i = 0;
for Ts = 1:7 %Sampling time has integer values: 1,2,...,6,7
    i = i+1; %Compute trajectories in discrete-time for six models
    
    T = Ts; % Sampling time
    ti = 0; % Initial time
    tf = 100; % Final time
    td = ti:T:tf; % Discrete time defined
    tf = td(end);
    sl = length(td); % Number of samples in discrete time frame
    wR = 0.3*cos(0.002*td); % Right wheel angular speed
    wL = sin(0.002*td); % Left wheel angular speed
    r = 5; % Radius of wheels [cm]
    L = 10; % Width of robot platform [cm]
    xd = zeros(6,sl); % Vector to store global x_n
    yd = zeros(6,sl); % Vector to store global y_n
    
    phd = zeros(1,sl); % Vector to store global \phi_n
    phd_exact = zeros(1,sl); % Vector to store global \phi_n for exact method
    
    for ii=2:sl-1
        % note that we need to run this loop for one less iteration because
        % the pose quantities have time index (ii+1)
        vR=r*wR(ii);
        vL=r*wL(ii);
        v=(vR+vL)/2;
        w=(vR-vL)/L;
        dph=T*w;
        xd(1,ii+1)=xd(1,ii)+v*T*cos(phd(ii)); % Rectangular model
        yd(1,ii+1)=yd(1,ii)+v*T*sin(phd(ii));

        xd(2,ii+1)=xd(2,ii)+(v*T./2)*(cos(phd(ii)) + cos(phd(ii) + T*w)); % Trapezoidal
        yd(2,ii+1)=yd(2,ii)+(v*T./2)*(sin(phd(ii)) + sin(phd(ii) + T*w));

        xd(3,ii)=xd(3,ii-1)+(v./w)*(sin(phd_exact(ii-1)+dph) - sin(phd_exact(ii-1))); % Exact
        yd(3,ii)=yd(3,ii-1)-(v./w)*(cos(phd_exact(ii-1)+dph) - cos(phd_exact(ii-1)));
        phd_exact(ii) = phd_exact(ii-1)+dph;

        xd(4,ii+1)=xd(4,ii)+v*T*(cos(phd(ii)+dph/2) ); % Geometry based I
        yd(4,ii+1)=yd(4,ii)+v*T*(sin(phd(ii)+dph/2) );

        xd(5,ii)=xd(5,ii-1)+((v*dph)./w)*(cos(phd_exact(ii-1)) - 0.5*sin(phd_exact(ii-1))*dph); % Taylor 
        yd(5,ii)=yd(5,ii-1)+((v*dph)./w)*(sin(phd_exact(ii-1)) + 0.5*cos(phd_exact(ii-1))*dph);
        phd_exact(ii) = phd_exact(ii-1)+dph;

        xd(6,ii+1) = xd(6,ii) + ((T*v)/6)*( cos(phd(ii)) + 2*cos(phd(ii)+dph/2) ... % RK4
            + 2*cos(phd(ii)+dph) + 2*cos(phd(ii)+(3*dph)/2) ...
            + cos(phd(ii)+4*dph) );
        yd(6,ii+1) = yd(6,ii) + ((T*v)/6)*( sin(phd(ii)) + 2*sin(phd(ii)+dph/2) ...
            + 2*sin(phd(ii)+dph) + 2*sin(phd(ii)+(3*dph)/2) ...
            + sin(phd(ii)+4*dph) );

        phd(ii+1)=phd(ii)+dph;
    end
    
    %Compare each trajectory above with that of continuous-time
    for j = 1:6
        error_x = xd(j,:) - xc(round(linspace(1,length(xc),length(xd(j,:)))));
        error_y = yd(j,:) - yc(round(linspace(1,length(yc),length(yd(j,:)))));
        sum_squared_error(j,i)= sum(error_x.^2)+sum(error_y.^2);
    end
end
all_Ts = [1,2,3,4,5,6,7];
for j =1:6
    hold on
    plot(all_Ts,sum_squared_error(j,:),'linewidth',2)
end
legend('Rectangular','Trapezoidal','Éxact','Geometry-Based 1','Taylor','RK4',...
    'fontsize',16,'location','best')
title('\Sigma error^2 of each model versus varying sampling time','fontsize', 18);
xlabel('Ts','fontsize', 16); ylabel('sum squared error','fontsize', 16);



% gravdragdemo.m
% Andy Ganse, Applied Physics Laboratory, Univ of WA, 2006.
% aganse@apl.washington.edu, http://staff.washington.edu/aganse
%
% A demo of basic nonlinear tracking examples 6.1-2 & 6.1-3 from the
% classic text, "Applied Optimal Estimation", editted by Gelb, 1974,
% for plotting and comparing the following recursive filters & smoothers:
% - Kalman Filter (KF),
% - Kalman Smoother (KS),
% - Linearized Kalman Filter (LKF),
% - Extended Kalman Filter (EKF), 
% - 2nd order EKF (EKF2)
%
% The user can change settings at the beginning of this script to choose
% which filters to run/display, set initial conds and noise values, etc.
% Please note the book's point about using Monte Carlo simulation to get
% the real cov estimates for EKF, and EKF2 - I don't do that here.

function gravdragdemo()
global nonlineardyn;
global nonlinearmeas;
global Qon;

% To understand this script note its five labeled sections:
% - USER CHANGEABLE STUFF,
% - PROBLEM SETUP,
% - FILTERS & SMOOTHERS,
% - PLOTTING,
% - FUNCTION DEFINITIONS.
% But the user can initially just concentrate on making changes in the USER
% CHANGEABLE STUFF section and then running the script to see their effects.


% -----------------------------------------
% USER CHANGEABLE STUFF...

nonlineardyn=0;  % problem dynamics: 1=nonlinear (both Gelb examples), 0=linear (to try KF)
nonlinearmeas=0;  % problem measurements: 1=nonlinear (Gelb ex6.1-3), 0=linear (Gelb ex6.1-2)

% Choose which filters to try:  1="on" and 0="off" on all these...
runKF=1;  % for when nonlineardyn=nonlinearmeas=0, or just to show it's bad otherwise
runKS=1;  % "      "       "
runLKF=0;
runEKF=0;
runEKF2=0;

% Initial conds and setup is the same for both 1D & 2D problems.
% (these default values here are taken from book):
x0true=[1e5+100,-6000+350,2500+500]';  % = elevation(ft), velocity(ft/sec), drag coeff
x0guess=[1e5,-6000,2500]';  % note this will cause LKF's xbar to cross zero a little early...
P0=diag([500,2e4,2.5e5]);
Qon=0;  % multiply this scalar times the Q matrix to "turn it on&off" (ie 1=on,0=off)
        % Gelb examples both have Q=0; you can change Q down in getQ() function, default is 1e2*I)
R=1e2;  % note meas noise cov matrix is a single scalar variance in this problem

% Note all filters in this script assume equally-spaced time points:
T=.01;  % Timestep for discretization of continuous functions of time.
        % While it's not specified in the book, the Gelb book examples very
        % likely used a timestep of order 0.01 and then only plotted every
        % 10th point or so.
        % Note if you want to use a larger T value then you may not want to
        % use the linearization approach to solve the nonlinear differential
        % equations below for xtrue(t) and xbar(t) (the latter in the LKF).
timerange=0:T:16;  % or might want to go till 17.7;

% narrow the plot boundaries?  (1=yes, 0=no)
axislimits=1;
newfigs=0;

% END OF USER CHANGEABLE STUFF.


% -----------------------------------------
% PROBLEM SETUP...
N=length(timerange);

% Calculate "true" state vector over time from dynamics.
% Note true model includes random element to dynamics (via "fplusw" func).
%    more accurate for nonlinear f :
%    [t,xcalc]=ode113(@fplusw,timerange,x0true);  % note t starts at zero
%    xtrue=xcalc';  % (because darned ODE routine transposes result)
% Getting away with linearization approach here with small enough T and
% weakly nonlinear function f:
xtrue=zeros(3,N);
xtrue(:,1)=x0true;
for i=2:N
    xtrue(:,i)=xtrue(:,i-1)+fplusw(xtrue(:,i-1))*T;
end;

% Create noisy scalar measurements:
%   For 2D problem (ex. 6.1-3), scalar nonlinear measurements are of form
%   z_k = h(x_k) + v_k = sqrt( r1^2 + ( x_k - r2 )^2 ) + v_k ,  v_k~N(0,R)
%   Note 1D problem (ex. 6.1-2), scalar linear measurements, is a special
%     case of 2D problem where r1=r2=0, becoming:
%   z(t) = h(x(t)) + v(t) = x1(t) + v(t) ,  v(t)~N(0,R)
randn('state',sum(100*clock));  % seed the rand # generator
z=h(xtrue)+sqrt(R)*randn(1,length(xtrue));  % add noise to scalar meas
% note here that first measurement was at t=0...

% We know the Q's in this synthetic problem because I chose them in true
% model.
Q=getQ();



% -----------------------------------------
% FILTERS & SMOOTHERS...

if runKF  % KF:  Kalman Filter
% use constant F=df/dx & H=dh/dx eval'd at x0:
F=jacf(x0guess);
H=jach(x0guess);
% Bring in first measurement at t=0:
K = P0*H'/(H*P0*H'+R);
P_kf(:,:,1) = (eye(size(P0))-K*H)*P0;
x_kf(:,1) = x0guess + K*(z(1)-H*x0guess);
% Now loop, iterating over dynamics-propagate / meas-update...
for i=2:N
    % Propagate to next measurement time:
     Pminus = P_kf(:,:,i-1) + (F*P_kf(:,:,i-1)+P_kf(:,:,i-1)*F'+Q)*T;
     xminus = x_kf(:,i-1) + F*x_kf(:,i-1)*T;
    % Measurement update:
    K = Pminus*H'/(H*Pminus*H'+R);
    P_kf(:,:,i) = (eye(size(Pminus))-K*H)*Pminus;
    x_kf(:,i) = xminus + K*(z(i)-H*xminus);
    drawnow; % lets us interrupt with Ctrl-C
end
end


if runKS  % KS:  Kalman Smoother, coded ala Gelb table 5.2-2, p164.
% Continuous time RTS fixed-interval smoother.
% use constant F=df/d x & H=dh/dx eval'd at x0:
F=jacf(x0guess);
H=jach(x0guess);
% Forward filter:
% bring in first measurement at t=0
K = P0*H'/(H*P0*H'+R);
P(:,:,1) = (eye(size(P0))-K*H)*P0;
x(:,1) = x0guess + K*(z(1)-H*x0guess);
for i=2:N
    % Propagate to next measurement time:
     Pminus = P(:,:,i-1) + (F*P(:,:,i-1)+P(:,:,i-1)*F'+Q)*T;
     xminus = x(:,i-1) + F*x(:,i-1)*T;
    % Measurement update:
    K = Pminus*H'/(H*Pminus*H'+R);
    P(:,:,i) = (eye(size(Pminus))-K*H)*Pminus;
    x(:,i) = xminus + K*(z(i)-H*xminus);
    drawnow; % lets us interrupt with Ctrl-C
end
% Smoothing (using fwd filter results):
P_ks(:,:,N)=P(:,:,N);
x_ks(:,N)=x(:,N);
for i=N-1:-1:1
    x_ks(:,i) = x_ks(:,i+1) - ...
        (F*x_ks(:,i+1)+Q*inv(P(:,:,i+1))*(x_ks(:,i+1)-x(:,i+1)))*T;
    P_ks(:,:,i) = ...
        P_ks(:,:,i+1) - ( (F+Q*inv(P(:,:,i+1)))*P_ks(:,:,i+1) + ...
        P_ks(:,:,i+1)*(F+Q*inv(P(:,:,i+1)))' - Q )*T;
    drawnow; % lets us interrupt with Ctrl-C
end
end


if runLKF   % LKF:  Linearized Kalman Filter.
% First we precompute xbar for prior linearization:
%    more accurate for nonlinear f :
%    [t,xcalc]=ode113(@fplusw,timerange,x0guess);  % note t starts at zero
%    xbar=xcalc';  % (because darned ODE routine transposes result)
% Getting away with linearization approach here with small enough T and
% weakly nonlinear function f:
xbar=zeros(3,N);
xbar(:,1)=x0guess;
for i=2:N
    xbar(:,i)=xbar(:,i-1)+fplusw(xbar(:,i-1))*T;
end;
% Bring in first measurement at t=0:
H=jach(x0guess);  % (note xbar(:,1)=x0guess so could write either way)
K = P0*H'*inv(H*P0*H'+R);
P_lkf(:,:,1) = (eye(size(P0))-K*H)*P0;
x_lkf(:,1) = x0guess + K*(z(1)-h(x0guess));  % -H*() term is zero since xbar(:,1)=x0guess
% Now loop, iterating over dynamics-propagate / meas-update...
for i=2:N
    % Propagate to next measurement time:
    F=jacf(xbar(:,i-1));
    Pminus = P_lkf(:,:,i-1) + (F*P_lkf(:,:,i-1)+P_lkf(:,:,i-1)*F'+Q)*T;
    xminus = x_lkf(:,i-1) + (f(xbar(:,i-1))+F*(x_lkf(:,i-1)-xbar(:,i-1)))*T;
    % Measurement update:
    H=jach(xbar(:,i));
    K = Pminus*H'*inv(H*Pminus*H'+R);
    P_lkf(:,:,i) = (eye(size(Pminus))-K*H)*Pminus;
    x_lkf(:,i) = xminus + K*(z(i)-h(xbar(:,i))-H*(xminus-xbar(:,i)));
    % Note I have H= and K= and P= inside the propagate/update loop here
    % unnecessarily (just for easier comparison to other filters' code) -
    % but notice that the sequences of H's and K's and P's can all be done
    % a-priori since we have the xbar sequence a-priori.  So the actual
    % "in-flight" computation loop is just the x's, greatly reducing 
    % computation time.
    drawnow; % lets us interrupt with Ctrl-C
end
end


if runEKF  % EKF:  Extended Kalman Filter
% Bring in first measurement at t=0:
H=jach(x0guess);
K = P0*H'*inv(H*P0*H'+R);
P_ekf(:,:,1) = (eye(size(P0))-K*H)*P0;
x_ekf(:,1) = x0guess + K*(z(1)-h(x0guess));
% Now loop, iterating over dynamics-propagate / meas-update...
for i=2:N
    % Propagate to next measurement time:
    F=jacf(x_ekf(:,i-1));
    Pminus = P_ekf(:,:,i-1) + (F*P_ekf(:,:,i-1)+P_ekf(:,:,i-1)*F'+Q)*T;
    xminus = x_ekf(:,i-1) + f(x_ekf(:,i-1))*T;
    % Measurement update:
    H=jach(xminus);
    K = Pminus*H'*inv(H*Pminus*H'+R);
    P_ekf(:,:,i) = (eye(size(Pminus))-K*H)*Pminus;
    x_ekf(:,i) = xminus + K*(z(i)-h(xminus));
    drawnow; % lets us interrupt with Ctrl-C
end
end


if runEKF2  % EKF2:  Gaussian 2nd order EKF
            % (the "Gaussian" part refers to the A matrix, ala Gelb 1974 p192)
% Bring in first measurement at t=0:
H=jach(x0guess);
A = calcA(x0guess,P0);
K = P0*H'*inv(H*P0*H'+R+A);
P_ekf2(:,:,1) = (eye(size(P0))-K*H)*P0;
x_ekf2(:,1) = x0guess + K*(z(1)-h(x0guess)-1/2*d2h(x0guess,P0));
% Now loop, iterating over dynamics-propagate / meas-update...
for i=2:N
    % Propagate to next measurement time:
    F=jacf(x_ekf2(:,i-1));
    Pminus = P_ekf2(:,:,i-1) + (F*P_ekf2(:,:,i-1)+P_ekf2(:,:,i-1)*F'+Q)*T;
    xminus = x_ekf2(:,i-1) + (f(x_ekf2(:,i-1))+1/2*d2f(x_ekf2(:,i-1),P_ekf2(:,:,i-1)))*T;
    % Measurement update:
    H=jach(xminus);
    % using Gaussian approximation to calc A, so this is a Gaussian 2nd order filter:
    A = calcA(xminus,Pminus);  % note A is a scalar in this problem since h and R scalar
    K = Pminus*H'*inv(H*Pminus*H'+R+A);
    P_ekf2(:,:,i) = (eye(size(Pminus))-K*H)*Pminus;
    x_ekf2(:,i) = xminus + K*(z(i)-h(xminus)-1/2*d2h(xminus,Pminus));
    drawnow; % lets us interrupt with Ctrl-C
end
end



% -----------------------------------------
% PLOTTING...

% Set corresponding plot title:
tstr='';
if runKF; tstr=strcat(tstr,' KF=magenta,'); end;
if runKS; tstr=strcat(tstr,' KS=green,'); end;
if runLKF; tstr=strcat(tstr,' LKF=cyan, '); end;
if runEKF; tstr=strcat(tstr,' EKF=blue, '); end;
if runEKF2; tstr=strcat(tstr,' EKF2=green,'); end;
tstr=tstr(1:end-1);  % <--  quick hack to remove final comma

% Initialize plot:
if newfigs
    figure;
else
    figure(1);
end
subplot(1,3,1); hold off;
subplot(1,3,2); hold off;
subplot(1,3,3); hold off;
t=timerange;  % just for brevity in plot lines below

% Note below that while we have N+1 time-points in the state estimates (one
% prior value plus N measurements), we skip the last time-point of the state
% estimate below since we subject the N-point "true" state sequence to look at
% estimatation errors.

% Plot KF results:
if runKF
s=[squeeze(sqrt(abs(P_kf(1,1,:)))),squeeze(sqrt(abs(P_kf(2,2,:)))),squeeze(sqrt(abs(P_kf(3,3,:))))];
subplot(1,3,1);
plot(t',x_kf(1,:)-xtrue(1,:),'m-',t',s(:,1),'m--',t',-s(:,1),'m--','linewidth',2)
hold on;
subplot(1,3,2);
plot(t',x_kf(2,:)-xtrue(2,:),'m-',t',s(:,2),'m--',t',-s(:,2),'m--','linewidth',2)
hold on;
subplot(1,3,3);
plot(t',x_kf(3,:)-xtrue(3,:),'m-',t',s(:,3),'m--',t',-s(:,3),'m--','linewidth',2)
hold on;
end

% Plot KS results:
if runKS
s=[squeeze(sqrt(abs(P_ks(1,1,:)))),squeeze(sqrt(abs(P_ks(2,2,:)))),squeeze(sqrt(abs(P_ks(3,3,:))))];
subplot(1,3,1);
plot(t',x_ks(1,:)-xtrue(1,:),'g-',t',s(:,1),'g--',t',-s(:,1),'g--','linewidth',2)
hold on;
subplot(1,3,2);
plot(t',x_ks(2,:)-xtrue(2,:),'g',t',s(:,2),'g--',t',-s(:,2),'g--','linewidth',2)
hold on;
subplot(1,3,3);
plot(t',x_ks(3,:)-xtrue(3,:),'g',t',s(:,3),'g--',t',-s(:,3),'g--','linewidth',2)
hold on;
end

% Plot LKF results:
if runLKF
s=[squeeze(sqrt(abs(P_lkf(1,1,:)))),squeeze(sqrt(abs(P_lkf(2,2,:)))),squeeze(sqrt(abs(P_lkf(3,3,:))))];
subplot(1,3,1);
plot(t',x_lkf(1,:)-xtrue(1,:),'c--',t',s(:,1),'c--',t',-s(:,1),'c--','linewidth',2)
hold on;
subplot(1,3,2);
plot(t',x_lkf(2,:)-xtrue(2,:),'c--',t',s(:,2),'c--',t',-s(:,2),'c--','linewidth',2)
hold on;
subplot(1,3,3);
plot(t',x_lkf(3,:)-xtrue(3,:),'c--',t',s(:,3),'c--',t',-s(:,3),'c--','linewidth',2)
hold on;
end

% Plot EKF results:
if runEKF
s=[squeeze(sqrt(abs(P_ekf(1,1,:)))),squeeze(sqrt(abs(P_ekf(2,2,:)))),squeeze(sqrt(abs(P_ekf(3,3,:))))];
subplot(1,3,1);
plot(t',x_ekf(1,:)-xtrue(1,:),'b-.',t',s(:,1),'b--',t',-s(:,1),'b--','linewidth',2)
hold on;
subplot(1,3,2);
plot(t',x_ekf(2,:)-xtrue(2,:),'b-.',t',s(:,2),'b--',t',-s(:,2),'b--','linewidth',2)
hold on;
subplot(1,3,3);
plot(t',x_ekf(3,:)-xtrue(3,:),'b-.',t',s(:,3),'b--',t',-s(:,3),'b--','linewidth',2)
hold on;
end

% Plot Gaussian 2nd order EKF results:
if runEKF2
s=[squeeze(sqrt(abs(P_ekf2(1,1,:)))),squeeze(sqrt(abs(P_ekf2(2,2,:)))),squeeze(sqrt(abs(P_ekf2(3,3,:))))];
subplot(1,3,1);
plot(t',x_ekf2(1,:)-xtrue(1,:),'g',t',s(:,1),'g--',t',-s(:,1),'g--','linewidth',2)
hold on;
subplot(1,3,2);
plot(t',x_ekf2(2,:)-xtrue(2,:),'g',t',s(:,2),'g--',t',-s(:,2),'g--','linewidth',2)
hold on;
subplot(1,3,3);
plot(t',x_ekf2(3,:)-xtrue(3,:),'g',t',s(:,3),'g--',t',-s(:,3),'g--','linewidth',2)
hold on;
end


% add labels, axes limits, etc...
subplot(1,3,1);
if axislimits
   axis([min(t) max(t) -5 5]);
end
ylabel('elev estim error (ft)','fontsize',12); 
xlabel('time (s)','fontsize',12);
grid on;
subplot(1,3,2);
if axislimits
   axis([min(t) max(t) -5 5]);
end
ylabel('vel estim error (ft/s)','fontsize',12); 
xlabel('time (s)','fontsize',12);
grid on;
title(tstr,'fontsize',16,'fontweight','bold');
subplot(1,3,3);
if axislimits
   axis([min(t) max(t) -600 600]);
end
ylabel('\beta estim error (lb/ft^2)','fontsize',12); 
xlabel('time (s)','fontsize',12);
grid on;


% Plot the true, measured, & estimated elevations over time
if newfigs
    figure;
else
    figure(2);
end
zlin=x0true(2)*t+x0true(1);
subplot(2,1,1); hold off;
plot(t,z,'k:',t,xtrue(1,:),'k-',t,zlin,'k--');  %,t,xbar(1,:),'k-.');
if axislimits; axis([1,max(t),-3e4,1e5]); end;
xlabel('time (s)'); ylabel('elevation (ft)');
title('Meas=dottedBlack, TruePos=solidBlack, x0line=dashedBlack');
if runKF; hold on; plot(t,x_kf(1,:),'y'); end;
if runKS; hold on; plot(t,x_ks(1,:),'g'); end;
if runLKF; hold on; plot(t,x_lkf(1,:),'c'); end;
if runEKF; hold on; plot(t,x_ekf(1,:),'m'); end;
if runEKF2; hold on; plot(t,x_ekf2(1,:),'g'); end;
subplot(2,1,2); hold off;
plot(t,z-zlin,'k:',t,xtrue(1,:)-zlin,'k-');  %,t,xbar(1,:)-zlin,'k-.');
if runKF; hold on; plot(t,x_kf(1,:)-zlin,'y'); end;
if runKS; hold on; plot(t,x_ks(1,:)-zlin,'g'); end;
if runLKF; hold on; plot(t,x_lkf(1,:)-zlin,'c'); end;
if runEKF; hold on; plot(t,x_ekf(1,:)-zlin,'m'); end;
if runEKF2; hold on; plot(t,x_ekf2(1,:)-zlin,'g'); end;
xlabel('time (s)'); ylabel('elevation - x0line');
if axislimits; axis([1,max(t),min(xtrue(1,:)-zlin)-150,150]); end;



% -----------------------------------------
% FUNCTION DEFINITIONS...

function out=h(x)
% measurement function
[r1,r2]=getgeom();
out=sqrt(r1^2+(x(1,:)-r2).^2);
% note if r1 and r2 are zero then this collapses to h(x)=x(1)

function out=jach(x)
% local jacobian of measurement function h(x)
[r1,r2]=getgeom();
out=[(x(1)-r2)/sqrt(r1+(x(1)-r2)^2),0,0];

function out=hessh(x)
% local hessian of measurement function h(x).
% used in EKF2.
% handily for this problem the output of h() is scalar, so hessh is one matrix
[r1,r2]=getgeom();
out(1,:)=[-(x(1)-r2)^2/(r1+(x(1)-r2)^2)^(3/2)+(r1+(x(1)-r2)^2)^(-1/2),0,0];%%%%
out(2,:)=[0,0,0];
out(3,:)=[0,0,0];

function out=d2h(x,P)
% used in EKF2.
out=trace(hessh(x)*P);  % (scalar for this problem since h and z are scalars)

function out=calcA(x,P)  % used in Gaussian 2nd order EKF
% Note out (=A) is a scalar in this problem since h and R scalar.
% Note also as per Gelb that this computation for A relies on a Gaussian
% approximation as described on p192 of Gelb.
hess=hessh(x);
out=0;
for p=1:3
    for q=1:3
        for m=1:3
            for n=1:3
                out=out+1/4*hess(p,q)*(P(p,m)*P(q,n)+P(p,n)*P(q,m))*hess(m,n);
            end
        end
    end
end

function out=hessf1(x)
% component 1 of hessian of dynamics function f(x)
out=zeros(3);

function out=hessf2(x)
% component 2 of hessian of dynamics function f(x)
rho0 = getrho0();
krho = 22000;
% formed from these which came from jacf() :
% df2/dx1 = -rho0/krho*exp(-x(1)/krho)*(x(2))^2/(2*x(3))
% df2/dx2 = rho0*exp(-x(1)/krho)*x(2)/x(3)
% df2/dx3 = -rho0*exp(-x(1)/krho)*(x(2))^2/(2*(x(3))^2)
out(1,:)=[ rho0/krho^2*exp(-x(1)/krho)*(x(2))^2/(2*x(3)), ...
           -rho0/krho*exp(-x(1)/krho)*x(2)/x(3), ...
           rho0/krho*exp(-x(1)/krho)*(x(2))^2/(2*(x(3))^2) ];
           % = [d2f2/dx1dx1, d2f2/dx1dx2, d2f2/dx1dx3]
out(2,:)=[ out(1,2), ...
           rho0/x(3)*exp(-x(1)/krho), ...
           -rho0/(x(3))^2*exp(-x(1)/krho)*x(2) ];
           % = [out(1,2), d2f2/dx2dx2, d2f2/dx2dx3]
out(3,:)=[ out(1,3), ...
           out(2,3), ...
           rho0*exp(-x(1)/krho)*(x(2))^2/(x(3))^3 ];
           % = [out(1,3), out(2,3), d2f2/dx3dx3]
% (note Hessian is symmetric so some elements here are multiply assigned)

function out=hessf3(x)
% component 3 of hessian of dynamics function f(x)
out=zeros(3);

function out=d2f(x,P)
% used in EKF2.
out=zeros(3,1);
out(1)=trace(hessf1(x)*P);
out(2)=trace(hessf2(x)*P);
out(3)=trace(hessf3(x)*P);

function out=fplusw(x)
% add on samples of the plant noise (theory error) to dynamics function f(x)
Q=getQ();
if Q~=zeros(3)
   L=chol(Q);  % since Q is a covariance we can use Cholesky decomp
   w=L*randn(3,size(x,2));
else
   w=zeros(3,size(x,2));
end;
out=f(x)+w;
   
function out=jacf(x)
% local jacobian of dynamics function f(x)
rho0 = getrho0();
krho = 22000;
out=zeros(3,3);
out(1,:)=[0,1,0];  % = df1/dx_{1-3}
out(2,:)=[-rho0/krho*exp(-x(1)/krho)*(x(2))^2/(2*x(3)),...
          rho0*exp(-x(1)/krho)*x(2)/x(3),...
          -rho0*exp(-x(1)/krho)*(x(2))^2/(2*(x(3))^2)];  % = df2/dx_{1-3}
out(3,:)=[0,0,0];  % = df3/dx_{1-3}

function xp=f(x)  % same nonlinear dynamics used in both Gelb examples 6.1-2 & 6.1-3
global nonlineardyn;
g = 32.2*nonlineardyn;  % nonlineardyn (=1 or 0) can zero out g and rho0 below,
                        % to make dynamics collapse to linear if desired.   
rho0 = getrho0();                        
krho = 22000;
xp=zeros(3,size(x,2));
xp(1,:) = x(2,:);
xp(2,:) = rho0*exp(-x(1,:)/krho).*(x(2,:)).^2./(2*x(3,:)) - g;
xp(3,:) = zeros(1,size(x,2));
% note the "nonlineardyn" factor (=1 or 0) allows us to zero g out here,
% making dynamics collapse to linear if desired.
% "nonlineardyn" also zeros out rho0 below.

function rho0=getrho0()
% air denisty.
global nonlineardyn;
rho0=3.4e-3*nonlineardyn;
% note the "nonlineardyn" factor (=1 or 0) allows us to zero this out from
% top of script, making dynamics collapse to linear if desired.
% "nonlineardyn" also zeros out the g up in f(x).

function Q=getQ()
% theory error covariance matrix Q (ie cov of the additive noise on f(x) )
global Qon;
%Q=Qon*diag([1e2 1e2 1e2]);
Q=Qon*diag([1e-2 1e-2 1e-2]);
% Note the "Qon" factor (=1 or 0) allows us to zero out the Q matrix from
% the top of the script, setting to no theory error as in Gelb examples.
% If Qon isn't zero, Q is specifying covariance of the theory error here.
% In this demo script we use Q to produce the "true" trajectory and then 
% use the same Q in the filtering since we know exactly what it is.
% However, for real life problems we must adaptively filter to estimate Q
% as we go; but that subject is outside the scope of this demonstation.

function [r1,r2]=getgeom()
% just a way to pass r1 & r2 to a bunch of other functions above.
global nonlinearmeas;
r1=1e3*nonlinearmeas;  % sensor horiz range from vertical target path.
r2=1e2*nonlinearmeas;  % sensor elevation                             
% note the "nonlinearmeas" factor (=1 or 0) allows us to zero these out
% from top of script, making the measurements collapse to linear if desired.

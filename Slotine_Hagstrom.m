clear all; %clear all variables 
close all; %close all open figures
%%%%%%%%%%%%%%%%%%%%%%%%%  GRID   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%N = 20; M = 20; % grid size
distFunc = ['euclidean']; %type of distance between nodes
%[X Y] = meshgrid(1:N,1:M);
%X = X(:); Y = Y(:);
%%%%%%%%%%%%%%%%%%%%%%%%%  GREY-SCALE-DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=20;
A=zeros(N);
%A(1:5,1:10)=0;
%A(6:15,1:5)=0;
%A(16:20,1:10)=0;
A(6:15,6:15)=255;
A(1:5,11:20)=127;
A(16:20,11:20)=127;
A(6:15,16:20)=127;
A=A+50.*(rand(size(A))-rand(size(A))); %noise
%A = imnoise(A,'gaussian');
%u=A(1:end)';
%minNewRange=0;
%maxNewRange=1;
%originalValue=u;
%minOriginalRange=min(u);
%maxOriginalRange=max(u);
%A = minNewRange + (((maxNewRange - minNewRange) * (originalValue - minOriginalRange))/(maxOriginalRange - minOriginalRange)); %normalize pixel potentials
A=reshape(A,[20,20]);
u=A(1:end)';
umin=min(u);
umax=max(u);

Imin=0.8;
Imax=2;
beta=20;
omega=5;

I=(u-umin).*((Imax-Imin)/(umax-umin))+Imin;
 %noise
dt = 0.01; %time step [ms]
t_end = 50; %total time of run [ms]

%%%%%%%%%%%%%%%%%%%%%%%   VIEW OF GREY-SCALE DATA
%%%%%%%%%%%%%%%%%%%%%%%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imshow(A, [0 1],'InitialMagnification','fit'); %plot grey-scale data
%%%%%%%%%%%%%%%%%%%%%%%% DISTANCE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xstart = 1;
ystart = 1;
xend = size(A,1);
yend = size(A,2);
N = size(A,1);
M = size(A,2);
Npts = N*M;
dist = zeros(Npts); %creates empty output matrix to be updated
xrange = linspace(xstart,xend,N);
yrange = linspace(ystart,yend,M);
[X, Y] = meshgrid(xrange, yrange);
X=X(:);
Y=Y(:);
pts = [X(:), Y(:)];
%distances = pdist2(pts, pts,'euclidean'); %Euclidean distance between points in grid.

%%%%%%%%1%%%%%%%%%%%%%%%%%  DEFINE PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha=12;
c=0.04;
rho=4;
h=dt;
%%%%%%%%%%%%%%%%%%%% CONNECTIVITY %%%%%%%%%%%%%%%%%%%%%%%%%%%

U=abs(u-u');
k=exp(-(U).^2/beta^2);

for i=1:400
    for j=1:400
        if abs(i-j)<omega
        k(i,j)=k(i,j);
        else
        k(i,j)=0;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%  DEFINE INITIAL VALUES AND VECTORS TO HOLD RESULTS %%%%%%%
t_vect = 0:dt:t_end; %will hold vector of times
filename='oscillators';
N = size(A,1);
M = size(A,2);
Npts = N*M;
v=zeros(N*M,t_end./dt);
w=zeros(M*N,t_end./dt);



%********************************************************************%
%Euler's method and Second Order Runge-Kutta
%********************************************************************%
%********************************************************************%
 v(:,1)=0.002*(rand(1,N*M)-rand(1,N*M))';
 w(:,1)=0.002*(rand(1,N*M)-rand(1,N*M))';

for i=1:t_end/dt
   
        
      v(:,i+1)=v(:,i)+h*dv_dt(v(:,i),w(:,i),I,k);
      w(:,i+1)=w(:,i)+h*dw_dt(v(:,i),w(:,i),c,alpha,rho,k);
  
   
     figure(2)
     surf(reshape(w(:,i),[N,M]))
    % view(0,90);
     zlim([-100 100])
end
%%

%********************************************************************%
%Functions to solve
%********************************************************************%
function[slope]=dv_dt(v,w,I,k)
slope=3.*v-v.^3-v.^7+2-w+I;
vv=v(:,ones(size(v,1),1));
slope=slope+sum(k.*(vv'-vv),2);
end

function[slope]=dw_dt(v,w,c,alpha,rho,k)
slope=c.*(alpha.*(tanh(rho*v))-w);
ww=w(:,ones(size(w,1),1));
slope=slope+sum(k.*(ww'-ww),2);
end









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
A=A+20.*(rand(size(A))-rand(size(A))); %noise
u=A(1:end)';
minNewRange=0;
maxNewRange=1;
originalValue=u;
minOriginalRange=min(u);
maxOriginalRange=max(u);
A = minNewRange + (((maxNewRange - minNewRange) * (originalValue - minOriginalRange))/(maxOriginalRange - minOriginalRange)); %normalize pixel potentials
u=A(1:end)';



 %noise
dt = 0.1; %time step [ms]
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




%%%%%%%%%%%%%%%%%%%%%%% 4-NEAREST NODE CONNECTIVITY  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
adj = squareform(pdist([X Y], distFunc) == 1 );
for a=1:N
  for b=1:M
      if adj(a,b)==true
         W(a,b)=1;
      else
         W(a,b)=0;
      end
   end
end
[i,j]=find(adj); %returns the row and column subscripts of each nonzero element in Adjacency matrix.
conn_result = accumarray(i,j,[size(adj,1), 1], @(x){sort(x).'}); %applies the function 'sort' to each group in data (j) specified by ind (i). Shows what nodes each node is connected to.

%Count the number of  of each node
Neighbors= cellfun(@numel, conn_result);
Neighbors=repmat(Neighbors,1,t_end./dt); 



%%%%%%%%1%%%%%%%%%%%%%%%%%  DEFINE PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
epsilon=0.02;
psi=3.0;
gamma=6.0;
beta=.10;
kappa=50;
thetax=-0.5;
thetazx=0.1;
thetaxz=0.2;
rho=0.02;
I=0.2;
Wik=2.5;
Wz=0.1;
%%%%%%%%%%%%%%%%%%%%%%%  DEFINE INITIAL VALUES AND VECTORS TO HOLD RESULTS %%%%%%%
t_vect = 0:dt:t_end; %will hold vector of times


filename='oscillators';
N = size(A,1);
M = size(A,2);
Npts = N*M;

v=zeros(N*M,t_end./dt);
w=zeros(M*N,t_end./dt);
v(:,1)=u;
w(:,1)=rand(400,1);
z(:,1)=rand(400,1);


h=dt;

%********************************************************************%
%Euler's method and Second Order Runge-Kutta
%********************************************************************%
%********************************************************************%

for i=1:t_end/dt
    for j=1:400
      Si(j,i)=sum(Wik*1./(1+exp(-kappa*v(Neighbors([conn_result{conn_result{j}}]),i)-thetax))-Wz*1./(1+exp(-kappa*z(Neighbors([conn_result{conn_result{j}}]),i)-thetazx)));
      v(j,i+1)=v(j,i)+h*dv_dt(v(j,i),w(j,i),I,Si(j,i),rho);
      w(j,i+1)=w(j,i)+h*dw_dt(v(j,i),w(j,i),beta,gamma,epsilon);
      z(j,i+1)=w(j,i)+h*dzdt(v(:,i),thetaxz,z(j,i),psi);
   
     figure(2)
     surf(reshape(v(:,i),[20,20]))
     view(0,90);

   end

end
%%

%********************************************************************%
%Functions to solve
%********************************************************************%
function[slope]=dv_dt(v,w,I,Si,rho)
slope=3.*v-v.^3+2-w+rho+I'+Si;
end


function[slope]=dw_dt(v,w,beta,gamma,epsilon)
slope=epsilon.*(gamma.*(tanh(v./beta))-w);
end


function result = dzdt(x,thetaxz,z,psi)
if all(x < thetaxz)
   sigmainf=0;
end
if any(x>=thetaxz)
    sigmainf=1;
end
    result=psi*(sigmainf-z);
end






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

A=A+105.*(rand(size(A))-rand(size(A))); %noise
%A = imnoise(A,'gaussian');
u=A(1:end);
black=find(u==0);
gray=find(u==127);
white=find(u==255);
minNewRange=0;
maxNewRange=1;
originalValue=u;
minOriginalRange=min(u);
maxOriginalRange=max(u);
A = minNewRange + (((maxNewRange - minNewRange) * (originalValue - minOriginalRange))/(maxOriginalRange - minOriginalRange)); %normalize pixel potentials
A=reshape(A,[20,20]);
 %noise
dt = 0.01; %time step [ms]
t_end = 10; %total time of run [ms]

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
for k=1:N*M
conn_result3{k}=[k,conn_result{k}];
end
conn_result3=conn_result3';
%Count the number of  of each node
Neighbors= cellfun(@numel, conn_result);
Neighbors=repmat(Neighbors,1,t_end./dt); 






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SLIDING BLOCK CONNECTIVITY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C = cell(20);
for i=1:20
    for j=1:20
C{i,1}=i;
C{i,j}=20*(j-1)+i;
    end
end
blockSize = 3;  % assumed square with sliding blocks so blockSize x blockSize
imageSize = 20; % assumed square
array_size=imageSize - blockSize +1;
blockArray = cell(imageSize - blockSize +1);
 for r=1:size(blockArray,1)
     for c=1:size(blockArray,2)
         blockArray{r,c} =C(r:r - 1 + blockSize, ...
                                   c:c - 1 + blockSize);
     end
 end
for r=size(blockArray,1)+1:imageSize
     for c=1:size(blockArray,2)
         blockArray{r,c}=C(r:end, c:c - 1 + blockSize);
     end
end

for r=1:array_size
for c=array_size+1:imageSize
      blockArray{r,c}=C(r:r - 1 + blockSize, c:end);
end
end

for r=array_size:imageSize
for c=array_size:imageSize
      blockArray{r,c}=C(r:end, c:end);
end
end
blockArray=reshape(blockArray,N*M,1);
for k=1:length(blockArray)
blockArray{k}=reshape(blockArray{k},size(blockArray{k},1)*size(blockArray{k},2),1);
end
for i=1:N*M
conn_result4{i,1}=cell2mat(blockArray{i,1});
end
for i=1:1:length(blockArray)-1
 blockArray{i,1}(1)=[];
end
for i=1:N*M
blockArray{i,1}=cell2mat(blockArray{i,1});
end
conn_result2=blockArray;




%%%%%%%%1%%%%%%%%%%%%%%%%%  DEFINE PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha=0.001;
V_th = 1; %spike threshold [mV]
V_reset = 0; %value to reset voltage to after a spike [mV] 
V_spike = 5; %value to draw a spike to, when cell spikes [mV]
%%%%%%%%%%%%%%%%%%%%%%%  DEFINE INITIAL VALUES AND VECTORS TO HOLD RESULTS %%%%%%%
t_vect = 0:dt:t_end; %will hold vector of times
Gamma=zeros(1,length(t_vect));
V_vect = zeros(N*M,length(t_vect)); %initialize the voltage vector
V_vect(:,1)=reshape(A,[N*M,1]); %first element of V, i.e. value of V at t=0
V_plot_vect = zeros(N*M,length(t_vect)); %pretty version of V_vect to be plotted displaying a spike
% whenever voltage reaches threshold





%%%%%%  INTEGRATE THE INTEGRATE-AND-FIRE NETWORK tau*dV/dt = -V + I+ Sum(J) WITH EULERS METHOD
h=dt;
I=zeros(N*M,t_end./dt); %initialize injected current
NumSpikes=0;       %initialize number of spikes so can count them
conn=zeros(400,t_end./dt);  %initialize connectivity which will change
for i=1:t_end./dt %loop through values of t in step sizes of dt. 
for j=1:N*M





%pixel test for injected current 
vv=V_vect(:,i*ones(size(V_vect,1),1)); % i'th time step for V_vect
diff2=abs(vv(conn_result4{j,1},1)-vv(conn_result4{j,1},1))' ; %column vector - row vector
uppertriang2=triu(diff2,1); %just look at elements in upper triangular portion of matrix
[row2,col2]=find(uppertriang2(1,:)<0.15); % just look at difference with current node (or first row)
NUMEL2=numel(col2);
if NUMEL2>=numel(conn_result4(j))/2 % half of nodes in the neighborhood
    I(j,i)=V_th+0.6;  %Injected current when half of nodes meet pixel test
elseif NUMEL2==0    %Injected current when none of nodes meet pixel test
     I(j,i)=0;
else
     I(j,i)=V_th-0.011;   %Injected current.
end




%pixel test for connectivity function
diff=abs(vv(conn_result3{j,1},1)-vv(conn_result3{j,1},1))'; %difference
uppertriang=triu(diff,1); % upper triangular portion of matrix
[row,col]=find(uppertriang(1,:)>0.15); %finds row indices where pixel test is not satisfied
NUMEL=numel(col); %number of neighbors not satisfying pixel test
Neighbors(j,i)=(Neighbors(j,i))-NUMEL; %subtract this number of neighbors from total number of neighbors
conn(j,i)=J(alpha,Neighbors([conn_result{conn_result{j}}],i));



V_vect(j,i+1) = V_vect(j,i)+dv_dt(V_vect(j,i),conn(j,1),I(j,i),Gamma(i))*h; %Euler's Method with h=0.01




%%%%%%%%%%%%%%%%%%%%%%% SPIKING NEURONS CODE %%%%%%%%%%%%%%%%%%%%%%
if V_vect(j,i+1) > V_th %cell spiked
      Gamma(1,i+1)=0.001; %update gamma, the amount of desynchrony
      V_vect(j,i+1)=V_reset;       %set voltage back to V_reset of nodes connected to the previous node 
      V_plot_vect(j,i+1) = V_spike; %set vector that will be plotted to show a spike here 
      V_vect(conn_result{j,:},i+1)=conn(j,i);
while V_vect(conn_result{j,:},i+1)>V_th %while greater than threshold, run this loop
      V_vect(conn_result{j,:},i+1)=V_vect(conn_result{j,:},i)+conn(j,i)-1; %update rule for neighboring nodes
      V_vect(conn_result{conn_result{j,:},:},i+1) =conn(j,i);   %update rule for neighbors of neighbors. They get an impulse.
      V_vect(j,i+1)=conn(j,i);
end
    NumSpikes=NumSpikes+1;  % count the number of spikes.
else %voltage didn't cross threshold so cell does not spike
    V_plot_vect(j,i+1) = V_vect(j,i+1); %plot the actual voltage
end
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%  PLOT THE RESULTS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
low=min(V_vect(:,i));
 high=max(V_vect(:,i));
 level=(high-low)/2;


surf(reshape(V_vect(:,i),[N,M])); %plot the pixels
view(0,90); %top-down view of oscillators
zlim([-10 10])
shading interp;
drawnow
%pause(0.5);

 
end




%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function slope1 =dv_dt(v,conn,I,Gamma) %standard integrate-and-fire network
slope1=-v+I-Gamma;
slope1=slope1+conn;
end

function conn=J(alpha,Neighbors) %connectivity
conn=sum(alpha./Neighbors,1);
end


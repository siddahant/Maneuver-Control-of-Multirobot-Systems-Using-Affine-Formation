clc; clear; close all
global nodenum neighborMat stressMat dim simTime leaderNum ifLeaderFollower kp kv samplePos_all
global rotateFlag tempp1 tempp2
%
ifLeaderFollower=1;
leaderNum=3; % number of Leaders
dim=2;
% configuration matrix
P=2*[2 0;1 1;1 -1;0 1;0 -1;-1 1;-1 -1];
nodenum=size(P,1); % number of nodes
% adjacent matrix:- tells about all the connections between the agents for
% information
neighborMat=zeros(nodenum,nodenum);
neighborMat(1,2)=1;neighborMat(1,3)=1;neighborMat(1,4)=1;neighborMat(1,5)=1;
neighborMat(2,4)=1;
neighborMat(3,5)=1;neighborMat(3,6)=1;
neighborMat(4,5)=1;neighborMat(4,6)=1;neighborMat(2,7)=1;
neighborMat(5,7)=1;
neighborMat(6,7)=1;
m=sum(sum(neighborMat)); % number of edges
neighborMat=neighborMat+neighborMat';
stressMat=StressMatrix(dim,nodenum,m,P,neighborMat); 

kv=2;kp=0.5; 
samplePos_all=P';
p=reshape(P',dim*nodenum,1);
x_init=p+[0 0 0 0 0 0 1 1 0 1 -1 1 -2 -1]'; 
for k=1:leaderNum
	x_init(dim*(k-1)+1:dim*k)=samplePos_all(:,k); % initial position of all leaders
end
v_leader_init=zeros(dim*leaderNum,1); % Initializing zero velocity for all the leaders and followers
v_init=zeros(dim*nodenum,1);
v_init(1:dim*leaderNum,1)=v_leader_init; 

rotateFlag=1; % Flag variable for rotation
tempp1=zeros(dim,1);
tempp2=zeros(dim,1);
%% SIMULINK
simTime=125;
stepsize=0.005 ;
option=simset('fixedstep',stepsize);
sim('affineManeuver.mdl', [0, simTime], option)


%% ANIMATION REQUEST PLOT RESULTS
simulation = 0;

if simulation==1
    Animate_Maneuver(p_all_time, error_all) % Animate the simulation
elseif simulation==0
    Plot_Results(p_all_time, v_all_time, a_all, error_all, aDiff_all) % Plot the results
end

%% FUNCTION STRESS MATRIX
function stressMat=StressMatrix(d,n,m,P,adjMat)
% incidence matrix
H=zeros(m,n);
pointer=1;
for i=1:n
    for j=i+1:n
        if adjMat(i,j)~=0
            H(pointer,i)=-1;
            H(pointer,j)=1;
            pointer=pointer+1;
        end
    end
end

% Single value decomposition of Configuration Matrix
Pbar=[ones(n,1),P];
[U,S,V]=svd(Pbar);
U2=U(:,d+1+1:end);
Q=U2';

% Plot the formation
hold on; axis equal
% plot edges
for i=1:n
    for j=1:n
        if adjMat(i,j)~=0
            pi=P(i,:);
            pj=P(j,:);
            line([pi(1),pj(1)], [pi(2),pj(2)], 'linewidth', 2, 'color', 'b');
        end
    end
end
% plot nodes
for i=1:n
    pi=P(i,:);
    plot(pi(1), pi(2), 'o', ...
        'MarkerSize', 15,...
        'linewidth', 2,...
        'MarkerEdgeColor', 'r',...
        'markerFaceColor', 'white');
    text(pi(1), pi(2), num2str(i),...
        'color', 'r', 'FontSize', 13, 'horizontalAlignment', 'center', 'FontName', 'times');
end
% numbering edges
pointer=1;
for i=1:n
    for j=i+1:n
        if adjMat(i,j)~=0
            pi=P(i,:);
            pj=P(j,:);
            text((pi(1)+pj(1))/2, (pi(2)+pj(2))/2, num2str(pointer),...
                'color', 'k', 'FontSize', 12, 'horizontalAlignment', 'center');
            pointer=pointer+1;
        end
    end
end

E=zeros(d*n,m);
for i=1:n
    hi=H(:,i);
    E((i-1)*d+1:(i-1)*d+d,:)=P'*H'*diag(hi);
end
% Single Value Decomposition of Energy Matrix
[U,S,V]=svd(E);
B=V(:,rank(S)+1:end); 
num=size(B,2);

% if num>1 -> LMI
M=zeros(n-d-1,n-d-1,num);
for i=1:num
    bi=B(:,i);
    M(:,:,i)=Q*H'*diag(bi)*H*Q';
end
setlmis([]);
LMI=newlmi;
if num==0
    return;
elseif num==1
    x1=lmivar(1,[1 0]);
    lmiterm([-LMI 1 1 x1],M(:,:,1),1); % right hand side
    LMIsys=getlmis;
    [tmin,x] = feasp(LMIsys,[0 0 1 0 0]);
elseif num==2
    x1=lmivar(1,[1 0]);
    x2=lmivar(1,[1 0]);
    lmiterm([-LMI 1 1 x1],M(:,:,1),1); % right hand side
    lmiterm([-LMI 1 1 x2],M(:,:,2),1); % right hand side
    LMIsys=getlmis;
    [tmin,x] = feasp(LMIsys,[0 0 1 0 0]); 
elseif num>2
    return;
end
% calculate omega and stress matrix
omega=B*x;
omega=omega/norm(omega);
stressMat=H'*diag(omega)*H;
if max(eig(stressMat))>0
    omega=-omega;
    stressMat=-stressMat;
end
end


%% FUNCTION PLOT RESULTS
function Plot_Results(p_all_time, v_all_time, a_all_time, error_all, aDiff_all)
close all

global nodenum neighborMat dim simTime leaderNum

delta=50;
time_all=v_all_time.time(1:delta:end);
p_all=p_all_time.signals.values(1:delta:end,:);
v_all=v_all_time.signals.values(1:delta:end,:);
a_all=a_all_time.signals.values(1:delta:end,:);
error_all=error_all.signals.values(1:delta:end,:);
aDiff_all=aDiff_all.signals.values(1:delta:end,:);

linewidth=0.5;
dotsize=7;
fontsize=7;%15;
formationColor=0*[0 0 1];
faceColorList=zeros(nodenum,3);
for i=1:nodenum
    if i<=leaderNum
        faceColorList(i,:)=[1,0,0];
    else
        faceColorList(i,:)=[0,0,1];
    end
end
for i=1:nodenum
    if i<=leaderNum
        markerList(i)={'o'};
    else
        markerList(i)={'o'};
    end
end
for i=1:nodenum
    if i<=leaderNum
        linestyleList(i)={'-'};
    else
        linestyleList(i)={'--'};
    end
end

% plot trajectory
figure; 
subplot(4,1,[1,2,3]); % to make the axis short
hold on; box on; 
set(gca, 'fontSize', fontsize)
set(get(gca, 'xlabel'), 'String', 'x (m)', 'fontSize', fontsize);
set(get(gca, 'ylabel'), 'String', 'y (m)', 'fontSize', fontsize);
axis equal
set(gca, 'fontsize', fontsize);

%plot trajectory
for i=1:nodenum
    if dim==2
        xi_all=p_all(:,2*i-1);
        yi_all=p_all(:,2*i);
        plot(xi_all, yi_all, ':', 'linewidth', linewidth, 'color', faceColorList(i,:));
    else
        xi_all=p_all(:,3*i-2);
        yi_all=p_all(:,3*i-1);
        zi_all=p_all(:,3*i);
        plot3(xi_all, yi_all, zi_all, ':', 'linewidth', linewidth, 'color', faceColorList(i,:));
    end
end

% plot obstacle
leftx=18;
width=5;
h=rectangle('Position', [leftx, 1.8, width, 15]);
set(h, 'faceColor', 0.4*ones(1,3))
h=rectangle('Position', [leftx, -16.5, width, 14.5]);
set(h, 'faceColor', 0.4*ones(1,3))
h=rectangle('Position', [leftx, -25, width, 6.5]);
set(h, 'faceColor', 0.4*ones(1,3))
text(20.5, -9, 'Obstacle','color', 'w', 'FontSize', 8, 'horizontalAlignment', 'center', 'FontWeight', 'normal', 'Rotation', 0, 'fontname', 'Times');

axis tight
xlim=get(gca,'xlim');
set(gca,'xlim', xlim+[-5,5]);
ylim=get(gca,'ylim');
set(gca,'ylim', ylim+[-4,3]);

% plot intermediate formations
dataNum=size(p_all,1);
hMarkerAll=zeros(nodenum,1);
index=[1,floor(dataNum/10*1.4),floor(dataNum/10*2.2),floor(dataNum/10*3.2),floor(dataNum/10*4.47),...
    floor(dataNum/10*5.3),floor(dataNum/10*6.45),floor(dataNum/10*7.2),floor(dataNum/10*8),floor(dataNum/10*9),floor(dataNum/10*10)]; % for 2D scale big example

for k=1:size(index,2)
    idx=index(k);
    for i=1:nodenum
        for j=1:nodenum
            if neighborMat(i,j)==1
                if dim==2
                    pi=p_all(idx,2*i-1:2*i)';
                    pj=p_all(idx,2*j-1:2*j)';
                    line([pi(1),pj(1)], [pi(2),pj(2)], 'linewidth', 0.5, 'color', formationColor);  
                else
                    pi=p_all(idx,3*i-2:3*i)';
                    pj=p_all(idx,3*j-2:3*j)';
                    line([pi(1),pj(1)], [pi(2),pj(2)], [pi(3),pj(3)], 'linewidth', 0.5, 'color', formationColor);
                end
            end
        end
    end    
    for i=1:nodenum
        if dim==2
            xi=p_all(idx,2*i-1);
            yi=p_all(idx,2*i);
            hMarkerAll(i)=plot(xi, yi, char(markerList(i)), 'MarkerEdgeColor', 'k', 'MarkerFaceColor', faceColorList(i,:), 'markersize', dotsize, 'linewidth', 0.5);
            text(xi, yi, num2str(i),'color', 'w', 'FontSize', 5, 'horizontalAlignment', 'center', 'FontWeight', 'bold');
        else
            xi=p_all(idx,3*i-2);
            yi=p_all(idx,3*i-1);
            zi=p_all(idx,3*i);
            hMarkerAll(i)=plot3(xi, yi, zi, char(markerList(i)), 'MarkerEdgeColor', 'k', 'MarkerFaceColor', faceColorList(i,:), 'markersize', dotsize, 'linewidth', 1);
        end
    end
end
hLegendTraj=legend([hMarkerAll(1),hMarkerAll(leaderNum+1)], 'Leader', 'Follower', 'location', 'West');

% plot time for corresponding formation
text(0, 5, strcat('t=',num2str(time_all(index(1)),'%.1f'),'s'),'color', 'k', 'FontSize', fontsize+1, 'horizontalAlignment', 'center', 'FontName', 'times');
text(10, 5, strcat('t=',num2str(time_all(index(2)),'%.1f'),'s'),'color', 'k', 'FontSize', fontsize+1, 'horizontalAlignment', 'center', 'FontName', 'times');
text(33, 5, strcat('t=',num2str(time_all(index(4)),'%.1f'),'s'),'color', 'k', 'FontSize', fontsize+1, 'horizontalAlignment', 'center', 'FontName', 'times');
text(45, 5, strcat('t=',num2str(time_all(index(5)),'%.1f'),'s'),'color', 'k', 'FontSize', fontsize+1, 'horizontalAlignment', 'center', 'FontName', 'times');
text(50, -7, strcat('t=',num2str(time_all(index(6)),'%.1f'),'s'),'color', 'k', 'FontSize', fontsize+1, 'horizontalAlignment', 'center', 'FontName', 'times');
text(50, -17, strcat('t=',num2str(time_all(index(7)),'%.1f'),'s'),'color', 'k', 'FontSize', fontsize+1, 'horizontalAlignment', 'center', 'FontName', 'times');
text(38, -22, strcat('t=',num2str(time_all(index(8)),'%.1f'),'s'),'color', 'k', 'FontSize', fontsize+1, 'horizontalAlignment', 'center', 'FontName', 'times');
text(28, -22, strcat('t=',num2str(time_all(index(9)),'%.1f'),'s'),'color', 'k', 'FontSize', fontsize+1, 'horizontalAlignment', 'center', 'FontName', 'times');
text(14, -22, strcat('t=',num2str(time_all(index(10)),'%.1f'),'s'),'color', 'k', 'FontSize', fontsize+1, 'horizontalAlignment', 'center', 'FontName', 'times');
text(3, -22, strcat('t=',num2str(time_all(index(11)),'%.1f'),'s'),'color', 'k', 'FontSize', fontsize+1, 'horizontalAlignment', 'center', 'FontName', 'times');

% plot error
figure
subplot(5,1,1) % to make the axis short
hold on; box on
% set the format of the axis
set(gca, 'fontSize', fontsize)
% set(get(gca, 'xlabel'), 'String', 'Time (s)', 'fontSize', fontsize);
set(get(gca, 'ylabel'), 'String', 'Tracking error', 'fontSize', fontsize);
set(gca, 'xlim', [0,simTime])
% set(gca, 'ylim', [0,10])
plot(time_all, error_all, 'color', 'm', 'linewidth', 1)
% set(get(gca, 'title'), 'String', strcat('Tracking error: kp=',num2str(kp),', kv=',num2str(kv)), 'fontsize', fontsize);

% plot velocity
hLineAll=zeros(nodenum,1);
if dim==2
    % x velocity
    subplot(5,1,2)
    hold on; box on;
    for i=1:nodenum
        hLineAll(i)=plot(time_all, v_all(:,2*i-1), 'color', faceColorList(i,:), 'linestyle', char(linestyleList(i)), 'linewidth', 1);
    end
    set(gca, 'fontSize', fontsize)
    set(get(gca, 'ylabel'), 'String', 'x-velocity (m/s)', 'fontSize', fontsize);
    set(gca, 'xlim', [0,simTime])
    set(gca, 'ylim', [-3,3])
    hLegend=legend([hLineAll(1),hLineAll(leaderNum+1)], 'Leader', 'Follower', 'location', 'southwest');
    set(hLegend, 'fontsize', fontsize);
    % y velocity
    subplot(5,1,3)
    hold on; box on;
    for i=1:nodenum
        hLineAll(i)=plot(time_all, v_all(:,2*i), 'color', faceColorList(i,:), 'linestyle', char(linestyleList(i)), 'linewidth', 1);
    end
    set(gca, 'fontSize', fontsize)
    set(get(gca, 'ylabel'), 'String', 'y-velocity (m/s)', 'fontSize', fontsize);
    set(gca, 'xlim', [0,simTime])
    set(gca, 'ylim', [-3,3])
    hLegend=legend([hLineAll(1),hLineAll(leaderNum+1)], 'Leader', 'Follower', 'location', 'southwest');
    set(hLegend, 'fontsize', fontsize);
elseif dim==3 % plot vx and vy
    % x velocity
    subplot(3,1,1)
    hold on; box on;
    for i=1:nodenum
        hLineAll(i)=plot(time_all, v_all(:,3*i-2), 'color', faceColorList(i,:), 'linestyle', char(linestyleList(i)), 'linewidth', 1);
    end
    set(gca, 'fontSize', fontsize)
    set(get(gca, 'ylabel'), 'String', 'x-velocity (m/s)', 'fontSize', fontsize);
    set(gca, 'xlim', [0,simTime])
    hLegend=legend([hLineAll(1),hLineAll(leaderNum+1)], 'Leader', 'Follower', 'location', 'southwest');
    set(hLegend, 'fontsize', fontsize);
    % title
    set(get(gca, 'title'), 'String', 'Velocity', 'fontsize', fontsize);
    % y velocity
    subplot(3,1,2)
    hold on; box on;
    for i=1:nodenum
        hLineAll(i)=plot(time_all, v_all(:,3*i-1), 'color', faceColorList(i,:), 'linestyle', char(linestyleList(i)), 'linewidth', 1);
    end
    set(gca, 'fontSize', fontsize)
    set(get(gca, 'ylabel'), 'String', 'y-velocity (m/s)', 'fontSize', fontsize);
    set(gca, 'xlim', [0,simTime])
    hLegend=legend([hLineAll(1),hLineAll(3)], 'Leader', 'Follower', 'location', 'southwest');
    set(hLegend, 'fontsize', fontsize);
    % z velocity
    subplot(3,1,3)
    hold on; box on;
    for i=1:nodenum
        hLineAll(i)=plot(time_all, v_all(:,3*i), 'color', faceColorList(i,:), 'linestyle', char(linestyleList(i)), 'linewidth', 1);
    end
    set(gca, 'fontSize', fontsize)
    set(get(gca, 'xlabel'), 'String', 'Time (sec)', 'fontSize', fontsize);
    set(get(gca, 'ylabel'), 'String', 'z-velocity (m/s)', 'fontSize', fontsize);
    set(gca, 'xlim', [0,simTime])
    hLegend=legend([hLineAll(1),hLineAll(3)], 'Leader', 'Follower', 'location', 'southwest');
    set(hLegend, 'fontsize', fontsize);
end

% plot acceleration
hLineAll=zeros(nodenum,1);
if dim==2
    % x velocity
    subplot(5,1,4)
    hold on; box on;
    for i=1:nodenum
        hLineAll(i)=plot(time_all, a_all(:,2*i-1), 'color', faceColorList(i,:), 'linestyle', char(linestyleList(i)), 'linewidth', 1);
    end
    set(gca, 'fontSize', fontsize)
    set(get(gca, 'ylabel'), 'String', 'x-acceleration (m^2/s)', 'fontSize', fontsize);
    set(gca, 'xlim', [0,simTime])
    hLegend=legend([hLineAll(1),hLineAll(leaderNum+1)], 'Leader', 'Follower', 'location', 'southwest');
    set(hLegend, 'fontsize', fontsize);
    % title
%     set(get(gca, 'title'), 'String', 'Acceleration', 'fontsize', fontsize);
    % y velocity
    subplot(5,1,5)
    hold on; box on;
    for i=1:nodenum
        hLineAll(i)=plot(time_all, a_all(:,2*i), 'color', faceColorList(i,:), 'linestyle', char(linestyleList(i)), 'linewidth', 1);
    end
    set(gca, 'fontSize', fontsize)
    set(get(gca, 'ylabel'), 'String', 'y-acceleration (m^2/s)', 'fontSize', fontsize);
    set(get(gca, 'xlabel'), 'String', 'Time (s)', 'fontSize', fontsize);
    set(gca, 'xlim', [0,simTime])
    hLegend=legend([hLineAll(1),hLineAll(leaderNum+1)], 'Leader', 'Follower', 'location', 'southwest');
    set(hLegend, 'fontsize', fontsize);
end

% plot acceleration error
figure
hLineAll=zeros(nodenum,1);
if dim==2
    % x velocity
    subplot(5,1,1)
    hold on; box on;
    for i=1:nodenum
        hLineAll(i)=plot(time_all, aDiff_all(:,2*i-1), 'color', faceColorList(i,:), 'linestyle', char(linestyleList(i)), 'linewidth', 1);
    end
    set(gca, 'fontSize', fontsize)
    set(get(gca, 'ylabel'), 'String', 'x-aDiff (m^2/s)', 'fontSize', fontsize);
    set(gca, 'xlim', [0,simTime])
    hLegend=legend([hLineAll(1),hLineAll(leaderNum+1)], 'Leader', 'Follower', 'location', 'southwest');
    set(hLegend, 'fontsize', fontsize);
    hLegend=legend([hLineAll(1),hLineAll(leaderNum+1)], 'Leader', 'Follower', 'location', 'southwest');
    % y velocity
    subplot(5,1,2)
    hold on; box on;
    for i=1:nodenum
        hLineAll(i)=plot(time_all, aDiff_all(:,2*i), 'color', faceColorList(i,:), 'linestyle', char(linestyleList(i)), 'linewidth', 1);
    end
    set(gca, 'fontSize', fontsize)
    set(get(gca, 'ylabel'), 'String', 'y-aDiff (m^2/s)', 'fontSize', fontsize);
    set(gca, 'xlim', [0,simTime])
    set(get(gca, 'xlabel'), 'String', 'Time (s)', 'fontSize', fontsize);
    hLegend=legend([hLineAll(1),hLineAll(leaderNum+1)], 'Leader', 'Follower', 'location', 'southwest');
    set(hLegend, 'fontsize', fontsize);
end
end

%% FUNCTION ANIMATE THE RESULTS
function Animate_Maneuver(p_all_time, error_all)
global nodenum edgenum neighborMat dim simTime leaderNum ki kp
close all;
figure;
set(gcf, 'unit', 'norm', 'pos', [0.1,0.1,0.72,0.8])
hAxisTraj=subplot(4,1,1:3);
hold on; box on; axis equal
xlim([-10,55])
ylim([-25,6])
delta=50;
time_all=p_all_time.time(1:delta:end);
p_all=p_all_time.signals.values(1:delta:end,:);
error_all=error_all.signals.values(1:delta:end,:);

linewidth=0.5;
fontsize=7;
dotsize=9;
formationColor=0*[0 0 1];

for i=1:nodenum
    if i<=leaderNum
        edgeColorList(i,:)=[1,0,0];
    else
        edgeColorList(i,:)=[0,0,1];
    end
end
faceColorList=edgeColorList;
set(gca, 'fontSize', fontsize)
set(get(gca, 'xlabel'), 'String', 'x (meter)', 'fontSize', fontsize);
set(get(gca, 'ylabel'), 'String', 'y (meter)', 'fontSize', fontsize);

for i=1:nodenum
    for j=1:nodenum
        if neighborMat(i,j)==1
            pi=p_all(1,2*i-1:2*i)';
            pj=p_all(1,2*j-1:2*j)';
            line([pi(1),pj(1)], [pi(2),pj(2)], 'linewidth', 1, 'color', formationColor);
        end
    end
    xi=p_all(1,2*i-1);
    yi=p_all(1,2*i);
    plot(xi, yi, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', faceColorList(i,:), 'markersize', dotsize-2, 'linewidth', 1)
end

% plot obstacle
leftx=18;
width=5;
h=rectangle('Position', [leftx, 1.8, width, 15]);
set(h, 'faceColor', 0.4*ones(1,3))
h=rectangle('Position', [leftx, -16.5, width, 14.5]);
set(h, 'faceColor', 0.4*ones(1,3))
h=rectangle('Position', [leftx, -25, width, 6.5]);
set(h, 'faceColor', 0.4*ones(1,3))
text(20.5, -9, 'Obstacle','color', 'w', 'FontSize', 8, 'horizontalAlignment', 'center', 'FontWeight', 'normal', 'Rotation', 0, 'fontname', 'Times');


% objects of the lines
hLine=zeros(nodenum,nodenum);
for i=1:nodenum
    for j=i+1:nodenum
        if neighborMat(i,j)==1
            pi=p_all(1,2*i-1:2*i)';
            pj=p_all(1,2*j-1:2*j)';
            hLine(i,j)=line([pi(1),pj(1)], [pi(2),pj(2)], 'linewidth', linewidth, 'color', formationColor);
        end
    end
end
% >>>objects of the dots
hMarker=zeros(1,nodenum);
hText=zeros(1,nodenum);
for i=1:nodenum
    xi=p_all(1,2*i-1);
    yi=p_all(1,2*i);
    hMarker(i) = plot(xi, yi, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', faceColorList(i,:), 'markersize', dotsize, 'linewidth', 1);
    hText(i)=text(xi, yi, num2str(i),'color', 'w', 'FontSize', fontsize, 'horizontalAlignment', 'center', 'FontWeight', 'bold');
end
hLegend=legend([hMarker(1),hMarker(leaderNum+1)], 'Leader', 'Follower', 'location', 'west');
set(hLegend, 'fontsize', fontsize);
% objects of the trajectories
hTraj=zeros(1,nodenum);
for i=1:nodenum
    xi_all=p_all(1,2*i-1);
    yi_all=p_all(1,2*i);
    hTraj(i) = plot(xi_all, yi_all, ':', 'linewidth', 1, 'color', edgeColorList(i,:));
end

% tracking error
hAxisError=subplot(4,1,4);
hold on; box on
% set the format of the axis
set(gca, 'fontSize', fontsize)
set(get(gca, 'xlabel'), 'String', 'Time (second)', 'fontSize', fontsize);
set(get(gca, 'ylabel'), 'String', 'Tracking error', 'fontSize', fontsize);
set(gca, 'xlim',[0 125])
set(gca, 'ylim', [0 4])
hError=plot(time_all(1), error_all(1),  'm', 'linewidth', 2);

% update the position of each object
for k=1:4:size(p_all, 1)
    for i=1:nodenum
        %
        xi_all=p_all(1:k,2*i-1);
        yi_all=p_all(1:k,2*i);
        set(hTraj(i), 'xdata', xi_all, 'ydata', yi_all);
        %
        xi=p_all(k,2*i-1);
        yi=p_all(k,2*i);
        set(hMarker(i), 'xdata', xi, 'ydata', yi);
        set(hText(i), 'Position', [xi,yi]);
        %
        for j=i+1:nodenum
            if neighborMat(i,j)==1
                pi=p_all(k,2*i-1:2*i)';
                pj=p_all(k,2*j-1:2*j)';
                set( hLine(i,j), 'xdata', [pi(1),pj(1)], 'ydata', [pi(2),pj(2)]);
            end
        end
        set(hError, 'xdata', time_all(1:k), 'ydata', error_all(1:k));
    end
    pause(0.1)
end
end


function output=LeaderAcceleration0(u)
global leaderNum dim nodenum 
global rotateFlag tempp2 tempp3

u_all = reshape(u(1:dim*nodenum*3),dim,nodenum*3);
v_all=u_all(:,1:nodenum);
p_all=u_all(:,nodenum+1:2*nodenum);
control_all=u_all(:,2*nodenum+1:3*nodenum);
currentTime=u(end);
%output
a_leader_all=zeros(dim,leaderNum);

% speed up to speed=1
time_start_speedup=0;
time_span_speedup=5;
if currentTime>time_start_speedup && currentTime<time_start_speedup+time_span_speedup
    ax=f_sin(time_start_speedup,time_span_speedup,currentTime)/pi;
    a_leader_all(:,1)=[ax,0]';
    a_leader_all(:,2)=[ax,0]';
    a_leader_all(:,3)=[ax,0]';
end

% slow down to speed=0
time_start_slowdown=10;
time_span_slowdown=5;
if currentTime>time_start_slowdown && currentTime<time_start_slowdown+time_span_slowdown
    %ax=fcn_normalDistribution(time_start_slowdown,time_span_slowdown, 2, currentTime);
    ax=-f_sin(time_start_slowdown,time_span_slowdown,currentTime)/pi;
    a_leader_all(:,1)=[ax,0]';
    a_leader_all(:,2)=[ax,0]';
    a_leader_all(:,3)=[ax,0]';
end

% rotate
v3=v_all(:,3);
v2=v_all(:,2);
v1=v_all(:,1);
p3=p_all(:,3);
p2=p_all(:,2);
p1=p_all(:,1);
R90=[cos(pi/2),-sin(pi/2);
     sin(pi/2),cos(pi/2)];
time_start_rotate=20;
time_span_rotate=10;
if currentTime>=time_start_rotate && currentTime<=time_start_rotate+time_span_rotate
    if rotateFlag==1
        tempp2=p2;
        tempp3=p3;
        rotateFlag=0;
    end
    % agent 2
    acce_mag=-Accel_slow_fast(time_start_rotate, time_span_rotate, currentTime, 4);
%     acce_vec=(p1+R90'*(tempp2-p1))-tempp2;
%     acce_vec=acce_vec/norm(acce_vec);
    acce_vec=[1 0]';
    a_leader_all(:,2)=acce_mag*acce_vec;
    % agent 3
    acce_mag=-Accel_slow_fast(time_start_rotate, time_span_rotate, currentTime, 4);
%     acce_vec=(p1+R90'*(tempp3-p1))-tempp3;
%     acce_vec=acce_vec/norm(acce_vec);
    acce_vec=[0 1]';
    a_leader_all(:,3)=acce_mag*acce_vec;
end

% move downward
time_start_downward=time_start_rotate;
time_span_downward=time_span_rotate;
if currentTime>=time_start_downward && currentTime<=time_start_downward+time_span_downward
    % agent 1
    ay=f_sin(time_start_downward,time_span_downward,currentTime)/pi;
    a_leader_all(:,1)=[0,-ay]';
    % agent 2
    a_leader_all(:,2)=a_leader_all(:,2)+[0,-ay]';
    % agent 3
    a_leader_all(:,3)=a_leader_all(:,3)+[0,-ay]';
end

% slow down to speed=0
time_start_slowdown2=time_start_downward+time_span_downward+5;
time_span_slowdown2=5;
if currentTime>time_start_slowdown2 && currentTime<time_start_slowdown2+time_span_slowdown2
    %ax=fcn_normalDistribution(time_start_slowdown,time_span_slowdown, 2, currentTime);
    ay=-f_sin(time_start_slowdown2,time_span_slowdown2,currentTime)/pi;
    a_leader_all(:,1)=[0,-ay]';
    a_leader_all(:,2)=[0,-ay]';
    a_leader_all(:,3)=[0,-ay]';
end

% output: replace the leaders' velocity with the one prescribed above
control_all(:,1:leaderNum) = a_leader_all;
output=reshape(control_all, dim*nodenum, 1);




% 模拟退火算法求解TSP
% 2-邻边法修改TSP回路，热平衡迭代次数L=10

clc,    clear all
tic      %计时开始

e=0.1^20;    L=10;    at=0.9995;    T=10^3;  %初始化退火算法因子

data0=load('TSP_Data.txt');    %加载100个目标的数据，数据按照表格中的位置保存在纯文本文件TSP_Data.txt中
x=data0(:,[1:2:8]);    x=x(:);               %目标数据经度坐标
y=data0(:,[2:2:8]);    y=y(:);               %目标数据纬度坐标
sj=[x y];
d1=[70,40];                                       %基地经纬度
sj=[d1;sj;d1];    sj=sj*pi/180;             %角度化成弧度
d=zeros(102);                                  %距离矩阵d初始化

%计算距离矩阵
for i=1:101
   for j=i+1:102
d(i,j)=6370*acos(cos(sj(i,1)-sj(j,1))*cos(sj(i,2))*cos(sj(j,2))+sin(sj(i,2))*sin(sj(j,2)));
   end
end
d=d+d';

%从10000组随机TSP回路中，选择最好的作为初始回路
path=[];    long=inf;                             %巡航路径及长度初始化
rand('state',    sum(clock));                 %初始化随机数发生器
for j=1:10000                                      %从10000组随机解中，求较好的初始解
    path0=[1 1+randperm(100),102];    temp=0;
    for i=1:101
        temp=temp+d(path0(i),path0(i+1));
    end
    if temp<long
        path=path0;    long=temp;
    end
end

% 模拟退火算法(SA)主循环迭代，2-邻边法修改回路
while T>e
    for k=1:L                                           % 退火过程在温度T时，进行L次迭代，趋向热平衡
        c=randi([2,101],1,2);                     % 随机确定2-邻边的点
        %c=2+floor(100*rand(1,2));
        c=sort(c);    c1=c(1);    c2=c(2);

        % 2-邻边法变换后，回路长度增量
        df=d(path(c1-1),path(c2))+d(path(c1),path(c2+1))...
             -d(path(c1-1),path(c1))-d(path(c2),path(c2+1));
        if df<0                              % df<0，回路更新
            path=[path(1:c1-1),path( c2:-1:c1),path(c2+1:102)]; long=long+df;
        elseif rand<=exp(-df/T)    % df>0，且符合概率条件，回路更新
            path=[path(1:c1-1),path(c2:-1:c1),path(c2+1:102)]; long=long+df;
        end
    end
    T=T*at;                                %降低温度
end

path,    long                             %输出巡航路径及路径长度
toc               %计时结束
%save('path_star.mat', 'path');
%save('path_long.txt','path','long','-ascii');

%画出巡航路径
xx=sj(path,1);    yy=sj(path,2);    plot(xx,yy,'-*')
% 模拟退火算法求解TSP
% 随机选择2-邻边法或3-邻边法修改回路
% 热平衡迭代次数L_T=1000

clc,    clear all
tic                                 %计时开始

T=10^3;    L_T=1000;    a=0.995;    T_f=0.1^12;  %初始化退火算法因子

data0=load('TSP_Data.txt');   %加载100个目标的经纬度数据，数据按照表格中的位置保存在纯文本文件TSP_Data.txt中
x=data0(:,[1:2:8]);    x=x(:);            %目标数据经度坐标
y=data0(:,[2:2:8]);    y=y(:);            %目标数据纬度坐标
data_bl=[x y];
d1=[70,40];                                       %基地经纬度
data_bl0=[d1;data_bl;d1];    data_bl=data_bl0*pi/180;           %角度化成弧度
D=zeros(102);                                  %距离矩阵D初始化

%计算距离矩阵
for i=1:101
   for j=i+1:102
D(i,j)=6370*acos(cos(data_bl(i,1)-data_bl(j,1))*cos(data_bl(i,2))*cos(data_bl(j,2))+sin(data_bl(i,2))*sin(data_bl(j,2)));
   end
end
D=D+D';

%从10000组随机TSP回路中，选择最好的作为初始回路
path=[];    long=inf;                          %巡航路径及长度初始化
rand('state',    sum(clock));              %初始化随机数发生器
for j=1:10000                                   %从10000组随机解中，求较好的初始解
    path0=[1 1+randperm(100),102];    temp=0;
    for i=1:101
        temp=temp+D(path0(i),path0(i+1));
    end
    if temp<long
        path=path0;    long=temp;
    end
end

%模拟退火主循环迭代，随机选择2-邻边法或3-邻边法修改回路
while T>T_f                                      %温度大于最低温度，继续冷却过程
    for k=1:L_T                                   %当前温度T，进行L_T次迭代，近似热平衡
        r=rand;
        if (r<0.5)
            %2-邻近点可行解
            c=2+floor(100*rand(1,2));     %产生两个2-邻近点，去构造新解
            c=sort(c);    c1=c(1);    c2=c(2);

            %计算代价函数值的增量
            df=D(path(c1-1),path(c2))+D(path(c1),path(c2+1))-D(path(c1-1),path(c1))-D(path(c2),path(c2+1));
        
            if df<0           %接受准则
                path=[path(1:c1-1),path( c2:-1:c1),path(c2+1:102)];    long=0;
                for i=1:101                        %计算当前TSP回路长
                    long=long+D(path(i),path(i+1));
                end
            elseif exp(-df/T)>rand           %虽然df>0，符合概率值要求，依然调整路径
                path=[path(1:c1-1),path(c2:-1:c1),path(c2+1:102)]; long=0;
                for i=1:101                        %计算当前TSP回路长
                    long=long+D(path(i),path(i+1));
                end
            end
        else            
            %3-邻近点可行解
            c=2+floor(100*rand(1,3));     %产生两个2-邻近点，去构造新解
            c=sort(c);    c1=c(1);    c2=c(2);    c3=c(3);
            
            %计算代价函数值的增量
            df=D(path(c1-1),path(c2))+D(path(c3-1),path(c1))+D(path(c2-1),path(c3))-D(path(c1-1),path(c1))-D(path(c2-1),path(c2))-D(path(c3-1),path(c3));
            
            if df<0           %接受准则
                path=[path(1:c1-1),path( c2:c3-1),path(c1:c2-1),path(c3:102)];    long=0;
                for i=1:101                        %计算当前TSP回路长
                    long=long+D(path(i),path(i+1));
                end
            elseif exp(-df/T)>rand          %虽然df>0，符合概率值要求，依然调整路径
                path=[path(1:c1-1),path( c2:c3-1),path(c1:c2-1),path(c3:102)]; long=0;
                for i=1:101                        %计算当前TSP回路长
                    long=long+D(path(i),path(i+1));
                end
            end
        end
    end
    T=T*a;                                          %降低温度
end

path,    long                                      %输出巡航路径及路径长度
toc                  %计时结束
%save('path-star.mat', 'path');
%save('path-long.txt','path','long','-ascii');

%画出巡航路径
xx=data_bl0(path,1);    yy=data_bl0(path,2);    plot(xx,yy,'-*')
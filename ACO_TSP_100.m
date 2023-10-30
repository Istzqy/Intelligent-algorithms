%蚁群算法求解TSP

%%1. 清空环境变量
clear all;      clc; 
tic

%%2. 导入数据
%load TSP_Data.txt
data0=load('TSP_Data.txt');   %加载100个目标的经纬度数据，数据按照表格中的位置保存在纯文本文件TSP_Data.txt中
x=data0(:,[1:2:8]);    x=x(:);                   %目标数据经度坐标
y=data0(:,[2:2:8]);    y=y(:);                   %目标数据纬度坐标
data_bl=[x y];
d1=[70,40];                                     %基地经纬度，基地作为顶点1，100个目标为顶点2-101
data_bl0=[d1;data_bl];    data_bl=data_bl0*pi/180;           %角度化成弧度
n = size(data_bl,1);
D=zeros(n)*1e-6;                                 %距离矩阵D初始化

%%3. 计算距离矩阵
for i=1:n
   for j=1:n
       if i ~= j
           D(i,j)=6370*acos(cos(data_bl(i,1)-data_bl(j,1))*cos(data_bl(i,2))*cos(data_bl(j,2))...
                     +sin(data_bl(i,2))*sin(data_bl(j,2)));
       end
   end
end

%%4. 初始化参数
m = 100;                            % 蚂蚁数量
alpha = 1;                           % 信息素重要程度因子
beta = 5;                             % 启发函数重要程度因子
rho = 0.1;                            % 信息素挥发因子
Q = 1;                                 % 常系数
Eta = 1./D;                          % 启发函数
Tau = ones(n,n);                 % 信息素矩阵
Table = zeros(m,n);             % 路径记录表(禁忌表)
iter = 1;                               % 迭代次数初值
iter_max = 800;                  % 最大迭代次数 
Route_best = zeros(iter_max,n);                 % 各代最佳路径       
Length_best = zeros(iter_max,1);               % 各代最佳路径的长度  
Length_ave = zeros(iter_max,1);                % 各代路径的平均长度 

%%5. 迭代寻找最佳路径
while iter <= iter_max
    %5.1 随机产生各个蚂蚁的起点城市
    start = zeros(m,1);
    for i = 1:m
        temp = randperm(n);
        start(i) = temp(1);
    end
    Table(:,1) = start; 
    
    % 构建解空间
    citys_index = 1:n;
      
    %5.2 蚂蚁随机选择回路
    for i = 1:m                     % 每个蚂蚁路径选择
       for j = 2:n                   % 逐个城市路径选择
           tabu = Table(i,1:(j - 1));                          % 已访问的城市集合(禁忌表)
           allow_index = ~ismember(citys_index,tabu);    % 不在禁忌表城市序号
           allow = citys_index(allow_index);           % 待访问的城市集合
           P = allow;
           % 计算城市间转移概率
           for k = 1:length(allow)
               P(k) = Tau(tabu(end),allow(k))^alpha * Eta(tabu(end),allow(k))^beta;
           end
           P = P/sum(P);
           % 用轮盘赌随机选择下一个访问城市
           Pc = cumsum(P);                             %与P同维概率向量，Pc是累计概率    
           target_index = find(Pc >= rand);      %大于随机数rand的Pc分量序号
           target = allow(target_index(1));       %蚂蚁随机选择的城市序号
           Table(i,j) = target;                            %选定城市序号加入路径记录表(禁忌表)
       end
    end
    
    %5.3 通过2-邻边算法优化路线
    K = 20;                              %20次2-邻边优化TSP路径            
    for k = 1:m                        %5.2得到的m个蚂蚁TSP路径
        c1 = Table(k,:);             %导出第k个蚂蚁路径
        for t = 1:K                     %2-邻边算法循环，假设进行K轮改进
            flag = 0;                    %2-邻边算法循环的退出标志
            for i = 1:98
                for j = i+2:100
                    if D(c1(i),c1(j)) + D(c1(i+1),c1(j+1)) < D(c1(i),c1(i+1)) + D(c1(j),c1(j+1))
                        c1(i+1:j) = c1(j:-1:i+1);    flag = 1;      %修改TSP路径                        
                    end
                end
            end
            if flag == 0;    break;    end           %2-邻边算法不能减少路径长度
        end
        Table(k,:) = c1;                                %经过2-邻边算法修正的第k个蚂蚁路径
    end
    
    %5.4 计算各个蚂蚁的路径距离
    Length = zeros(m,1);
    for i = 1:m
        Route = Table(i,:);
        for j = 1:(n - 1)
            Length(i) = Length(i) + D(Route(j),Route(j + 1));
        end
        Length(i) = Length(i) + D(Route(n),Route(1));
    end
      
    %5.5 计算最短路径距离及平均距离
    if iter == 1
        [min_Length,min_index] = min(Length);
        Length_best(iter) = min_Length;  
        Length_ave(iter) = mean(Length);
        Route_best(iter,:) = Table(min_index,:);
    else
        [min_Length,min_index] = min(Length);
        Length_best(iter) = min(Length_best(iter - 1),min_Length);     %迭代iter次后，最短路径长
        Length_ave(iter) = mean(Length);
        if Length_best(iter) == min_Length
            Route_best(iter,:) = Table(min_index,:);                               %本次迭代路径为最短
        else
            Route_best(iter,:) = Route_best((iter-1),:);                           %前面路径最短
        end
    end
    
    %5.6 更新信息素
    Delta_Tau = zeros(n,n);
    %5.6.1 逐个蚂蚁计算
    for i = 1:m
        % 逐个城市计算，路径上每条边信息素更新
        for j = 1:(n - 1)
            Delta_Tau(Table(i,j),Table(i,j+1)) = Delta_Tau(Table(i,j),Table(i,j+1)) + Q/Length(i);
        end
        %终点->起点信息素更新
        Delta_Tau(Table(i,n),Table(i,1)) = Delta_Tau(Table(i,n),Table(i,1)) + Q/Length(i);
    end
    Tau = (1-rho) * Tau + Delta_Tau;
    
    %5.6.2 迭代次数加1，清空路径记录表
    iter = iter + 1;
    Table = zeros(m,n);
end
toc

%%6. 结果显示
[Shortest_Length,index] = min(Length_best);
Shortest_Route = Route_best(index,:);
disp(['TSP最短长度:' num2str(Shortest_Length)]);
disp(['TSP最短路径:' num2str([Shortest_Route Shortest_Route(1)])]); 

%%7. 绘图
%7.1 图1：最短TSP路径
figure(1)
plot([data_bl0(Shortest_Route,1);data_bl0(Shortest_Route(1),1)],...
     [data_bl0(Shortest_Route,2);data_bl0(Shortest_Route(1),2)],'o-');
grid on
for i = 1:size(data_bl0,1)
    text(data_bl0(i,1),data_bl0(i,2),['   ' num2str(i)]);
end
%text(citys(Shortest_Route(1),1),citys(Shortest_Route(1),2),'       起点');
%text(citys(Shortest_Route(end),1),citys(Shortest_Route(end),2),'       终点');
xlabel('目标经度')
ylabel('目标纬度')
title(['蚁群算法-TSP回路(最短长度:' num2str(Shortest_Length) ')'])

%7.2 图2：每轮迭代TSP平均路长
figure(2)
plot(1:iter_max,Length_best,'b',1:iter_max,Length_ave,'r:')
legend('当前最短TSP长度','每轮平均TSP长度')
xlabel('迭代次数')
ylabel('TSP长度')
title('TSP最短长度与每轮迭代TSP平均长度图')

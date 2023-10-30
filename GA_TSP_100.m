% 遗传算法求解TSP
% 根据《用 Python 学数学》第12章，对 TSP 进行整数编码、遗传的交叉、变异和选择操作

clc;    clear all;
tic;      % 计时开始

% 1.数据预处理与参数设置
data0=load('TSP_Data.txt');           % 加载100个目标的经纬度数据
x=data0(:,1:2:8);    x=x(:);
y=data0(:,2:2:8);    y=y(:);
sj=[x y];    d1=[70,40];
sj=[d1;sj;d1];    sj=sj*pi/180;      % 单位化成弧度
d=zeros(102);                            % 距离矩阵d的初始值
for i=1:101
    for j=i+1:102
        d(i,j)=6370*acos(cos(sj(i,1)-sj(j,1))*cos(sj(i,2))*cos(sj(j,2))+sin(sj(i,2))*sin(sj(j,2)));
    end
end
d=d+d';    

% 设置遗传算法参数
w=1600;            % w种群规模
g=50;                 % g进化代数      
rand('state',sum(clock));         % 初始化随机数发生器
C_best = zeros(g,102);          % 每代种群的最优路径数组
long_best = zeros(1,g);          % 每代种群的最小长度数组
long_ave = zeros(1,g);           % 每代种群的平均长度数组

% 2.构造初始种群J，J有w行，每行表示一个染色体
% 随机构造TSP顺序(染色体)，经3轮2-邻边法改进染色体
for k=1:w                                % 构造w个初始解，作为初始种群
    c=randperm(100);              % 产生1，...，100的一个全排列
    c1=[1,c+1,102];                  % 构造一个待改进初始解
    for t=1:3                             % 每个初始解进行3轮2-邻边改进
        flag=0;
        for m=1:100
            for n=m+2:101
                if d(c1(m),c1(n))+d(c1(m+1),c1(n+1))<d(c1(m),c1(m+1))+d(c1(n),c1(n+1))
                    c1(m+1:n)=c1(n:-1:m+1);
                    J(k,:)=c1;            % 改进第k条染色体的TSP路径顺序
                    flag=1;
                end
            end
        end
        if flag==0                        % 不能改进当前初始解，则退出当前初始解
            J(k,:)=c1;     break
        end
    end
end

% 3.遗传算法，主循环层
for k=1:g
    % 3.1.调用函数crossover.m，进行染色体单点交叉操作，得到染色体数组A
    for i = 1:w
        parent = randsample(1:w,2);       % 选择种群中两个染色体序号，准备进行交叉操作
        pA = parent(1);
        pB = parent(2);
        A(i,:) = crossover(J,pA,pB);
    end

    % 3.2.调用函数mutate.m，执行染色体变异操作，得到染色体数组B
    for j = 1:w
        B(j,:)=mutate(d,J(j,:));                 % 多次3-邻边，...，10-邻边算法变异，选择TSP回路最小染色体
    end

    % 3.3.从J、A、B中，选择TSP回路最短的w个染色体，构成下一代种群
    G=[J;A;B];
    [N1,N2]=size(G);    
    T_long=zeros(N1,1);
    for ii = 1:N1                                    % 计算每个染色体的TSP回路长度
        for jj = 1:N2-1
            T_long(ii) = T_long(ii)+d(G(ii,jj),G(ii,jj+1));
        end
    end
    [T_min_long,ind]=sort(T_long);
    J=G(ind(1:w),:);                             % w个长度最短染色体，构成下一代种群J
    C_best(k,:)=J(1,:);    
    long_best(k)=T_min_long(1);
    long_ave(k)=mean(T_min_long(1:w));
end

% 4.遗传算法求解结果输出
disp('每一代最小TSP回路长度:');     disp(long_best);
c1=G(ind(1),:);
disp('GA算法最优TSP回路:');     disp(c1); 
flong=T_min_long(1);
disp('GA算法最优TSP回路长度:');     disp(flong); 

toc

% 图形绘制
xx=sj(c1,1);    yy=sj(c1,2);
figure(1)
plot(xx,yy,'-o')
title('图1-GA算法的TSP回路路径')

figure(2)
plot(1:g,long_best,'b',1:g,long_ave,'r')
legend('当前最短TSP长度','每轮平均TSP长度')
xlabel('种群代数')
ylabel('TSP长度')
title('TSP每轮迭代最短长度与平均长度图')

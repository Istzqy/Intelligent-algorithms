%��Ⱥ�㷨���TSP

%%1. ��ջ�������
clear all;      clc; 
tic

%%2. ��������
%load TSP_Data.txt
data0=load('TSP_Data.txt');   %����100��Ŀ��ľ�γ�����ݣ����ݰ��ձ���е�λ�ñ����ڴ��ı��ļ�TSP_Data.txt��
x=data0(:,[1:2:8]);    x=x(:);                   %Ŀ�����ݾ�������
y=data0(:,[2:2:8]);    y=y(:);                   %Ŀ������γ������
data_bl=[x y];
d1=[70,40];                                     %���ؾ�γ�ȣ�������Ϊ����1��100��Ŀ��Ϊ����2-101
data_bl0=[d1;data_bl];    data_bl=data_bl0*pi/180;           %�ǶȻ��ɻ���
n = size(data_bl,1);
D=zeros(n)*1e-6;                                 %�������D��ʼ��

%%3. ����������
for i=1:n
   for j=1:n
       if i ~= j
           D(i,j)=6370*acos(cos(data_bl(i,1)-data_bl(j,1))*cos(data_bl(i,2))*cos(data_bl(j,2))...
                     +sin(data_bl(i,2))*sin(data_bl(j,2)));
       end
   end
end

%%4. ��ʼ������
m = 100;                            % ��������
alpha = 1;                           % ��Ϣ����Ҫ�̶�����
beta = 5;                             % ����������Ҫ�̶�����
rho = 0.1;                            % ��Ϣ�ػӷ�����
Q = 1;                                 % ��ϵ��
Eta = 1./D;                          % ��������
Tau = ones(n,n);                 % ��Ϣ�ؾ���
Table = zeros(m,n);             % ·����¼��(���ɱ�)
iter = 1;                               % ����������ֵ
iter_max = 800;                  % ���������� 
Route_best = zeros(iter_max,n);                 % �������·��       
Length_best = zeros(iter_max,1);               % �������·���ĳ���  
Length_ave = zeros(iter_max,1);                % ����·����ƽ������ 

%%5. ����Ѱ�����·��
while iter <= iter_max
    %5.1 ��������������ϵ�������
    start = zeros(m,1);
    for i = 1:m
        temp = randperm(n);
        start(i) = temp(1);
    end
    Table(:,1) = start; 
    
    % ������ռ�
    citys_index = 1:n;
      
    %5.2 �������ѡ���·
    for i = 1:m                     % ÿ������·��ѡ��
       for j = 2:n                   % �������·��ѡ��
           tabu = Table(i,1:(j - 1));                          % �ѷ��ʵĳ��м���(���ɱ�)
           allow_index = ~ismember(citys_index,tabu);    % ���ڽ��ɱ�������
           allow = citys_index(allow_index);           % �����ʵĳ��м���
           P = allow;
           % ������м�ת�Ƹ���
           for k = 1:length(allow)
               P(k) = Tau(tabu(end),allow(k))^alpha * Eta(tabu(end),allow(k))^beta;
           end
           P = P/sum(P);
           % �����̶����ѡ����һ�����ʳ���
           Pc = cumsum(P);                             %��Pͬά����������Pc���ۼƸ���    
           target_index = find(Pc >= rand);      %���������rand��Pc�������
           target = allow(target_index(1));       %�������ѡ��ĳ������
           Table(i,j) = target;                            %ѡ��������ż���·����¼��(���ɱ�)
       end
    end
    
    %5.3 ͨ��2-�ڱ��㷨�Ż�·��
    K = 20;                              %20��2-�ڱ��Ż�TSP·��            
    for k = 1:m                        %5.2�õ���m������TSP·��
        c1 = Table(k,:);             %������k������·��
        for t = 1:K                     %2-�ڱ��㷨ѭ�����������K�ָĽ�
            flag = 0;                    %2-�ڱ��㷨ѭ�����˳���־
            for i = 1:98
                for j = i+2:100
                    if D(c1(i),c1(j)) + D(c1(i+1),c1(j+1)) < D(c1(i),c1(i+1)) + D(c1(j),c1(j+1))
                        c1(i+1:j) = c1(j:-1:i+1);    flag = 1;      %�޸�TSP·��                        
                    end
                end
            end
            if flag == 0;    break;    end           %2-�ڱ��㷨���ܼ���·������
        end
        Table(k,:) = c1;                                %����2-�ڱ��㷨�����ĵ�k������·��
    end
    
    %5.4 ����������ϵ�·������
    Length = zeros(m,1);
    for i = 1:m
        Route = Table(i,:);
        for j = 1:(n - 1)
            Length(i) = Length(i) + D(Route(j),Route(j + 1));
        end
        Length(i) = Length(i) + D(Route(n),Route(1));
    end
      
    %5.5 �������·�����뼰ƽ������
    if iter == 1
        [min_Length,min_index] = min(Length);
        Length_best(iter) = min_Length;  
        Length_ave(iter) = mean(Length);
        Route_best(iter,:) = Table(min_index,:);
    else
        [min_Length,min_index] = min(Length);
        Length_best(iter) = min(Length_best(iter - 1),min_Length);     %����iter�κ����·����
        Length_ave(iter) = mean(Length);
        if Length_best(iter) == min_Length
            Route_best(iter,:) = Table(min_index,:);                               %���ε���·��Ϊ���
        else
            Route_best(iter,:) = Route_best((iter-1),:);                           %ǰ��·�����
        end
    end
    
    %5.6 ������Ϣ��
    Delta_Tau = zeros(n,n);
    %5.6.1 ������ϼ���
    for i = 1:m
        % ������м��㣬·����ÿ������Ϣ�ظ���
        for j = 1:(n - 1)
            Delta_Tau(Table(i,j),Table(i,j+1)) = Delta_Tau(Table(i,j),Table(i,j+1)) + Q/Length(i);
        end
        %�յ�->�����Ϣ�ظ���
        Delta_Tau(Table(i,n),Table(i,1)) = Delta_Tau(Table(i,n),Table(i,1)) + Q/Length(i);
    end
    Tau = (1-rho) * Tau + Delta_Tau;
    
    %5.6.2 ����������1�����·����¼��
    iter = iter + 1;
    Table = zeros(m,n);
end
toc

%%6. �����ʾ
[Shortest_Length,index] = min(Length_best);
Shortest_Route = Route_best(index,:);
disp(['TSP��̳���:' num2str(Shortest_Length)]);
disp(['TSP���·��:' num2str([Shortest_Route Shortest_Route(1)])]); 

%%7. ��ͼ
%7.1 ͼ1�����TSP·��
figure(1)
plot([data_bl0(Shortest_Route,1);data_bl0(Shortest_Route(1),1)],...
     [data_bl0(Shortest_Route,2);data_bl0(Shortest_Route(1),2)],'o-');
grid on
for i = 1:size(data_bl0,1)
    text(data_bl0(i,1),data_bl0(i,2),['   ' num2str(i)]);
end
%text(citys(Shortest_Route(1),1),citys(Shortest_Route(1),2),'       ���');
%text(citys(Shortest_Route(end),1),citys(Shortest_Route(end),2),'       �յ�');
xlabel('Ŀ�꾭��')
ylabel('Ŀ��γ��')
title(['��Ⱥ�㷨-TSP��·(��̳���:' num2str(Shortest_Length) ')'])

%7.2 ͼ2��ÿ�ֵ���TSPƽ��·��
figure(2)
plot(1:iter_max,Length_best,'b',1:iter_max,Length_ave,'r:')
legend('��ǰ���TSP����','ÿ��ƽ��TSP����')
xlabel('��������')
ylabel('TSP����')
title('TSP��̳�����ÿ�ֵ���TSPƽ������ͼ')

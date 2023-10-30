% �Ŵ��㷨���TSP
% ���ݡ��� Python ѧ��ѧ����12�£��� TSP �����������롢�Ŵ��Ľ��桢�����ѡ�����

clc;    clear all;
tic;      % ��ʱ��ʼ

% 1.����Ԥ�������������
data0=load('TSP_Data.txt');           % ����100��Ŀ��ľ�γ������
x=data0(:,1:2:8);    x=x(:);
y=data0(:,2:2:8);    y=y(:);
sj=[x y];    d1=[70,40];
sj=[d1;sj;d1];    sj=sj*pi/180;      % ��λ���ɻ���
d=zeros(102);                            % �������d�ĳ�ʼֵ
for i=1:101
    for j=i+1:102
        d(i,j)=6370*acos(cos(sj(i,1)-sj(j,1))*cos(sj(i,2))*cos(sj(j,2))+sin(sj(i,2))*sin(sj(j,2)));
    end
end
d=d+d';    

% �����Ŵ��㷨����
w=1600;            % w��Ⱥ��ģ
g=50;                 % g��������      
rand('state',sum(clock));         % ��ʼ�������������
C_best = zeros(g,102);          % ÿ����Ⱥ������·������
long_best = zeros(1,g);          % ÿ����Ⱥ����С��������
long_ave = zeros(1,g);           % ÿ����Ⱥ��ƽ����������

% 2.�����ʼ��ȺJ��J��w�У�ÿ�б�ʾһ��Ⱦɫ��
% �������TSP˳��(Ⱦɫ��)����3��2-�ڱ߷��Ľ�Ⱦɫ��
for k=1:w                                % ����w����ʼ�⣬��Ϊ��ʼ��Ⱥ
    c=randperm(100);              % ����1��...��100��һ��ȫ����
    c1=[1,c+1,102];                  % ����һ�����Ľ���ʼ��
    for t=1:3                             % ÿ����ʼ�����3��2-�ڱ߸Ľ�
        flag=0;
        for m=1:100
            for n=m+2:101
                if d(c1(m),c1(n))+d(c1(m+1),c1(n+1))<d(c1(m),c1(m+1))+d(c1(n),c1(n+1))
                    c1(m+1:n)=c1(n:-1:m+1);
                    J(k,:)=c1;            % �Ľ���k��Ⱦɫ���TSP·��˳��
                    flag=1;
                end
            end
        end
        if flag==0                        % ���ܸĽ���ǰ��ʼ�⣬���˳���ǰ��ʼ��
            J(k,:)=c1;     break
        end
    end
end

% 3.�Ŵ��㷨����ѭ����
for k=1:g
    % 3.1.���ú���crossover.m������Ⱦɫ�嵥�㽻��������õ�Ⱦɫ������A
    for i = 1:w
        parent = randsample(1:w,2);       % ѡ����Ⱥ������Ⱦɫ����ţ�׼�����н������
        pA = parent(1);
        pB = parent(2);
        A(i,:) = crossover(J,pA,pB);
    end

    % 3.2.���ú���mutate.m��ִ��Ⱦɫ�����������õ�Ⱦɫ������B
    for j = 1:w
        B(j,:)=mutate(d,J(j,:));                 % ���3-�ڱߣ�...��10-�ڱ��㷨���죬ѡ��TSP��·��СȾɫ��
    end

    % 3.3.��J��A��B�У�ѡ��TSP��·��̵�w��Ⱦɫ�壬������һ����Ⱥ
    G=[J;A;B];
    [N1,N2]=size(G);    
    T_long=zeros(N1,1);
    for ii = 1:N1                                    % ����ÿ��Ⱦɫ���TSP��·����
        for jj = 1:N2-1
            T_long(ii) = T_long(ii)+d(G(ii,jj),G(ii,jj+1));
        end
    end
    [T_min_long,ind]=sort(T_long);
    J=G(ind(1:w),:);                             % w���������Ⱦɫ�壬������һ����ȺJ
    C_best(k,:)=J(1,:);    
    long_best(k)=T_min_long(1);
    long_ave(k)=mean(T_min_long(1:w));
end

% 4.�Ŵ��㷨��������
disp('ÿһ����СTSP��·����:');     disp(long_best);
c1=G(ind(1),:);
disp('GA�㷨����TSP��·:');     disp(c1); 
flong=T_min_long(1);
disp('GA�㷨����TSP��·����:');     disp(flong); 

toc

% ͼ�λ���
xx=sj(c1,1);    yy=sj(c1,2);
figure(1)
plot(xx,yy,'-o')
title('ͼ1-GA�㷨��TSP��··��')

figure(2)
plot(1:g,long_best,'b',1:g,long_ave,'r')
legend('��ǰ���TSP����','ÿ��ƽ��TSP����')
xlabel('��Ⱥ����')
ylabel('TSP����')
title('TSPÿ�ֵ�����̳�����ƽ������ͼ')

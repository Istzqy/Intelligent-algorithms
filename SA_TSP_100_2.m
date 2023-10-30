% ģ���˻��㷨���TSP
% 2-�ڱ߷��޸�TSP��·����ƽ���������L=10

clc,    clear all
tic      %��ʱ��ʼ

e=0.1^20;    L=10;    at=0.9995;    T=10^3;  %��ʼ���˻��㷨����

data0=load('TSP_Data.txt');    %����100��Ŀ������ݣ����ݰ��ձ���е�λ�ñ����ڴ��ı��ļ�TSP_Data.txt��
x=data0(:,[1:2:8]);    x=x(:);               %Ŀ�����ݾ�������
y=data0(:,[2:2:8]);    y=y(:);               %Ŀ������γ������
sj=[x y];
d1=[70,40];                                       %���ؾ�γ��
sj=[d1;sj;d1];    sj=sj*pi/180;             %�ǶȻ��ɻ���
d=zeros(102);                                  %�������d��ʼ��

%����������
for i=1:101
   for j=i+1:102
d(i,j)=6370*acos(cos(sj(i,1)-sj(j,1))*cos(sj(i,2))*cos(sj(j,2))+sin(sj(i,2))*sin(sj(j,2)));
   end
end
d=d+d';

%��10000�����TSP��·�У�ѡ����õ���Ϊ��ʼ��·
path=[];    long=inf;                             %Ѳ��·�������ȳ�ʼ��
rand('state',    sum(clock));                 %��ʼ�������������
for j=1:10000                                      %��10000��������У���Ϻõĳ�ʼ��
    path0=[1 1+randperm(100),102];    temp=0;
    for i=1:101
        temp=temp+d(path0(i),path0(i+1));
    end
    if temp<long
        path=path0;    long=temp;
    end
end

% ģ���˻��㷨(SA)��ѭ��������2-�ڱ߷��޸Ļ�·
while T>e
    for k=1:L                                           % �˻�������¶�Tʱ������L�ε�����������ƽ��
        c=randi([2,101],1,2);                     % ���ȷ��2-�ڱߵĵ�
        %c=2+floor(100*rand(1,2));
        c=sort(c);    c1=c(1);    c2=c(2);

        % 2-�ڱ߷��任�󣬻�·��������
        df=d(path(c1-1),path(c2))+d(path(c1),path(c2+1))...
             -d(path(c1-1),path(c1))-d(path(c2),path(c2+1));
        if df<0                              % df<0����·����
            path=[path(1:c1-1),path( c2:-1:c1),path(c2+1:102)]; long=long+df;
        elseif rand<=exp(-df/T)    % df>0���ҷ��ϸ�����������·����
            path=[path(1:c1-1),path(c2:-1:c1),path(c2+1:102)]; long=long+df;
        end
    end
    T=T*at;                                %�����¶�
end

path,    long                             %���Ѳ��·����·������
toc               %��ʱ����
%save('path_star.mat', 'path');
%save('path_long.txt','path','long','-ascii');

%����Ѳ��·��
xx=sj(path,1);    yy=sj(path,2);    plot(xx,yy,'-*')
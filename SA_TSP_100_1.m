% ģ���˻��㷨���TSP
% ���ѡ��2-�ڱ߷���3-�ڱ߷��޸Ļ�·
% ��ƽ���������L_T=1000

clc,    clear all
tic                                 %��ʱ��ʼ

T=10^3;    L_T=1000;    a=0.995;    T_f=0.1^12;  %��ʼ���˻��㷨����

data0=load('TSP_Data.txt');   %����100��Ŀ��ľ�γ�����ݣ����ݰ��ձ���е�λ�ñ����ڴ��ı��ļ�TSP_Data.txt��
x=data0(:,[1:2:8]);    x=x(:);            %Ŀ�����ݾ�������
y=data0(:,[2:2:8]);    y=y(:);            %Ŀ������γ������
data_bl=[x y];
d1=[70,40];                                       %���ؾ�γ��
data_bl0=[d1;data_bl;d1];    data_bl=data_bl0*pi/180;           %�ǶȻ��ɻ���
D=zeros(102);                                  %�������D��ʼ��

%����������
for i=1:101
   for j=i+1:102
D(i,j)=6370*acos(cos(data_bl(i,1)-data_bl(j,1))*cos(data_bl(i,2))*cos(data_bl(j,2))+sin(data_bl(i,2))*sin(data_bl(j,2)));
   end
end
D=D+D';

%��10000�����TSP��·�У�ѡ����õ���Ϊ��ʼ��·
path=[];    long=inf;                          %Ѳ��·�������ȳ�ʼ��
rand('state',    sum(clock));              %��ʼ�������������
for j=1:10000                                   %��10000��������У���Ϻõĳ�ʼ��
    path0=[1 1+randperm(100),102];    temp=0;
    for i=1:101
        temp=temp+D(path0(i),path0(i+1));
    end
    if temp<long
        path=path0;    long=temp;
    end
end

%ģ���˻���ѭ�����������ѡ��2-�ڱ߷���3-�ڱ߷��޸Ļ�·
while T>T_f                                      %�¶ȴ�������¶ȣ�������ȴ����
    for k=1:L_T                                   %��ǰ�¶�T������L_T�ε�����������ƽ��
        r=rand;
        if (r<0.5)
            %2-�ڽ�����н�
            c=2+floor(100*rand(1,2));     %��������2-�ڽ��㣬ȥ�����½�
            c=sort(c);    c1=c(1);    c2=c(2);

            %������ۺ���ֵ������
            df=D(path(c1-1),path(c2))+D(path(c1),path(c2+1))-D(path(c1-1),path(c1))-D(path(c2),path(c2+1));
        
            if df<0           %����׼��
                path=[path(1:c1-1),path( c2:-1:c1),path(c2+1:102)];    long=0;
                for i=1:101                        %���㵱ǰTSP��·��
                    long=long+D(path(i),path(i+1));
                end
            elseif exp(-df/T)>rand           %��Ȼdf>0�����ϸ���ֵҪ����Ȼ����·��
                path=[path(1:c1-1),path(c2:-1:c1),path(c2+1:102)]; long=0;
                for i=1:101                        %���㵱ǰTSP��·��
                    long=long+D(path(i),path(i+1));
                end
            end
        else            
            %3-�ڽ�����н�
            c=2+floor(100*rand(1,3));     %��������2-�ڽ��㣬ȥ�����½�
            c=sort(c);    c1=c(1);    c2=c(2);    c3=c(3);
            
            %������ۺ���ֵ������
            df=D(path(c1-1),path(c2))+D(path(c3-1),path(c1))+D(path(c2-1),path(c3))-D(path(c1-1),path(c1))-D(path(c2-1),path(c2))-D(path(c3-1),path(c3));
            
            if df<0           %����׼��
                path=[path(1:c1-1),path( c2:c3-1),path(c1:c2-1),path(c3:102)];    long=0;
                for i=1:101                        %���㵱ǰTSP��·��
                    long=long+D(path(i),path(i+1));
                end
            elseif exp(-df/T)>rand          %��Ȼdf>0�����ϸ���ֵҪ����Ȼ����·��
                path=[path(1:c1-1),path( c2:c3-1),path(c1:c2-1),path(c3:102)]; long=0;
                for i=1:101                        %���㵱ǰTSP��·��
                    long=long+D(path(i),path(i+1));
                end
            end
        end
    end
    T=T*a;                                          %�����¶�
end

path,    long                                      %���Ѳ��·����·������
toc                  %��ʱ����
%save('path-star.mat', 'path');
%save('path-long.txt','path','long','-ascii');

%����Ѳ��·��
xx=data_bl0(path,1);    yy=data_bl0(path,2);    plot(xx,yy,'-*')
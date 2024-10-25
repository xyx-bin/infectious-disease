clc;
clear;
global beta1 beta2 gamma N % 定义全局变量
gamma=0.15;beta2=0.0625;N=2000;
load('facet1=poissrnd_facet2=Zipf.mat');
Beta10=[];
%%
Beta1=[0:0.0001:0.005];%给定beta1要取的范围
Infect1=zeros(1,length(Beta1));%每个beta1都会对应一个I，所以需要给I申请一个空间来存放
%通过循环将每个beta1对应的I算出来并放在Infect中
for i=1:length(Beta1)
    beta1=Beta1(i);
    M1=length(uniqueDegree_1_ult);%一维面的长度
    N1=length(uniqueDegree_2_ult);%二维面的长度
    X0=zeros(1,2*M1*N1);
    S0=1950; I0=50;%初始时刻所有易感者和所有感染者的数量
    for k=1:M1 %i代表一维面
        for j=1:N1 %j代表二维面
            X0(2*(k-1)*N1+2*j-1)=S0*Com_num(j,k);
            X0(2*(k-1)*N1+2*j)=I0*Com_num(j,k);
        end
    end
    [t,x]=ode45(@distribution_fun,[0 200],X0,[]);
    sum=0;
    for j=2:2:length(x(1,:))
        sum = sum+x(end,j);
    end
    Infect1(i)=sum;
    fprintf('beta1=%f时完成！\n',Beta1(i));   
end
plot(Beta1,Infect1);
hold on;
Beta1=[0:0.0001:0.005];%给定beta1要取的范围
Infect2=zeros(1,length(Beta1));%每个beta1都会对应一个I，所以需要给I申请一个空间来存放
%通过循环将每个beta1对应的I算出来并放在Infect中
for i=1:length(Beta1)
    beta1=Beta1(i);
    M1=length(uniqueDegree_1_ult);%一维面的长度
    N1=length(uniqueDegree_2_ult);%二维面的长度
    X0=zeros(1,2*M1*N1);
    S0=1500; I0=500;%初始时刻所有易感者和所有感染者的数量
    for k=1:M1 %i代表一维面
        for j=1:N1 %j代表二维面
            X0(2*(k-1)*N1+2*j-1)=S0*Com_num(j,k);
            X0(2*(k-1)*N1+2*j)=I0*Com_num(j,k);
        end
    end
    [t,x]=ode45(@distribution_fun,[0 200],X0,[]);
    sum=0;
    for j=2:2:length(x(1,:))
        sum = sum+x(end,j);
    end
    Infect2(i)=sum;
    fprintf('beta1=%f时完成！\n',Beta1(i));   
end
plot(Beta1,Infect2);
hold on;
if Infect1(1) > 100
    beta10 = 0;
else
    for i=2:length(Infect1)
        if Infect1(i)-Infect1(i-1) > 20
            beta10 = Beta1(i-1);
            break;
        end
    end
end
if Infect2(1) > 100
    beta11 = 0;
else
    for i=2:length(Infect2)
        if Infect2(i)-Infect2(i-1) > 20
            beta11 = Beta1(i-1);
            break;
        end
    end
end
Beta10(1)=min(beta10,beta11);%%%%%%%调
disp(Beta10);

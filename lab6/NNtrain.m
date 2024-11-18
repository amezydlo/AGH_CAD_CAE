function NNtrain()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Preparation of dataset

%collection of dataset
A = [1/5 1/10 1/30; 1/10 2/15 1/10;  1/30 1/10 1/5];
i=1;
for n=0.01:0.01:0.5
rhs= [  (pi*pi*n*n+2*cos(pi*n)-2)/(pi*pi*pi*n*n*n);
(-2*pi*n*sin(pi*n)-4*cos(pi*n)+4)/(pi*pi*pi*n*n*n); 
((2-pi*pi*n*n)*cos(pi*n)+2*pi*n*sin(pi*n)-2)/(pi*pi*pi*n*n*n) ];
A;
rhs;
u=A\rhs;
dataset_in(i)=n; 
dataset_u1(i)=u(1); 
dataset_u2(i)=u(2); 
dataset_u3(i)=u(3); 
i=i+1;
end
ndataset=i-1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%training
a1 = 1.0; b1 = 1.0; c1 = 1.0; d1 = 1.0;
a2 = 1.0; b2 = 1.0; c2 = 1.0; d2 = 1.0;
a3 = 1.0; b3 = 1.0; c3 = 1.0; d3 = 1.0;

% a1 = 0.27; a2=2.78; a3=3.11;
% b1=-0.5072; b2=  -0.612; b3= -0.63353;
% c1= 0.10709; c2= 2.541; c3= 2.8908;
% d1=-0.047701; d2= -0.94432; d3= -0.89337;


eta1=0.1;
eta2=0.1;
eta3=0.1;

r = 0 + (1-0).*rand(ndataset,1);
r=r.*ndataset;


num_epochs = 100; % Liczba epok
for epoch = 1:num_epochs
    disp("Epoch: " + epoch);
    % Istniejąca pętla dla zbioru danych
    for j = 1:ndataset
        i = floor(r(j));
        if i == 0
            i = 1;
        end
     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %1
  a1_=a1; b1_=b1; c1_=c1; d1_=d1;
  %evaluation of the first NN modeling a1
  eval1 = c1*1.0/(1.0+exp(-(a1*dataset_in(i)+b1)))+d1;
  error1 = 0.5*(eval1-dataset_u1(i))^2;
  %a1=a1-eta*derivative of error with respect to a
  %d error/da = d/da [[1/2*(c sigma (ax+b)+d-y(x))^2]=cxe^(-ax-b)(eval1-y(x))/(e^{-ax-b)+1)^2
  derrorda = c1*dataset_in(i)*exp(-a1*dataset_in(i)-b1)*(eval1-dataset_u1(i))/(exp(-a1*dataset_in(i)-b1)+1)^2;
  a1=a1-eta1*  derrorda;
  %b1=b1-eta*derivative of error with respect to b
  %d error/db = d/db [1/2*(c sigma (ax+b)+d-y(x))^2]=ce^(-ax-b)(eval1-y(x))/(e^{-ax-b)+1)^2
  derrordb = c1*exp(-a1_*dataset_in(i)-b1)*(eval1-dataset_u1(i))/(exp(-a1_*dataset_in(i)-b1)+1)^2;
  b1=b1-eta1*  derrordb;
  %c1=c1-eta*derivative of error with respect to c
  %d error/dc = d/dc [1/2*(c sigma (ax+b)+d-y(x))^2]=(eval1-y(x))/(e^{-ax-b)+1)
  derrordc = (eval1-dataset_u1(i))/(exp(-a1*dataset_in(i)-b1)+1);
  c1=c1-eta1*  derrordc;
  %d1=d1-eta*derivative of error with respect to d
  %d error/dd = d/dd [1/2*(c sigma (ax+b)+d-y(x))^2]=eval1-y
  derrordd = (eval1-dataset_u1(i));
  d1=d1-eta1*  derrordd;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %2
  a2_=a1; b2_=b2; c2_=c2; d2_=d2;
  %evaluation of the first NN modeling a2  
  eval2 = c2*1.0/(1.0+exp(-(a2*dataset_in(i)+b2)))+d2;
  error2 = 0.5*(eval2-dataset_u2(i))^2;
  %a2=a2-eta*derivative of error with respect to a
  %d error/da = d/da [[1/2*(c sigma (ax+b)+d-y(x))^2]=cxe^(-ax-b)(eval2-y(x))/(e^{-ax-b)+1)^2
  derrorda = c2*dataset_in(i)*exp(-a2*dataset_in(i)-b2)*(eval2-dataset_u2(i))/(exp(-a2*dataset_in(i)-b2)+1)^2;
  a2=a2-eta2*  derrorda;
  %b2=b2-eta*derivative of error with respect to b
  %d error/db = d/db [1/2*(c sigma (ax+b)+d-y(x))^2]=ce^(-ax-b)(eval2-y(x))/(e^{-ax-b)+1)^2
  derrordb = c2*exp(-a2_*dataset_in(i)-b2)*(eval2-dataset_u2(i))/(exp(-a2_*dataset_in(i)-b2)+1)^2;
  b2=b2-eta2*  derrordb;
  %c2=c2-eta*derivative of error with respect to c
  %d error/dc = d/dc [1/2*(c sigma (ax+b)+d-y(x))^2]=(eval2-y(x))/(e^{-ax-b)+1)
  derrordc = (eval2-dataset_u2(i))/(exp(-a2*dataset_in(i)-b2)+1);
  c2=c2-eta2*  derrordc;
  %d2=d2-eta*derivative of error with respect to d
  %d error/dd = d/dd [1/2*(c sigma (ax+b)+d-y(x))^2]=
  derrordd = (eval2-dataset_u2(i));
  d2=d2-eta2*  derrordd;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  %3
  a3_=a1; b3_=b2; c3_=c2; d3_=d2;
  %evaluation of the first NN modeling a3 
  eval3 = c3*1.0/(1.0+exp(-(a3*dataset_in(i)+b3)))+d3;
  error3 = 0.5*(eval3-dataset_u3(i))^2;
  %a3=a3-eta*derivative of error with respect to a
  %d error/da = d/da [[1/2*(c sigma (ax+b)+d-y(x))^2]=cxe^(-ax-b)(eval3-y(x))/(e^{-ax-b)+1)^2
  derrorda = c3*dataset_in(i)*exp(-a3*dataset_in(i)-b3)*(eval3-dataset_u3(i))/(exp(-a3*dataset_in(i)-b3)+1)^2;
  a3=a3-eta3*  derrorda;
  %b3=b3-eta*derivative of error with respect to b
  %d error/db = d/db [1/2*(c sigma (ax+b)+d-y(x))^2]=ce^(-ax-b)(eval3-y(x))/(e^{-ax-b)+1)^2
  derrordb = c3*exp(-a3_*dataset_in(i)-b3)*(eval3-dataset_u3(i))/(exp(-a3_*dataset_in(i)-b3)+1)^2;
  b3=b3-eta3*  derrordb;
  %c3=c3-eta*derivative of error with respect to c
  %d error/dc = d/dc [1/2*(c sigma (ax+b)+d-y(x))^2]=(eval3-y(x))/(e^{-ax-b)+1)
  derrordc = (eval3-dataset_u3(i))/(exp(-a3*dataset_in(i)-b3)+1);
  c3=c3-eta3*  derrordc;
  %d3=d3-eta*derivative of error with respect to d
  %d error/dd = d/dd [1/2*(c sigma (ax+b)+d-y(x))^2]=eval3-y
  derrordd = (eval3-dataset_u3(i));
  d3=d3-eta3*  derrordd;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%store for making plots later
  y1a(j)=a1;
  y1b(j)=b1;
  y1c(j)=c1;
  y1d(j)=d1;
  y2a(j)=a2;
  y2b(j)=b2;
  y2c(j)=c2;
  y2d(j)=d2;
  y3a(j)=a3;
  y3b(j)=b3;
  y3c(j)=c3;
  y3d(j)=d3;
  yeval1(j)=eval1;
  yeval2(j)=eval2;
  yeval3(j)=eval3;
  yu1(j)=dataset_u1(i);
  yu2(j)=dataset_u2(i);
  yu3(j)=dataset_u3(i);
  ye1(j)=error1;
  ye2(j)=error2;
  ye3(j)=error3;
%
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot convergence
x=1:1:ndataset;
hold on

%loglog(x,ye1,x,ye2,x,ye3);
%h=legend('error1','error2','error3');

plot(x,y2a,x,y2b,x,y2c,x,y2d);
h=legend('a2','b2','c2','d2');
set(h,'FontSize',20);

figure

plot(x,yeval1,x,yeval2,x,yeval3,x,yu1,x,yu2,x,yu3);
h=legend('ANN approx of u1','ANN approx of u2','ANN approx of u3','u1','u2','u3');
set(h,'FontSize',20);
set(h,'Location','northeast');
set(gca,'FontSize',20);
%a2
%b2
%c2
%d2


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot error 
%z(x)=sum_{i=1,2,3}ANN(i)*B_{i,2}
%z(x)=ANN1(n)*(1-x)^2 +ANN2(n)*2*x*(1-x)+ANN3(n)*x^2

figure

i=1;
for n=0.01:0.01:0.5
  %evaluated by ANN
  ANN1 = c1*1.0/(1.0+exp(-(a1*n+b1)))+d1;
  ANN2 = c2*1.0/(1.0+exp(-(a2*n+b2)))+d2;
  ANN3 = c3*1.0/(1.0+exp(-(a3*n+b3)))+d3;
%solved by projection  
rhs= [  (pi*pi*n*n+2*cos(pi*n)-2)/(pi*pi*pi*n*n*n);
(-2*pi*n*sin(pi*n)-4*cos(pi*n)+4)/(pi*pi*pi*n*n*n); 
((2-pi*pi*n*n)*cos(pi*n)+2*pi*n*sin(pi*n)-2)/(pi*pi*pi*n*n*n) ];
u=A\rhs;
  %compute error between projection and ANN
  error1(i)=0.0;
  error2(i)=0.0;
  h=0.01;
  for x=0.0:h:1.0
    z=ANN1*(1-x)^2+ANN2*2*x*(1-x)+ANN3*x^2;
    us=u(1)*(1-x)^2+u(2)*2*x*(1-x)+u(3)*x^2;
    s=sin(n*pi*x);
    error1(i)=error1(i)+(z-s)^2*h;
    error2(i)=error2(i)+(z-us)^2*h;
%    error1(i)=error1(i)+(us-s)^2*h;
  end
  i=i+1;
end
n=0.01:0.01:0.5;
%plot(n,error1);
%h=legend('||sum_i=1,2,3 u(i)B_{i,2}(x)- sin(n*pi*x)');
plot(n,error1,n,error2);
h=legend('||sum_i=1,2,3 u(i)B_{i,2}(x)- ANN(i)B_{i,2}(x)||^2_L2','||sin(n*pi*x)-sum_i=1,2,3 ANN(i)B_{i,2}(x)||^2_L2');
set(h,'FontSize',20);
set(h,'Location','northeast');
set(gca,'FontSize',20);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test single sin(n*pi*x)

figure
%test n
n=0.123;
%evaluation of ANN
eval1 = c1*1.0/(1.0+exp(-(a1*n+b1)))+d1;
eval2 = c2*1.0/(1.0+exp(-(a2*n+b2)))+d2;
eval3 = c3*1.0/(1.0+exp(-(a3*n+b3)))+d3;

%computing B-spline values
x=0:0.01:0.5;
y=sin(n*pi.*x);
z=eval1*(1-x).^2+eval2*2*x.*(1-x)+eval3*x.^2;
plot(x,y,x,z);
h=legend('sin(n*pi*x)','u1*B(1,2)(x)+u2*B(2,2)(x)+u3*B(3,3)(x)');
set(h,'FontSize',20);
set(h,'Location','northwest');
set(gca,'FontSize',20);


disp("a1 = " + round(a1, 2) + "; b1 = " + round(b1, 2) + "; c1 = " + round(c1, 2) + "; d1 = " + round(d1, 2))
disp("a2 = " + round(a2, 2) + "; b2 = " + round(b2, 2) + "; c2 = " + round(c2, 2) + "; d2 = " + round(d2, 2))
disp("a3 = " + round(a3, 2) + "; b3 = " + round(b3, 2) + "; c3 = " + round(c3, 2) + "; d3 = " + round(d3, 2))
end

function y=sigma(x)
  y=1.0/(1.0+exp(-x));
  return
end
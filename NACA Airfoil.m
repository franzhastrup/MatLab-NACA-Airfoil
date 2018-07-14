clear all
close all
%Flow for 2D NACA section with camber, thickness and angle of attack
%https://en.wikipedia.org/wiki/NACA_airfoil
%%% User-Defined values:
a = 1; 
tc = 0.06 %thickness ratio, standard = 0.06
fc=0.06651 %camber ratio, standard = 0.06651
alpha=1.4 %Angle of attack
delta_x=1.25 % Percent between each x-point
beta=0 %Flap angle
b=0 %Length of flap, note that c=2a. Please note that b/delta_x has to give an integer!
%%%

U = 1;
C_Li=1.0;
alpha_i=1.4*pi/180;
rho=1.0;%t/m3
alpha=alpha*pi/180;
beta=beta*pi/180;

%xm is taken with the same interval over the whole distance
xm1=0.01*2*a*(-50:delta_x:50);
xm1=transpose(xm1);

%Tau is found from a polynomium 
%The polynomium is not well defined at the edges, so tau is set to 0 at the edges:
Tau1=zeros(length(xm1),1);
for i=1:length(xm1)-2
    Tau1(i+1,1)=a*(-0.0399036788*xm1(i+1)^6 - 0.0311936266*xm1(i+1)^5 + 0.0595367329*xm1(i+1)^4 + 0.0136473772*xm1(i+1)^3 - 0.0750686865*xm1(i+1)^2 + 0.0208886069*xm1(i+1) + 0.0581841773);
end
%Scaling of tau
Tau1=Tau1*tc/0.06; %The thicknes ratio t/c in the table is 0.06 which is scsled with the user defined thickness ratio

%yf from parabel in Excel - more smooth CAMBER
yf1=zeros(length(xm1),1);
for i=1:length(xm1)-2
    yf1(i+1,1)=a*(-0.0427209421*xm1(i+1)^5 - 0.0303017360*xm1(i+1)^4 + 0.0572310577*xm1(i+1)^3 - 0.1025677303*xm1(i+1)^2 - 0.0130122840*xm1(i+1) + 0.1333747330);
end
%Scaling of camber height
yf1=yf1/0.06651*fc; %Where the scaled original f/c value is 0.06651

%Slope of y_f and xm %The slope of xm and yf is calculated from the differentiated polynomium
for i=1:length(xm1)-1
   dyf(i,1)=a*(-5*0.0427209421*xm1(i)^4 - 4*0.0303017360*xm1(i)^3 + 3*0.0572310577*xm1(i)^2 - 2*0.1025677303*xm1(i) - 0.0130122840);
end
dyf=dyf+tan(alpha);
dyf=dyf*a*fc/0.06651; %scaling


%Slope of Tau is found from the integrated polynomium. 
for i=1:length(xm1)
    Tau_d(i,1)=a*(-6*0.0399036788*xm1(i)^5 - 5*0.0311936266*xm1(i)^4 + 4*0.0595367329*xm1(i)^3 + 3*0.0136473772*xm1(i)^2 - 2*0.0750686865*xm1(i) + 0.0208886069);
end

%The angle of attack is added to the xm and yf + Tau vectors
%The thickness Tau has two sides and they both needs to get their
%coordinates moved
Tau_minus=-Tau1;
for i=1:length(xm1)
xm(i,1)=cos(alpha)*xm1(i)-sin(alpha)*yf1(i);
yf(i,1)=sin(alpha)*xm1(i)+cos(alpha)*yf1(i);
%Tau_new(i,1)=sin(alpha)*xm1(i)+cos(alpha)*Tau_new1(i);
%Tau_new_minus(i,1)=sin(alpha)*xm1(i)+cos(alpha)*Tau_new_minus1(i);
end

%Flaps
b_leng=b/(delta_x*2*a*0.01);%No of x-coordinates for flaps
for i=1:b_leng
    xm(i,1)=(cos(beta)*(xm(i)-xm(b_leng))-sin(beta)*(yf(i)-yf(b_leng)))+xm(b_leng);
    yf(i,1)=(sin(beta)*(xm(i)-xm(b_leng))+cos(beta)*(yf(i)-yf(b_leng)))+yf(b_leng);
    dyf(i,1)=dyf(i,1) + tan(beta);
    Tau_d(i,1)=Tau_d(i,1) + tan(beta);
end

%Velocity vektor V_i is defined from the boundary conditions
vi=U*dyf;%BC - finds Gamma

%The thickness is supposed to have the camber as the mean line
Tau_new=Tau1+yf;
Tau_new_minus=yf-Tau1;

x_c=xm(1:end-1); % Control Points are similar with the xm values, removing the front because a vortex has to be at the leading edge
y_v=zeros(length(x_c),1);
%Defines the vortex points as the point in the middle of each pair of
%control points
%Define vortex coordinates
for i=1:(length(x_c))
   x_v(i,1)=xm(i)+(xm(i+1)-xm(i))/2 ;
end
%y-koordinater vortix
for i=1:(length(y_v))
   y_v(i,1)=yf(i)+(yf(i+1)-yf(i))/2 ;
end

%Definining the matrix from the boundary conditions
%see notes
for i=1:length(x_c)
    for j=1:length(x_v)
   coef(i,j)= -1/(2*pi*(x_c(i)-x_v(j)));
    end
end
Gamma=coef\vi;
%xlswrite('tempdata1.xls', Gamma); %Take out Gamme vector into Excel

for i=1:length(xm)-1
gamma(i,1)=Gamma(i)/(xm(i+1)-xm(i)); %strength at each vortex
end

%Vertical velocity u=+-gamma
for i=1:length(gamma)
   u_plus(i,1)=-gamma(i)/2; %upper-side of wing
   u_minus(i,1)=gamma(i)/2; %Lower side of wing
end

Lift=U*rho*sum(Gamma)
Lift_flat=2*pi*alpha
C_L=2*pi*(alpha-alpha_i)+C_Li

%Meshgrid for velocity plot
x=0.01*2*a*(-100:delta_x:100);
x=transpose(x);
y=zeros(length(x),1);
for i = 0:length(x)-1
    y(i+1,1) = -a+(i*2*a/length(x));
end

u_c=zeros(length(y),length(x));
v_c=zeros(length(y),length(x));
%Velocity fields Camber + Angle of attack
for j = 1:length(y)
    for i = 1:length(x)
        for k=1:length(Gamma)
        u_c(j,i)= u_c(j,i)+ (-1)/2*Gamma(k)*(y(j)-y_v(k))/((x(i)-x_v(k))^2+(y(j)-y_v(k))^2*pi);
        v_c(j,i)=v_c(j,i)+1/2*Gamma(k)/((x(i)-x_v(k))*(1+(y(j)-y_v(k))^2/(x(i)-x_v(k))^2)*pi);
        end
        end
end
u_c_plot = u_c;
u_c_plot=u_c_plot-U; %The free stream is subtracted from the u velocity for both the thickness and the camber but when superpositioning the velocities only one free stream velocity should be taken. 
%Because of this the u_c_plot here is still used to plot the camber and
%angle of attack separately from the thickness. u_c will be used for the
%final superposition because it is without the free stream

figure
plot(xm,yf,'red')
xlim([-1.5 1.5])
ylim([-1 1])
hold on
streamslice(x,y,u_c_plot,v_c)
set(gca,'LooseInset',get(gca,'TightInset'))
title('Camber with angle of attack, before superposition','Fontsize',14)
xlabel('x','Fontsize',14)
ylabel('y_f','Fontsize',14)
hold off

figure
plot(xm,yf,'red')
xlim([-1.5 1.5])
ylim([-1 1])
hold on
quiver(x,y,u_c_plot,v_c)
set(gca,'LooseInset',get(gca,'TightInset'))
title('Camber with angle of attack, before superposition','Fontsize',14)
xlabel('x','Fontsize',14)
ylabel('y_f','Fontsize',14)
hold off

%Velocity fields thickness before superposition
for j = 1:length(y)
    for i = 1:length(x)
        u_t_plot(j,i) = -U/pi*trapz(xm1,Tau_d.*(x(i)-xm1)./((x(i)-xm1).^2+y(j)^2))-U;
        v_t_plot(j,i) = -U/pi*trapz(xm1,Tau_d.*y(j)./((x(i)-xm1).^2+y(j).^2));
    end
end



figure
 plot(xm1,Tau1,'red',xm1,-Tau1,'red')
 hold on
streamslice(x,y,u_t_plot,v_t_plot);
set(gca,'LooseInset',get(gca,'TightInset'))
hold off
xlim([-2*a 2*a])
ylim([-a a])
title('Thickness velocity field, before superposition','Fontsize',14)
xlabel('x','Fontsize',14)
ylabel('\tau','Fontsize',14)

figure
 plot(xm1,Tau1,'red',xm1,-Tau1,'red')
 hold on
quiver(x,y,u_t_plot,v_t_plot);
set(gca,'LooseInset',get(gca,'TightInset'))
hold off
xlim([-1.5*a 1.5*a])
ylim([-0.3 0.3])
title('Thickness velocity field, before superposition','Fontsize',14)
xlabel('x','Fontsize',14)
ylabel('\tau','Fontsize',14)
 
%The thickness velocities used for the superposition
dx=delta_x*0.01*2*a;
%Both velocity fields are multiplied with dx, "the distance between each xm
%point" because I am summing over the integral which is over a distance
u_t=zeros(length(x),length(y));
v_t=zeros(length(x),length(y));
for j = 1:length(y)
    for i = 1:length(x)
       for k=1:length(Tau_d)
           u_t(j,i) =u_t(j,i) + (-U*Tau_d(k)*(x(i)-xm(k))./(pi*((x(i)-xm(k)).^2+(y(j)-yf(k)).^2))).*dx;
           v_t(j,i) =v_t(j,i) +  (-U*Tau_d(k).*(y(j)-yf(k))/(pi*((x(i)-xm(k)).^2+(y(j)-yf(k)).^2)))*dx;
        end
    end
end
%Summation of the velocities:
u=u_c+u_t-U;
v=v_c+v_t;

%Test of tangential velocities at thickness surface
x_qt=xm;
y_qt=Tau_new;
u_qt=zeros(length(x_qt),length(y_qt));
v_qt=zeros(length(x_qt),length(y_qt));
for j = 1:length(y_qt)
    for i = 1:length(x_qt)
       for k=1:length(Tau_d)
           u_qt(j,i) =u_qt(j,i) + (-U*Tau_d(k)*(x_qt(i)-xm(k))./(pi*((x_qt(i)-xm(k)).^2+(y_qt(j)-yf(k)).^2))).*dx;
           v_qt(j,i) =v_qt(j,i) +  (-U*Tau_d(k).*(y_qt(j)-yf(k))/(pi*((x_qt(i)-xm(k)).^2+(y_qt(j)-yf(k)).^2)))*dx;
        end
    end
end
u_qt=u_qt-U;
beta2 = atan(Tau_d);
for i=1:length(x_qt)
q_t(i,1)= cos(beta2(i))*u_qt(i,i) - sin(beta2(i))*v_qt(i,i);    
end
%xlswrite('q_t.xls', [x_qt,q_t]); %The tangential velocities are sent to an Excel spreadsheet


figure
 plot(xm,Tau_new,'red',xm,Tau_new_minus,'red',xm,yf,'black')
 hold on
streamslice(x,y,u,v,3);
set(gca,'LooseInset',get(gca,'TightInset'))
hold off
xlim([-1.5 1.5])
ylim([-1 1])
title('Thickness with camber and angle of attack','Fontsize',14)
xlabel('x','Fontsize',14)
ylabel('y_f and \tau','Fontsize',14)

figure
 plot(xm,Tau_new,'red',xm,Tau_new_minus,'red',xm,yf,'black')
 hold on
quiver(x,y,u,v,1);
set(gca,'LooseInset',get(gca,'TightInset'))
hold off
xlim([-1.5*a 1.5*a])
ylim([-1 1])
title('Thickness with camber and angle of attack','Fontsize',14)
xlabel('x','Fontsize',14)
ylabel('y_f and \tau','Fontsize',14)

for i=1:length(gamma)
Delta_C_p(i,1)=2*(u_minus(i)-u_plus(i))/U+(u_plus(i)^2-u_minus(i)^2)/U^2;
end
%xlswrite('Cp_Matlab.xls', Delta_C_p);

figure
plot(x_v,Delta_C_p)
set(gca,'LooseInset',get(gca,'TightInset'))
xlim([-a a])
%ylim([0 2])
title('Pressure distribution for camber zero angle of attack','Fontsize',14)
xlabel('x','Fontsize',14)
ylabel('Delta C_p','Fontsize',14)
     
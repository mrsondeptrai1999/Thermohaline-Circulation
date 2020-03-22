%%% Stommel two-box model
clear 
close all
clf
x = linspace(-0.5*10^6,1.5*10^6,100);
alpha = 2*10^-4;
beta = 0.8*10^-3;
k = 0.5 * 10^10;
delta_t = 20;
s_0 = 35;

% Hysteresis
fimplicit(@(x,y) k*beta*x*s_0 + y.^2 - k*alpha*delta_t*y,...
    [0 1000000 -10*10^6 30*10^6],'k','LineWidth',2)
hold on
xlabel('E') 
ylabel('q') 
fimplicit(@(x,y) k*beta*x*s_0 - y.^2 + k*alpha*delta_t*y,...
    [0 1000000 -10*10^6 0],'LineWidth',2)
line([-5,20],k*alpha*delta_t/2)
hold off
disp('Press enter to continue')
pause

% Evolution
tspan = [0 9.5*10^10];
y0 = [0.5*10^10*2*10^-4*20-0.5*10^10*0.8*10^-3*(0.35-0.1)];
[t_6,y_6] = ode15s(@(t,y) odefun7(t,y,2*10^5,35,0.5*10^10,2*10^-4,...
    0.8*10^-3,20,1,3.2*10^10), tspan, y0);
[t_7,y_7] = ode15s(@(t,y) odefun7(t,y,2*10^5,35,0.5*10^10,2*10^-4,...
    0.8*10^-3,20,2.5,3.2*10^10), tspan, y0);
plot(t_6/31536000,y_6(:)/y_6(1))
hold on
grid on
plot(t_7/31536000,y_7(:)/y_7(1),'LineStyle','--')
ylabel('q(t)/q(0)')
xlabel('time (in years)')
legend('t0 = 1000 y; \Delta H_1 =  1 H_1(0)',...
     't0 = 1000 y;  \Delta H_1 = 2.5 H_1(0)',...
     'Location','east');
hold off
disp('Press enter to continue')
pause


%%% Welander three-box model
clear 
close all
clf
Delta_T = 20;
S_0 = 35;
F = 10^6;
k = 0.5 * 10^4;
alpha = 2*10^-4;
beta = 0.8*10^-3;
D = 0.5;

tspan = [0 9.5*10^10];
y0 = [34.7,34.1];
[t_8,y_8] = ode15s(@(t,y) odefun8(t,y,F,Delta_T,S_0,k,alpha,beta,...
    0.4,3.2*10^10), tspan, y0);
[t_9,y_9] = ode15s(@(t,y) odefun8(t,y,F,Delta_T,S_0,k,alpha,beta,...
    5,3.2*10^10), tspan, y0);
[t_10,y_10] = ode15s(@(t,y) odefun8(t,y,F,Delta_T,S_0,k,alpha,beta,...
    10,3.2*10^10), tspan, y0);
plot(t_8,k*alpha*Delta_T-k*beta*((S_0-y_8(:,1)-y_8(:,2))/2-y_8(:,2)))
hold on
grid on
plot(t_9,k*alpha*Delta_T-k*beta*((S_0-y_9(:,1)-y_9(:,2))/2-y_9(:,1)),...
    'LineStyle','--')
plot(t_10,k*alpha*Delta_T-k*beta*((S_0-y_10(:,1)-y_10(:,2))/2-y_10(:,1)),...
    'LineStyle','--')  
lgd = legend('q_1; for t0 = 1000 y; smallest forcing',...
     'q_1; for t0 = 1000 y;	second highest forcing',...
     'q_1; for t0 = 1000 y;	highest forcing',...
     'Location','east');
hold off
disp('Press enter to continue')
pause 


fimplicit(@(x,y) y^2 - y*k*alpha*Delta_T+2*x)
axis([-100 100 -40 40])
hold on
grid on
fimplicit(@(x,y) -y^2 + y*k*alpha*Delta_T+2*x,[0 100 -40 0],...
    'LineStyle','--')
fimplicit(@(x,y) y^2 - y*k*alpha*Delta_T+2*x)
hold off
disp('Press enter to continue')
pause 


%%% Rooth three-box model
clear 
close all
clf
k = 1.5 * 10^(-6);
alpha = 1.5 * 10^(-4);
beta = 8*10^(-4);
q_0 = 15.5 * 10^6;
tau_3 = 0;
T_3 = 0.3;
lambda = 1.29*10^-9;
H_0 = 0.16*10^(-12);
H_3 = lambda * (tau_3 - T_3);
tau_1 = 0;
T_1 = 2.9;
H_1 = lambda * (tau_1 - T_1);
T_2 = 28.4;
S_1 = 34.7;
S_2 = 35.6;
S_3 = 34.1;
F_3_init = 9*10^-11;
F_1_init = 13.5*10^-11;
q_init = sqrt(k*(alpha*H_3+beta*F_3_init));
q_init_2 = -alpha*k*(T_3-T_1)/2 + sqrt(k^2*alpha^2*(T_3-T_1)^2+4*k*beta...
    *0.3*F_1_init)/2;
scale2 = -(alpha*k*(2.9-0.3)/2 - sqrt(k^2*alpha^2*(2.9-0.3)^2+4*k*beta...
    *(F_3_init+(1-1)*F_1_init/(10/3)))/2);

%Hysteresis Freshwater Forcing
temp_sol_3 = zeros(3,11);
temp_sol_1 = zeros(3,11);

for c = 2:12
    x0 = [0 30 0 34 35 34];
    ratio = 0.3;
    real_x = c/2;
    F_1_eq = real_x*F_1_init;
    F_3_eq = ratio*F_1_eq + F_3_init -ratio*F_1_init;
    x = fsolve(@(x) myfun(x,F_3_eq,F_1_eq),x0);
    temp_sol_3(1,c-1)=x(1); 
    temp_sol_3(2,c-1)=x(2);
    temp_sol_3(3,c-1)=x(3);
end

for c = 2:12
    x0 = [0 30 0 34 35 34];
    ratio = 0.1;
    real_x = c/2;
    F_1_eq = real_x*F_1_init;
    F_3_eq = ratio*F_1_eq + F_3_init -ratio*F_1_init;
    x = fsolve(@(x) myfun(x,F_3_eq,F_1_eq),x0);
    temp_sol_1(1,c-1)=x(1); 
    temp_sol_1(2,c-1)=x(2);
    temp_sol_1(3,c-1)=x(3);
end

T_n_3 = temp_sol_3(2,:)-temp_sol_3(1,:);
T_s_3 = temp_sol_3(2,:)-temp_sol_3(3,:);
T_n_1 = temp_sol_1(2,:)-temp_sol_1(1,:);
T_s_1 = temp_sol_1(2,:)-temp_sol_1(3,:);
q_eq_3 = -alpha*k*(T_n_3-T_s_3)/2 + sqrt(k^2*alpha^2*(T_n_3-T_s_3).^2+...
    4*k*beta*0.3*F_1_init*[2:12]/2)/2;
q_eq_1 = -alpha*k*(T_n_1-T_s_1)/2 + sqrt(k^2*alpha^2*(T_n_1-T_s_1).^2+...
    4*k*beta*0.3*F_1_init*[2:12]/2)/2;
plot([2:12]/2,q_eq_3/q_eq_3(1))

fimplicit(@(x,y) y - sqrt(k*(alpha*H_3+beta*F_3_init))/q_init,'-.b',...
    [0 5.8 -10 5 ])
axis([0 18 -10 5 ])
grid on
hold on
fimplicit(@(x,y) y - sqrt(k*(alpha*H_3+beta*(F_3_init+(x-1)*F_1_init/5)))...
    /q_init,':r','LineWidth',2)
fimplicit(@(x,y) y - sqrt(k*(alpha*H_3+beta*F_3_init))/q_init,...
    [0 2.2 -10 5 ],'-k')
fimplicit(@(x,y) y*scale2 + alpha*k*(2.9-0.3)/2 - sqrt(k^2*alpha^2*...
    (2.9-0.3)^2+4*k*beta*(F_3_init+(x-1)*F_1_init*0.3))/2,[0 16 -10 5 ],...
    'Color','Magenta')
plot([2:12]/2,q_eq_3/q_eq_3(1),'--g')
fimplicit(@(x,y) y*scale2 - alpha*k*(2.9-0.3)/2 + sqrt(k^2*alpha^2*...
    (2.9-0.3)^2+4*k*beta*x*F_1_init)/2,[0 6.2 -10 5 ],'--g')
fimplicit(@(x,y) y*scale2 - alpha*k*(2.9-0.3)/2 + sqrt(k^2*alpha^2*...
    (2.9-0.3)^2+4*k*beta*x*F_1_init)/2,[0 2.2 -10 5 ],'-k')
fimplicit(@(x,y) y*q_init + sqrt(k*(alpha*H_1+beta*x*F_1_init)),...
    ':r','LineWidth',2)
fimplicit(@(x,y) y*q_init + sqrt(k*(alpha*H_1+beta*x*F_1_init)),...
    '-.b',[0 5.8 -10 5 ])
line([2.2 2.2],[1 -2.6],'Color','black')
line([16.5 16.6],[5 -9.5],'Color','red','LineStyle',':','LineWidth',2)
line([6 6.1],[1.7 -5.1],'Color','green','LineStyle','--')
line([5.8 5.8],[1 -2.2],'Color','blue','LineStyle','-.')
line([16 16.1],[4.5 -9.2],'Color','Magenta')
plot(2.2,1,'r*','linewidth',2,'Color','Black')
plot(16.5,5,'r*','linewidth',2,'Color','Black')
plot(6,1.7,'r*','linewidth',2,'Color','Black')
plot(5.8,1,'r*','linewidth',2,'Color','Black')
plot(16,4.5,'r*','linewidth',2,'Color','Black')
ylabel('q^{eq}/q(0)') 
xlabel('F_1^{eq}/F_1(0)')
lgd = legend('\Delta F_3/\Delta F_1 = 0 (Flux BCs)',...
    '\Delta F_3/\Delta F_1 = 0.2 (Flux BCs)',...
    '\Delta F_3/\Delta F_1 = 0 (Mixed BCs)',...
    '\Delta F_3/\Delta F_1 = 0.3 (Mixed BCs non numerical)',...
    '\Delta F_3/\Delta F_1 = 0.3 (Mixed BCs numerical)','Location',...
    'southwest');
lgd.FontSize = 12;
hold off

disp('Press enter to continue')
pause

% Evolution Freshwater Forcing (Mixed BC)
tspan = [0 9.5*10^10];
y0 = [2.9 28.4 0.3 34.7 35.6 34.1 1.45*10^-10];
[t,y] = ode15s(@(t,y) odefun22(t,y,lambda,k,alpha,beta,0.3,4.4,3.2*10^10,...
    F_1_init,F_3_init), tspan, y0);
[t_1,y_1] = ode15s(@(t,y) odefun22(t,y,lambda,k,alpha,beta,0.3,4.8,...
    3.2*10^10,F_1_init,F_3_init), tspan, y0);
plot(t/31536000,y(:,7)/y(1,7))
ylim([-9 3])
xlim([0 3000])
grid on
hold on
plot(t_1/31536000,y_1(:,7)/y(1,7),'LineStyle','--')
ylabel('q(t)/q(0)') 
xlabel('time (in years)')
legend('t0 = 1000 y;\Delta F_3/\Delta F_1 = 0.3;\Delta F_1 = 4.4 F_1(0)',...
     't0 = 1000 y;	\Delta F_3/\Delta F_1 = 0.3;   \Delta F_1 = 4.8 F_1(0)',...
     'Location','southeast');
hold off
disp('Press enter to continue')
pause

% Evolution Freshwater Forcing 2 (ratio > 0.5) % (Mixed BC)
tspan = [0 9.5*10^10];
y0 = [2.9 28.4 0.3 34.7 35.6 34.1 1.45*10^-10];
[t,y] = ode15s(@(t,y) odefun22(t,y,lambda,k,alpha,beta,0.5,4.4,3.2*10^10,...
    F_1_init,F_3_init), tspan, y0);
[t_1,y_1] = ode15s(@(t,y) odefun22(t,y,lambda,k,alpha,beta,0.8,4.4,...
    3.2*10^10,F_1_init,F_3_init), tspan, y0);
[t_a,y_a] = ode15s(@(t,y) odefun22(t,y,lambda,k,alpha,beta,0.8,8,...
    3.2*10^10,F_1_init,F_3_init), tspan, y0);
[t_1_a,y_1_a] = ode15s(@(t,y) odefun22(t,y,lambda,k,alpha,beta,1,4.4,...
    3.2*10^10,F_1_init,F_3_init), tspan, y0);
plot(t/31536000,y(:,7)/y(1,7))
ylim([-2 5])
xlim([0 3000])
grid on
hold on
plot(t_1/31536000,y_1(:,7)/y(1,7),'LineStyle','--')
plot(t_a/31536000,y_a(:,7)/y(1,7),'LineStyle','-.')
plot(t_1_a/31536000,y_1_a(:,7)/y(1,7),'LineStyle',':')
ylabel('q(t)/q(0)') 
xlabel('time (in years)')
legend('t0 = 1000 y;\Delta F_3/\Delta F_1 = 0.5;   \Delta F_1 = 4.4 F_1(0)',...
     't0 = 1000 y;	\Delta F_3/\Delta F_1 = 0.5;   \Delta F_1 = 4.8 F_1(0)',...
     't0 = 1000 y;	\Delta F_3/\Delta F_1 = 0.8;   \Delta F_1 = 4.4 F_1(0)',...
     't0 = 1000 y;	\Delta F_3/\Delta F_1 = 0.8;   \Delta F_1 = 4.8 F_1(0)',...
     'Location','southeast');
hold off 
disp('Press enter to continue')
pause

% Evolution thermal forcing (Mixed BC)
tspan = [0 9.5*10^10];
y0 = [2.9 28.4 0.3 34.7 35.6 34.1 1.47*10^-10];
[t_2,y_2] = ode15s(@(t,y) odefun5(t,y,lambda,k,alpha,beta,0.3,6.35,...
    3.2*10^10,F_1_init,F_3_init), tspan, y0);
[t_3,y_3] = ode15s(@(t,y) odefun5(t,y,lambda,k,alpha,beta,0.3,6.6,...
    3.2*10^10,F_1_init,F_3_init), tspan, y0);
plot(t_2/31536000,y_2(:,7)/y_2(1,7))
hold on 
grid on
plot(t_3/31536000,y_3(:,7)/y_3(1,7),'LineStyle','--')
ylim([-4 1])
xlim([0 3000])
ylabel('q(t)/q(0)')
xlabel('time (in years)')
legend('t0 = 1000 y;\Delta \tau_3/\Delta \tau_1 = 0.3;\Delta \tau_1 = 6.35',...
     't0 = 1000 y;\Delta \tau_3/\Delta \tau_1 = 0.3;   \Delta \tau_1 = 6.6',...
     'Location','southwest');
hold off
disp('Press enter to continue')
pause

% Evolution freshwater flux (flux BC)
tspan = [0 9.5*10^10];
y0 = [2.9 28.4 0.3 34.7 35.6 34.1 1.47*10^-10];
[t_4,y_4] = ode15s(@(t,y) odefun4(t,y,lambda,k,alpha,beta,0.3,13.5,...
    3.2*10^10,F_1_init,F_3_init), tspan, y0);
[t_5,y_5] = ode15s(@(t,y) odefun4(t,y,lambda,k,alpha,beta,0.3,14.5,...
    3.2*10^10,F_1_init,F_3_init), tspan, y0);
plot(t_4/31536000,y_4(:,7)/y_4(1,7))
hold on 
grid on
plot(t_5/31536000,y_5(:,7)/y_5(1,7),'LineStyle','--')
ylim([-9 7])
xlim([0 3000])
ylabel('q(t)/q(0)')
xlabel('time (in years)')
legend('t0 = 1000 y;\Delta F_3/\Delta F_1 = 0.3;\Delta F_1 = 13.5 F_1(0)',...
     't0 = 1000 y;\Delta F_3/\Delta F_1 = 0.3;\Delta F_1 = 14.5 F_1(0)',...
     'Location','southeast');
hold off
disp('Press enter to continue')
pause

% Evolution thermal flux (flux BC)
tspan = [0 9.5*10^10];
y0 = [2.9 28.4 0.3 34.7 35.6 34.1 1.47*10^-10];
[t_4,y_4] = ode15s(@(t,y) odefun6(t,y,lambda,k,alpha,beta,0.3,16,...
    3.2*10^10,F_1_init,F_3_init), tspan, y0);
[t_5,y_5] = ode15s(@(t,y) odefun6(t,y,lambda,k,alpha,beta,0.3,2,...
    3.2*10^10,F_1_init,F_3_init), tspan, y0);
[t_4_a,y_4_a] = ode15s(@(t,y) odefun6(t,y,lambda,k,alpha,beta,0.03,...
    16,3.2*10^10,F_1_init,F_3_init), tspan, y0);
[t_5_a,y_5_a] = ode15s(@(t,y) odefun6(t,y,lambda,k,alpha,beta,0.03,2,...
    3.2*10^10,F_1_init,F_3_init), tspan, y0);
plot(t_4/31536000,y_4(:,7)/y_4(1,7))
hold on
grid on
plot(t_5/31536000,y_5(:,7)/y_5(1,7),'LineStyle','--')
plot(t_4_a/31536000,y_4_a(:,7)/y_4_a(1,7),'LineStyle','-.')
plot(t_5_a/31536000,y_5_a(:,7)/y_5_a(1,7),'LineStyle',':')
ylim([-2 2])
xlim([0 3000])
ylabel('q(t)/q(0)')
xlabel('time (in years)')
legend('t0 = 1000 y;\Delta F_3/\Delta F_1 = 0.3;   \Delta H_1 =  16 H_1(0)',...
     't0 = 1000 y;	\Delta F_3/\Delta F_1 = 0.3;   \Delta H_1 = 3 H_1(0)',...
     't0 = 1000 y;	\Delta F_3/\Delta F_1 = 0.03;   \Delta H_1 = 16 H_1(0)',...
     't0 = 1000 y;	\Delta F_3/\Delta F_1 = 0.03;   \Delta H_1 = 2 H_1(0)',...
     'Location','southeast');
hold off
disp('Press enter to continue')
pause

 
 


function F = myfun(y,F_3_init,F_1_init)
lambda = 1.29*10^-9;
k = 1.5 * 10^(-6);
alpha = 1.5 * 10^(-4);
beta = 8*10^(-4);
F(1) = k*(alpha*(y(3)-y(1))+beta*(y(4)-y(6)))*(y(2)-y(1))*10^9 +...
    lambda*(0-y(1))*10^9;
F(2) = k*(alpha*(y(3)-y(1))+beta*(y(4)-y(6)))*(y(3)-y(2))*10^9/2 +...
    lambda*(30-y(2))*10^9;
F(3) = k*(alpha*(y(3)-y(1))+beta*(y(4)-y(6)))*(y(1)-y(3))*10^9 +...
    lambda*(0-y(3))*10^9;
F(4) = k*(alpha*(y(3)-y(1))+beta*(y(4)-y(6)))*(y(5)-y(4))*10^12 -...
    F_1_init*10^12;
F(5) = k*(alpha*(y(3)-y(1))+beta*(y(4)-y(6)))*(y(6)-y(5))*10^12/2 +...
    (F_1_init+F_3_init)*10^12/2;
F(6) = k*(alpha*(y(3)-y(1))+beta*(y(4)-y(6)))*(y(4)-y(6))*10^12 -...
    F_3_init*10^12 ;
end
function dydt = odefun(t,y,lambda,k,alpha,beta,ratio,t0,F_1_init,F_3_init)
dydt = zeros(6,1);
q = k*(alpha*(y(3)-y(1))+beta*(y(4)-y(6)));
f1 = F_1_init + 5.2*F_1_init*t/t0;
f3 = ratio*f1 + F_3_init - ratio*F_1_init;
dydt(1) = q*(y(2)-y(1)) + lambda*(0-y(1));
dydt(2) = q*(y(3)-y(2))/2 + lambda*(30-y(2));
dydt(3) = q*(y(1)-y(3)) + lambda*(0-y(3));
dydt(4) = q*(y(5)-y(4)) - f1;
dydt(5) = q*(y(6)-y(5))/2 + (f1+f3)/2;
dydt(6) = q*(y(4)-y(6)) - f3;
end
function dydt = odefun2(t,y,lambda,k,alpha,beta,ratio,t0,F_1_init,F_3_init)
dydt = zeros(6,1);
dydt(1) = k*(alpha*(y(3)-y(1))+beta*(y(4)-y(6)))*(y(2)-y(1)) +...
    lambda*(0-y(1));
dydt(2) = k*(alpha*(y(3)-y(1))+beta*(y(4)-y(6)))*(y(3)-y(2))/2 +...
    lambda*(30-y(2));
dydt(3) = k*(alpha*(y(3)-y(1))+beta*(y(4)-y(6)))*(y(1)-y(3)) + ...
    lambda*(0-y(3));
dydt(4) = k*(alpha*(y(3)-y(1))+beta*(y(4)-y(6)))*(y(5)-y(4)) - ...
    F_1_init + 5.8*F_1_init*min(t,t0)/t0;
dydt(5) = k*(alpha*(y(3)-y(1))+beta*(y(4)-y(6)))*(y(6)-y(5))/2 + ...
    (F_1_init + 5.8*F_1_init*min(t,t0)/t0+F_3_init +...
    5.8*F_1_init*min(t,t0)*ratio/t0)/2;
dydt(6) = k*(alpha*(y(3)-y(1))+beta*(y(4)-y(6)))*(y(4)-y(6)) +...
    F_3_init - 5.8*F_1_init*min(t,t0)*ratio/t0;
end
function dydt = odefun3(t,y,lambda,k,alpha,beta,ratio,t0,F_1_init,F_3_init)
dydt = zeros(7,1);
dydt(1) = y(7)*(y(2)-y(1)) + lambda*(0-y(1));
dydt(2) = y(7)*(y(3)-y(2))/2 + lambda*(30-y(2));
dydt(3) = y(7)*(y(1)-y(3)) + lambda*(0-y(3));
dydt(4) = y(7)*(y(5)-y(4)) - F_1_init + 5.2*F_1_init*t/t0;
dydt(5) = y(7)*(y(6)-y(5))/2 + (F_1_init + 5.2*F_1_init*t/t0+...
    F_3_init + 5.2*F_1_init*t*ratio/t0)/2;
dydt(6) = y(7)*(y(4)-y(6)) - F_3_init + 5.2*F_1_init*t*ratio/t0;
dydt(7) = -y(7)^2+k*y(7)*alpha*(y(1)-y(2))-k*y(7)*beta*(y(4)-y(5))+...
    k*alpha*(lambda*(0-y(1))-lambda*(0-y(3)))+k*beta*...
    (F_1_init + 5.2*F_1_init*t/t0-F_3_init + 5.2*F_1_init*t*ratio/t0);
end
function dydt = odefun4(t,y,lambda,k,alpha,beta,ratio,grad,t0,F_1_init,...
    F_3_init)
dydt = zeros(7,1);
dydt(1) = y(7)*(y(2)-y(1)) + lambda*(0-2.9);
dydt(2) = y(7)*(y(3)-y(2))/2 + lambda*(30-28.4);
dydt(3) = y(7)*(y(1)-y(3)) + lambda*(0-0.3);
dydt(4) = y(7)*(y(5)-y(4)) - F_1_init - grad*F_1_init*min(t,t0)/t0;
dydt(5) = y(7)*(y(6)-y(5))/2 + (F_1_init + grad*F_1_init*min(t,t0)/t0 +...
    F_3_init + ratio*grad*F_1_init*min(t,t0)/t0)/2;
dydt(6) = y(7)*(y(4)-y(6)) - F_3_init - ratio*grad*F_1_init*min(t,t0)/t0;
dydt(7) = -y(7)^2+k*y(7)*alpha*(y(1)-y(2))-k*y(7)*beta*(y(4)-y(5))+...
    k*alpha*(lambda*(0-0.3)-lambda*(0-2.9))+...
    k*beta*(F_3_init + ratio*grad*F_1_init*min(t,t0)/t0 ...
    - F_1_init - grad*F_1_init*min(t,t0)/t0);
end
function dydt = odefun5(t,y,lambda,k,alpha,beta,ratio,grad,t0,F_1_init,...
    F_3_init)
dydt = zeros(7,1);
dydt(1) = y(7)*(y(2)-y(1)) + lambda*((0+grad*min(t,t0)/t0)-y(1));
dydt(2) = y(7)*(y(3)-y(2))/2 + lambda*(30-y(2));
dydt(3) = y(7)*(y(1)-y(3)) + lambda*((0+ratio*grad*min(t,t0)/t0)-y(3));
dydt(4) = y(7)*(y(5)-y(4)) - F_1_init ;
dydt(5) = y(7)*(y(6)-y(5))/2 + (F_1_init+F_3_init)/2;
dydt(6) = y(7)*(y(4)-y(6)) - F_3_init ;
dydt(7) = -y(7)^2+k*y(7)*alpha*(y(1)-y(2))-k*y(7)*beta*(y(4)-y(5))+...
    k*alpha*(lambda*((0+ratio*grad*min(t,t0)/t0)-y(3))-...
    lambda*((0+grad*min(t,t0)/t0)-y(1)))+k*beta*(F_3_init-F_1_init);
end
function dydt = odefun6(t,y,lambda,k,alpha,beta,ratio,grad,t0,F_1_init,...
    F_3_init)
dydt = zeros(7,1);
dydt(1) = y(7)*(y(2)-y(1)) + lambda*((0+grad*min(t,t0)/t0)-2.9);
dydt(2) = y(7)*(y(3)-y(2))/2 + lambda*(30-28.4);
dydt(3) = y(7)*(y(1)-y(3)) + lambda*((0+ratio*grad*min(t,t0)/t0)-0.3);
dydt(4) = y(7)*(y(5)-y(4)) - F_1_init ;
dydt(5) = y(7)*(y(6)-y(5))/2 + (F_1_init+F_3_init)/2;
dydt(6) = y(7)*(y(4)-y(6)) - F_3_init ;
dydt(7) = -y(7)^2+k*y(7)*alpha*(y(1)-y(2))-k*y(7)*beta*(y(4)-y(5))+...
    k*alpha*(lambda*((0+ratio*grad*min(t,t0)/t0)-0.3)-...
    lambda*((0+grad*min(t,t0)/t0)-2.9))+k*beta*(F_3_init-F_1_init);
end
function dydt = odefun7(t,y,E,S_0,k,alpha,beta,Delta_T,grad,t0)
dydt = zeros(1,1);
dydt(1) = -2*(y(1).^2 - y(1)*k*alpha*Delta_T + ...
    k*beta*(E+E*grad*min(t,t0)/t0)*S_0);
end

function F = myfun2(x)
f1 = (F_1_init - F_3_init/0.3 + x/0.3);
c = ((3*alpha.^2*k*2.6^2)/(32*beta) -35*f1/4);
F(1) = ((3*alpha.^2*k*2.6^2)/(32*beta) -...
    35*(F_1_init - F_3_init/0.3 + x/0.3)/4)/35 +...
    sqrt(((3*alpha.^2*k*2.6^2)/(32*beta) -...
    35*(F_1_init - F_3_init/0.3 + x/0.3)/4).^2+...
    35*(F_1_init - F_3_init/0.3 + x/0.3)*(3*(alpha.^2/beta)*k*2.6^2-...
    35*(F_1_init - F_3_init/0.3 + x/0.3))/16);
end
function dydt = odefun22(t,y,lambda,k,alpha,beta,ratio,grad,t0,F_1_init,...
    F_3_init)
dydt = zeros(7,1);
dydt(1) = y(7)*(y(2)-y(1)) + lambda*(0-y(1));
dydt(2) = y(7)*(y(3)-y(2))/2 + lambda*(30-y(2));
dydt(3) = y(7)*(y(1)-y(3)) + lambda*(0-y(3));
dydt(4) = y(7)*(y(5)-y(4)) - F_1_init - grad*F_1_init*min(t,t0)/t0;
dydt(5) = y(7)*(y(6)-y(5))/2 + (F_1_init + grad*F_1_init*min(t,t0)/t0 +...
    F_3_init + ratio*grad*F_1_init*min(t,t0)/t0)/2;
dydt(6) = y(7)*(y(4)-y(6)) - F_3_init - ratio*grad*F_1_init*min(t,t0)/t0;
dydt(7) = -y(7)^2+k*y(7)*alpha*(y(1)-y(2))-k*y(7)*beta*(y(4)-y(5))+...
    k*alpha*(lambda*(0-y(3))-lambda*(0-y(1)))+k*beta*(F_3_init +...
    ratio*grad*F_1_init*min(t,t0)/t0 - F_1_init -...
    grad*F_1_init*min(t,t0)/t0);
end
function dydt = odefun8(t,y,F,Delta_T,S_0,k,alpha,beta,grad,t0)
dydt = zeros(2,1);
dydt(1) = -F*(1+grad*min(t,t0)) + abs(k*alpha*Delta_T-k*beta*...
    ((S_0-y(1)-y(2))/2-y(1)))*((S_0-y(1)-y(2))/2-y(1));
dydt(2) = -F*(1+grad*min(t,t0)) + abs(k*alpha*Delta_T-k*beta*...
    ((S_0-y(1)-y(2))/2-y(2)))*((S_0-y(1)-y(2))/2-y(2));
end
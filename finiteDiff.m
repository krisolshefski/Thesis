clear all; close all; clc

%%
%---------------------------%
%     dh/dx [6 Eqns]        %
%---------------------------%

syms h a3 a_3 a2 a_2 a0 a1 a_1
eqn1= a_3                +a_2               +a_1             +a1            +a2              +a3              ==  0;
eqn2= a_3*(-5*h/2)       +a_2*(-3*h/2)      +a_1*(-h/2)      +a1*(h/2)      +a2*(3*h/2)      +a3*(5*h/2)      == -1;
eqn3= a_3*(-5*h/2)^2/2   +a_2*(-3*h/2)^2/2  +a_1*(-h/2)^2/2  +a1*(h/2)^2/2  +a2*(3*h/2)^2/2  +a3*(5*h/2)^2/2  ==  0;
eqn4= a_3*(-5*h/2)^3/6   +a_2*(-3*h/2)^3/6  +a_1*(-h/2)^3/6  +a1*(h/2)^3/6  +a2*(3*h/2)^3/6  +a3*(5*h/2)^3/6  ==  0;
eqn5= a_3*(-5*h/2)^4/24  +a_2*(-3*h/2)^4/24 +a_1*(-h/2)^4/24 +a1*(h/2)^4/24 +a2*(3*h/2)^4/24 +a3*(5*h/2)^4/24 ==  0;
eqn6= a_3*(-5*h/2)^5/120 +a_2*(-3*h/2)^5/120+a_1*(-h/2)^5/120+a1*(h/2)^5/120+a2*(3*h/2)^5/120+a3*(5*h/2)^5/120==  0;

sol = solve( [eqn1, eqn2, eqn3, eqn4, eqn5, eqn6] , [a1, a2,a3,a_1,a_2,a_3]);
h  = 1;
format long
a_3x = eval(sol.a_3);
a_2x = eval(sol.a_2);
a_1x = eval(sol.a_1);
a1x  = eval(sol.a1);
a2x  = eval(sol.a2);
a3x  = eval(sol.a3);
hx5   = -[a_3x; a_2x; a_1x; a1x; a2x; a3x];

%%
%---------------------------%
%    d2h/dx2 [6 Eqns]       %
%---------------------------%

syms a3 a_3 a2 a_2 a0 a1 a_1 h
eqn1= a_3                +a_2               +a_1             +a1            +a2              +a3              ==  0;
eqn2= a_3*(-5*h/2)       +a_2*(-3*h/2)      +a_1*(-h/2)      +a1*(h/2)      +a2*(3*h/2)      +a3*(5*h/2)      ==  0;
eqn3= a_3*(+5*h/2)^2/2   +a_2*(+3*h/2)^2/2  +a_1*(+h/2)^2/2  +a1*(h/2)^2/2  +a2*(3*h/2)^2/2  +a3*(5*h/2)^2/2  == -1;
eqn4= a_3*(-5*h/2)^3/6   +a_2*(-3*h/2)^3/6  +a_1*(-h/2)^3/6  +a1*(h/2)^3/6  +a2*(3*h/2)^3/6  +a3*(5*h/2)^3/6  ==  0;
eqn5= a_3*(+5*h/2)^4/24  +a_2*(+3*h/2)^4/24 +a_1*(+h/2)^4/24 +a1*(h/2)^4/24 +a2*(3*h/2)^4/24 +a3*(5*h/2)^4/24 ==  0;
eqn6= a_3*(-5*h/2)^5/120 +a_2*(-3*h/2)^5/120+a_1*(-h/2)^5/120+a1*(h/2)^5/120+a2*(3*h/2)^5/120+a3*(5*h/2)^5/120==  0;

sol = solve([eqn1, eqn2, eqn3,eqn4, eqn5, eqn6], [a1, a2,a3,a_1,a_2,a_3]);
h  = 1;
format long
a_3xx = eval(sol.a_3);
a_2xx = eval(sol.a_2);
a_1xx = eval(sol.a_1);
a1xx  = eval(sol.a1);
a2xx  = eval(sol.a2);
a3xx  = eval(sol.a3);
hxx5   = -[a_3xx; a_2xx; a_1xx; a1xx; a2xx; a3xx];


x=linspace(0,5,500);
dx=linspace(1,1e-3, length(x));
ana = exp(x);
f   = [ exp(x-5.*dx/2);
        exp(x-3.*dx/2);
        exp(x-   dx/2);
        exp(x+   dx/2);
        exp(x+3.*dx/2);
        exp(x+5.*dx/2) ];
hxxw = [-1/8;7/8;-3/4;-3/4;7/8;-1/8];

Hxx5 =  sum(hxx5.*f)./dx.^2;
Hxxw =  sum(hxxw.*f)./dx.^2;

er5xx = sqrt((Hxx5-ana).^2) ./ sqrt(ana.^2);


%%
% ===============
%     PLOTS
% ===============

figure(1); clf(1)
plot (x,ana,'ko')
hold on
plot (x,Hxx5,'b--')
xlabel('Values of "x"')
ylabel('Values of e^x','interpreter','latex')
legend('analytic','findiff-5th','findiff-2nd','Location','NorthWest')
title('$\frac{d^2h}{dx^2}$','interpreter','latex','fontsize',24)

Hx5 =  sum(hx5.*f)./dx;

er5 = sqrt((Hx5-ana).^2) ./ sqrt(ana.^2);

figure(2); clf(2)
plot (x,ana,'ko')
hold on
plot (x,Hx5,'b--')
xlabel('Values of "x"')
ylabel('Values of exp(x)')
legend('analytic','findiff-5th','findiff-2nd','Location','NorthWest')
title('$\frac{dh}{dx}$','interpreter','latex','fontsize',24)

figure(3); clf(3)
loglog(dx,er5,'-^');
hold on ;
loglog(dx,dx.^(2)*1e-2);
hold on ;
loglog(dx,dx.^(3)*1e2);
hold on ;
loglog(dx,dx.^(4)*1e6);
hold on ;
loglog(dx,dx.^(5)*1e10);
xlabel_h=xlabel({'','dx',''},'fontsize',20,'FontWeight','bold');
ylabel_h=ylabel({'','Error',''},'fontsize',20,'FontWeight','bold');
title('Convergence of $\frac{dh}{dx}$ Error With Mesh Refinement','interpreter','latex','fontsize',18);
legend('Calculated-2nd','2nd Order','3rd Order','4th Order','5th Order', 'Location', 'NorthWest')
set(gca, 'Color', 'none'); % Sets axes background


figure(4); clf(4)
loglog(dx,er5xx,'-^');
hold on ;
loglog(dx,dx.^(2)*1e-2);
hold on ;
loglog(dx,dx.^(3)*1e2);
hold on ;
loglog(dx,dx.^(4)*1e6);
hold on ;
loglog(dx,dx.^(5)*1e10);
xlabel_h=xlabel({'','dx',''},'fontsize',20,'FontWeight','bold');
ylabel_h=ylabel({'','Error',''},'fontsize',20,'FontWeight','bold');
title('Convergence of $\frac{d^2h}{dx^2}$ Error With Mesh Refinement','interpreter','latex','fontsize',18);
legend('Calculated-2nd','2nd Order','3rd Order','4th Order','5th Order', 'Location', 'NorthWest')
set(gca, 'Color', 'none'); % Sets axes background



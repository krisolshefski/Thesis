%function [ hx,hxx ] = weights
% Compute weights for 4th order height function method
clear; clc

syms x z xik zik dx_ dz

% Order of expansion
M=5; N=0;

% Taylor series incluing upto max(M,N) order terms
h=sym(zeros(M+1,N+1));
for n=0:N
    for m=0:M
        h(m+1,n+1)=(x-xik)^m*(z-zik)^n/(factorial(m)*factorial(n));
    end
end

% Allocate Matrices
Z=sym(zeros((M+1)*(N+1),(M+1)*(N+1)));

B=sym(eye  ((M+1)*(N+1),(M+1)*(N+1)));

% Loop over cells in stencil
c=0;
for kp=0;%-3:2
    for ip=-3:2
        c=c+1; % Column number
        %fprintf('%i/6 \n',c)
        % Integrate h over this cell to form height
        H=int(h,x,xik+(ip)*dx_,xik+(ip+1)*dx_)/(dx_);
        % Get coefficents from F
        Z(:,c)=H(:);
    end
end

% Constrained least squares that enforces 2nd-order 
% finite difference opperators and pushes magnitude
%%Z(6,

% Solve for A's (in a vector)
Av=Z\B;

%%   Weights for needed derivatives    %
% ------------------------------------ %
dx_=1;
c=0;
for n=0;%:N
    for m=0:M
        c=c+1;
        if     m==1 && n==0;  hx =eval(reshape(Av(:,c),6,1)); % h_x
        elseif m==0 && n==1;  hz =eval(reshape(Av(:,c),5,5)); % h_z
        elseif m==2 && n==0;  hxx=eval(reshape(Av(:,c),6,1)); % h_xx
        elseif m==0 && n==2;  hzz=eval(reshape(Av(:,c),5,5)); % h_zz
        elseif m==1 && n==1;  hxz=eval(reshape(Av(:,c),5,5)); % h_xz
        end
    end
end
% Test accuracy
figure(1); clf(1)
syms x z
% Function
f=x^5+2*x*cos(x);
% Location to compute derivative
xo=1;
zo=1;

% Compute integral of f for construction of heights
syms z1 z2 x1 x2
intf=matlabFunction(int(int(f,z,z1,z2),x,x1,x2),'Vars',[x1,x2,z1,z2]);

% Loop over derivatives to test
c=0;
for n=0:N % d/dy^n
    for m=0:M % d/dx^m
        c=c+1;
        % Get correct A's for this derivative
        A=reshape(Av(:,c),6,1);
        
        % Grid size to compute heights on
        dxs=[1,1/2,1/4,1/8,1/16,1/32,1/64];
        error=zeros(1,length(dxs));
        for ndx=1:length(dxs);
            dx=dxs(ndx);
            dz=1; %dxs(ndx);
            exact=subs(diff(diff(f,z,n),x,m),[x,z],[xo,zo]);
            % Compute heights
            ht=zeros(6,1);
            for j=0;%-2:2
                for i=-3:2
                    ht(i+4,j+1)=intf(xo+(i)*dx,xo+(i+1)*dx,zo+(j-1/2)*dz,zo+(j+1/2)*dz)/(dx*dz);
                end
            end
            
            % Perform sum(A*heights) (Evaluates A with this dx/dy)
            dx_=dx;
            computed=sum(sum(eval(A).*ht));
            
            % Analysis
            error(ndx)=abs(exact-computed);
            %fprintf('Exact=%5.5f, Computed=%5.5f, Error=%5.5e\n',exact,computed,error(ndx))
        end
%         subplot(N+1,M+1,c)
%         loglog(dxs,error,'-o')
%         hold on
%         loglog(dxs,dxs.^5)
%         title(['d/dx^',num2str(m),' d/dy^',num2str(n)])
%         drawnow 
 if c==2 || c ==3 
        h=figure(c); clf(c)  
        plot(N+1,M+1)
        loglog(dxs,error,'-.','LineWidth',3)
        hold on
        loglog(dxs,dxs.^5,'LineWidth',3)
        title(['d/dx^',num2str(m),' d/dy^',num2str(n)],'interpreter','latex','fontsize',16)
        legend('5th Order Method','5th Order', 'Location', 'NorthWest')
        set(gca, 'Color', 'none','Fontsize',16);
        drawnow 
        basefile = 'figs';
        saveas(h,fullfile(basefile,['WeightsDer',num2str(c),'.png']));
 end
    end
end

%% Convert to Fortran notation (Removing
% hx = [-1/15; -1/9; -1/3;  1/3; 1/9; 1/15];
% hx = [-1/8; -1/8;    0;    0; 1/8; 1/8];
% hxx= [ 1/8;  1/8; -1/4; -1/4; 1/8; 1/8];
dhs={'hx','hxx'};
dx_=1;
dz_=1;
for n=1:length(dhs)
    dh=eval(dhs{n});
    for kp=0;%-2:2
%         for ip=-3:2
%            strs{ip+4}=eval(dh(ip+4,kp+1));
%         end
        fprintf('data dhd%s(:,%+i) / %+5.15e_WP, %+5.15e_WP, %+5.15e_WP, %+5.15e_WP, %+5.15e_WP, %+5.15e_WP / \n',dhs{n}(2:end),kp,dh);
    end
end
%end


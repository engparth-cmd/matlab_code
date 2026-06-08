 
clc;
clear all;
L=1;
error=1;
nume=1;
denom=1;
nx=282;  %  x domain
ny=282;  % y domain
[X,Y] = meshgrid(1:nx,1:ny);
% frame_count=0;
% figh=figure;
% frame_count=0;
t=500000000;  % time steps
max_steps=60000000; % maximum time steps
p=zeros(ny,nx);   % presuure
so=zeros(ny,nx); 
f=zeros(ny,nx,9); % distribution function for phase field
h =zeros(ny,nx,9); % distribution function for electric field
ht=zeros(ny,nx,9); % streaming
h_eq=zeros(ny,nx,9);
g=zeros(ny,nx,9); % distribution function for the fluid flow
gt=zeros(ny,nx,9);
g_eq=zeros(ny,nx,9);
ft=zeros(ny,nx,9);
f_eq=zeros(ny,nx,9);
Fs_x=zeros(ny,nx); 
Fs_y=zeros(ny,nx);
FE_x=zeros(ny,nx);
FE_y=zeros(ny,nx);
FE1=zeros(ny,nx);
FE2=zeros(ny,nx);
FE3=zeros(ny,nx);
FE4=zeros(ny,nx);
uy=zeros(ny,nx);
ux=zeros(ny,nx);
E2=zeros(ny,nx);
phi_old=zeros(ny,nx);
phi_new=zeros(ny,nx);
phi_now=zeros(ny,nx);
ux_old=zeros(ny,nx);
uy_old=zeros(ny,nx);
ux_new=zeros(ny,nx);
uy_new=zeros(ny,nx);
u2=zeros(ny,nx);
rho=zeros(ny,nx);
Ri=zeros(ny,nx,9);
nu=zeros(ny,nx);
mu=zeros(ny,nx);
sigma=zeros(ny,nx);
epsi=zeros(ny,nx);
g_phix=zeros(ny,nx);
g_phiy=zeros(ny,nx);
E_x=zeros(ny,nx);
E_y=zeros(ny,nx);
rho_e=zeros(ny,nx);
z=zeros(ny,nx);
L_phi=zeros(ny,nx);
F=zeros(ny,nx,9);
Fx=zeros(ny,nx);
Fy=zeros(ny,nx);
Fxx=zeros(ny,nx);
Fyy=zeros(ny,nx);
dCu_x=zeros(ny,nx);
dCu_y=zeros(ny,nx);
si=zeros(ny,nx,9);
ex=[0,1,0,-1,0,1,-1,-1,1];
ey=[0,0,1,0,-1,1,1,-1,-1];
exx=[0,1,0,-1,0];
eyy=[0,0,1,0,-1];
wt=[4/9 1/9 1/9 1/9 1/9 1/36 1/36 1/36 1/36];
wwt=[1/3 1/6 1/6 1/6 1/6];
phi=zeros(ny,nx);
phi_h=1;
phi_l=0;
rho_H=1;
rho_L=1;
nu_H=0.1;
nu_L=0.1;
sigma_H=0.5;
sigma_L=0.1;
epsi_H=0.06;
epsi_L=0.001;
W=4.7;
R1=47;
sig=0.001;
kappa=3*sig*W/8;
bita=0.75*sig/W;
max_iter=100000000;
tol=1e-7;
phi_ot=zeros(ny,nx);
%% Electric potential field
phi_e=zeros(ny,nx);  % Electric potential field
gradphi_ex=zeros(ny,nx);
gradphi_ey=zeros(ny,nx);
phi_eold=zeros(ny,nx);
phi_enew=zeros(ny,nx);
%Set Dirichlet boundary conditions
phi_e(1, :)=0;  % top boundary (y = 1)
phi_e(282, :)=28.2;  % bottom boundary (y = ny)
PH=zeros(ny,nx);
%%%%%% New Definition %%%%%
dG=zeros(ny,nx,9);   %%% new
dGt=zeros(ny,nx,9);%% new


% for i=1:ny
% for j=1:nx   
%    phi(i,j)=((phi_h+phi_l)/2)+((phi_h-phi_l)/2)*tanh(2*(R1-sqrt((j-141)^2+(i-141)^2))/W);
% end
% end
for i=1:ny
for j=1:nx   
 phi(i,j)=0.5+0.5*tanh(2*(R1-sqrt((j-141)^2+(i-141)^2))/W);
end
end
for j=1:nx   
phi(1,j)=phi(2,j) ; 
phi(ny,j)=phi(ny-1,j) ;
mu(1,j)   = mu(2,j);
mu(ny,j)  = mu(ny-1,j);
end

for i =1:ny
for j =1:nx
      %rho(i,j)=((rho_H-rho_L)/2)*phi(i,j)+(rho_H+rho_L)/2;
      rho(i,j)=((rho_H-rho_L))*phi(i,j)+rho_L;
      p(i,j)=rho(i,j)/3;
     % nu(i,j)=((nu_H-nu_L)/2)*phi(i,j)+(nu_H+nu_L)/2;
     nu(i,j)=((nu_H-nu_L))*phi(i,j)+nu_L;
     % sigma(i,j)=((sigma_H-sigma_L)/2)*phi(i,j)+(sigma_H+sigma_L)/2;
     % epsi(i,j)=((epsi_H-epsi_L)/2)*phi(i,j)+(epsi_H+epsi_L)/2;
     sigma(i,j)=((sigma_H-sigma_L))*phi(i,j)+sigma_L;
     epsi(i,j)=(epsi_L-epsi_H)*phi(i,j)*phi(i,j)*(2*phi(i,j)-3)+epsi_L;
end
end
% initialization of phase field distribution function
% for i=1:ny
% for j=1:nx
% for k=1:9
% f(i,j,k)=wt(k)*phi(i,j)*(1+3*(ex(k)*ux(i,j)+ey(k)*uy(i,j)));
% end
% end
% end
% %% STEP 2 %%%%%%%% initialization  for phase field  distribution function
for i=1:ny
for j=1:nx
for k=1:9
    if k==1
f(i,j,k)=phi(i,j)-(1-wt(k))*(4*bita*phi(i,j)*(phi(i,j)-0.5)*(phi(i,j)-1)-kappa*L_phi(i,j)-3*phi(i,j)*(phi(i,j)-1)*(epsi_L-epsi_H)*(E_x(i,j)*E_x(i,j)+E_y(i,j)*E_y(i,j)));
    else
f(i,j,k)=wt(k)*(4*bita*phi(i,j)*(phi(i,j)-0.5)*(phi(i,j)-1)-kappa*L_phi(i,j)-3*phi(i,j)*(phi(i,j)-1)*(epsi_L-epsi_H)*(E_x(i,j)*E_x(i,j)+E_y(i,j)*E_y(i,j)));
   end
end
end
end

%% initialization  for potential distribution function
for i=1:ny
    for j=1:nx
        for k=1:9
           if k==1
                h(i,j,k)=(wt(k)-1)*phi_e(i,j);
            else
                h(i,j,k)=wt(k)*phi_e(i,j);   
           end
        end
     end
end
% %% initialization  for hydrodynamic distribution function 

for i=2:ny-1
 for j=1:nx
     for k=1:9
    cu= ex(k)*ux(i,j) + ey(k)*uy(i,j);   
    u2= ux(i,j)^2 + uy(i,j)^2;          
    si(i,j,k)=wt(k)*(3*cu +4.5*(cu^2)-1.5*u2);
      end
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%distribution function for fluid flow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=2:ny-1
    for j=1:nx
        for k=1:9
           if k==1
                g(i,j,k)=(wt(k)-1)*p(i,j)*3+rho(i,j)*si(i,j,k);
            else
                g(i,j,k)=3*wt(k)*p(i,j)+rho(i,j)*si(i,j,k);   
           end
        end
     end
end

for i=1:ny
for j=1:nx
phi_eold(i,j) = phi_e(i,j);
phi_new(i,j) =phi(i,j);
phi_old(i,j) =phi(i,j);
phi_ot(i,j) =phi(i,j);
ux_new(i,j)=ux(i,j);
ux_old(i,j)=ux(i,j);
end
end

% %%%%%calculation of hydrodynamic pressure %%%%%%%%%%
for i =2:ny-1
 for j =1:nx
     for k = 1
    cu= ex(k)*ux(i,j) + ey(k)*uy(i,j);   
    u2= ux(i,j)^2 + uy(i,j)^2;          
    so(i,j) = wt(k)*(3*cu +4.5*(cu^2)-1.5*u2);
      end
  end
end

%% defining the gradient phi and laplace phi

for i=2:ny-1
    for j=1:nx
        tempx=0;
        tempy=0;
         for k=2:9 
             ia=i+ey(k);
             ja=j+ex(k);

            if ja>nx
              ja=1;
            elseif ja<1
              ja=nx;
            end
             tempx= tempx+ex(k)*wt(k)*phi(ia,ja);
             tempy=tempy+ey(k)*wt(k)*phi(ia,ja);
         end
g_phix(i,j)=tempx*3;
g_phiy(i,j)=tempy*3;

    end
end

for i=2:ny-1
    for j=1:nx
      p(i,j)=0;
        for k=2:9
    p(i,j)=p(i,j) + g(i,j,k);
        end
      p(i,j)= 0.6*((p(i,j)+0.5*0.5*(rho_H-rho_L)*(ux(i,j)*g_phix(i,j)+uy(i,j)*g_phiy(i,j))+rho(i,j)*so(i,j)));
   end 
end

% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%loop starts from here
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for st=1:10000000000
    fprintf('st = %d\n', st);
    % pause(0.2)
if (error > 0.0000001)
for step=1:1000
    % pause(0.5)
   fprintf('step = %d\n', step);
   fprintf('error = %d\n', error);
% for i = 2:ny-1
%     for j =1:nx
%         dCu_x(i,j) = (phi_new(i,j)*ux_new(i,j)-phi_old(i,j)*ux_old(i,j));
%         dCu_y(i,j) = (phi_new(i,j)*uy_new(i,j)-phi_old(i,j)*uy_old(i,j));
% 
%     end
% end
% %%%%%%%%%%%%neumann boundary condition for phase field
% for i=2:ny-1
% for j=1:nx   
% phi(1,j)=phi(2,j) ; 
% phi(ny,j)=phi(ny-1,j) ;
% mu(1,j)=mu(2,j);
% mu(ny,j)= mu(ny-1,j);
% end
% end
for i = 2:ny-1
    for j =1:nx
        for k=1:9
        dG(i,j,k) = wt(k)*(g_phix(i,j)*ux_new(i,j)+g_phiy(i,j)*uy_new(i,j))*(1.5*(ex(k)*ex(k)+ey(k)*ey(k))-2);
        dGt(i,j,k) = (wt(k)*(g_phix(i,j)*ux_new(i,j)+g_phiy(i,j)*uy_new(i,j))*(1.5*(ex(k)*ex(k)+ey(k)*ey(k))-2))-(wt(k)*(g_phix(i,j)*ux_old(i,j)+g_phiy(i,j)*uy_old(i,j))*(1.5*(ex(k)*ex(k)+ey(k)*ey(k))-2));
        end
    end
end
for i = 2:ny-1
    for j =1:nx
        phi_old(i,j)=phi(i,j);
        ux_old(i,j)=ux(i,j);
        uy_old(i,j)=uy(i,j);
    end
end
%%% defining the gradient phi and laplace phi

for i=2:ny-1
    for j=1:nx
        tempx=0;
        tempy=0;
         for k=2:9 
             ia=i+ey(k);
             ja=j+ex(k);

            if ja>nx
              ja=1;
            elseif ja<1
              ja=nx;
            end
             tempx= tempx+ex(k)*wt(k)*phi(ia,ja);
             tempy=tempy+ey(k)*wt(k)*phi(ia,ja);
         end
g_phix(i,j)=tempx*3;
g_phiy(i,j)=tempy*3;

    end
end
% % calculation of normal

for i=2:ny-1
    for j=1:nx

           z(i,j)=sqrt(g_phix(i,j)^2+g_phiy(i,j)^2)+0.00000000001;

    end
end
% % % % calculation of laplace phi
for i=2:ny-1
    for j=1:nx
        tempx1=0;
         for k=2:9 
             ia=i+ey(k);
             ja=j+ex(k);

            if ja>nx
              ja=1;
            elseif ja<1
              ja=nx;
            end
             tempx1= tempx1 + wt(k)*(phi(ia,ja)-phi(i,j));

         end
L_phi(i,j)=tempx1*6;
    end
end
%%%%% calculation of force for phase field source term %%%%

% for i=2:ny-1
%     for j=1:nx
% 
%         Fx(i,j) =((1/(3*W))*(g_phix(i,j)/z(i,j))*2*(1-phi(i,j)^2))+dCu_x(i,j);
%         Fy(i,j) =((1/(3*W))*(g_phiy(i,j)/z(i,j))*2*(1-phi(i,j)^2))+dCu_y(i,j);
% 
%         for k = 1:9
%             F(i,j,k) =(3*wt(k))*(ex(k)*Fx(i,j)+ey(k)*Fy(i,j));
%         end
%     end
% end
% for i=2:ny-1
%     for j=1:nx
%         for k = 1:9
%           f_eq(i,j,k) =wt(k)*phi(i,j)*(1+3*(ex(k)*ux(i,j)+ey(k)*uy(i,j)));
%           ft(i,j,k)=-1.25*((f(i,j,k)-f_eq(i,j,k)))+f(i,j,k)+0.375*F(i,j,k);
%         end
%     end
% end

% %%% calculation of the equilibrium distribution function for Cahn-Hilliard
for i=1:ny
for j=1:nx
for k=1:9
    if k==1
f_eq(i,j,k)=phi(i,j)-(1-wt(k))*(4*bita*phi(i,j)*(phi(i,j)-0.5)*(phi(i,j)-1)-kappa*L_phi(i,j)-3*phi(i,j)*(phi(i,j)-1)*(epsi_L-epsi_H)*(E_x(i,j)*E_x(i,j)+E_y(i,j)*E_y(i,j)));
    else
f_eq(i,j,k)=wt(k)*(4*bita*phi(i,j)*(phi(i,j)-0.5)*(phi(i,j)-1)-kappa*L_phi(i,j)-3*phi(i,j)*(phi(i,j)-1)*(epsi_L-epsi_H)*(E_x(i,j)*E_x(i,j)+E_y(i,j)*E_y(i,j)));
   end
end
end
end
% %%%%% calculation of force for phase field source term %%%%
for i=2:ny-1
    for j=1:nx
        for k = 1:9
          ft(i,j,k)=-1.25*((f(i,j,k)-f_eq(i,j,k)))+f(i,j,k)+dG(i,j,k)+0.5*dGt(i,j,k);
        end
    end
end

% streaming of post collision particle distribution

for i=2:ny-1
    for j=1:nx
       for k=1:9
        ia=i+ey(k);
        ja=j+ex(k);

        if ja>nx
          ja=1;
        elseif ja<1
          ja=nx;
        end
         f(ia,ja,k)=ft(i,j,k);  

        end  
    end
end 


% % Boundary condition
for j=1:nx
for i=2
        f(i,j,6)=ft(i,j,8);
         f(i,j,3)=ft(i,j,5);
          f(i,j,7)=ft(i,j,9);

end
for i=ny-1

        f(i,j,9)=ft(i,j,7);
        f(i,j,5)=ft(i,j,3);
        f(i,j,8)=ft(i,j,6);

end
end

%     % calculation of phi
for i=2:ny-1
    for j=1:nx
        PH(i,j)=0;
        for k=1:9
            PH(i,j)=PH(i,j)+f(i,j,k);

        end
        phi(i,j)= PH(i,j);
    end
end
%%%%%%%%%%%neumann boundary condition for phase field
% for i=1:ny
% for j=1:nx   
% phi(1,j)=phi(2,j) ; 
% phi(ny,j)=phi(ny-1,j) ;
% mu(1,j)= mu(2,j);
% mu(ny,j)= mu(ny-1,j);
% end
% end
% %%% calculation of equation 19
% 
% for i = 1:ny
% for j = 1:nx
% 
%     rho(i,j)=((rho_H-rho_L)/2)*phi(i,j)+(rho_H+rho_L)/2;
%     nu(i,j)=((nu_H-nu_L)/2)*phi(i,j)+(nu_H+nu_L)/2;
%     sigma(i,j)=((sigma_H-sigma_L)/2)*phi(i,j)+(sigma_H+sigma_L)/2;
%     epsi(i,j)=((epsi_H-epsi_L)/2)*phi(i,j)+(epsi_H+epsi_L)/2;
% end
% end

for i =1:ny
for j =1:nx
     rho(i,j)=((rho_H-rho_L))*phi(i,j)+rho_L;
     nu(i,j)=((nu_H-nu_L))*phi(i,j)+nu_L;
     sigma(i,j)=((sigma_H-sigma_L))*phi(i,j)+sigma_L;
     epsi(i,j)=(epsi_L-epsi_H)*phi(i,j)*phi(i,j)*(2*phi(i,j)-3)+epsi_L;
end
end

%%% defining the gradient  new phi and new laplace phi

for i=2:ny-1
    for j=1:nx
        tempx=0;
        tempy=0;
         for k=2:9 
             ia=i+ey(k);
             ja=j+ex(k);

            if ja>nx
              ja=1;
            elseif ja<1
              ja=nx;
            end
             tempx= tempx+ex(k)*wt(k)*phi(ia,ja);
             tempy=tempy+ey(k)*wt(k)*phi(ia,ja);
         end
g_phix(i,j)=tempx*3;
g_phiy(i,j)=tempy*3;

    end
end
% % calculation of new normal

for i=2:ny-1
    for j=1:nx

           z(i,j)=sqrt(g_phix(i,j)^2+g_phiy(i,j)^2)+0.00000000001;

    end
end
% % % % calculation of new laplace phi
for i=2:ny-1
    for j=1:nx
        tempx1=0;
         for k=2:9 
             ia=i+ey(k);
             ja=j+ex(k);

            if ja>nx
              ja=1;
            elseif ja<1
              ja=nx;
            end
             tempx1= tempx1 + wt(k)*(phi(ia,ja)-phi(i,j));

         end
L_phi(i,j)=tempx1*6;
    end
end

% % calculation of mu

% for i=2:ny-1
%     for j=1:nx
%  mu(i,j)=4*bita*phi(i,j)*(phi(i,j)-1)*(phi(i,j)+1)-kappa*L_phi(i,j);
%     end
% end
% % % calculation of mu according to Cahn-Hilliard equation
for i=2:ny-1
    for j=1:nx
 mu(i,j)=4*bita*phi(i,j)*(phi(i,j)-0.5)*(phi(i,j)-1)-kappa*L_phi(i,j)-3*phi(i,j)*(phi(i,j)-1)*(epsi_L-epsi_H)*(E_x(i,j)*E_x(i,j)+E_y(i,j)*E_y(i,j));
    end
end
% %calculation of Fs
for i=2:ny-1
    for j=1:nx
            Fs_x(i,j)=mu(i,j)*g_phix(i,j);
            Fs_y(i,j)=mu(i,j)*g_phiy(i,j);
    end
end


% for i=1:ny
% for j=1:nx   
% mu(1,j)=mu(2,j) ; 
% mu(ny,j)=mu(ny-1,j) ;
% end
% end


%calculation of Fs
% for i=2:ny-1
%     for j=1:nx
% 
%             Fs_x(i,j)=mu(i,j)*g_phix(i,j);
%             Fs_y(i,j)=mu(i,j)*g_phiy(i,j);
%     end
% end
%%%%%%%%%%%%%%%%%%%%% ELECTRIC SOLVER  STARTS

% %%%% implementing the equation of electric
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% electric part of model time step loop%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  for l=1:t
    if (L > 0.0000001)
        % pause(0.1)
     fprintf('l = %d\n', l);
     fprintf('L= %d\n', L);
 for i=1:ny
    for j=1:nx
        for k=1:9
           if k==1
                h_eq(i,j,k)=(wt(k)-1)*phi_e(i,j);
            else
                h_eq(i,j,k)=wt(k)*phi_e(i,j);   
           end
        end
     end
end
% %% calculation of distribution function for electric potential

for i=1:ny
    for j=1:nx
        for k=1:9
          tau_h=0.5+3*sigma(i,j);                                   %%%% leaky dielectric model
          ht(i,j,k)=-((h(i,j,k)-h_eq(i,j,k))/tau_h)+h(i,j,k);
        end
    end
end

% %%% streaming of electric potential
% 
% % streaming of post collision particle distribution
% 
for i=2:ny-1
    for j=1:nx
       for k=1:9
        ia=i+ey(k);
        ja=j+ex(k);

        if ja>nx
          ja=1;
        elseif ja<1
          ja=nx;
        end
           h(ia,ja,k)=ht(i,j,k);  

        end  
    end
end 


% % Boundary condition

for j=1:nx
for i=2

         % h(i,j,3)=-ht(1,j,5);

         h(i,j,3)=-ht(1,j,5);

         h(i,j,6)=-ht(1,j,8);

         h(i,j,7)=-ht(1,j,9);
end
for i=ny-1

        % h(i,j,5)=-ht(ny,j,3)+2*wwt(3)*28.2;   
        h(i,j,5)=-ht(ny,j,3)+2*wt(3)*28.2;     
        h(i,j,8)=-ht(ny,j,6)+2*wt(6)*28.2;
        h(i,j,9)=-ht(ny,j,7)+2*wt(7)*28.2;

end
end



% %%% convergence criteria for electric potential
% % calculation of phi_e
for i=2:ny-1
    for j=1:nx
        phi_e(i,j)=0;
        for k=2:9
            phi_e(i,j)=phi_e(i,j)+1.8*h(i,j,k);

        end

    end
end
for i=1:ny
    for j=1:nx
phi_enew(i,j) = phi_e(i,j);
    end
end

% % --- Compute error using loops ---
    numerator=0;
    denominator=0;
    for i=1:ny
        for j=1:nx
            numerator=numerator+abs(phi_enew(i,j)-phi_eold(i,j));
            denominator=denominator+abs(phi_enew(i,j));
        end
    end

    L=numerator/(denominator);
     for    i=2:ny-1
        for j = 1:nx
          phi_eold(i,j) = phi_enew(i,j);
        end
     end


    end 
  end


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          L=1;
%%%%%%Compute E_x and E_y (Electric field):%%%%%%%%%%%%%%%%%%

for i =2:ny-1
    for j =1:nx
        dx=0;
        dy=0;
        for k=2:9
            ia=i+ey(k);
            ja=j+ex(k);
            if ja>nx
                ja=1;
            elseif ja<1
                ja=nx;
            end
            dx=dx+ex(k)*wt(k)*phi_e(ia,ja);
            dy=dy+ey(k)*wt(k)*phi_e(ia,ja);

        end
        E_x(i,j)=-3*dx;
        E_y(i,j)=-3*dy;
    end
end

% %%%%%%%%%%%%%%%%%%%%%% Compute E^2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 2:ny-1
    for j = 1:nx
   E2(i,j)=E_x(i,j)*E_x(i,j)+E_y(i,j)*E_y(i,j);
    end
end


% %%%%%%%%%%%%%%%%%%%%% for dielectric model %%%%%%%%%%%%%%%%%%%%%%


for i=2:ny-1
    for j=1:nx
        tempx=0;
        tempy=0;
         for k=2:9
             ia=i+ey(k);
             ja=j+ex(k);

            if ja>nx
              ja=1;
            elseif ja<1
              ja=nx;
            end
             tempx= tempx+ex(k)*wt(k)*phi_e(ia,ja);
             tempy=tempy+ey(k)*wt(k)*phi_e(ia,ja);
         end
gradphi_ex(i,j)=tempx*3;
gradphi_ey(i,j)=tempy*3;

    end
end

for i=2:ny-1
    for j=1:nx
        div=0;
        for k=2:9
            ia=i+ey(k);
            ja=j+ex(k);
            if ja > nx
                ja=1;
            elseif ja < 1
                ja=nx;
            end
            % central difference approximation for divergence
            % div=div+wt(k)*(epsi(ia,ja)*(phi_e(ia,ja)-phi_e(i,j)));     %%%changes 1       
            div=div+wt(k)*(phi_e(ia,ja)-phi_e(i,j));

        end
        rho_e(i,j)=-((6*epsi(i,j)*div+0.5*(epsi_H-epsi_L)*(gradphi_ex(i,j)*g_phix(i,j)+gradphi_ey(i,j)*g_phiy(i,j)))); % factor depends on your stencil (3*2)
    end
end

% 
% %%%%%%%%%%%%%%%%%%Electric force calculation %%%%%%%%%%%%%%%%%
% 
% 
for i = 2:ny-1
    for j = 1:nx
        FE_x(i,j)=-0.5*0.5*(epsi_H-epsi_L)*g_phix(i,j)*E2(i,j)+rho_e(i,j)*E_x(i,j);
        FE_y(i,j)=-0.5*0.5*(epsi_H-epsi_L)*g_phiy(i,j)*E2(i,j)+rho_e(i,j)*E_y(i,j);
       FE1(i,j)=-0.5*0.5*(epsi_H-epsi_L)*g_phix(i,j)*E2(i,j);
       FE2(i,j)=-0.5*0.5*(epsi_H-epsi_L)*g_phiy(i,j)*E2(i,j);
       FE3(i,j)=rho_e(i,j)*E_x(i,j);     
       FE4(i,j)=rho_e(i,j)*E_y(i,j); 

    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%ELECTRIC SOLVER ENDS

% %%%%%%%%%%%%%%%%%%%%%%%%% calculation of force %%%%%%%%%%%%%
% 
for i = 2:ny-1
    for j = 1:nx
     Fxx(i,j)=Fs_x(i,j)+FE_x(i,j);
     Fyy(i,j)=Fs_y(i,j)+FE_y(i,j);

  
    end
end
% 
% 
% %%%%%%%%%%%%%%%%%%% calculation of Ri %%%%%%%%%%%%%%%%%%%%%%%%
% 
for i = 2:ny-1
    for j = 1:nx
        for k = 1:9
            % Dot product u ⋅ ∇ρ
            A1 = 0.5*(rho_H-rho_L)*(ux(i,j)*g_phix(i,j)+uy(i,j)*g_phiy(i,j));

            % Force term:  ci⋅F
            A2 = ex(k)*Fxx(i,j)+ey(k)*Fyy(i,j);

            % Tensor contraction term

            Q =0.5*(rho_H-rho_L)*((ex(k)*ex(k)-(1/3))*ux(i,j)*g_phix(i,j)+(ey(k)*ey(k)-(1/3))*uy(i,j)*g_phiy(i,j)+ex(k)*ey(k)*(ux(i,j)*g_phiy(i,j)+uy(i,j)*g_phix(i,j)));  % cross terms

            % Final R_i
            Ri(i,j,k) = wt(k)*(A1+3*A2+3*Q);
        end
    end
end
% %%%%%%%%%%%%%%%%%%%calculation of si %%%%%%%%%%%%%%%%
% 
% 
for i = 2:ny-1
 for j = 1:nx
     for k = 1:9
    DD= ex(k)*ux(i,j) + ey(k)*uy(i,j);   
    SS= ux(i,j)^2 + uy(i,j)^2;          
    si(i,j,k) = wt(k)*(3*DD +4.5*((DD)^2)-1.5*SS);
      end
  end
end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%distribution function for fluid flow
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
for i=2:ny-1
    for j=1:nx
        for k=1:9
           if k==1
                g_eq(i,j,k)=(wt(k)-1)*p(i,j)*3+rho(i,j)*si(i,j,k);
            else
                g_eq(i,j,k)=3*wt(k)*p(i,j)+rho(i,j)*si(i,j,k);   
           end
        end
     end
end

for i=2:ny-1
    for j=1:nx
        for k=1:9
          tau_g=0.5+3*(nu(i,j)/rho(i,j));
          gt(i,j,k)=-((g(i,j,k)-g_eq(i,j,k))/tau_g)+g(i,j,k)+(1-0.5/tau_g)*Ri(i,j,k);
        end
    end
end

for i=2:ny-1
    for j=1:nx
       for k=1:9
        ia=i+ey(k);
        ja=j+ex(k);

        if ja>nx
          ja=1;
        elseif ja<1
          ja=nx;
        end
                g(ia,ja,k)=gt(i,j,k);  

        end  
    end
end 
% %%% calculating the order parameter

% % Boundary condition
for j=1:nx
for i=2
        g(i,j,6)=gt(i,j,8);
         g(i,j,3)=gt(i,j,5);
          g(i,j,7)=gt(i,j,9);

end
for i=ny-1

        g(i,j,9)=gt(i,j,7);
        g(i,j,5)=gt(i,j,3);
        g(i,j,8)=gt(i,j,6);

end
end
% %%%%%%%%%%%% calculation of velocity  
% %%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%
for i=2:ny-1
    for j=1:nx
      ux(i,j)=0;
      uy(i,j)=0;
for k = 1:9
    ux(i,j) = ux(i,j) + ex(k)*g(i,j,k);
    uy(i,j) = uy(i,j)+ ey(k)*g(i,j,k);
end
ux(i,j)=(ux(i,j)+0.5*Fxx(i,j))/rho(i,j);
uy(i,j)=(uy(i,j)+0.5*Fyy(i,j))/rho(i,j);
   end 
end
% %%%%%%%%%%%%%%%%%%%%%%%%%% calculation of hydrodynamic pressure%%%%%%%%%%%%%%%%55
for i =2:ny-1
 for j =1:nx
     for k = 1
    DD= ex(k)*ux(i,j) + ey(k)*uy(i,j);   
    SS= ux(i,j)^2 + uy(i,j)^2;          
    so(i,j) = wt(k)*(3*DD +4.5*(DD^2)-1.5*SS);
      end
  end
end


for i=2:ny-1
    for j=1:nx
      p(i,j)=0;
        for k=2:9
    p(i,j)=p(i,j) + g(i,j,k);
        end
      p(i,j)= 0.6*((p(i,j)+0.5*0.5*(rho_H-rho_L)*(ux(i,j)*g_phix(i,j)+uy(i,j)*g_phiy(i,j))+rho(i,j)*so(i,j)));
   end 
end
for i=2:ny-1
    for j=1:nx
     phi_new(i,j)=phi(i,j);
     ux_new(i,j) =ux(i,j);
     uy_new(i,j) =uy(i,j);
     phi_now(i,j) = phi(i,j);
   end 
end

 %% Every 1000 steps, check convergence
    if mod(step, 1000)== 0
        nume = 0;
        denom= 0;
         for i= 2:ny-1
            for j = 1:nx
                nume = nume + abs(phi_now(i,j)-phi_ot(i,j));
                denom = denom+abs(phi_now(i,j));
            end
        end
        error = nume/(denom);  % avoid divide-by-zero
        for i=2:ny-1
        for j=1:nx
        phi_ot(i,j) = phi(i,j);
       end 
        end
    end
%%%%%%%%%%%%%%%%%%% Streamlines  %%%%%%%%%%%%
% i=1:ny;
% j=1:nx;
% [I,J]= meshgrid(i,j);
% % contour(I,J,psi,50)
% u11=ux(i,j);
% vy11=uy(i,j);
% streamslice(I,J,u11,vy11,3)
% xlabel('x')
% ylabel('y')
% grid on


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Uncomment it for observing the droplet
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%

% figure(302); clf; set(gcf,'Color','w');
% 
% imagesc([1 size(X,2)], [1 size(Y,1)], phi);
% axis image tight;
% set(gca,'YDir','normal');           % origin at bottom-left
% colormap('jet'); colorbar;
% 
% xlabel('x'); ylabel('y');
% title(sprintf('Phase field  \\phi  (step %d)', step), ...
%       'Interpreter','tex','FontSize',12);               
% 
% Optional: emphasize droplet interface (φ=0 contour)
% hold on;
% contour(X, Y, phi, [0 0], 'w', 'LineWidth', 1.2);
% 
% drawnow;


% 
% 
% 
% 
% 
% 
% 
% 
end


end

end


save('results.mat');

load('results.mat');





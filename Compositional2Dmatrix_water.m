clc
clear all

% length=22;
% height=6;
length=1;
height=1;

% ngx=218;
% ngy=58;
ngx=50;
ngy=50;

dx=length/(ngx+2);
dy=height/(ngy+2);

h=500*dx;

pInit=80*10^5;
pInj=120*10^5;
pPro=80*10^5;

% hetro_l=10;
% hetro_h=10;
% hetro_ij=[20 2];
 K=ones(ngy+2,ngx+2);
 K=K*2*10^(-12);
% K(hetro_ij(1):hetro_ij(1)+hetro_l,hetro_ij(2):hetro_ij(2)+hetro_h)=0.2*10^-12;
%KKK(1,:,:)=K;
% %KKK=smooth3(smooth3(KKK));
% KKK=5*exp(3*(smooth3(smooth3(randn(1,ngy+2,ngx+2)))))*10^-13;
% 
% for i=1:ngy+2
%     for j=1:ngx+2
%         K(i,j)=KKK(1,i,j);
%     end
% end
% K=load('Permeability.txt');
% K=K*10^-15;
n=100
%K=3*K*10^-3*10^-12
%K(1:200,1:60)=2*10^(-12);
%K=(log(K))*10^-12;
Sor=0;
Scg=0;
Swc=0;

vo=0.00001;
vg=0.000005;

p=ones(ngy+2,ngx+2);

p(1:ngy+2,1:ngx+2)=pInit;
% p(1:ngy+2,1)=pInj;
% p(2:ngy+1,ngx+2)=pPro;

% p(2:ngy+1,1)=pInj;
% p(2:ngy+1,ngx+2)=pPro;
p(2,2)=pInj;
p(ngy+2,ngx+2)=pPro;

sg(1:ngy+2,1:ngx+2)=0;
%sg(2:ngy+1,1)=1;
sg(2,2)=1;

kro=((1-sg-Sor)/(1-Scg-Sor)).^2;
krg=((sg-Scg)/(1-Scg-Sor)).^2;

Kx=(2*dx)./(dx./K(2:ngy+1,1:ngx+1)+dx./K(2:ngy+1,2:ngx+2));
Ky=(2*dy)./(dy./K(1:ngy+1,2:ngx+1)+dy./K(2:ngy+2,2:ngx+1));

Tox=Kx.*(kro(2:ngy+1,1:ngx+1))/(vo);
Tgx=Kx.*(krg(2:ngy+1,1:ngx+1))/(vg);
Toy=Ky.*(kro(1:ngy+1,2:ngx+1))/(vo);
Tgy=Ky.*(krg(1:ngy+1,2:ngx+1))/(vg);

% Tox=(2*dx)./(dx./(K(2:ngy+1,1:ngx+1).*kro(2:ngy+1,1:ngx+1)/vo)+dx./(K(2:ngy+1,2:ngx+2).*kro(2:ngy+1,2:ngx+2)/vo));
% Tgx=(2*dx)./(dx./(K(2:ngy+1,1:ngx+1).*krg(2:ngy+1,1:ngx+1)/vg)+dx./(K(2:ngy+1,2:ngx+2).*krg(2:ngy+1,2:ngx+2)/vg));
% Toy=(2*dy)./(dy./(K(1:ngy+1,2:ngx+1).*kro(1:ngy+1,2:ngx+1)/vo)+dy./(K(2:ngy+2,2:ngx+1).*kro(2:ngy+2,2:ngx+1)/vo));
% Tgy=(2*dy)./(dy./(K(1:ngy+1,2:ngx+1).*krg(1:ngy+1,2:ngx+1)/vg)+dy./(K(2:ngy+2,2:ngx+1).*krg(2:ngy+2,2:ngx+1)/vg));

% while 1
%      
%     for i=1:ngy
%         for j=1:ngx
%             
%             a(ngx*(i-1)+j)= (Tgy(i,j)  +Toy(i,j))/dy^2;
%             b(ngx*(i-1)+j)= (Tgx(i,j)  +Tox(i,j))/dx^2;
%             c(ngx*(i-1)+j)=(-Tgx(i,j+1)-Tgx(i,j))/dx^2+(-Tgy(i+1,j)-Tgy(i,j))/dy^2+(-Tox(i,j+1)-Tox(i,j))/dx^2+(-Toy(i+1,j)-Toy(i,j))/dy^2;
%             d(ngx*(i-1)+j)= (Tgx(i,j+1)+Tox(i,j+1))/dx^2;
%             e(ngx*(i-1)+j)= (Tgy(i+1,j)+Toy(i+1,j))/dy^2;
%             
%             F(ngx*(i-1)+j)=0;
%             
%         end
%     end
%     
%     %build ceof matrix
%     mn=ngx*ngy;
%     c(1:ngx)=c(1:ngx)+a(1:ngx);
%     c(mn-ngx+1:mn)=c(mn-ngx+1:mn)+e(mn-ngx+1:mn);
%       c(ngx+1:ngx:ngx*ngy-ngx+1)=c(ngx+1:ngx:ngx*ngy-ngx+1)+b(ngx+1:ngx:ngx*ngy-ngx+1);
%       c(ngx:ngx:ngx*ngy-ngx)=c(ngx:ngx:ngx*ngy-ngx)+d(ngx:ngx:ngx*ngy-ngx);
%     
%     aa=sparse(1:mn,1:mn,c(1:mn),mn,mn);
%     bb=sparse(1:mn-ngx,ngx+1:mn,e(1:mn-ngx),mn,mn);
%     cc=sparse(ngx+1:mn,1:mn-ngx,a(ngx+1:mn),mn,mn);
%     dd=sparse(1:mn-1,2:mn,d(1:mn-1),mn,mn);
%     ee=sparse(2:mn,1:mn-1,b(2:mn),mn,mn);
%     
%     A=aa+bb+cc+dd+ee;
%     
%     A(ngx:ngx:mn-ngx,ngx+1:ngx:mn-ngx+1)=0;
%     A(ngx+1:ngx:mn-ngx+1,ngx:ngx:mn-ngx)=0;
%     
%     
% %     F(ngx:ngx:mn)=F(ngx:ngx:mn)-d(ngx:ngx:mn).*(p(2:ngy+1,ngx+2)');
% %     F(1:ngx:mn-ngx+1)=F(1:ngx:mn-ngx+1)-b(1:ngx:mn-ngx+1).*(p(2:ngy+1,1)');
%     
%     
%       F(1)=F(1)-b(1)*pInj;
%       F(ngx*ngy)=F(ngx*ngy)-d(ngx*ngy)*pPro;
%     
%     
%     u=A\F';
%     for i=1:ngy
%         for j=1:ngx
%             p(i+1,j+1)=u(ngx*(i-1)+j);
%         end
%     end
% %     p(1,2:ngx+1)=p(2,2:ngx+1);
% %     p(ngy+2,2:ngx+1)=p(ngy+1,2:ngx+1);
% %     p(3:ngy+1,1)=p(3:ngy+1,2);
% %     p(2:ngy,ngx+2)=p(2:ngy,ngx+1);
% %     p(1,2:ngx+1)=p(2,2:ngx+1);
% %     p(ngy+2,2:ngx+1)=p(ngy+1,2:ngx+1);
% 
%     
%     for i=1:ngy
%         for j=1:ngx
%             transibility=((Tgx(i,j+1)*((p(i+1,j+1+1)-p(i+1,j+1))/dx)-Tgx(i,j)*((p(i+1,j+1)-p(i+1,j))/dx))/dx+(Tgy(i+1,j)*((p(i+1+1,j+1)-p(i+1,j+1))/dy)-Tgy(i,j)*((p(i+1,j+1)-p(i,j+1))/dy))/dy);
%             sg(i+1,j+1)=sg(i+1,j+1)+transibility*h;
%          end
%     end
%      
% %     sg(1,2:ngx+1)=sg(2,2:ngx+1);
% %     sg(ngy+2,2:ngx+1)=sg(ngy+1,2:ngx+1);
% %     sg(3:ngy+1,1)=sg(3:ngy+1,2);
% %     sg(2:ngy,ngx+2)=sg(2:ngy,ngx+1);
% %     sg(
%     
% 
%     kro=((1-sg-Sor)/(1-Scg-Sor)).^2;
%     krg=((sg-Scg)/(1-Scg-Sor)).^2;
% % % 
% for i=1:ngy+1
%     for j=1:ngx
%         if p(i,j+1)>=p(i+1,j+1)
%         Toy(i,j)=Ky(i,j)*kro(i,j+1)/vo;
%         Tgy(i,j)=Ky(i,j)*krg(i,j+1)/vg;
%         end
%         if p(i,j+1)<p(i+1,j+1)
%         Toy(i,j)=Ky(i,j)*kro(i+1,j+1)/vo;
%         Tgy(i,j)=Ky(i,j)*krg(i+1,j+1)/vg;
%         end
%     end
% end
% for i=1:ngy
%     for j=1:ngx+1
%         if p(i+1,j)>=p(i+1,j+1)
%         Tox(i,j)=Kx(i,j)*kro(i+1,j)/vo;
%         Tgx(i,j)=Kx(i,j)*krg(i+1,j)/vg;
%         end
%         if p(i+1,j)<p(i+1,j+1)
%         Tox(i,j)=Kx(i,j)*kro(i+1,j+1)/vo;
%         Tgx(i,j)=Kx(i,j)*krg(i+1,j+1)/vg;
%         end
%     end
% end
%         
% %     Tox=Kx.*kro(2:ngy+1,1:ngx+1)/(vo);
% %     Tgx=Kx.*krg(2:ngy+1,1:ngx+1)/(vg);
% %     Toy=Ky.*kro(1:ngy+1,2:ngx+1)/(vo);
% %     Tgy=Ky.*krg(1:ngy+1,2:ngx+1)/(vg);
% 
% % Tox=(2*dx)./(dx./(K(2:ngy+1,1:ngx+1).*kro(2:ngy+1,1:ngx+1)/vo)+dx./(K(2:ngy+1,2:ngx+2).*kro(2:ngy+1,2:ngx+2)/vo));
% % Tgx=(2*dx)./(dx./(K(2:ngy+1,1:ngx+1).*krg(2:ngy+1,1:ngx+1)/vg)+dx./(K(2:ngy+1,2:ngx+2).*krg(2:ngy+1,2:ngx+2)/vg));
% % Toy=(2*dy)./(dy./(K(1:ngy+1,2:ngx+1).*kro(1:ngy+1,2:ngx+1)/vo)+dy./(K(2:ngy+2,2:ngx+1).*kro(2:ngy+2,2:ngx+1)/vo));
% % Tgy=(2*dy)./(dy./(K(1:ngy+1,2:ngx+1).*krg(1:ngy+1,2:ngx+1)/vg)+dy./(K(2:ngy+2,2:ngx+1).*krg(2:ngy+2,2:ngx+1)/vg));
% 
%     imagesc(sg);
%     
%     pause(0.01)
% end


pcolor(K)
grid off
shading flat
title('Permeability')
xlabel('X direction 50X0.2(m)')
ylabel('Y direction 50X0.2(m)')

figure
pcolor(sg)

grid off
shading flat
title('Water saturation')
xlabel('X direction 220X0.1(m)')
ylabel('Y direction 60X0.1(m)')

figure
 pcolor(p)

grid off
shading flat
title('Vapor pressure')
xlabel('X direction 50X0.2(m)')
ylabel('Y direction 50X0.2(m)')

figure
 pcolor(z(:,:,2))

grid off
shading flat
title('CO2 concentration')
xlabel('X direction 50X0.2(m)')
ylabel('Y direction 50X0.2(m)')
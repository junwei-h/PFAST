% This is a matlab script to import and visualize 3D Field models for 
% Transmission Tomography
% last updated on Sep 22 2010 by Junwei Huang @ NRCan, March 28 2011

% Note: this matlab script must be used together with the program manual.
% Note: authors strongly recommend users to understand the script before
% making any changes. 
clear all

% or load from mat file
load T_Field3D

a=size(modvp0);

% dx=25;dy=25;dz=30;
dx=20;dy=20;dz=10;datum=600;
[x,y,z]=meshgrid(dx:dx:a(2)*dx,dy:dy:a(1)*dy,dz:dz:a(3)*dz);
%%
figure(1)

slice(x,y,z,modvp0,dx,dy,180);hold on;shading interp;colorbar;
box on;title('Initial Model');zlabel('Depth (m)');
xlabel('X (m)');ylabel('Y (m)');
plot3(rec(1:1:end,1),rec(1:1:end,2),rec(1:1:end,3),'.');
plot3(src(1:1:end,1),src(1:1:end,2),src(1:1:end,3),'r*');
set(gca,'zdir','reverse')
% daspect([1 1 0.3])
axis equal tight
view([-23 50]);
figure(2)
slice(x,y,z,modvp,dx,dy,180);shading interp;colorbar;
hold on;box on;title('Final Model');zlabel('Depth (m)');
xlabel('X (m)');ylabel('Y (m)');
set(gca,'zdir','reverse')
% daspect([1 1 0.3])
axis equal tight
view([-23 50]);

figure(3)
slice(x,y,z,modvp,75*dx,75*dy,180);shading interp;colorbar;
hold on;box on;title('Final Model');zlabel('Depth (m)');
xlabel('X (m)');ylabel('Y (m)');
set(gca,'zdir','reverse')
% daspect([1 1 0.3])
axis([0 inf 0 inf 0 500]);
view([0 0])
daspect([1 1 0.3])
%%
pathdir='..\models\';
fid=fopen([pathdir,'Field_Tr3D.bin'],'r');tmp=zeros(776146,1);tmp(:)=fread(fid,inf,'double');
fclose(fid);Tr=tmp;clear tmp

tri=284350:285000;
plot(tri,1000*Tr(tri,1),'k.-');grid on;
axis tight; legend('Transmission');
xlabel('Trace Number');ylabel('Travel Time (ms)');set(gca,'ydir','reverse');
title('Sample First Arrivals');
%%
%pick the depth of the refractor
[xx,yy]=meshgrid(dx:dx:a(2)*dx,dy:dy:a(1)*dy);
modvp(modvp==330)=NaN;
modvp0(modvp0==330)=NaN;
a=size(modvp);
zz=zeros(size(xx));C=zeros(size(xx));
for i=1:a(1)
    for j=1:a(2)
        f6=permute(modvp(i,j,:),[3,2,1]);
        [spmax, ~, ~]= extrema(f6);spmax([1 end],:)=[];           
        if ~isempty(spmax)
            idx=find(spmax(:,2)==max(spmax(:,2)));
            C(i,j)=spmax(idx,2);
            zz(i,j)=spmax(idx,1)*dz;
        else
            fprintf('No maxima found at (%d,%d)\n',i,j);
            zz(i,j)=NaN;
        end            

    end
end

figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,4,[1 2 3])
pcolor(xx/1e3,yy/1e3,datum-zz);shading interp
colorbar;box on;caxis([350 450]);hold on;
ix=ceil(2.0316*1e3/dx);iy=ceil(1.7009*1e3/dy);
plot(ix*dx/1e3,iy*dy/1e3,'rs')
axis equal tight
xlabel('Xline (km)');ylabel('Inline (km)');
title('Refractor Elevation (m)');

subplot(2,4,[5 6 7])
pcolor(xx/1e3,yy/1e3,C);shading interp
colorbar;box on;hold on;
plot(ix*dx/1e3,iy*dy/1e3,'rs')
axis equal tight
xlabel('Xline (km)');ylabel('Inline (km)');
title('Refractor Velocity (m/s)');
caxis([5000 6500])

subplot(2,4,[4 8])
f0=permute(modvp0(iy,ix,:),[3,2,1]);
f6=permute(modvp(iy,ix,:),[3,2,1]);
plot(f0,[10:10:200*10]','.-',f6,[dz:dz:200*dz]','r.-');
set(gca,'ydir','reverse')
grid on;xlabel('Vp (m/s)');ylabel('Depth (m)')


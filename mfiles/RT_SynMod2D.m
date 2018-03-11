% This is a matlab script to import and visualize 2D synthetic models for 
% Transmission/Reflection Tomography
% last updated on Sep 22 2010 by Junwei Huang @ NRCan, March 28 2011

% Note: this matlab script must be used together with the program manual.
% Note: authors strongly recommend users to understand the script before
% making any changes. 
clear all
load model2d
x=dx:dx:a(2)*dx;y=dy:dy:a(1)*dy;

figure(1)
subplot(311)
imagesc(x,y,modvp);colorbar;
hold on;box on;title('True Model');ylabel('Depth (m)');
plot(rec(1:end,2),rec(1:end,1),'g^');
%plot3(src(25:50:end,2),src(25:50:end,1),src(25:50:end,3),'r*');
plot(src(:,2),src(:,1),'r*');
for i=1:size(RR,1),plot(dx:dx:dx*a(2),RR(i,:),'.-');hold on;end
axis equal tight;hold off;
daspect([1 0.5 1])
subplot(312)
imagesc(x,y,modvp0);colorbar;
box on;title('Background Model');ylabel('Depth (m)');
axis equal tight;
daspect([1 0.5 1])
subplot(313)
imagesc(x,y,modvp-modvp0);colorbar;
hold on;plot(dx:dx:dx*a(2),z,'r:');hold off
box on;title('Anomaly');xlabel('Distance (m)');ylabel('Depth (m)');
axis equal tight;
daspect([1 0.5 1])

%%
pathdir='..\models';
fid=fopen([pathdir,'\RT_SynMod2D_FW.rec'],'wb');
fwrite(fid,permute(rec,[2 1]),'double');
fclose(fid);
fid=fopen([pathdir,'\RT_SynMod2D.src'],'wb');
fwrite(fid,permute(src,[2 1]),'double');
fclose(fid);
fid=fopen([pathdir,'\RT_SynMod2D.ref'],'wb');
fwrite(fid,permute(Z,[2 1]),'double');
fclose(fid);
fid=fopen([pathdir,'\RT_SynMod2D_FW.vp'],'wb');
fwrite(fid,permute(modvp,[2 1]),'double');
fclose(fid);
fid=fopen([pathdir,'\RT_SynMod2D_Inv.vp'],'wb');
fwrite(fid,permute(modvp0,[2 1]),'double');
fclose(fid);
fid=fopen([pathdir,'\RT_SynMod2D_Inv.wr'],'wb');
fwrite(fid,permute(WR,[2 1]),'double');
fclose(fid);

pathdir='..\par';
fid=fopen([pathdir,'\RT_SynMod2D_FW.inp'],'w');
fprintf(fid,'#The Input file for Travel time Tomography using Fast Sweeping Method \n');
fprintf(fid,'#software developed by Junwei Huang @ NRCan, Ottawa, starting from Mar 08, 2010\n');
fprintf(fid,'#====================\n');
fprintf(fid,'#The file name of the initial velocity\n');
fprintf(fid,'../models/RT_SynMod2D_FW.vp\n');
fprintf(fid,'#The file name of the source locations\n');
fprintf(fid,'../models/RT_SynMod2D.src\n');
fprintf(fid,'%d %d\n',size(src,1),size(src,2));
fprintf(fid,'#The file name of the receiver locations for Tran and Reflection\n');
fprintf(fid,'../models/RT_SynMod2D_FW.rec\n');
fprintf(fid,'%d %d\n',size(rec,1),size(rec,2));
fprintf(fid,'#The file name of the reflector locations\n');
fprintf(fid,'../models/RT_SynMod2D.ref\n');
fprintf(fid,'%d %d\n',size(Z,1),size(Z,2));
fprintf(fid,'#Size of the model, rows, columns, sample interval\n');
fprintf(fid,'%f %f %f %f\n',a(1),a(2),dx,dy);
fprintf(fid,'#===Output the inversed model===\n');
fprintf(fid,'../models/RT_SynMod2D.finalvp\n');
fprintf(fid,'#===Regularization Parameter: nux, nuz\n');
fprintf(fid,'60 10\n');
fprintf(fid,'#===Weighting factor of Reflection (WR)\n');
fprintf(fid,'../models/RT_SynMod2D_Inv.wr\n');
fprintf(fid,'#===Inversion Scheme tag (schemetag) and Line search tag (lsrc)\n');
fprintf(fid,'#recommend: 3 2, 3 3, 1 2, 1 3, 3 1, 2 2, 2 3, 2 1, 0 0\n');
fprintf(fid,'#schemetag:\n');
fprintf(fid,'#	1-HS Nonlinear Conjugate Method\n');
fprintf(fid,'#	2-CGDESCENT Nonlinear Conjugate Method\n');
fprintf(fid,'#	3-L-BFGS quasi-newton method\n');
fprintf(fid,'#	0-Steepest Decent method\n');
fprintf(fid,'#lsr:\n');
fprintf(fid,'#	1-Secant method with exact Conditions\n');
fprintf(fid,'#	2-Secant method Line Search with Strong Wolfe Conditions\n');
fprintf(fid,'#	3-Cubic Interpolation method with Strong Wolfe Conditions\n');
fprintf(fid,'#	0-No Line Search for Steepest Decent only\n');
fprintf(fid,'3 2\n');
fclose(fid);

%%
% run after travel time is calculated
% to create inversion input files
pathdir='..\models\';
fid=fopen([pathdir,'RT_SynMod2D_FW_Tr2D.bin'],'r');tmp=zeros(size(Z,1)/a(2)+1,size(rec,1));tmp(:)=fread(fid,inf,'double');
tmp=permute(tmp,[2 1]);fclose(fid);Tr=tmp;clear tmp

tri=22801:24900;
plot(tri,1000*Tr(tri,1),'k.-',tri,1000*Tr(tri,2),'r.-',tri,1000*Tr(tri,3),'g.-');grid on;
axis tight; legend('Transmission','1^s^t Reflection','2^n^d Reflection');
xlabel('Trace Number');ylabel('Travel Time (ms)');set(gca,'ydir','reverse');hold off

nTr=zeros(size(Tr));
for j=1:length(Tr(1,:))
    for i=1:length(Tr(:,j))
        nTr(i,j)=Tr(i,j)+(rand(1)-0.5)*sqrt(12)*Tr(i,1)*2/100;
    end
end

plot(tri,1000*nTr(tri,1),'k.-',tri,1000*nTr(tri,2),'r.-',tri,1000*nTr(tri,3),'g.-');grid on;
axis([-inf inf -50 2050]); legend('Transmission','1^s^t Reflection','2^n^d Reflection');
title('Travel Time with 2% Noise')
xlabel('Trace Number');ylabel('Noisy Travel Time (ms)');set(gca,'ydir','reverse');hold off

%%
rec2=[rec(:,1:4) nTr];
rec2(:,6)=0;
pathdir='..\models';
fid=fopen([pathdir,'\RT_SynMod2D_Inv_rand2p.rec'],'wb');
fwrite(fid,permute(rec2,[2 1]),'double');
fclose(fid);

rec2=[rec(:,1:4) Tr];
rec2(:,6)=0;
pathdir='..\models';
fid=fopen([pathdir,'\RT_SynMod2D_Inv.rec'],'wb');
fwrite(fid,permute(rec2,[2 1]),'double');
fclose(fid);
        
%%
pathdir='..\par';
fid=fopen([pathdir,'\RT_SynMod2D_Inv.inp'],'w');
fprintf(fid,'#The Input file for Travel time Tomography using Fast Sweeping Method \n');
fprintf(fid,'#software developed by Junwei Huang @ NRCan, Ottawa, starting from Mar 08, 2010\n');
fprintf(fid,'#====================\n');
fprintf(fid,'#The file name of the initial velocity\n');
fprintf(fid,'../models/RT_SynMod2D_Inv.vp\n');
fprintf(fid,'#The file name of the source locations\n');
fprintf(fid,'../models/RT_SynMod2D.src\n');
fprintf(fid,'%d %d\n',size(src,1),size(src,2));
fprintf(fid,'#The file name of the receiver locations for Tran and Reflection\n');
fprintf(fid,'../models/RT_SynMod2D_Inv_rand2p.rec\n');
fprintf(fid,'%d %d\n',size(rec2,1),size(rec2,2));
fprintf(fid,'#The file name of the reflector locations\n');
fprintf(fid,'../models/RT_SynMod2D.ref\n');
fprintf(fid,'%d %d\n',size(Z,1),size(Z,2));
fprintf(fid,'#Size of the model, rows, columns, sample interval\n');
fprintf(fid,'%f %f %f %f\n',a(1),a(2),dx,dy);
fprintf(fid,'#===Output the inversed model===\n');
fprintf(fid,'../models/RT_SynMod2D.finalvp\n');
fprintf(fid,'#===Regularization Parameter: nux, nuz\n');
fprintf(fid,'60 20\n');
fprintf(fid,'#===Weighting factor of Reflection (WR)\n');
fprintf(fid,'../models/RT_SynMod2D_Inv.wr\n');
fprintf(fid,'#===Inversion Scheme tag (schemetag) and Line search tag (lsrc)\n');
fprintf(fid,'#recommend: 3 2, 3 3, 1 2, 1 3, 3 1, 2 2, 2 3, 2 1, 0 0\n');
fprintf(fid,'#schemetag:\n');
fprintf(fid,'#	1-HS Nonlinear Conjugate Method\n');
fprintf(fid,'#	2-CGDESCENT Nonlinear Conjugate Method\n');
fprintf(fid,'#	3-L-BFGS quasi-newton method\n');
fprintf(fid,'#	0-Steepest Decent method\n');
fprintf(fid,'#lsr:\n');
fprintf(fid,'#	1-Secant method with exact Conditions\n');
fprintf(fid,'#	2-Secant method Line Search with Strong Wolfe Conditions\n');
fprintf(fid,'#	3-Cubic Interpolation method with Strong Wolfe Conditions\n');
fprintf(fid,'#	0-No Line Search for Steepest Decent only\n');
fprintf(fid,'3 2\n');
fclose(fid);

%%
%to read downT and upT
pathdir='..\models';
airidx=find(modvp<500);
for srcid=1:length(js)
    fpname=[pathdir,'\RT_SynMod2D_FW_Src',num2str(srcid-1),'_TTd.2d'];
    fid=fopen(fpname,'r');tmp1=zeros(a(2),a(1),size(RR,1));
    tmp1(:)=fread(fid,inf,'double');
    tmp1=permute(tmp1,[2 1 3]);
    fclose(fid);
    tmp1(tmp1>1000)=NaN;
%     eval(['dT',num2str(ir-1),'=tmp;']);
%     eval(['dT',num2str(ir-1),'(find(dT',num2str(ir-1),'>1000))=NaN;']);
    
    fpname=[pathdir,'\RT_SynMod2D_FW_Src',num2str(srcid-1),'_TTu.2d'];
    fid=fopen(fpname,'r');tmp2=zeros(a(2),a(1),size(RR,1));
    tmp2(:)=fread(fid,inf,'double');
    tmp2=permute(tmp2,[2 1 3]);
    fclose(fid);
    tmp2(tmp2>1000)=NaN;
%     eval(['uT',num2str(ir-1),'=tmp;']);
%     eval(['uT',num2str(ir-1),'(find(uT',num2str(ir-1),'>1000))=NaN;']);
%     clear tmp 
    screen_size = get(0, 'ScreenSize');
    f1 = figure(1);
    set(f1, 'Position', [0 0 screen_size(3) screen_size(4) ] );
    set(gcf,'PaperPositionMode','auto')   
    for ir=1:size(RR,1)
        subplot(size(RR,1),1,ir)
        imagesc(x,y,modvp);colorbar;hold on;
        for i=1:size(RR,1),plot(dx:dx:dx*a(2),RR(i,:),'k');hold on;end
        dT(:,:)=tmp1(:,:,ir);dT(airidx)=NaN;
        uT(:,:)=tmp2(:,:,ir);uT(airidx)=NaN;
        contour(x,y,dT,30,'w','linewidth',2);
        contour(x,y,uT,30,'r','linewidth',2);
%         eval(['contour(x2,y2,dT',num2str(ir-1),'(:,:),30,''w'')']);
%         eval(['contour(x2,y2,uT',num2str(ir-1),'(:,:),30,''r'')']);
        axis equal tight
%         set(gca,'dataaspectratio',[1 0.2 1])
        caxis([2000 4000])
        plot(rec(1:10:end,2),rec(1:10:end,1),'^','MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',12);
        plot(src(srcid,2),src(srcid,1),'*','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',12);
        title(['Shot ',num2str(srcid)],'fontsize',18,'fontweight','b');
        xlabel('Horizon X (m)','fontsize',18,'fontweight','b');
        ylabel('Depth Z (m)','fontsize',18,'fontweight','b');
        set(gca,'fontsize',18,'fontweight','b');
        

    end
    pause(0.1);
%         filename=['..\models\RT2layer\snapshot\shot',num2str(srcid),'.tif'];
%      print('-dtiff','-r100',filename);
end

%%
pathdir='..\models';
% pathdir='.\';
for iter=1
    fid=fopen([pathdir,'/RT_SynMod2D_Inv_bf_regu_RT',num2str(iter)],'r');tmp=zeros(a(2),a(1));tmp(:)=fread(fid,inf,'double');
    tmp=permute(tmp,[2 1 ]);fclose(fid);RTbf=tmp;
    fid=fopen([pathdir,'/RT_SynMod2D_Inv_af_regu_RT',num2str(iter)],'r');tmp=zeros(a(2),a(1));tmp(:)=fread(fid,inf,'double');
    tmp=permute(tmp,[2 1 ]);fclose(fid);RTaf=tmp;

    fid=fopen([pathdir,'/RT_SynMod2D_Inv_bf_regu_R',num2str(iter)],'r');tmp=zeros(a(2),a(1));tmp(:)=fread(fid,inf,'double');
    tmp=permute(tmp,[2 1 ]);fclose(fid);Rbf=tmp;
    fid=fopen([pathdir,'/RT_SynMod2D_Inv_af_regu_R',num2str(iter)],'r');tmp=zeros(a(2),a(1));tmp(:)=fread(fid,inf,'double');
    tmp=permute(tmp,[2 1 ]);fclose(fid);Raf=tmp;

    fid=fopen([pathdir,'/RT_SynMod2D_Inv_bf_regu_T',num2str(iter)],'r');tmp=zeros(a(2),a(1));tmp(:)=fread(fid,inf,'double');
    tmp=permute(tmp,[2 1 ]);fclose(fid);Tbf=tmp;
    fid=fopen([pathdir,'/RT_SynMod2D_Inv_af_regu_T',num2str(iter)],'r');tmp=zeros(a(2),a(1));tmp(:)=fread(fid,inf,'double');
    tmp=permute(tmp,[2 1 ]);fclose(fid);Taf=tmp;

    subplot(311)
    imagesc([1:a(2)]*dx,[1:a(1)]*dy,Taf);colorbar;axis tight
    % hold on;for i=1:size(RR,2),plot(dx:dx:dx*a(2),RR(:,i),'w.-');end
    title('Gradient of Transmission','fontsize',12,'fontweight','b');
    ylabel('Depth Z (m)','fontsize',12,'fontweight','b');
    set(gca,'fontsize',12,'fontweight','b');
    subplot(312)
    imagesc([1:a(2)]*dx,[1:a(1)]*dy,Raf);colorbar;axis tight
    % hold on;for i=1:size(RR,2),plot(dx:dx:dx*a(2),RR(:,i),'w.-');end
    title('Gradient of Reflection','fontsize',12,'fontweight','b');
    ylabel('Depth Z (m)','fontsize',12,'fontweight','b');
    set(gca,'fontsize',12,'fontweight','b');
    subplot(313)
    imagesc([1:a(2)]*dx,[1:a(1)]*dy,RTaf);colorbar;axis tight
    % hold on;for i=1:size(RR,2),plot(dx:dx:dx*a(2),RR(:,i),'w.-');end
    title('Gradient in use','fontsize',12,'fontweight','b');
    xlabel('Distance X (m)','fontsize',12,'fontweight','b');
    ylabel('Depth Z (m)','fontsize',12,'fontweight','b');
    set(gca,'fontsize',12,'fontweight','b');
end

%%
pathdir='..\models';
fid=fopen([pathdir,'\RT_SynMod2D.finalvp'],'r');tmp=zeros(a(2),a(1));tmp(:)=fread(fid,inf,'double');
tmp=permute(tmp,[2 1]);fclose(fid);

subplot(211)
imagesc([1:a(2)]*dx,[1:a(1)]*dy,tmp-modvp0);
hold on;plot(rec(1:5:end,2),rec(1:5:end,1),'k^');plot(src(:,2),src(:,1),'r*');
for i=1:size(RR,1),plot(dx:dx:dx*a(2),RR(i,:),'.-');hold on;end;hold off
colorbar;caxis([ -400 400])
title('Reflection Tomography');ylabel('Depth Z (m)');
subplot(212)
imagesc([1:a(2)]*dx,[1:a(1)]*dy,modvp-modvp0);
hold on;plot(rec(1:5:end,2),rec(1:5:end,1),'k^');plot(src(:,2),src(:,1),'r*');
for i=1:size(RR,1),plot(dx:dx:dx*a(2),RR(i,:),'.-');hold on;end;hold off
colorbar;caxis([ -400 400])
title('True Perturbation');ylabel('Depth Z (m)');xlabel('Distance X (km)');





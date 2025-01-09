function plot_comparison(BINNED,profID,filename,folder_main,folder_out)

FieldList = fieldnames(BINNED);
for iField = 3:numel(FieldList)
    Field    = FieldList{iField};
    BINNEDplot.(Field) = BINNED.(Field)(:,profID);
end


idx_good1=find(BINNEDplot.flag_T1==1 & BINNEDplot.flag_S1==1);
idx_good2=find(BINNEDplot.flag_T2==1 & BINNEDplot.flag_S1==1);

figure
subplot(2,2,1)
loglog(BINNEDplot.eps_S1(:),BINNEDplot.eps_T1(:),'ob','MarkerSize',2); hold on
loglog(BINNEDplot.eps_S1(:),BINNEDplot.eps_T2(:),'or','MarkerSize',2);
loglog(BINNEDplot.eps_S1(idx_good1),BINNEDplot.eps_T1(idx_good1),'ob','MarkerFaceColor','b','MarkerSize',2); hold on
loglog(BINNEDplot.eps_S1(idx_good2),BINNEDplot.eps_T2(idx_good2),'or','MarkerFaceColor','r','MarkerSize',2);
loglog([1E-12 1E-5],[1E-12 1E-5],'-k')
loglog([1E-12 1E-5],2.8*[1E-12 1E-5],'--k')
loglog([1E-12 1E-5],2.8^(-1)*[1E-12 1E-5],'--k')
loglog([1E-12 1E-5],[1E-11 1E-4],':k')
loglog([1E-12 1E-5],[1E-13 1E-6],':k')
% count n element that differ more than a factor 2.8 (propagation error) and a factor 10
nel_perr1 = sum(BINNEDplot.eps_T1(idx_good1)./BINNEDplot.eps_S1(idx_good1) < 2.8 & BINNEDplot.eps_T1(idx_good1)./BINNEDplot.eps_S1(idx_good1) >1/2.8);
nel_perr2 = sum(BINNEDplot.eps_T2(idx_good2)./BINNEDplot.eps_S1(idx_good2) < 2.8 & BINNEDplot.eps_T2(idx_good2)./BINNEDplot.eps_S1(idx_good2) >1/2.8);
nel_1ord1 = sum(BINNEDplot.eps_T1(idx_good1)./BINNEDplot.eps_S1(idx_good1) < 10 & BINNEDplot.eps_T1(idx_good1)./BINNEDplot.eps_S1(idx_good1) >0.1);
nel_1ord2 = sum(BINNEDplot.eps_T2(idx_good2)./BINNEDplot.eps_S1(idx_good2) < 10 & BINNEDplot.eps_T2(idx_good2)./BINNEDplot.eps_S1(idx_good2) >0.1);
nel_tot = sum(~isnan(BINNEDplot.eps_S1(:)));
nel_good1 = length(idx_good1);
nel_good2 = length(idx_good2);
% some statistics
text(1.5E-12,5E-6,['stats good fit'],'fontsize',5);
text(1.5E-12,2E-6,['2.8^{\pm 1}=',num2str(nel_perr1/nel_good1*100,'%3.0f\n'),'%, ',num2str(nel_perr2/nel_good2*100,'%3.0f\n'),'%'],'fontsize',5);
text(1.5E-12,8E-7,['10^{\pm 1}=',num2str(nel_1ord1/nel_good1*100,'%3.0f\n'),'%, ',num2str(nel_1ord2/nel_good2*100,'%3.0f\n'),'%'],'fontsize',5);
[r1, p1] = corrcoef(log10(BINNEDplot.eps_S1(idx_good1)),log10(BINNEDplot.eps_T1(idx_good1)),'rows','complete');
[r2, p2] = corrcoef(log10(BINNEDplot.eps_S1(idx_good2)),log10(BINNEDplot.eps_T2(idx_good2)),'rows','complete');
if sum(size(r1))~=2 && sum(size(r2))~=2
    text(1.5E-12,3E-7,['R^2=',num2str(r1(2,1)^2,1),', ',num2str(r2(2,1)^2,1)],'fontsize',5);
end

set(gca,'FontSize',7);
l=legend('Tf1','Tf2','location','southeast','fontsize',7);
set(l,'position',[0.52,0.63,0,0])
grid on
axis square; xlim([1E-12 1E-5]);ylim([1E-12 1E-5])
xlabel('\epsilon_{sh} (m^2 s^{-3})','fontsize',7)
ylabel('\epsilon_{T} (m^2 s^{-3})','fontsize',7)

subplot(2,2,2)
loglog(BINNEDplot.KOsborn_S1(:),BINNEDplot.KOsbornCox_T1(:),'ob','MarkerSize',2); hold on
loglog(BINNEDplot.KOsborn_S1(:),BINNEDplot.KOsbornCox_T2(:),'or','MarkerSize',2);
loglog(BINNEDplot.KOsborn_S1(idx_good1),BINNEDplot.KOsbornCox_T1(idx_good1),'ob','MarkerFaceColor','b','MarkerSize',2); hold on
loglog(BINNEDplot.KOsborn_S1(idx_good2),BINNEDplot.KOsbornCox_T2(idx_good2),'or','MarkerFaceColor','r','MarkerSize',2);
loglog([1E-8 1E-2],[1E-8 1E-2],'-k')
loglog([1E-8 1E-2],2.8*[1E-8 1E-2],'--k')
loglog([1E-8 1E-2],2.8^(-1)*[1E-8 1E-2],'--k')
loglog([1E-8 1E-2],[1E-7 1E-1],':k')
loglog([1E-8 1E-2],[1E-9 1E-3],':k')
% count n element that differ more than a factor 2.8 (propagation error) and a factor 10
nel_perr1 = sum(BINNEDplot.KOsbornCox_T1(idx_good1)./BINNEDplot.KOsborn_S1(idx_good1) < 2.8 & BINNEDplot.KOsbornCox_T1(idx_good1)./BINNEDplot.KOsborn_S1(idx_good1) >1/2.8);
nel_perr2 = sum(BINNEDplot.KOsbornCox_T2(idx_good2)./BINNEDplot.KOsborn_S1(idx_good2) < 2.8 & BINNEDplot.KOsbornCox_T2(idx_good2)./BINNEDplot.KOsborn_S1(idx_good2) >1/2.8);
nel_1ord1 = sum(BINNEDplot.KOsbornCox_T1(idx_good1)./BINNEDplot.KOsborn_S1(idx_good1) < 10 & BINNEDplot.KOsbornCox_T1(idx_good1)./BINNEDplot.KOsborn_S1(idx_good1) >0.1);
nel_1ord2 = sum(BINNEDplot.KOsbornCox_T2(idx_good2)./BINNEDplot.KOsborn_S1(idx_good2) < 10 & BINNEDplot.KOsbornCox_T2(idx_good2)./BINNEDplot.KOsborn_S1(idx_good2) >0.1);
% some statistics
text(1.5E-8,5E-3,['stats good fit'],'fontsize',5);
text(1.5E-8,2E-3,['2.8^{\pm 1}=',num2str(nel_perr1/nel_good1*100,'%3.0f\n'),'%, ',num2str(nel_perr2/nel_good2*100,'%3.0f\n'),'%'],'fontsize',5);
text(1.5E-8,8E-4,['10^{\pm 1}=',num2str(nel_1ord1/nel_good1*100,'%3.0f\n'),'%, ',num2str(nel_1ord2/nel_good2*100,'%3.0f\n'),'%'],'fontsize',5);
[r1, p1] = corrcoef(log10(BINNEDplot.KOsborn_S1(idx_good1)),log10(BINNEDplot.KOsbornCox_T1(idx_good1)),'rows','complete');
[r2, p2] = corrcoef(log10(BINNEDplot.KOsborn_S1(idx_good2)),log10(BINNEDplot.KOsbornCox_T2(idx_good2)),'rows','complete');
if sum(size(r1))~=2 && sum(size(r2))~=2
    text(1.5E-8,3E-4,['R^2=',num2str(r1(2,1)^2,1),', ',num2str(r2(2,1)^2,1)],'fontsize',5);
end

set(gca,'FontSize',7);
% legend('Tf1','Tf2','location','southeast','fontsize',7);
grid on
axis square; xlim([1E-8 1E-2]);ylim([1E-8 1E-2])
xlabel('K_{sh} (m^2 s^{-1})','fontsize',7)
ylabel('K_{T} (m^2 s^{-1})','fontsize',7)

ax3=subplot(2,2,3);
clear points
points(:,1)=[BINNEDplot.eps_S1(idx_good1); BINNEDplot.eps_S1(idx_good2) ];
points(:,2)=[BINNEDplot.eps_T1(idx_good1); BINNEDplot.eps_T2(idx_good2) ];
idx=find(~isnan(points(:,1)));points=points(idx,:);
minvals = -12; maxvals = -5;
pts = logspace(minvals, maxvals, 41);
N = histcounts2(points(:,2), points(:,1), pts, pts); N=N/sum(N(:)); N(N==0)=NaN;
b=imagesc(pts, pts, N); hold on
% b=imagesc(pts, pts, N,'AlphaData',N/max(N(:)),'AlphaDataMapping','scaled'); hold on
colormap(ax3,parula(13)); set(b,'AlphaData',~isnan(N)); caxis([-0.0005, 0.0125]); cb=colorbar;
cb.Ticks = linspace(0, 0.012,7);
cb.TickLabels = num2cell(cb.Ticks); cb.TickLabels{end}='>0.012';
pos = get(cb,'Position');
set(cb,'Position',pos + [0.05,0,0,0])
ylabel(cb, 'Probability')
set(gca, 'XLim', pts([1 end]), 'YLim', pts([1 end]), 'YDir', 'normal','xscale','log','yscale','log','FontSize',7);
loglog([1E-12 1E-5],[1E-12 1E-5],'-k')
loglog([1E-12 1E-5],2.8*[1E-12 1E-5],'--k')
loglog([1E-12 1E-5],2.8^(-1)*[1E-12 1E-5],'--k')
loglog([1E-12 1E-5],[1E-11 1E-4],':k')
loglog([1E-12 1E-5],[1E-13 1E-6],':k')
grid on
axis square;
text(1.5E-12,5E-6,['N_{good}=',num2str(nel_good1),', ',num2str(nel_good2)],'fontsize',7);
text(1.5E-12,2E-6,['%_{good}=',num2str(nel_good1/nel_tot*100,'%3.0f\n'),'%, ',num2str(nel_good2/nel_tot*100,'%3.0f\n'),'%'],'fontsize',7);

xlabel('\epsilon_{sh} (m^2 s^{-3})','fontsize',7)
ylabel('\epsilon_{T} (m^2 s^{-3})','fontsize',7)

ax4=subplot(2,2,4);
% loglog(BINNEDplot.eps_S1(:),BINNEDplot.eps_T1(:)./BINNEDplot.eps_S1(:),'ob','MarkerSize',2); hold on
% loglog(BINNEDplot.eps_S1(:),BINNEDplot.eps_T2(:)./BINNEDplot.eps_S1(:),'or','MarkerSize',2);
% loglog(BINNEDplot.eps_S1(idx_good1),BINNEDplot.eps_T1(idx_good1)./BINNEDplot.eps_S1(idx_good1),'ob','MarkerFaceColor','b','MarkerSize',2); hold on
% loglog(BINNEDplot.eps_S1(idx_good2),BINNEDplot.eps_T2(idx_good2)./BINNEDplot.eps_S1(idx_good2),'or','MarkerFaceColor','r','MarkerSize',2);
clear points
points(:,1)=[BINNEDplot.eps_S1(idx_good1); BINNEDplot.eps_S1(idx_good2) ];
points(:,2)=[BINNEDplot.eps_T1(idx_good1)./BINNEDplot.eps_S1(idx_good1); BINNEDplot.eps_T2(idx_good2)./BINNEDplot.eps_S1(idx_good2) ];
idx=find(~isnan(points(:,1))); points=points(idx,:);
minvals = -12; maxvals = 4;
pts = logspace(minvals, maxvals, 101);
N = histcounts2(points(:,2), points(:,1), pts, pts); N=N/sum(N(:)); N(N==0)=NaN;
b=imagesc(pts, pts, N); hold on
colormap(ax3,parula(13)); set(b,'AlphaData',~isnan(N)); caxis([-0.0005, 0.0125]); 

loglog([1E-12 1E-5],[1 1],'-k')
loglog([1E-12 1E-5],[2.8 2.8],'--k')
loglog([1E-12 1E-5],[2.8^(-1) 2.8^(-1)],'--k')
loglog([1E-12 1E-5],[10 10],':k')
loglog([1E-12 1E-5],[10^(-1) 10^(-1)],':k')
set(gca, 'xlim',[1e-12 1E-5], 'YLim', [1e-4 1E4], 'YDir', 'normal','xscale','log','yscale','log','FontSize',7);
% set(gca,'xscale','log','yscale','log','xlim',[1e-12 1E-5],'ylim',[1e-4 1E4])
set(gca,'FontSize',7);
% legend('Tf1','Tf2','location','northeast','fontsize',7)
% cb=colorbar; colormap(ax4,parula(Nprf));
% caxis([-6,-2])
% pos = get(cb,'Position');
% set(cb,'Position',pos + [0.05,0,0,0])
% ylabel(cb, 'prof ID')
grid on
axis square;
xlabel('\epsilon_{sh} (m^2 s^{-3})','fontsize',7)
ylabel('\epsilon_{T}/\epsilon_{sh}','fontsize',7)

saveas(gcf,[folder_main,'/comparison_',filename,'.png'])
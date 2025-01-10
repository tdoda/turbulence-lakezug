function [] = plot_profile(DATA,SLOW,FAST,param,pmaxplot)
%PLOT_PROFILE Make plots of the profiles.
%   Detailed explanation goes here

info=param.info;
pplot=[info.pmin,info.pmax];


% 1. Pressure time series
figure('Units','centimeters','Position',[1 1 18 15])

clf
subplot(4,4,[1 4])
plot(DATA.t_slow/60,DATA.P_slow) % Complete pressure time series in the file
ylabel('p (db)')
hold on
% plot(DATA.t_slow(iips)/60,SLOW.pressure);
plot(SLOW.dtsec/60,SLOW.pressure); % Pressure time series of the selected profile
ylim(pplot)
xlabel('time [min]')

% 2. Shear profiles
subplot(4,4,[5,9,13])
if param.config.S1
    plot(FAST.fast_S1,FAST.pressure); hold on;
end
%if isfield(DATA,'sh2')
if param.config.S2
    plot(FAST.fast_S2+10,FAST.pressure); hold on;
end
plot(get(gca,'xlim'),[pmaxplot pmaxplot],'--k')
set(gca,'ydir','reverse');
ylabel('p(db)'); xlabel('sh (1/s)'); ylim(pplot);
%     subplot(4,4,[6,10,14])
%     plot(AA(:,1),Pf); hold on
%     plot(AA(:,2),Pf); set(gca,'ydir','reverse');
%     xlabel('Acc (1/s)')

% 3. FP07 data
subplot(4,4,[7,11,15]);
% if param.config.S2
%     plot(gradT1f,Pf)
%     hold on
%     plot(gradT1f+10,Pf);
% else
%     plot(T1f,Pf)
%     hold on
%     plot(T2f,Pf);
% end
% xlabel('T fast (°C)');

% T. Doda: modified to plot gradT only
if param.config.T1
    plot(FAST.grad_T1,FAST.pressure); hold on;
end
if param.config.T2
    plot(FAST.grad_T2+10,FAST.pressure); hold on;
end
xlabel('gradT fast (°C/m)');

plot(get(gca,'xlim'),[pmaxplot pmaxplot],'--k')
set(gca,'ydir','reverse');
ylim(pplot)
set(gca,'yticklabel',[])

% 4. Velocity and inclination
ax1=subplot(4,4,[8,12,16]);
plot(FAST.velocity,FAST.pressure,'-k'); hold on;
plot(get(gca,'xlim'),[pmaxplot pmaxplot],'--k')
set(gca,'ydir','reverse');
xlabel('W (m/s)');yticklabels([]);ylim(pplot);
if ~isempty(find(~isnan(SLOW.Incl_x),1)) % At least one non NaN value
    ax2=axes('Position',ax1.Position,'XAxisLocation','top',...
        'YAxisLocation','right','color','none',...
        'xColor','r','yColor','k');
    set(ax1,'box','off')
    line(SLOW.Incl_x,SLOW.pressure,'color','r'); set(gca,'ydir','reverse');
    xlabel('Inclination [°]')
    yticklabels([]); xlim([-5 5]);
end


% 5. add despiked and HP signals to the plot
subplot(4,4,[5,9,13])
if param.config.S1
    plot(FAST.fast_S1_hp,FAST.pressure)
end

if param.config.S2
    plot(FAST.fast_S2_hp+10,FAST.pressure)
end

subplot(4,4,[6,10,14])
plot(FAST.AA_filt,FAST.pressure); hold on
set(gca,'ydir','reverse');
plot([-200 -200],[min(FAST.pressure) max(FAST.pressure)],'--k'); plot([200 200],[min(FAST.pressure) max(FAST.pressure)],'--k');
xlim([-500 500]);  xlabel('Acc (1/s)'); ylim(pplot);
plot(get(gca,'xlim'),[pmaxplot pmaxplot],'--k')




end
function [hp]=plot_prof_ql(DATA,filenames_all,filenames,col,select_profile,varname,indprobe,smoothing_window,absolute)
% Plot all the profiles of a binned variable

hp=[];

if nargin<9
    absolute=false;
end

indprof=0;
for kf=1:length(filenames)
    ind_file=find(strcmp(filenames_all,filenames{kf}),1);
    prof_data=DATA(ind_file).DISS_QL;
    if strcmp(select_profile,'all')
        for kprof=1:length(prof_data)
            indprof=indprof+1;
            if absolute
                %scatter(abs(prof_data{kprof}.(varname)),prof_data{kprof}.depth,'o','SizeData',3,'MarkerFaceColor',col(indprof,:),'MarkerEdgeColor','none','MarkerFaceAlpha',0.5);
                %hold on
                hp(end+1)=plot(movmean(abs(prof_data{kprof}.(varname)(indprobe,:)),smoothing_window),prof_data{kprof}.P,'Color',col(indprof,:),'LineWidth',1);
                hold on
            else
                hp(end+1)=plot(movmean(prof_data{kprof}.(varname)(indprobe,:),smoothing_window),prof_data{kprof}.P,'Color',col(indprof,:),'LineWidth',1);
                hold on
            end
        end
    else  % Plot the longest profile only
        depth_range=NaN(1,length(prof_data)); % [m]
        for kprof=1:length(prof_data)
            depth_range(kprof)=max(prof_data{kprof}.P)-min(prof_data{kprof}.P);
        end
        [~,indmax]=max(depth_range);
        if absolute
            %scatter(abs(prof_data{indmax}.(varname)),prof_data{indmax}.depth,'o','SizeData',3,'MarkerFaceColor',col(kf,:),'MarkerEdgeColor','none','MarkerFaceAlpha',0.1);
            %hold on
            hp(end+1)=plot(movmean(abs(prof_data{indmax}.(varname)(indprobe,:)),smoothing_window),prof_data{indmax}.P,'Color',col(kf,:),'LineWidth',1);
            hold on
        else
            hp(end+1)=plot(movmean(prof_data{indmax}.(varname)(indprobe,:),smoothing_window),prof_data{indmax}.P,'Color',col(kf,:),'LineWidth',1);
            hold on
        end
    end
end
end
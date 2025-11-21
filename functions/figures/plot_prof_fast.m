function [hp]=plot_prof_fast(DATA,filenames_all,filenames,col,select_profile,varname,smoothing_window,absolute)
% Plot all the profiles of a fast variable

hp=[];

if nargin<8
    absolute=false;
end

indprof=0;
for kf=1:length(filenames)
    ind_file=find(strcmp(filenames_all,filenames{kf}),1);
    prof_data=DATA(ind_file).FAST;
    if strcmp(select_profile,'all')
        for kprof=1:length(prof_data)
            indprof=indprof+1;
            if absolute
                %scatter(abs(prof_data{kprof}.(varname)),prof_data{kprof}.depth,'o','SizeData',3,'MarkerFaceColor',col(indprof,:),'MarkerEdgeColor','none','MarkerFaceAlpha',0.5);
                %hold on
                hp(end+1)=plot(movmean(abs(prof_data{kprof}.(varname)),smoothing_window),prof_data{kprof}.depth,'Color',col(indprof,:),'LineWidth',1);
                hold on
            else
                hp(end+1)=plot(movmean(prof_data{kprof}.(varname),smoothing_window),prof_data{kprof}.depth,'Color',col(indprof,:),'LineWidth',1);
                hold on
            end
        end
    else  % Plot the longest profile only
        depth_range=NaN(1,length(prof_data)); % [m]
        for kprof=1:length(prof_data)
            depth_range(kprof)=max(prof_data{kprof}.depth)-min(prof_data{kprof}.depth);
        end
        [~,indmax]=max(depth_range);
        if absolute
            %scatter(abs(prof_data{indmax}.(varname)),prof_data{indmax}.depth,'o','SizeData',3,'MarkerFaceColor',col(kf,:),'MarkerEdgeColor','none','MarkerFaceAlpha',0.1);
            %hold on
            hp(end+1)=plot(movmean(abs(prof_data{indmax}.(varname)),smoothing_window),prof_data{indmax}.depth,'Color',col(kf,:),'LineWidth',1);
            hold on
        else
            hp(end+1)=plot(movmean(prof_data{indmax}.(varname),smoothing_window),prof_data{indmax}.depth,'Color',col(kf,:),'LineWidth',1);
            hold on
        end
    end
end
end
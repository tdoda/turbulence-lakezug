function [hp]=plot_prof_bin(DATA,filenames_all,filenames,col,select_profile,varname,smoothing_window,flagname,absolute)
% Plot all the profiles of a binned variable

hp=[];

if nargin<8
    flagname='';
end

if nargin<9
    absolute=false;
end

indprof=0;
for kf=1:length(filenames)
    ind_file=find(strcmp(filenames_all,filenames{kf}),1);
    prof_data=DATA(ind_file).BINNED;
    if strcmp(select_profile,'all')
        for kprof=1:length(prof_data)
            % datavar=prof_data{kprof}.(varname);
            if ~isempty(flagname)
                %datavar(prof_data{kprof}.(flagname)==0)=NaN;
                datavar=prof_data{kprof}.(varname)(prof_data{kprof}.(flagname)==1);
                datadepth=prof_data{kprof}.depth(prof_data{kprof}.(flagname)==1);
            else
                datavar=prof_data{kprof}.(varname);
                datadepth=prof_data{kprof}.depth;
            end
            indprof=indprof+1;
            if absolute
                %scatter(abs(prof_data{kprof}.(varname)),prof_data{kprof}.depth,'o','SizeData',3,'MarkerFaceColor',col(indprof,:),'MarkerEdgeColor','none','MarkerFaceAlpha',0.5);
                %hold on
                hp(end+1)=plot(movmean(abs(datavar),smoothing_window),datadepth,'-o','Color',col(indprof,:),'MarkerSize',2,'MarkerFaceColor',col(indprof,:),'LineWidth',1);
                hold on
            else
                hp(end+1)=plot(movmean(datavar,smoothing_window,"omitnan"),datadepth,'-o','Color',col(indprof,:),'MarkerSize',2,'MarkerFaceColor',col(indprof,:),'LineWidth',1);
                hold on
            end
        end
    else  % Plot the longest profile only
        depth_range=NaN(1,length(prof_data)); % [m]
        for kprof=1:length(prof_data)
            depth_range(kprof)=max(prof_data{kprof}.depth)-min(prof_data{kprof}.depth);
        end
        [~,indmax]=max(depth_range);
        % datavar=prof_data{indmax}.(varname);
        if ~isempty(flagname)
            % datavar(prof_data{indmax}.(flagname)==0)=NaN;
            datavar=prof_data{indmax}.(varname)(prof_data{indmax}.(flagname)==1);
            datadepth=prof_data{indmax}.depth(prof_data{indmax}.(flagname)==1);
        else
            datavar=prof_data{indmax}.(varname);
            datadepth=prof_data{indmax}.depth;
        end
        if absolute
            %scatter(abs(prof_data{indmax}.(varname)),prof_data{indmax}.depth,'o','SizeData',3,'MarkerFaceColor',col(kf,:),'MarkerEdgeColor','none','MarkerFaceAlpha',0.1);
            %hold on
            hp(end+1)=plot(movmean(abs(datavar),smoothing_window),datadepth,'-o','Color',col(kf,:),'MarkerSize',2,'MarkerFaceColor',col(kf,:),'LineWidth',1);
            hold on
        else
            hp(end+1)=plot(movmean(datavar,smoothing_window,"omitnan"),datadepth,'-o','Color',col(kf,:),'MarkerSize',2,'MarkerFaceColor',col(kf,:),'LineWidth',1);
            hold on
        end
    end
end
end

function [depth_allprof,var_allprof,depth_cells,logmed_cells] = compute_logmedian(DATA,cell_size,varname,ind_file)

% ind_file is optional

% Returns median of the log of the variable between cell depth (i-1) and cell
% depth i
%%

if nargin<4
    ind_file=1:length(DATA);
end

%% Combined the data
var_allprof=[];
depth_allprof=[];

for k=1:length(ind_file)
    for kprof=1:length(DATA(ind_file(k)).BINNED)
        datavar=DATA(ind_file(k)).BINNED{kprof}.(varname);
        var_allprof=[var_allprof,datavar(~isnan(datavar))];
        depth_allprof=[depth_allprof,DATA(ind_file(k)).BINNED{kprof}.depth(~isnan(datavar))];
    end
end

%% Compute the median of the log
depth_start=0:cell_size:max(depth_allprof)-cell_size;
depth_end=cell_size:cell_size:max(depth_allprof);
depth_cells=sort([depth_start,depth_end]);
logmed_cells=NaN(size(depth_cells));
for k=1:length(depth_start)
    indcell=find(depth_cells==depth_start(k),1,'last');
    logmed_cells(indcell:indcell+1)=median(log10(var_allprof(depth_allprof>=depth_start(k)&depth_allprof<depth_end(k))));
end


function [axlim1_scaled,axlim2_scaled] = scale_limits(axlim1,axlim2)
%SCALE_LIMITS Change the axis limits to have the same scale between two
%axes.
%   Detailed explanation goes here
    axlim=[axlim1;axlim2];
    [~,indmax]=max([diff(axlim(1,:)),diff(axlim(2,:))]); % Get the axis with the largest limits difference = limits to keep
    axis_middle=axlim(3-indmax,1)+0.5*diff(axlim(3-indmax,:)); % Get the middle of the axis to modify
    axlim(3-indmax,:)=[axis_middle-0.5*diff(axlim(indmax,:)),axis_middle+0.5*diff(axlim(indmax,:))]; % Modify the limits for the axis with the smallest limits difference
    axlim1_scaled=axlim(1,:);
    axlim2_scaled=axlim(2,:);
    
end
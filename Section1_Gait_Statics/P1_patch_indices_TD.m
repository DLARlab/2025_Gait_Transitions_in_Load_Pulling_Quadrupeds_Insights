function [TD_XData,TD_YData] = P1_patch_indices_TD(touchdown_sequence)

% Extract patch vertices for touchdown plots
TD_index  = find(touchdown_sequence==1);
TD_segments  = find(diff(TD_index)~=1);

% Create XData for patching
TD_XData = [];
if ~isempty(TD_segments)
    for i = 1:length(TD_segments)+1
    
        if i==1
            segment = [TD_index(1); TD_index(TD_segments(i)); TD_index(TD_segments(i)); TD_index(1)];
        elseif  i>length(TD_segments)
            segment = [TD_index(TD_segments(end)+1); TD_index(end); TD_index(end); TD_index(TD_segments(end)+1)];
        else
            segment = [TD_index(TD_segments(i-1)+1); TD_index(TD_segments(i)); TD_index(TD_segments(i)); TD_index(TD_segments(i-1)+1)];
        end
    
        TD_XData = [TD_XData segment];
       
    end
else
    TD_XData = [TD_index(1); TD_index(end); TD_index(end); TD_index(1)];
end

% Create YData for patching
width = 0.10;

YData = [1-width;1-width;1+width;1+width];

TD_YData = ones(size(TD_XData)).*YData;

end


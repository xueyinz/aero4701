%% 440305585
% AERO4701
% Assignment 3
%
% get_dominant_rotation_axis.m

function dominant_rotation = get_dominant_rotation_axis(I, w_array)

    % check which principal axis there is dominant rotation
    if (w_array(1) > w_array(2)) && (w_array(1) > w_array(3))
        dominant_axis = 1;    % rotation about x-axis
        dominant_rotation.axis = 'x';
    elseif (w_array(2) > w_array(1)) && (w_array(2) > w_array(3))
        dominant_axis = 2;    % rotation about y-axis
        dominant_rotation.axis = 'y';
    elseif (w_array(3) > w_array(1)) && (w_array(3) > w_array(2))
        dominant_axis = 3;    % rotation about z-axis     
        dominant_rotation.axis = 'z';
    end
    
    % get ordered indices where idx = [min_axis, inter_axis, max_axis]
    [~, idx_moi] = sort([I.xx I.yy I.zz]);
    
    % get ordered indices where idx = [x_axis_rank, y_axis_rank, z_axis_rank]
    [~, idx_axes] = sort(idx_moi);
    moi_rank = idx_axes(dominant_axis);
    
    if moi_rank == 1
        dominant_rotation.I = 'minimum axis';
    elseif moi_rank == 2
        dominant_rotation.I = 'intermediate axis';
    elseif moi_rank == 3
        dominant_rotation.I = 'maximum axis';
    else
        error('Something wrong with the indexing logic in get_dominant_rotation_axis');
    end
    
end
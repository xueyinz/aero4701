function UAV_truth = load_UAV_position(text_file)

    % open text file
    fid = fopen(text_file);
    UAV_truth = struct;
    UAV_truth.time = [];
    UAV_truth.pos_LGCV = [];

    % loop through text file line by line
    while ~feof(fid)
        
        tline = fgetl(fid);
        tline = str2double(strsplit(tline, ' '));
        tline = tline(2:end);           % get rid of empty column
        
        % store uav data in struct
        UAV_truth.time = [UAV_truth.time tline(1)];
        UAV_truth.pos_LGCV = [UAV_truth.pos_LGCV tline(2:end)'];

    end
    
end
%% 440305585
% AERO4701
%
% input: [time (s), sat. #, a (m), pseudorange (m)]

function GPS_pseudorange = load_GPS_pseudorange(text_file)

    TIME = 1;       % time stamp
    SAT_NUM = 2;    % satellite number
    RANGE = 3;      % pseudorange

    fid = fopen(text_file);
    GPS_pseudorange = [];
    
    tline = fgetl(fid);
    tline = str2double(strsplit(tline, '  '));
    tline = tline(:, 2:end);                % get rid of empty column
    observation.time = tline(TIME);
    observation.sat_num = tline(SAT_NUM);
    observation.range = tline(RANGE);
    
    while ~feof(fid)
        
        tline = fgetl(fid);
        tline = str2double(strsplit(tline, '  '));
        tline = tline(:, 2:end);            % get rid of empty column
        
        if tline(TIME) ~= observation.time
            GPS_pseudorange = [GPS_pseudorange; observation];   % add the readings for the current timestamp
            % reset for the new observation
            observation.time = tline(TIME);
            observation.sat_num = tline(SAT_NUM);
            observation.range = tline(RANGE);
        else
            observation.sat_num = [observation.sat_num; tline(SAT_NUM)];
            observation.range = [observation.range; tline(RANGE)];
        end
        
    end
    
    GPS_pseudorange = [GPS_pseudorange; observation];   % add the readings for the current timestamp
end
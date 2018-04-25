%% 440305585
% AERO4701
%
% Load GPS pseudorange data from text file. Group the data into timestamp
% groups so that each group is the pseudorange observations seen at the
% unique timestamp.
%
% input text file columns: [time (s), sat. #, pseudorange (m)]
% output structure: struct
%   - each row is one timestamp group
%   - each row has field time, satellite number, and pseudorange
%   - time is a scalar
%   - satellite number and pseudorange are vectors; length = # observations

function GPS_pseudorange = load_GPS_pseudorange(text_file)

    % column constants
    TIME = 1;       % time stamp
    SAT_NUM = 2;    % satellite number
    PSEUDORANGE = 3;      % pseudorange

    % open text file
    fid = fopen(text_file);
    GPS_pseudorange = [];
    
    % get the first line of pseudorange data
    tline = fgetl(fid);
    tline = str2double(strsplit(tline, '  '));
    tline = tline(:, 2:end);                % get rid of the empty column
    
    % store the first line of pseudorange data
    observation.time = tline(TIME);
    observation.sat_num = tline(SAT_NUM);
    observation.range = tline(PSEUDORANGE);
    
    % loop through text file line by line
    while ~feof(fid)
        
        % get the next line of pseudorange data
        tline = fgetl(fid);
        tline = str2double(strsplit(tline, '  '));
        tline = tline(:, 2:end);            % get rid of the empty column
        
        % group the pseudorange data by time
        if tline(TIME) ~= observation.time
            
            % current line has a new timestamp; save the previous timestamp
            % group data
            GPS_pseudorange = [GPS_pseudorange; observation];
            
            % reset for the new timestamp group
            observation.time = tline(TIME);
            observation.sat_num = tline(SAT_NUM);
            observation.range = tline(PSEUDORANGE);
            
        else
            
            % same timestamp group; save data into temporary vectors
            observation.sat_num = [observation.sat_num; tline(SAT_NUM)];
            observation.range = [observation.range; tline(PSEUDORANGE)];
            
        end
        
    end
    
    % save the last timestamp group data
    GPS_pseudorange = [GPS_pseudorange; observation];   
    
end
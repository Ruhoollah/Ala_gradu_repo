function [match, num_match] = find_filematch_bynum(srch_num, srch_dir, exts)
    
    if nargin < 3
        exts = ''; 
    else
        exts = unique(exts, 'stable');
    end

    if ~isempty(exts) && ~iscell(exts) && ischar(exts)
        exts = {exts};
    end
    
    % find all files in srch_dir with the given extension(s)
    if isempty(exts)
        fs = dir(srch_dir);
    else
        fs = cellfun(@dir, fullfile(srch_dir, strcat('*', exts)), 'un', 0);
    end
    
    % concat all matching files to vertical struct array, with non-empty rows
    fs = vertcat(fs{~cellfun(@isempty, fs)});
    
    % use only the filename part of the found names
    [~, fnm, ~] = cellfun(@fileparts, {fs.name}, 'Un', 0);
    
    % Try completely numeric way
    fnums = cellfun(@(x) regexprep(x, '\D', ''), fnm, 'Un', 0);
    if ~isnumeric(srch_num)
        num_num = str2double(regexprep(srch_num, '\D', ''));
    end
    matches = str2double(fnums) == num_num;
    
    % if no unique match, try string match
    if sum(matches) ~= 1
        matches = contains(fnm, srch_num, 'IgnoreCase', 1);
    end
    
    if ~any(matches)
        match = '';
        num_match = 0;
        return
    end
    match = fullfile(fs(find(matches, 1)).folder, fs(find(matches, 1)).name);
    num_match = sum(matches);

end
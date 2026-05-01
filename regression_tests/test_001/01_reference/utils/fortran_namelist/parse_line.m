function [key, value] = parse_line(line)
    % Remove comments from the line
    %line = regexprep(line, '!!.*', '');

    % Remove comments that start with single '!' or '!!' anywhere in the line
    line = regexprep(line, '^\s*!.*|!!.*', '');

    % Split the line into key and value
    parts = strsplit(line, '=');
    key = strtrim(parts{1});
    value_str = strtrim(parts{2});

    % Use regular expression to capture everything up to the first '!' or to the end if '!' is absent
    matches = regexp(value_str, '^(.*?)\s*(?:!|$)', 'tokens');
    value_str = matches{1}{:};

    % Convert the value string to a MATLAB data type
    if startsWith(value_str, '''') && endsWith(value_str, '''')
        % String value
        value = strrep(value_str, '''', '');
    elseif contains(value_str, ',')
        % Array value
        value = str2num(value_str); %#ok<ST2NM>
    else
        % Single numeric value
        value = str2double(value_str);
        if isnan(value)
            value = strrep(value_str, '"', ''); % Keep as string if not a number
        end
    end
end
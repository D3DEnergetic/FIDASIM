function data = read_fortran_namelist(file_path)
    % Initialize an empty structure
    data = struct();

    % Open the file for reading
    fid = fopen(file_path, 'r');
    if fid == -1
        error('Cannot open the file.');
    end

    % Read the file line by line
    while ~feof(fid)
        line = fgetl(fid);

        % Skip comments, empty lines, and lines with whitespace before comments
        if isempty(line) || ~isempty(regexp(line, '^\s*[!&/]', 'once'))
            continue;
        end

        % Parse the line
        [key, value] = parse_line(line);

        % Assign the value to the structure
        if ~isempty(key)
            % Handle array elements
            array_expr = regexp(key, '(.*)\((\d+)\)', 'tokens');
            if ~isempty(array_expr)
                base_key = array_expr{1}{1};
                index = str2double(array_expr{1}{2});
                if ~isfield(data, base_key)
                    data.(base_key) = [];
                end
                data.(base_key)(index) = value;
            else
                data.(key) = value;
            end
        end
    end

    % Close the file
    fclose(fid);
end
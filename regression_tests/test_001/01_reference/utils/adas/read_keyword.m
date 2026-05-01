function val = read_keyword(line, key, fmt, default)
%READ_KEYWORD: Read value following KEY= from a line.
%   val = read_keyword(line, key, fmt, default)
%   - key     : string, e.g. 'SCREF'
%   - fmt     : sscanf format, e.g. '%f' or '%s'
%   - default : value to return if key not found ([] if required)

    token = [key '='];
    if contains(line, token)
        val = sscanf(extractAfter(line, token), fmt, 1);
    else
        if nargin < 4
            error('Missing %s=', key);
        else
            val = default;
        end
    end
end
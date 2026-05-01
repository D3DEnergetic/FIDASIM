function save_figure_to_file(output_dir, file_name, formats, pdfContentType)
%SAVE_FIGURE_TO_FILE Save the current figure to disk in one or more formats.
%
%   SAVE_FIGURE_TO_FILE(OUTPUT_DIR, FILE_NAME) saves the current figure as both PDF and
%   PNG into OUTPUT_DIR (created if it does not exist).
%
%   SAVE_FIGURE_TO_FILE(OUTPUT_DIR, FILE_NAME, FORMATS) saves only the requested
%   format(s). FORMATS may be:
%       - char/string scalar:          'pdf' or 'png'
%       - string array:               ["pdf","png"]
%       - cell array (mixed ok):       {'pdf','png'}
%
%   SAVE_FIGURE_TO_FILE(..., PDFCONTENTTYPE) controls how PDF content is written:
%       'auto'   -> uses 'image' internally (robust for dense plots)
%       'vector' -> vector PDF (best for simple line/marker plots)
%       'image'  -> raster image embedded in PDF (best for surf/contourf)
%
%   Behavior:
%     - If the figure contains exactly one non-legend, non-colorbar axes,
%       the axes are exported (tight framing, avoids clipping).
%     - Otherwise, the entire figure is exported (supports subplots/layouts).
%
%   Notes:
%     - Resolution is fixed at 600 DPI.
%     - Uses exportgraphics; PDF vectorization can be slow for dense graphics.
%
%   Example:
%     save_figure_to_file('output','lambda_map')                 % PDF + PNG
%     save_figure_to_file('output','line_plot','pdf','vector')   % vector PDF only
%
%   Notes:
%   - Next logical upgrades could include:
%       • Adding other formats like 'svg' or 'eps'
%       • Option to include timestamps in file names
%       • Option to specify figure handle instead of assuming gcf

    % Choose the target Handle:
    fig = gcf;
    ax = findall(fig,'Type','axes');
    ax = setdiff(ax, findall(fig,'Type','axes','Tag','legend'));
    ax = setdiff(ax, findall(fig,'Type','axes','Tag','Colorbar'));

    if numel(ax) == 1
        target = ax;
    else
        target = fig;
    end

    % Normalize inputs and set defaults:
    output_dir = char(output_dir);
    file_name  = char(file_name);
    if nargin < 3; formats = []; end
    formats = normalize_formats(formats);
    formats = unique(formats,'stable');
    if nargin < 4 || isempty(pdfContentType)
        pdfContentType = 'auto';  % 'auto' | 'vector' | 'image'
    end
    pdfContentType = char(pdfContentType);

    % Validate pdfContentType:
    validTypes = {'auto','vector','image'};
    if ~any(strcmpi(pdfContentType, validTypes))
        error('pdfContentType must be ''auto'', ''vector'', or ''image''');
    end

    % Validate formats:
    valid = {'pdf','png'};
    unknown = setdiff(formats, valid);
    if ~isempty(unknown)
        error('Unsupported format(s): %s. Valid: %s', ...
            strjoin(unknown, ', '), strjoin(valid, ', '));
    end

    % Ensure output directory exists:
    if ~exist(output_dir,'dir')
        mkdir(output_dir);
    end

    % Construct the base filename string
    figurePath = fullfile(output_dir, file_name);

    % Complete rendering figure:
    drawnow

    % PDF export:
    if any(strcmpi(formats, 'pdf'))
        if strcmpi(pdfContentType,'auto')
            ct = 'image';    % safest default for dense plots
        else
            ct = pdfContentType;
        end

        exportgraphics(target, [figurePath '.pdf'], ...
            'Resolution', 600, ...
            'ContentType', ct);
    end

    % PNG export:
    if any(strcmpi(formats,'png'))
        exportgraphics(target, [figurePath '.png'], ...
            'Resolution', 600);
    end
end

function formats = normalize_formats(formats)

    if nargin == 0 || isempty(formats)
        formats = {'pdf','png'};
        return
    end

    % If user passed a single string/char (e.g. "pdf" or 'pdf')
    if isstring(formats) || ischar(formats)
        formats = cellstr(formats);          % -> {'pdf'} or {'pdf','png'} if string array
    elseif iscell(formats)
        % If user passed a cell array, make sure every element is char
        formats = cellfun(@char, formats, 'UniformOutput', false);
    else
        error('formats must be char, string, string array, or cell array.');
    end

    % Normalize case + whitespace
    formats = cellfun(@(s) lower(strtrim(s)), formats, 'UniformOutput', false);

    % Remove empties (optional)
    formats = formats(~cellfun('isempty', formats));
end
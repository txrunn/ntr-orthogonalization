classdef visualizeDistanceMatrix

    methods (Static)

        function visualizeAndSave(allDistanceMatrices, wordlist)
            % This function is the main entry point. It generates color-coded matrices
            % from the input data, creates a heatmap for each of them and saves
            % the resulting heatmaps into an Excel file.

            % Generate the color-coded matrices
            colorCodedMatrices = visualizeDistanceMatrix.generateColorCodedMatrices(allDistanceMatrices);

            % Iterate through each color-coded matrix and create a heatmap
            % The heatmap will not be visible if you simply run this script,
            % because it will be directly saved to the Excel file.
            for i = 1:length(colorCodedMatrices)
                hFig = visualizeDistanceMatrix.createAndSaveHeatmaps(colorCodedMatrices(i), wordlist);

                % Define the filename, sheetname, and cell location to save the heatmap in Excel
                filename = 'distanceMatrix.xlsx';
                sheetname = ['Sheet' num2str(i)];
                xlcell = 'A1';

                visualizeDistanceMatrix.xlswritefig(hFig, filename, sheetname, xlcell);
            end

        end

        function colorCodedMatrices = generateColorCodedMatrices(filename)
            % This function loads distance matrices from a file and transforms them
            % into a format that can be visualized.

            % Load the distance matrices
            load(filename, 'allMatrices');
            matrix = struct2cell(allMatrices);
            matrixNames = fieldnames(allMatrices);

            % Initialize the output structure
            colorCodedMatrices(length(matrix)).matrix = [];
            colorCodedMatrices(length(matrix)).name = '';

            % Transform the matrices for visualization
            for k = 1:length(matrix)
                letterTriTrans = transpose(matrix{k});
                colorCodedMatrices(k).matrix = letterTriTrans;
                colorCodedMatrices(k).name = matrixNames{k};
            end

        end

        function hFig = createAndSaveHeatmaps(colorCodedMatrix, wordlist)
            % This function creates a heatmap for each color-coded matrix.
            % Please note that the heatmap will not be visible if you simply run this script,
            % because it will be directly saved to the Excel file.
            % If you want to visually inspect the heatmaps, you can add a line to display the figure (like figure(hFig)),
            % or save the figure to an image file.

            % Create a new figure
            hFig = figure;

            % Create the heatmap
            heatmap(colorCodedMatrix.matrix, 'Colormap', parula);

            % Set the x and y data for the axes
            set(gca, 'XData', wordlist, 'YData', wordlist);

            % Set the fontsize
            set(gca, 'FontSize', 8);

            % Set the title
            title(colorCodedMatrix.name);
        end

        function xlswritefig(hFig, filename, sheetname, xlcell)
            % XLSWRITEFIG  Write a MATLAB figure to an Excel spreadsheet
            %
            % xlswritefig(hFig,filename,sheetname,xlcell)
            %
            % All inputs are optional:
            %
            %    hFig:      Handle to MATLAB figure.  If empty, current figure is
            %                   exported
            %    filename   (string) Name of Excel file, including extension.  If not specified, contents will
            %                  be opened in a new Excel spreadsheet.
            %    sheetname:  Name of sheet to write data to. The default is 'Sheet1'
            %                       If specified, a sheet with the specified name must
            %                       exist
            %    xlcell:     Designation of cell to indicate the upper-left corner of
            %                  the figure (e.g. 'D2').  Default = 'A1'
            %
            % Requirements: Must have Microsoft Excel installed.  Microsoft Windows
            % only.
            %
            % Ex:
            % Paste the current figure into a new Excel spreadsheet which is left open.
            %         plot(rand(10,1))
            %         drawnow    % Maybe overkill, but ensures plot is drawn first
            %         xlswritefig
            %
            % Specify all options.
            %         hFig = figure;
            %         surf(peaks)
            %         xlswritefig(hFig,'MyNewFile.xlsx','Sheet2','D4')
            %         winopen('MyNewFile.xlsx')
            % Michelle Hirsch
            % The MathWorks
            % mhirsch@mathworks.com
            %
            % Is this function useful?  Drop me a line to let me know!
            if nargin == 0 || isempty(hFig)
                hFig = gcf;
            end

            if nargin < 2 || isempty(filename)
                filename = '';
                dontsave = true;
            else
                dontsave = false;

                % Create full file name with path
                filename = fullfilename(filename);
            end

            if nargin < 3 || isempty(sheetname)
                sheetname = 'Sheet1';
            end;

            if nargin < 4
                xlcell = 'A1';
            end

            % Put figure in clipboard
            if ~verLessThan('matlab', '9.8')
                copygraphics(hFig)
            else
                % For older releases, use hgexport. Set renderer to painters to make
                % sure it looks right.
                r = get(hFig, 'Renderer');
                set(hFig, 'Renderer', 'Painters')
                drawnow
                hgexport(hFig, '-clipboard')
                set(hFig, 'Renderer', r)
            end

            % Open Excel, add workbook, change active worksheet,
            % get/put array, save.
            % First, open an Excel Server.
            Excel = actxserver('Excel.Application');
            % Two cases:
            % * Open a new workbook, save with given file name
            % * Open an existing workbook
            if exist(filename, 'file') == 0
                % The following case if file does not exist (Creating New File)
                op = invoke(Excel.Workbooks, 'Add');
                %     invoke(op, 'SaveAs', [pwd filesep filename]);
                new = 1;
            else
                % The following case if file does exist (Opening File)
                %     disp(['Opening Excel File ...(' filename ')']);
                op = invoke(Excel.Workbooks, 'open', filename);
                new = 0;
            end

            % set(Excel, 'Visible', 0);
            % Make the specified sheet active.
            try
                Sheets = Excel.ActiveWorkBook.Sheets;
                target_sheet = get(Sheets, 'Item', sheetname);
            catch %#ok<CTCH>   Suppress so that this function works in releases without MException
                % Add the sheet if it doesn't exist
                target_sheet = Excel.ActiveWorkBook.Worksheets.Add();
                target_sheet.Name = sheetname;
            end;

            invoke(target_sheet, 'Activate');
            Activesheet = Excel.Activesheet;
            % Paste to specified cell
            Paste(Activesheet, get(Activesheet, 'Range', xlcell, xlcell))
            % Save and clean up
            if new && ~dontsave
                invoke(op, 'SaveAs', filename);
            elseif ~new
                invoke(op, 'Save');
            else % New, but don't save
                set(Excel, 'Visible', 1);
                return % Bail out before quitting Excel
            end

            invoke(Excel, 'Quit');
            delete(Excel)
        end

        function filename = fullfilename(filename)
            [filepath, filename, fileext] = fileparts(filename);

            if isempty(filepath)
                filepath = pwd;
            end

            if isempty(fileext)
                fileext = '.xlsx';
            end

            filename = fullfile(filepath, [filename fileext]);
        end

    end

end

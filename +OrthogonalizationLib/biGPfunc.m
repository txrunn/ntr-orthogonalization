%When function is included in orthog code, set up input so it eliminates
%the same words as parameters and jpglove inputs
function bigpTable = biGPfunc(biGPorig)
    biGP = table2array(biGPorig);
    biGP = regexprep(biGP, '^\$', ''); %If a cell starts with $, replace it with nothing
    biGP = regexprep(biGP, '\$$', ''); %If a cell ends with $, replace it with nothing
    [numrows, numcol] = size(biGP); %eliminate +/ or /+ in monGP syllables

    for i = 1:numrows

        for j = 2:numcol - 1

            if not(strcmp(biGP{i, j}, '+/+'))

                if strcmp(biGP{i, j - 1}, '+/+')

                    if strcmp(biGP{i, j + 1}, '+/+')
                        biGP{i, j} = regexprep(biGP{i, j}, '/\+$', '');
                        biGP{i, j} = regexprep(biGP{i, j}, '^\+/', '');
                    end

                end

            end

        end

    end

    biGP = regexprep(biGP, '^\+/\S*', '+/+'); %Replace all non-white character space strings starting with + with +/+
    biGP = regexprep(biGP, '\S*\+$', '+/+'); %Replace all non-white character space strings starting with + with +/+
    biGP = regexprep(biGP, '^#', ''); %If a cell starts with #, replace it with nothing
    bigpTable = regexprep(biGP, '+/\+', '+'); %Replace all +/+ with + to match gpTable
    bigpTable = array2table(bigpTable);
end

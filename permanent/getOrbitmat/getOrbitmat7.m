% This file makes 'orbitmat7.dat' from 'orbitmatRedundant7part1.dat' .. 'orbitmatRedundant7part40.dat'.
for i = 1:40
    % Concatenate parts
    if i == 1
        orbitmatRedundant = readmatrix(strcat('orbitmatRedundant7part', int2str(i), '.dat'));
    else
        orbitmatRedundant = cat(1, orbitmatRedundant, readmatrix(strcat('orbitmatRedundant7part', int2str(i), '.dat')));
    end

    % Check if the file 'orbitmatRedundant7part{i}.dat' is complete.
    [numrow,~] = size(orbitmatRedundant);
    if i < 40
        if numrow ~= 500*i
            error(strcat('orbitmatRedundant7part', int2str(i), '.dat is incomplete.'));
        end
    elseif i == 40
        if numrow ~= 19940
            error(strcat('orbitmatRedundant7part', int2str(i), '.dat is incomplete.'));
        end
    end
    disp(i);
end

[numrow,numcol] = size(orbitmatRedundant);
orbitmatRedundant = orbitmatRedundant(:,2:numcol); % eliminate the first column which are just row indices

% Eliminate redundant rows
orbitmat = orbitmatRedundant(1,:);
firstNumListInOrbits = orbitmatRedundant(1,1);
for i=2:numrow
    if ~ismember(orbitmatRedundant(i,1),firstNumListInOrbits)
        orbitmat = cat(1, orbitmat, orbitmatRedundant(i,:));
        firstNumListInOrbits = cat(1, firstNumListInOrbits, orbitmatRedundant(i,1));
    end
    disp(i);
end

% Write to a file
writematrix(orbitmat,'orbitmat7.dat','Delimiter','space');
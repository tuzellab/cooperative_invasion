% Mean-squared displacement calculation code
% Written by Tüzel Lab
% Temple University, 2022

directory = uigetdir();
maxSteps = input('What is the maximum number of time steps in your data?');
numberOfColumns = input('What is the number of columns in the data files?');
trajectoryFiles = dir(directory);
trajectoryFiles([trajectoryFiles.isdir]) = []; % remove directories

msd = zeros(maxSteps, numberOfColumns);
counts = zeros(maxSteps, 1);

for trajectoryFile = trajectoryFiles'
    trajectory = dlmread(strcat(directory, '/', trajectoryFile.name), '\t');
    tasd = zeros(size(trajectory, 1), size(trajectory, 2));
    % For every possible pair of data points, calculate the displacement
    % between those points and add the squared displacement to the element
    % of tasd that represents the separation between those data points in
    % time.
    for i = 1:size(trajectory, 1)
        tasd(i, 1) = trajectory(i, 1) - trajectory(1, 1);
        for j = i:size(trajectory, 1)
            for k = 2:size(trajectory, 2)
                tasd(j - i + 1, k) = tasd(j - i + 1, k) + (trajectory(j, k) - trajectory(i, k)) * (trajectory(j, k) - trajectory(i, k));
            end
        end
    end
    for i = 1:size(trajectory, 1)
        msd(i, 1) = tasd(i, 1);
        for k = 2:size(trajectory, 2)
            msd(i, k) = msd(i, k) + tasd(i, k) / (size(trajectory, 1) + 1 - i);
        end
        counts(i) = counts(i) + 1;
    end
end

for i = 1:size(msd, 1)
    for k = 2:size(msd, 2)
        msd(i, k) = msd(i, k) / counts(i);
    end
end

dlmwrite(strcat(directory, '/msd.txt'), msd, '\t');
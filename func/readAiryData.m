function imgAiry = readAiryData(filePath,set)

dataPath = [filePath,'Raw data of each detector\E1\'];
dataOutput = dir(fullfile(dataPath,'*.mat'));
dataNames = {dataOutput.name};
dataName = dataNames{contains(dataNames,['set',num2str(set-1)])};
load([dataPath,dataName]);

for ii = 1:19
    eval(sprintf('imgAiry(:,:,ii) = IN0_n0_set%d_d%d;',set-1,ii-1));
end

end
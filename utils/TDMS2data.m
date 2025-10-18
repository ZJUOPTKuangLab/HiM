function [img,pixelsize] = TDMS2data(filePath) 

% img(imageRes,imageRes,19,stack)

[finalOutput,metaStruct] = TDMS_readTDMSFile([filePath,'.tdms']);
IN = finalOutput.propValues{1, 5}{5};
n = finalOutput.propValues{1, 5}{2};
imgRes = finalOutput.propValues{1, 4}{1};
pixelsize = finalOutput.propValues{1, 4}{4};

img = zeros(imgRes,imgRes,19,IN*n);
num = 0;
for kr = 1:IN*n
    DataPos = contains(metaStruct.chanNames,'set0')&contains(metaStruct.groupNames,['IN',num2str(kr-1),'-E1']);
    AiryPos = metaStruct.chanNames(DataPos);
    [~,AiryIndex] = sort_nat(AiryPos,'ascend');
    DataNow = finalOutput.data(DataPos);
    DataNow = DataNow(AiryIndex);
    num = num+1;
    for i = 1:length(DataNow)
        name = AiryPos{AiryIndex(i)};
%         nn = regexp(name, 'n(\d+)', 'tokens');
        dd = regexp(name, 'd(\d+)', 'tokens');
%         nn_value = 1 + str2double(nn{1});
        dd_value = 1 + str2double(dd{1});
        temp = DataNow{i};
        imgtemp = reshape(temp,imgRes,imgRes);
        img(:,:,dd_value,num) = imgtemp;
    end
end
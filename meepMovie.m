function meepMovie(fileName,component)

component = strcat('/',component);

data = readMeepH5Component(fileName,component);
maxValue = max(data(:));
minValue = min(data(:));


figure;
for cnt = 1:size(data,1)
    temp = squeeze(data(cnt,:,:));
    plotFieldComponent(temp,cnt,size(data,1),minValue,maxValue)
    pause(0.1);
end


end


function plotFieldComponent(data,cnt,len,minValue,maxValue)
imagesc(data)
colormap(jet(256))
colorbar;
tstring = sprintf('%d of %d',cnt,len);
title(tstring);
grid on
% caxis([minValue maxValue])

end
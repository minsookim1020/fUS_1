tic
[B,A]=butter(4,0.3,'high');    %coefficients for the high pass filter

    % sustraction of the first image
    % the signal stat at 0 and minimises filter oscilatons
    for i=1:size(frame_007,1)
        for j=1:size(frame_007,2)
            frame_007(i,j,:)=frame_007(i,j,:)-frame_007(i,j,1); 
        end
    end

  sb=filter(B,A,frame_007,[],3);    % blood signal (filtering in the time dimension)
    sb=sb(:,:,5:end);           % the first 4 temporal samples are eliminates (filter oscilations)
    PDI=mean(abs(sb).^2,3);     % computing the intensity of the blood signal the 
toc 
    
    
    
    % subplot(1,2,2);
imagesc(10*log10(PDI./max(PDI(:))));
% imagesc(10*log10(PDI./max(max(PDI(:,100:175))))); 
% caxis([-35 0]);
caxis auto
colormap gray;
title('PDI image')
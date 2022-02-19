%%

IQR = squeeze(Data.D.IQR);

%%
frame_030 = IQR (:,:,:,030);

[B,A]=butter(5,0.6,'high');    %coefficients for the high pass filter

    % sustraction of the first image
    % the signal stat at 0 and minimises filter oscilatons
    for i=1:size(frame_030,1)
        for j=1:size(frame_030,2)
            frame_030(i,j,:)=frame_030(i,j,:)-frame_030(i,j,1); 
        end
    end

  sb=filter(B,A,frame_030,[],3);    % blood signal (filtering in the time dimension)
    sb=sb(:,:,5:end);           % the first 4 temporal samples are eliminates (filter oscilations)
    PDI=mean(abs(sb).^2,3);     % computing the intensity of the blood signal the 
 
    
    
%%    
    % subplot(1,2,2);
imagesc(10*log10(PDI./max(PDI(:))));
% imagesc(10*log10(PDI./max(max(PDI(:,100:175))))); 
% caxis([-35 0]);
caxis auto
colormap gray;
title('Please')
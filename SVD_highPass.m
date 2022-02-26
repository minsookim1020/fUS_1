%% load and squeeze IQR data from .mat file
load (file_name, 'Data') ; %file_name is that what you want to do
IQR = squeeze(Data.D.IQR) ; % IQR is 4-D, raw data form I saved. (depth x width x # of CUIs x # of frames)

for fr = 1:size(IQR,4); 
    frame = IQR (:,:,:,fr);
    [nz,nx,nt]=size(frame);
    [nz,nx,nt]=size(frame);
    m=zeros(nx*nz,nt);
    i=0;
    for ix=1:nx
        for iz=1:nz
            i=i+1;
            m(i,:)=squeeze(frame(iz,ix,:)-mean(frame(iz,ix,:)));
        end
    end

    [U,S,V]=svd(m,'econ'); SVD
    
    nfilt=50;
    Sf=S;
    Sf(1:nfilt,1:nfilt)=0;
    mf=U*Sf*V';

    i=0;
    IQR1=zeros(nz,nx,nt);
    for ix=1:nx
        for iz=1:nz
            i=i+1;
            IQR1(iz,ix,:)=squeeze(mf(i,:));
        end
    end
    
    [B,A]=butter(4,0.6,'high'); % 4th Butter, 0.6=30Hz/100*2 
    IQR2=filter(B,A,IQR1,[],3);
    IQR2=IQR2(:,:,5:end);           % the first 4 temporal samples are eliminates (filter oscilations)
    PDI=mean(abs(IQR2).^2,3);     % computing the intensity of the blood signal 
    
    subplot (5/fr, 5, fr)
    imagesc(10*log10(PDI./max(PDI(:))));
    caxis([-15 0]); % or caxis auto
    colormap gray;
title(fr)

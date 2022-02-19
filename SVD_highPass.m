tic
[nz,nx,nt]=size(frame_007);

    m=zeros(nx*nz,nt);
    i=0;
    for ix=1:nx
        for iz=1:nz
            i=i+1;
            m(i,:)=squeeze(frame_007(iz,ix,:)-mean(frame_007(iz,ix,:)));
        end
    end

    [U,S,V]=svd(m,'econ');
    
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
    
    [B,A]=butter(4,0.3,'high');
    IQR2=filter(B,A,IQR1,[],3);
    IQR2=IQR2(:,:,5:end);           % the first 4 temporal samples are eliminates (filter oscilations)
    PDI=mean(abs(IQR2).^2,3);     % computing the intensity of the blood signal the 
toc    
        % subplot(1,2,2);
imagesc(10*log10(PDI./max(PDI(:))));
% imagesc(10*log10(PDI./max(max(PDI(:,100:175))))); 
% caxis([-35 0]);
caxis auto
colormap gray;
title('PDI image')

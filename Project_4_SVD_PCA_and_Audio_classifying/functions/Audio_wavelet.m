function Audio_Data = Audio_wavelet(audio_csv_file)

    [m,n] = size(audio_csv_file)
    pxl = (sqrt(m))
    nw = round(m/4) % wavelet resolution
    Audio_Data = zeros(nw,n);
    
    for k = 1:n
       X = im2double(reshape(audio_csv_file(:,k),pxl,pxl)); 
       [~,cH,cV,~]=dwt2(X,'haar');
       cod_cH1 = rescale(abs(cH))
       cod_cV1 = rescale(abs(cV))
       cod_edge = cod_cH1+cod_cV1
       Audio_Data(:,k) = reshape(cod_edge,nw,1);
    end

end


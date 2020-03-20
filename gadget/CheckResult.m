function [ SNR1,SNR2 ] = CheckResult( HS,MS,X,R,b )
%CHECKRESULT 
RX = R*img2mat(X);
XBS = DownsampleIMG(BlurIMG(X,b),length(MS)/length(HS));
XBS = img2mat(XBS);
Nh = XBS-img2mat(HS);
Nm = RX-img2mat(MS);
SNR1 = snr(XBS,Nh);
SNR2 = snr(RX,Nm);
end


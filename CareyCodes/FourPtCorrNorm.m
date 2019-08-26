function [FourCorr,FourCorrCombo,FourCorrNorm,FourCorrComboNorm] = FourPtCorrNorm(f,tau1range,tau2,tau3range,res)
tic
FourCorr = zeros(length(tau1range),length(tau3range));
FourCorrCombo = zeros(length(tau1range),length(tau3range));

fsiz = size(f);

if fsiz(2) >= 2
f=f';
end

for a = 1:length(tau1range)
    tau1 = tau1range(a);
    
parfor b = 1:length(tau3range)
    tau3 = tau3range(b);
    
    f2 = circshift(f,-round(tau1/res));
    f3 = circshift(f2,-round(tau2/res));
    f4 = circshift(f3,-round(tau3/res));
    
    ft = f(1:length(f)-round((tau3+tau2+tau1)/res));
    f2 = f2(1:length(ft));
    f3 = f3(1:length(ft));
    f4 = f4(1:length(ft));
    
% sep = (tau1+tau2+tau3)/res;
% 
% FourProduct = zeros(length(f)-sep,1);
% for k = 1:length(f)-sep
%     FourProduct(k) = f(k)*f(k+tau1/res)*f(k+tau1/res+tau2/res)*f(k+tau1/res+tau2/res+tau3/res);
% end

    FourCorrCombo(a,b) = mean(ft.*f2.*f3.*f4);
    FourCorr(a,b) = FourCorrCombo(a,b) - (mean(ft.*f2))*mean((f3.*f4));

end
a
end
fnew = circshift(f,-round(tau2/res));
FourCorrNorm = FourCorr/(mean(f.^2.*fnew.^2));
FourCorrComboNorm = FourCorrCombo/(mean(f.^2.*fnew.^2));
toc
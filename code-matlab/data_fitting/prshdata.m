%plot SF versus Hb

figure(1)
subplot(1,2,1)
[hb2,ix]=sort(hb);
plot(hb,sf,'.', hb2, smooth(sf(ix),5),'-')
xlabel('hemoglobin (g/dL)')
ylabel('serum ferritin (ug/L)')

%plot stfr versus Hb
subplot(1,2,2)
plot(hb,stfr,'.',hb2, smooth(stfr(ix),5),'-')
xlabel('hemoglobin (g/dL)')
ylabel('serum transferritin receptor (mg/L)')

%plot Hb versus absorp
figure(2)
subplot(1,3,1)
plot(hb,absorp,'.', hb2, smooth(absorp(ix),5),'-')
xlabel('hemoglobin (g/dL)')
ylabel('iron absorption (%)')

%plot SF versus absorp
subplot(1,3,2)
[sf2,ix]=sort(sf);
plot(sf,absorp,'.', sf2, smooth(absorp(ix),5),'-')
xlabel('serum ferritin (ug/L)')
ylabel('iron absorption (%)')

%plot stfr versus absorp
subplot(1,3,3)
[stfr2,ix]=sort(stfr);
plot(stfr,absorp,'.',stfr2, smooth(absorp(ix),5),'-')
xlabel('serum transferritin receptor (mg/L)')
ylabel('iron absorption (%)')
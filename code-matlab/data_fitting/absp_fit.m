

figure(4)
bins=[6 7 8 9 9.5 10 10.5 11 11.5 12 12.5 13 14 15];
nhb=histc(hb,bins);
bar(bins, nhb)

%get regression for Hb versus absorption


[hb2,ix]=sort(hb);
absorp2=absorp(ix);

ind=find(10.5<hb2 & hb2<12.5);
hb3=hb2(ind);
absorp3=absorp2(ind);

%choose subset of data to fit line to

whichstats = {'beta','rsquare','r','tstat','fstat','yhat'};

stats=regstats(absorp3,hb3,'linear',whichstats);
df=stats.tstat.dfe;
tci=tinv(0.975,df);
ha_b=stats.beta
ha_berr=stats.tstat.se
ha_bp=stats.tstat.pval
ha_r=stats.r;  %residuals (I think)
ha_f=stats.fstat.pval
ha_ypred=stats.yhat;
rsquare=stats.rsquare
ha_SSE=sum(ha_r.^2)


%plot Hb versus absorp

x=10.5:0.1:12.5;

figure(10)
hold on
plot(x,ha_b(1)+ha_b(2)*x,'-r')
%{
plot(hb,absorp,'.', hb2, smooth(absorp2,5),'-',x,ha_b(1)+ha_b(2)*x,'-r')
xlabel('hemoglobin (g/dL)')
ylabel('iron absorption (%)')
legend('data','smoothed data',sprintf('slope=%4.2f intercept=%5.2f \n R^2=%4.2f p=%7.5f',...
    ha_b(2),ha_b(1),rsquare,ha_bp(2)))
%}
size(ind)

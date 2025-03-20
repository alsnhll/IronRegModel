
fix=find(10.5<=hb);
hb2=hb(fix);
absorp2=absorp(fix);

ind=find(0<absorp2 & absorp2<5);
hb5=mean(hb2(ind));
hb5err=std(hb2(ind));

ind=find(5<=absorp2 & absorp2<10);
hb10=mean(hb2(ind));
hb10err=std(hb2(ind));

ind=find(10<=absorp2 & absorp2<15);
hb15=mean(hb2(ind));
hb15err=std(hb2(ind));

ind=find(15<=absorp2 & absorp2<20);
hb20=mean(hb2(ind));
hb20err=std(hb2(ind));

ind=find(20<=absorp2 & absorp2<25);
hb25=mean(hb2(ind));
hb25err=std(hb2(ind));

hb_vec=[hb5 hb10 hb15 hb20 hb25];
hberr_vec=[hb5err hb10err hb15err hb20err hb25err];

absorp2_vec=[2.5:5:22.5];

figure(10)
herrorbar(hb_vec,absorp2_vec,hberr_vec,hberr_vec)


function a=absp(hb,conv)
        
        %s=0.15; %slope of absorption versus fractional HB count. baseline 0.15
        %a=1/30-s*(hb-2.5)/2.5; %regular absp
        
        %sigmoid absp
        B=0.1; %"fermi temperature" of absoprtion versus fractional HB count. baseline 0.1
        %mina=0.01;maxa=0.2;norma=0.06; MALE VALUES
        
        %to define based on three absoprtion values: maximum, minimum and a
        %value of at normal Hb of 14
        mina=0.01;maxa=0.20;norma=0.03;a=mina+(maxa-mina)./(1+exp((hb-13/conv+B*log((maxa-norma)/(norma-mina)))/B));
        
        %to define based on three absorption values: maximum, minimum and
        %Hb value where half max occurs
         %mina=0.01;maxa=0.25;hb_mid=11;a=mina+(maxa-mina)./(1+exp((hb-hb_mid/5.6)/B));
        
        %inverse absp (following Lao paper)
        %responding to Hb(y)
        %mina=0.01;maxa=0.2;norma=0.06;a=-(maxa-mina)*hb./((2.5*(norma-mina)/(maxa-norma))+hb)+maxa;
        
        %responding to body Fe (x)
        %mina=0.01;maxa=0.2;norma=0.06;a=-(maxa-mina)*hb./((1.5*(norma-mina)/(maxa-norma))+1.5)+maxa;
        
        %a=0.06-0.14*(hb-2.5)/2.5; %Anura's absp
        %a=0.06; %constant
    end
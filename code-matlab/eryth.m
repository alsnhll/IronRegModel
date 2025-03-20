function h=eryth(hb,h0,conv)
        
        %h=h0; %constant
        %h=(3/4*(2.5-hb)+1)*h0; %linearly increases as Hb fals
        
        %dose response curve
        %B=0.2;
        B=0.062;
        mina=0.7*h0;
        %maxa=4*h0;
        maxa=5*h0;
        norma=h0;h=mina+(maxa-mina)./(1+exp((hb-13/conv+B*log((maxa-norma)/(norma-mina)))/B));

        %piecewise linear
        %{
        h=zeros(size(hb));
        hb_min=6; %hb value where minimum eryth occurs
        h_slope=4; %factor by which h increases between h_min and h_max
        for i=1:1:length(hb)
            if hb(i)*5.6<hb_min
                h(i)=h_slope*h0;
            elseif hb(i)*5.6>14
                h(i)=h0;
            else
                h(i)=((h_slope-1)/((14-hb_min)/5.6)*(2.5-hb(i))+1)*h0; %linearly increases as Hb falls
            end
        end
        %}
        
    end
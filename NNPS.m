function [PAIR_I,PAIR_J,NIAC_o,W_o, DWDX_o]=NNPS(x,N,h)
    NIAC=0;
    COUNTIAC(1:N)=0;
    for i=1:(N-1)
        for j=(i+1):N
            if abs(x(i)-x(j))<=2.0*(h(i)+h(j))/2
                NIAC=NIAC+1;
                PAIR_i(NIAC)=i;
                PAIR_j(NIAC)=j;
                hij=(h(i)+h(j))/2;
                R=abs(x(i)-x(j))/hij;
                if R<1
                    W(NIAC)=(2/3-R^2+R^3/2)/hij;                                                      
                    DWDX(NIAC)=(-2+3*R/2)*(x(i)-x(j))/hij^3;
                elseif R<2
                    W(NIAC)=(2-R)^3/6/hij;
                    DWDX(NIAC)=-(2-R)^2*(x(i)-x(j))/hij^2/abs(x(i)-x(j))/2;
                elseif R>=2
                    W(NIAC)=0;
                    DWDX(NIAC)=0;
                end
            end
        end
    end
    PAIR_I=PAIR_i;
    PAIR_J=PAIR_j;
    NIAC_o=NIAC;
    W_o=W;
    DWDX_o=DWDX;
end
function [Tcc,Tuu] = chris_mode_to_mode_transfer(u,w,c)

[nz,nx] = size(c);

n=0;
    for k = 0.02:0.02:0.4
        n=n+1;
        m=0;
        for q = 0.02:0.02:0.4
            m=m+1;
            for i = 1:nz
                uk(i,:) = fftf([1:nx], squeeze(u(i,:)),k,k-0.02);
                wk(i,:) = fftf([1:nx], squeeze(w(i,:)),k,k-0.02);
                ck(i,:) = fftf([1:nx], squeeze(c(i,:)),k,k-0.02);
            
                uq(i,:) = fftf([1:nx], squeeze(u(i,:)),q,q-0.02);
                wq(i,:) = fftf([1:nx], squeeze(w(i,:)),q,q-0.02);
                cq(i,:) = fftf([1:nx], squeeze(c(i,:)),q,q-0.02);
            end
        
            [dcqz,dcqx] = chris_gradient(cq,1,1);
            [duqz,duqx] = chris_gradient(uq,1,1);
            [dwqz,dwqx] = chris_gradient(wq,1,1);
        
            Tcc1 = ck.*u.*dcqx +ck.*w.*dcqz;
            Tuu1 = uk.*u.*duqx+uk.*w.*duqz+wk.*u.*dwqx + wk.*w.*dwqz;
        
            Tcc(n,m) = trapz([1:nz],trapz([1:nx],Tcc1,2));
            Tuu(n,m) = trapz([1:nz],trapz([1:nx],Tuu1,2));
        
        end
    end

end


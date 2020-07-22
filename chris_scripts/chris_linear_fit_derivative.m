function [dfdx] = chris_linear_fit_derivative(f,x,n)

[nz,nx] = size(f);



for j = 1:nz
    
    a = smooth(squeeze(f(j,:)),n);
    
    for i = 1:nx-n
        linefit=polyfit(x(i:i+n),smooth(a(i:i+n))',1);
        dfdx(j,i)=linefit(1);
    end
end

end


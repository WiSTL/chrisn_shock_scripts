function [P,x,moments] = chris_pdf_moments(U,interval)

[P,x] = chris_pdf(U,interval);

mean = trapz(x,x.*P);
variance = trapz(x,((x-mean).^2).*P);
skewness = trapz(x,((x-mean).^3).*P)/variance^(3/2);
kurtosis = trapz(x,((x-mean).^4).*P)/variance^(2);

moments(1) = mean;
moments(2) = variance;
moments(3) = skewness;
moments(4) = kurtosis;

end


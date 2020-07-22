function [C,Ec,Ch,phi] = IC_reconstruct(Ec,cp2,cm)

nz = length(cp2);
nx = length(Ec);

Ec = smooth(Ec./trapz([1:nx],Ec))/2;

for i = 1:nx
    E(:,i) = cp2.*Ec(i);
    Cm(:,i) = cm;
    
    if i<nx/2
        phi(:,i) = 0*cm+0.666*(1 - i*(4/nx));
    else
        phi(:,i) = 0*cm+0.666*(1 - (i-nx/2)*(4/nx));
    end
    
end

%phi = 0*Cm;%pi*(2*rand(nz,nx)-1);

Ch = sqrt(E.*exp(-1j*phi));

cp = 0.5*(fft(Ch')' + conj(fft(Ch')'));

C = Cm+cp;

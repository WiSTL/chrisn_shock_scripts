function [] = chris_pdf_matrixplot(varargin)

figure

for i = 1:nargin
    for j = 1:i
        subplot(nargin,nargin,(i-1)*nargin+1 + (j-1))
        
        if j==i
             yyaxis right
             [A,x] = chris_pdf(varargin{i},(nanmean(smooth(abs(varargin{i}))))/100,nan,nan);
             plot(x,smooth(A,5))
        else
             [P,x,y] = chris_jpdf_moments(varargin{j},varargin{i},100,[nanmin(smooth(varargin{j})) nanmin(smooth(varargin{i}))],[nanmax(smooth(varargin{j})) nanmax(smooth(varargin{i}))]);
             pcolor(x,y,medfilt2(squeeze(P),[3 3]))
             shading interp
             cmocean('thermal-')
             yyaxis right
             set(gca,'ytick',[])             
        end      
        
        if i~=nargin
           set(gca,'xtick',[]) 
        end
        
        if i==j && i~=1
            yyaxis left
            set(gca,'ytick',[])
        end
        
        if j~=1 || (i==j && i==1)
                yyaxis left
                set(gca,'ytick',[]) 
        end
    end
end

end


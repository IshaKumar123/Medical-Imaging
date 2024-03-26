function gfs=glcmfeature(im,angle,NG,NL)

global Idy Idx;

norient=orgoffsets(angle);
nor=size(norient,1); 
gfs=zeros(length(Idy),nor*4);
R = ceil(NG/2);
imp=padarray(im,[R R],'replicate');
Idy2= Idy+R;Idx2=Idx+R;
 for k=1 : length(Idy)
   m=Idy2(k);n=Idx2(k);     
   imf=imp(m-R:m+R,n-R:n+R);     
   gfs(k,:)=texturefea(imf,norient,NL);   
 end






function norient=orgoffsets(angle)

ng=length(angle);
norient=[];    
k=1;

for n=1 : ng
    
 switch angle(n)
       case 0
        norient(k,1) = 0;
        norient(k,2) = 1;
        k=k+1;
       case 45
        norient(k,1) = -1;
        norient(k,2) = 1;   
        k=k+1;
      case 90
        norient(k,1) = -1;
        norient(k,2) = 0;
        k=k+1;
      case 135
        norient(k,1) = -1;
        norient(k,2) = -1;   
        k=k+1;
 end     
   
end



function tfs=texturefea(imf,norient,NL)

glcm= graycomatrix(imf,'NumLevels',NL,'Offset',norient);
stats = graycoprops(glcm, {'contrast','Correlation','Energy','Homogeneity'});
tfea =zeros(size(norient,1),4);
for k=1 : size(norient,1)
tfea(k,1)=stats.Contrast(k);
tfea(k,2)=stats.Correlation(k) ;
tfea(k,3)=stats.Energy(k);
tfea(k,4)=stats.Homogeneity(k);
end

tfs=[tfea(1,:) tfea(2,:)];


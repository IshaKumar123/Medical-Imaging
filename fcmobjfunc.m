function J=fcmobjfunc(V,fm)

global X;

[N,~]=size(X);
nc=size(V,1);

U=zeros(nc,1);
fs = 2/(fm-1);
J=0;

 for n=1 : N               
                                     
     Dis=  sum( (repmat(X(n,:),nc,1)   -V).^2,2 )   ;                                     
     for c=1 : nc       
      U(c) = sum((Dis(c)./Dis).^fs)  ;      
     end              
     U= 1./U;     
     J=J + sum( (U.^fm).*Dis) ;
     
 end 



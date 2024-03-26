function mlbp=mulreslbp(im,R,P)

%global Id;
global Idy Idx;


mlbp=zeros(length(Idy),length(R));

for r=1 : length(R)
mlbp(:,r)=imagelbp(im,R(r),P(r));
S=sum(2.^(P(r)-1:-1:0));
mlbp(:,r) = mlbp(:,r)./S;
end

%[M,N]=size(im);
%mflb=zeros(M,N,length(R));
%for r=1 : length(R)
%mflb(Id,r) = mlbp(:,r);
%end

%Idy2=Idy+R; Idx2=Idx+R;
%for k=1 : length(Idy)    
%   m=Idy2(k);n=Idx2(k);  
%   sim=mflb(m-R:m+R,n-R:n+R,:);  
%end



function lbp=imagelbp(im,R,P)

global Idy Idx;
lbp=zeros(length(Idy),1);
nyx=neighboryx([0 0],R,P);
yx0= nyx+(R+1);

imp=padarray(im,[R R],'replicate');
Idy2= Idy+R;Idx2=Idx+R;
for k=1 : length(Idy)    
   m=Idy2(k);n=Idx2(k);  
   sim=imp(m-R:m+R,n-R:n+R);
   cx=sim(R+1,R+1);
   lbp(k)=lbpattern(sim,yx0,cx,0.01);   
end



function test()

yx=[0 0]; R=1;P=8;
nyx=neighboryx(yx,R,P)
nyx+(R+1)

yx=[0 0]; R=2;P=12;
nyx2=neighboryx(yx,R,P)
nyx2+(R+1)

yx=[0 0]; R=3;P=16;
nyx3=neighboryx(yx,R,P)
nyx3+(R+1)

%nyx3(1,:) = nyx3(1,:) +4;
%nyx3(2,:) = nyx3(2,:) +4;

figure;plot(nyx(2,:),nyx(1,:),'b*-',nyx2(2,:),nyx2(1,:),'r*-',nyx3(2,:),nyx3(1,:),'m*-');
axis([-10 10 -10 10]);


function lbp=lbpattern(sim,yx0,cx,incr)
  
  P=size(yx0,2);  
  bt=zeros(1,P);
  for m=1 : P      
     bt(m)= sim(yx0(1,m),yx0(2,m)) > (cx+incr); 
  end
 lbp=sum( bt.* (2.^(P-1:-1:0)) );     


function nyx=neighboryx(yx,R,P)

 angle = ( 0 : (2*pi)/P : 2*pi - (2*pi)/P ) -  pi/2; 
 xd = R.*cos(angle)+ yx(2) ;
 yd = R.*sin(angle)+ yx(1) ;    
 nyx=[round(yd) ; round(xd)];

function crowoptsegment


clc;clear all; close all;
global X; global Idx Idy;
global Id; global fm;

fname='liver/1.png';
sz=0.5;  % image size

[im,map]=readimage(fname,sz);
[Row,Col]=size(im);
disp(sprintf('Image Size : %d  x %d \n',Row,Col));
Id=im>(5/255);
[Idy,Idx]=find(Id); 


% Features based 
%X=feauresbased(im);
%X=X[1 2 3 4 5 6 7 8]

% gaussian and median filter based 
X=nonfeatures(im);

AP=0.1; % Awareness probability
fl=0.5; % Flight length (fl)
nP=10;nx=size(X,2);nc=3;fm=3.0;T=10;
[gcs,gfs,tbfs]=crowfcm(fm,nP,AP,fl,nc,nx,T);
figure;plot(1:length(tbfs),tbfs,'*-'); title('CROW - Optimization');
xlabel('Iteration');ylabel('Fitness value');

disp('Initial Cluster Centers');
gcs

[IDX,V,objf,Um]=membership(gcs,fm);   %%% Changes 
sim=coloring(IDX,V);
sim2=zeros(Row*Col,3);
sim2(Id,:)= sim;
sim3=reshape(sim2,Row,Col,3);
disp( sprintf('Best Objective values %.2f \n',objf(end)) );



[CS,U,objv]=fcm(X,nc);
[~,IDX]=max(U);
sim=coloring(IDX,CS);
sim2=zeros(Row*Col,3);
sim2(Id,:)= sim;
sim4=reshape(sim2,Row,Col,3);
disp( sprintf('Best Objective values %.2f',objv(end)) );


figure;imshow(im,map); title('Input-Image');
figure;subplot(1,2,1);imshow(sim4,map); title('FCM -Clustering');
subplot(1,2,2);imshow(sim3,map); title('CROW-FCM -Clustering');
figure;plot(1:length(objv),objv,'-*r',1:length(objf),objf,'-ob'); title('Performance : FCM/CROW-FCM');
xlabel('Iteration'); ylabel('Objective value'); 
legend('FCM','CROW-FCM');


%%% Changes  start

Img=im(Id);
cent=V;
cNum=nc;
[Vpc,Vpe,Vfs,Vxb]=fcmclsvalidate(Img,Um,cent,fm);
disp(sprintf('\n CROW-FCM Segmentation Validation'));
disp(sprintf('Num-clusters : %d',cNum));
disp(sprintf( 'Vpc  : %.4f',Vpc) );
disp(sprintf( 'Vpe  : %.4f',Vpe) );
disp(sprintf( 'Vxb  : %.4f',Vxb) );
disp(sprintf( 'Vfs  : %e',Vfs) );

%%% Changes end



function sim=coloring(IDX,gcs)

colr =[1 0 0; 0 1 0; 0 0 1;  
       1 0.5 0.5; 0.5 1 0; 0 0.5 1;  
       1 0.5 0; 0.5 1 0.5; 0.5 0.5 1; 
       1 0.25 0.25; 0.25 1 0.25; 0.25 0.25 1;
       0.25 0.5 0.5; 0.5 0.25 0.5; 0.5 0.5 0.25 ];
   
   
sim =zeros(length(IDX),3);
[nc,nd]=size(gcs);
if (nd>1) 
ps=sum(gcs,2);    
else
ps=gcs;        
end
[~,si]=sort(ps);

for c=1 : nc    
  S=find( IDX==si(c));  
  sim(S,:) = repmat( colr(c,:),length(S),1);       
end



function mfs=nonfeatures(im)

global Id;

[Row,Col]=size(im);

H=fspecial('gaussian',3,1);
img=conv2(im,H,'same');
imc2=reshape(img,Row*Col,1);

%H=ones(3,3)./9;
%iml=conv2(im,H,'same');
%imc3=reshape(iml,Row*Col,1);

imm=medfilt2(im,[3 3]);
imc4=reshape(imm,Row*Col,1);

%mfs =[imc2(Id) imc4(Id)];

mfs =[imc2(Id) ];   %%% changes 

%figure;subplot(1,2,1);imshow(img);
%subplot(1,2,2);imshow(iml);



function mfs=feauresbased(im)

gfs=glcmfeatures(im);
mlbp=lbpfeatures(im);

%Concat features GLCM & LBP
mfs=cat(2,gfs,mlbp);
n=isnan(mfs);
mfs(n)=0;

str1 ={'GLCM DEG :0 Contrast','GLCM DEG :0 Correlation','GLCM DEG :0 Energy','GLCM DEG :0 Homogenity'};
dispfeatures(mfs,[1 2 3 4],Row,Col,str1);
str2 ={'GLCM DEG :45 Contrast','GLCM DEG :45 Correlation','GLCM DEG :45 Energy','GLCM DEG :45 Homogenity'};
dispfeatures(mfs,[5 6 7 8],Row,Col,str2);
str2 ={'MLB R=1','MLB R=2','MLB R=3'};
dispfeatures(mfs,[9 10 11],Row,Col,str2);


save('features.mat', 'mfs','Id','Idx','Idy');



function dispfeatures(mfs,fs,Row,Col,str)
global Id;
figure;
for n=1: length(fs)    
 mf=zeros(Row,Col);
 mf(Id) = mfs(:,fs(n));
 subplot(2,2,n); imshow(mf); title(str{n});
end


function mlbp=lbpfeatures(im)
% LBP features
R=[1,2,3];   % radius 
P=[8,10,12];  % number of points in radius
mlbp=mulreslbp(im,R,P);
disp(sprintf('\nMulti Resolution LBD Features '));
disp(sprintf('\tRadius : %s', sprintf('%d ',R)));
disp(sprintf('\tLBP Features Size : %d %d',size(mlbp)));

function gfs=glcmfeatures(im)

%GLCM features
angle=[0 45]; NG=11;  % neighbor size
NL=12; % number of bin
gfs=glcmfeature(im,angle,NG,NL);
for n=1 : size(gfs,2)    
mx=max(gfs(:,n));
mn=min(gfs(:,n));    
gfs(:,n)=(gfs(:,n)- mn)./(mx-mn);
end
disp(sprintf('GLCM Features : \t Contrast \t Correlation \t Energy \t Homogenity'))
disp(sprintf('\tOrientation Angle : %d %d',angle(1),angle(2)));
disp(sprintf('\tTexture Features Size : %d %d',size(gfs)));


function [IDX,V,objf,U]=membership(V,fm)   %%% changes 

global X;
[N,dim]=size(X);
nc=size(V,1);
U=zeros(N,nc);
fs=2/(fm-1);
pV=V;
objf=[];
for t=1 : 30
    J=0;
 for n=1 : N                    
     Dis=  sum( (repmat(X(n,:),nc,1)   -V).^2,2 )   ;         
     for c=1 : nc       
      U(n,c) = sum((Dis(c)./Dis).^fs)  ;      
     end              
     U(n,:)= 1./U(n,:);   
     J=J + sum( (U(n,:).^fm).*Dis') ;
 end 
 
 
Um= U.^fm;  
for c =1 : nc            
    nm=sum(X.*repmat(Um(:,c),1,dim) ) ;  
    dm=sum(Um(:,c));
    V(c,:) = nm./dm;
end
 
 disp( sprintf( 'Iteration count = %d, obj. fcn = %.4f',round(t),J) );
 objf=[objf J];
 if (  max( abs(V(:)-pV(:)) )<1e-3 )     
     break;
 end
 
 
end
[~,IDX]=max(U,[],2);

 
function [im,map]=readimage(fname,sz)

[im,map]= imread(fname);
if ( size(im,3) ~=1 ) 
im=rgb2gray(im);    
end
im=im2double(im);
if ( sz<1 )
 im=imresize(im,sz);
end

function [gcs,gfs,tbfs]=crowfcm(fm,nP,AP,fl,nc,nvars,T)


range=[0 1];
Xc=population(nP,nc,nvars,range);
tbfs =zeros(T,1);   % history of fitness iteration
fs=fitness(Xc,fm);
[gfs,gbi]=min(fs);  % global best particle index
gcs = Xc{gbi};

Xmem=Xc; % Memory initialization
fs_mem=fs; % Fitness of memory positions



for t=1: T

    num=ceil(nP*rand(1,nP)); % Generation of random candidate crows for following (chasing)
    Xnew=cell(size(Xc));
    
    for i=1:nP
        if rand>AP
          tmp=fl*rand*(Xmem{num(i)}-Xc{i});  
          Xnew{i}= Xc{i}+ tmp ; % Generation of a new position for crow i (state 1)         
        else        
          Xnew{i}= newcrow(nc,nvars,range); % Generation of a new position for crow i (state 2)        
        end             
        Xnew{i}=lulimit(Xnew{i});
    end

    fs=fitness(Xnew,fm);
       
    for i=1:nP % Update position and memory
        cs=Xnew{i};
        if ( all(cs(:)>=range(1)) && all(cs(:)<=range(2)) )        
            Xc{i}=Xnew{i}; % Update position
            if fs(i)<fs_mem(i)
                Xmem{i}=Xnew{i}; % Update memory
                fs_mem(i)=fs(i);
            end
        end
        
    end

    tbfs(t)=min(fs_mem); % Best found value until iteration t    
    
    % Output the results  to screen    
    ngbest=find(fs_mem== min(fs_mem));
    gcs=Xmem{ngbest(1)}; % Solution of the problem
    gfs= fs_mem(ngbest(1));       
    disp(sprintf('Iteration Count :%d Objective value : %.2f',t,gfs));      
        
end



function fs=fitness(Xc,fm)

nP=length(Xc);
fs=zeros(1,nP);

for i=1:nP           
   cs= Xc{i};
   if ( all(cs(:)>=0) && all(cs(:)<=1) )
      fs(i)=fcmobjfunc(Xc{i},fm);
   else
      fs(i)=1e6;
   end
end



function ps=lulimit(ps)

I1 =(ps<0);   
I2= (ps>1);   
ps(I1) = 0;
ps(I2) = 1;


function  Xc=population(np,nc,nd,const)

Xc=cell(np,1); 
for p=1 : np     
 Xc{p}=newcrow(nc,nd,const);
end


function cs=newcrow(nc,nd,const)

cs=zeros(nc,nd);     
for c=1 :  nc     
  cs(c,:) = const(1)  + rand(1,nd).*(const(2)- const(1));    
end 


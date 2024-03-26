function [Vpc,Vpe,Vfs,Vxb]=fcmclsvalidate(ab1,U,center,m)

    %%%%%%%%%%%%Vpc%%%%%%%%%%%%%
  vpc=sum(sum(U.*U));
  div=size(ab1,1);
  Vpc=vpc/div;
  
  
  %%%%%%%%%%%%%Vpe%%%%%%%%%%%%
  vpe=-(sum(sum(U.* log(U))));
  Vpe=vpe/div;
  
  
%%%%%%%%%%%%Vxb%%%%%%%%%%%%%%%%  
  center1=center;
for j=1: size(center,1)
  for i=1: size(center,1)
      if(center(j)~=center1(i))
      cen_f(i)=((center1(i)-center(j))^2) ;  
      end
  end
end



for j=1: size(center,1)
   dis(:,j)=(ab1-center(j)).^2; 
end


cen=min(cen_f);
vxb=(sum(sum((U.^m).*dis)));
Vxb=-vxb/(cen*div);

%%%%%%%%Vfs%%%%%%
V_bar=(sum(center))/size(center,1);
for j=1: size(center,1)
      dis(:,j)= dis(:,j) - (center(j)-V_bar)^2;
end


Vfs=(sum(sum((U.^m).*dis)));

%This function is used to perform ANM calculation, return B-factor(from 7th modes to all modes)
%input: median structure
%output: Bfactor

function [b]=anm(median_str)
 anm_ca=median_str;
 [r,c]=size(anm_ca);
 
 i=1;
 m=3*r;
 aa=zeros(m,m);
 cutoff=15;
 gama=1.0;
 eigen=0.00001;
 
 while(i<=r)
    j=1;
   while(j<=r)
     dx=anm_ca(i,1)-anm_ca(j,1);
     dy=anm_ca(i,2)-anm_ca(j,2);
     dz=anm_ca(i,3)-anm_ca(j,3);
     distance=dx*dx+dy*dy+dz*dz;

     if i~=j&&distance<=cutoff*cutoff
	aa(3*i-2,3*i-2)=aa(3*i-2,3*i-2)+gama*dx*dx/distance;
	aa(3*i-1,3*i-1)=aa(3*i-1,3*i-1)+gama*dy*dy/distance;
	aa(3*i,3*i)=aa(3*i,3*i) +gama*dz*dz/distance;
	aa(3*i-2,3*i-1)=aa(3*i-2,3*i-1)+gama*dx*dy/distance;
	aa(3*i-2,3*i)=aa(3*i-2,3*i) +gama*dx*dz/distance;
	aa(3*i-1,3*i-2)=aa(3*i-1,3*i-2)+gama*dy*dx/distance;
	aa(3*i-1,3*i)=aa(3*i-1,3*i) +gama*dy*dz/distance;
	aa(3*i,3*i-2)=aa(3*i,3*i-2) +gama*dx*dz/distance;
	aa(3*i,3*i-1)=aa(3*i,3*i-1) +gama*dy*dz/distance;

	aa(3*i-2,3*j-2)= -gama*dx*dx/distance;
	aa(3*i-1,3*j-1)= -gama*dy*dy/distance;
	aa(3*i,3*j)= -gama*dz*dz/distance;
	aa(3*i-2,3*j-1)= -gama*dx*dy/distance;
	aa(3*i-2,3*j)= -gama*dx*dz/distance;
	aa(3*i-1,3*j-2)= -gama*dy*dx/distance;
	aa(3*i-1,3*j)= -gama*dy*dz/distance;
	aa(3*i,3*j-2)= -gama*dx*dz/distance;
	aa(3*i,3*j-1)= -gama*dy*dz/distance;
     end 
       j=j+1;
    end
    i=i+1;
 end
        
       [V_eigs,D_eigs]=eig(aa);       
       EigenVector=V_eigs;
       EigenValue=diag(D_eigs);

       [anm_m,anm_n]=size(EigenVector);
       
%calculate B-factor from the slowest mode
 for i=1:r
      TBx(i)=0;
      TBy(i)=0;
      TBz(i)=0;
      for j=1:anm_n
        if EigenValue(j)>eigen
      TBx(i)=TBx(i)+EigenVector(3*i-2,j)^2/EigenValue(j);
      TBy(i)=TBy(i)+EigenVector(3*i-1,j)^2/EigenValue(j);
      TBz(i)=TBz(i)+EigenVector(3*i,j)^2/EigenValue(j);
        end
      end
      TB(i)=TBx(i)+TBy(i)+TBz(i);
 end
 b=TB;

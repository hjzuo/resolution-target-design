% clear all
% n_aSi=3.5;
% n_air=1;
% lambda=4;%um
% f=100;%um
% k0=2*pi/lambda;
% height=2;%height of silicon pillar 
% %% boundry of zero order
% %basic equation: phi(0)+f*k0=phi(r)+sqrt(f^2+r^2)*k0
% % p=phase_delay_height;
% % 
% % p=p-p(1);
% N_zone=1:5;R_zone=zeros(length(N_zone),1);
% for i=1:length(N_zone)
% R_zone(i)=sqrt(((2*pi*N_zone(i))/k0+f)^2-f^2);
% end
% 
% 
% 
% 
% 
% phi=[];
% rPillar=[];
% phi_whole_len=[];
% cutoff=lambda/n_aSi;%struture cutoff
% pitch=2;%periodity
% 
% %% hexogonal lattice vector
% 
% vector_a=[pitch,0];
% vector_b=[pitch*cos(2*pi/6),pitch*sin(2*pi/6)];
% 
% %% phi and radius
% disp('calculating grid position......');
% len_size=max(R_zone);
% k=1;
% for i=0:1:2*round(len_size/pitch)
%     for j=0:1:2*round(len_size/pitch)
% %         temp=vector_a*i+vector_b*j;
% %         r=sqrt(temp(1)^2+temp(2)^2);
% %         if r<len_size
% %             
%         temp=vector_a*i+vector_b*j;
%         Grid_x(k)=temp(1);
%         Grid_y(k)=temp(2);
%         k=k+1;
%         temp=-vector_a*i+vector_b*j;
%         Grid_x(k)=temp(1);
%         Grid_y(k)=temp(2);
%         k=k+1;
%         temp=vector_a*i-vector_b*j;
%         Grid_x(k)=temp(1);
%         Grid_y(k)=temp(2);
%         k=k+1;
%         temp=-vector_a*i-vector_b*j;
%         Grid_x(k)=temp(1);
%         Grid_y(k)=temp(2);
%         k=k+1;
% %         end
%         
%     end
% end
% 
% %% phi and radius
% disp('calculating phi and radius......');
% k=1;
% len_size=max(R_zone);
% for i=1:length(Grid_x)
%     r=sqrt(Grid_x(i)^2+Grid_y(i)^2);
%     if (r<=len_size)&&(i>1)&&((Grid_x(i)~=Grid_x(i-1))||(Grid_y(i)~=Grid_y(i-1)))
%         EBL_x(k)=Grid_x(i);
%         EBL_y(k)=Grid_y(i);
%         phi(k)=(f-sqrt(f^2+r^2))*k0+2*pi;       %center phi=2*pi  phi(0)+f*k0=phi(r)+sqrt(f^2+r^2)*k0
%         phi(k)=mod(phi(k),2*pi)/pi;
%         rPillar(k)=phi2radius(phi(k));
%         
%         k=k+1;
%     end
% end
% phi(k)=2;rPillar(k)=phi2radius(phi(k));EBL_x(k)=0;EBL_y(k)=0;  % center pillar
% % len_size=max(R_zone)+pitch;
% % for i=1:1:round(len_size/pitch)
% %     for j=1:1:round(len_size/pitch)
% %         r=sqrt(i^2+j^2)*pitch;
% %         if r<len_size
% %             
% %             phi(i,j)=(f-sqrt(f^2+r^2))*k0+2*pi;       %center phi=2*pi  phi(0)+f*k0=phi(r)+sqrt(f^2+r^2)*k0
% %             temp_phi=mod(phi(i,j),2*pi)/(pi);
% %             rPillar(i,j)=phi2radius(temp_phi);
% %         else
% %             phi(i,j)=0;
% %             rPillar(i,j)=0;
% %         end
% %     end
% % end
% 
% % wrap_phi=mod(phi,2*pi);
% % x_grid=linspace(0,len_size,i);
% % y_grid=linspace(0,len_size,i);
% % imagesc(x_grid,y_grid,wrap_phi)
% % figure
% % imagesc(x_grid,y_grid,rPillar)
% disp('generating GDS');
% %% generate GDS file
% name=['RCWA_metalen_' 'pitch' num2str(pitch) 'um_' 'focal' num2str(f) 'um' 'span' num2str(round(len_size*2)) 'um' ];
% l=length(rPillar);
% k=1;
% for i=1:l
% %     for j=1:l
%         if rPillar(i)~=0
%         E(k)=Raith_element('circle',0,[EBL_x(i) EBL_y(i)],rPillar(i),[],60,1);k=k+1;
% %         E(k)=Raith_element('circle',0,[(-i+0.5)*pitch (j-0.5)*pitch],rPillar(i,j),[],60,1);k=k+1;
% %         E(k)=Raith_element('circle',0,[(-i+0.5)*pitch (-j+0.5)*pitch],rPillar(i,j),[],60,1);k=k+1;
% %         E(k)=Raith_element('circle',0,[(i-0.5)*pitch (-j+0.5)*pitch],rPillar(i,j),[],60,1);k=k+1;
%         end
% axis equal;
% % end
% end
% E(k)=Raith_element('text',0,[0 -80],5,0,[1 1],name,1.5);
% S=Raith_structure(name,E);
% clf;
% axis equal;
% S.plot;
% span=Raith_library(name,S);
% span.writegds;
name='resolution target';
k=1;
for k=1:5
E(k)=Raith_element('path',0,[0 0 1 1 2;1 0 0 1 1],0.5,1.3);
end
S=Raith_structure(name,E);
clf;
axis equal;
S.plot;
span=Raith_library(name,S);
span.writegds;
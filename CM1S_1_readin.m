% #########################################################################
% Rochelle Worsnop
%Program:CM1S_1_readin.m Based on CM1_readin_1.m & calc_hurricane_nonavg_1.m
%Purpose: reads in variables from netcdf files
%Last updated: 10/06/2014,09/09/2014, 01/06/2015

% Code used for the "simple hurricane-like profile" case
%dx=dy=62.5 m dz= 31.25 m dt = 0.375s
% This data set is used to compare with the CM1 data of the same temporal
% and spatial resolution to determine if the simplistic model can be used
% to represent the winds in a hurricane.
% #########################################################################

clear % Clear the workspace
close all % Closes all plots (figures)
clc % Clears the command line

dir='E:\Matlab programs_v3 data\!!!!' ; %Change to match directory of data from each simulation
prefix='simple_gp';
radius='_r0';
azimuth='_a0';
coord='_s' ;
suffix='.nc' ;

n=0 ; %used later to start the loop

% Change below to match what the new simulations have!!!
radii_num=2 ; %I actually have 7 radii, but the other two are far away from the center.
%radius 1 =4 radius 2 = 6
altitude_num =33 ;
time_num=19201 ; %28801 ; % 2 hours at 0.25 second time step. before 19200 for 0.375s ;
location_num =1;
grid=7 ; %I have a 7x7 grid.

p0=100000.0 ; %initial pressure, sea level? Pa
rair=287.04 ; %Universal gas constant for dry air J/kgK
cp=1005. ; %specific heat of dry air J/kg K
cnst=rair/cp ;

disp('before loop')
% for i=1 ; %loop over number of radii
%for j=1:location_num ; %loop over number of locations
for i=1:radii_num ; %loop over number of radii
    
    if i ==1
        i_it=4 ; %4th radius. I was only given two radii to work with. 
    end;
    if i==2
        i_it=6 ; %6th radius
    end
    % i=1 ; %only one radii in this data set
    j=location_num ; %only one location used for the grid data
    location=strcat(radius,num2str(i_it),azimuth,num2str(j)) ; % radii_num used to just be '1'
    file=strcat(dir,prefix,location,suffix) ;
    fileid=netcdf.open(file,'nowrite') ;
    
    % Identify variables in the file
    x = read_ncdfvar(fileid,'xlocation') ;
    y = read_ncdfvar(fileid,'ylocation') ;
    u = read_ncdfvar(fileid,'u') ;
    v = read_ncdfvar(fileid,'v') ;
    w = read_ncdfvar(fileid,'w') ;
    %     qv = read_ncdfvar(fileid,'qv') ;
    ust = read_ncdfvar(fileid,'ust') ;
    %     u10 = read_ncdfvar(fileid,'u10') ;
    %     v10 = read_ncdfvar(fileid,'v10') ;
    %     prs = read_ncdfvar(fileid,'prs') ;
    %     time = read_ncdfvar(fileid,'time') ;
    %     t= read_ncdfvar(fileid,'t') ; %temperature
    disp('read var')
    
    % Shift dimensions of the variables to be largest to smallest
    %(time, height,nj,ni)
    u=permute(u,[4,3,2,1]);
    v=permute(v,[4,3,2,1]);
    w=permute(w,[4,3,2,1]);
    %     t=permute(t,[4,3,2,1]);
    %     qv=permute(qv,[4,3,2,1]);
    %     prs=permute(prs,[4,3,2,1]);
    
    ust=permute(ust,[3,2,1]); %(time,nj,ni)
    %     u10=permute(u10,[3,2,1]);
    %     v10=permute(v10,[3,2,1]);
    disp('done permute')
    %Initialize the variables with the dimensions of the 7 radii
    %and 1 location, 7x7 grid. Fill the arrays with zeros first.
    %suffix 'a' means 'all'
    if n==0 ; %only do this one time. 
        s=size(u) ;
        nh=s(2) ;
        nt=s(1) ;
        gridj=s(3) ; %south-north location of scalar grid points
        gridi=s(4) ; %west-eat location of scalar grid points
        %          xa = zeros(gridj,gridi,radii_num) ;
        %          ya = zeros(gridj,gridi) ; %ya is the same for every radii, so only need two dimensions
        ua = zeros(nt,nh,gridj,gridi,radii_num); %time,height,j,i,location, radii
        va = zeros(nt,nh,gridj,gridi,radii_num);
        wa = zeros(nt,nh,gridj,gridi,radii_num);
        %         ta = zeros(nt,nh,gridj,gridi,radii_num);
        %         qva = zeros(nt,nh,gridj,gridi,radii_num);
        %         prsa = zeros(nt,nh,gridj,gridi,radii_num);
        usta = zeros(nt,gridj,gridi,radii_num);
        %         u10a = zeros(nt,gridj,gridi,radii_num);
        %         v10a = zeros(nt,gridj,gridi,radii_num);
    end
    
    disp('done creating arrays')
    % What do I do with nj and ni? Do I make it another dimension?
    % for each time, height,position north and east, there is a location and
    % radius for each point.
    ind=1 ;
    %      xa(:,:,i)       =round(x)           ; %i=radius, j=location I took out j, because there is only one location with this data set
    %      ya(:,:)         =round(y)           ; %ya is the same for all radii
    ua(:,:,:,:,i)   = u(:,:,:,:)     ;%time,alt,north,east,location,radius
    va(:,:,:,:,i)   = v(:,:,:,:)     ;
    wa(:,:,:,:,i)   = w(:,:,:,:)     ;
    %     ta(:,:,:,:,i)   = t(:,:,:,:)     ;
    %     qva(:,:,:,:,i)  = qv(:,:,:,:)    ;
    %     prsa(:,:,:,:,i) = prs(:,:,:,:)   ;
    usta(:,:,:,i)   = ust(:,:,:)     ;
    %     u10a(:,:,:,i)   = u10(:,:,:)     ;
    %     v10a(:,:,:,i)   = v10(:,:,:)     ;
    
    n=n+1;
end %end loop over radii
disp('done filling array')

% Add on specific attachments for each simulaiton, such as "V45_R5_dx31_dt0_375"
save('E:\Matlab code for CM1 simulations output\save_files\NewHurricaneLES_grid_nonavg_wind_v3.mat!!!!',...
    'ua','va','wa','usta','-v7.3');
            % 'xa','ya','ua','va','wa','ta','qva','prsa','usta','u10a','v10a','-v7.3');





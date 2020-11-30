
%% This file is called by 'SetUp_Initiate' to initiate Arp2/3 actin nodes and nucleation points

disp('initializing Npfs')

%%-----Following lines of code will create a circular random area
%%containing nucleation points

for i=1:length(cellCenters)
% %% for circular
cellCenTemp = cellCenters(i,:);
cellCenter_x= cellCenTemp(1,2); %% col
cellCenter_y= - cellCenTemp(1,1); %% row
for x =cellCenter_x-CapR1 : cellCenter_x+CapR1
    y1= cellCenter_y - sqrt(CapR1^2- (x-cellCenter_x)^2);
    y2= cellCenter_y + sqrt((CapR1)^2- (x-cellCenter_x)^2);
    for y = y1:y2
        row =round(-y);
        col = round(x);
        rno=rand(1);
        if rno< Den_c1             % Density of inner Actin Nodes in cap
            npfChart(row,col)= Vnpf;    % Initializing Nucleation points
        end
    end
end
end

for i=1:length(cellCenters)
% %% for circular
cellCenTemp = cellCenters(i,:);
cellCenter_x= cellCenTemp(1,2); %% col
cellCenter_y= - cellCenTemp(1,1); %% row
for x =cellCenter_x-CapR : cellCenter_x+CapR
    y1= cellCenter_y - sqrt(CapR^2- (x-cellCenter_x)^2);
    y2= cellCenter_y + sqrt((CapR)^2- (x-cellCenter_x)^2);
    for y = y1:y2
        row =round(-y);
        col = round(x);
        rno=rand(1);
        if rno< Den_c             % Density of outer Actin Nodes in cap
            npfChart(row,col)= Vnpf;    % Initializing Nucleation points
        end
    end
end
end

% Count and initialize Nucleation Points
% %
disp('Making NFPIds')
NPFCount=1;

for colNo=1:bound
    for rowNo=1:bound
        if npfChart(rowNo,colNo)==Vnpf
            NPF(NPFCount,1)=rowNo;%NPF stores nucleation points location data
            NPF(NPFCount,2)=colNo;
            NPFCount=NPFCount+1;
            
        end
    end
end
NPFCount=NPFCount-1;
%%
imagesc(npfChart)
pause(1)

%%---Following lines of code will initialize actin nodes using the
%%nucleation points 

Ac_NodeCount=1;
for i=1:NPFCount
    row= NPF(i,1);
    col= NPF(i,2);
    dir=randi(8);
    theta = directions(dir,3);
    dir_c=  cos(deg2rad(theta));
    dir_r=  -sin(deg2rad(theta));
    
    len=randi(8);
    Fnpf_r= Knpf*dir_r*len;%polymerization force
    Fnpf_c= Knpf*dir_c*len;
    
    rowNo=row+Fnpf_r; % actin position
    colNo=col+Fnpf_c;  % actin position
    
    Ac_Node(Ac_NodeCount,ROW)=rowNo; %Update Actin array with rowNo
    Ac_Node(Ac_NodeCount,COL)=colNo; 
    Ac_Node(Ac_NodeCount,BRANCHORNOT)= 0; % A zero at this position represents unbranched
    Ac_Node(Ac_NodeCount,THETA)=theta ;%direction of actin filamnet
    Ac_Node(Ac_NodeCount,LEN)= sqrt((rowNo-row)^2+(colNo-col)^2); % will store length of filament
    Ac_Node(Ac_NodeCount,ASSCNPF)= i;  % associated NPF
    Ac_Node(Ac_NodeCount,MCSATBIRTH)= 0;  % Simulation step of creation like a date of birth
    Ac_Node(Ac_NodeCount,AGE)= len;  % age of actin filament is made different initially for variability
    
  Ac_NodeNo(round(rowNo),round(colNo))= Ac_NodeCount;  % this stores Actin node ID in an array
   Combo(round(rowNo),round(colNo))= Ac;  % Update Combo(visualization lattice) with Actin
    
    Ac_NodeCount=Ac_NodeCount+1;
end
Ac_NodeCount=Ac_NodeCount-1;
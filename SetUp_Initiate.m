%% This file is used for initialization of the two materials i.e., defining initial position of actomyosin and actin nodes.
%% After defining the node properties, it connects the nodes with neighors to create a connected network.
%% This is the first file to be run followed by 'Simulate_CombinedModel'

clc
% clear all
close all

% Following module is called to define different parameters that are used
% in this simulation

ParameterDefination

%% --------Actomyosin material Initialization-------------

disp('initializing Nodes')


cellCenTemp = cellCenters(1,:);
cellCenter_x = 500;% cellCenTemp(1,2); %% col
cellCenter_y =-500;% - cellCenTemp(1,1); %% row


% fill whole area with passive myosin
for x =cellCenter_x-(bound/2-6) : cellCenter_x+ (bound/2-6)  %% 6 px are deducted to prevent lattice boundary colision
    for y = cellCenter_y-(bound/2-6) : cellCenter_y+ (bound/2-6)
        row =round(-y);
        col = round(x);
        rno=rand(1);
        if rno< Den_m             % Density of Myosin Nodes
            Combo(row,col)= Dm;    %  Myosin in combo
        end
    end
end

% In the desired area convert passive myosin to active
for x =cellCenter_x-MyoR : cellCenter_x+ MyoR
    y1= cellCenter_y - sqrt((MyoR)^2- (x-cellCenter_x)^2);
    y2= cellCenter_y + sqrt((MyoR)^2- (x-cellCenter_x)^2);
    for y = y1:y2
        row =round(-y);
        col = round(x);
        % Density of Myosin Nodes
        if Combo(row,col)== Dm
            Combo(row,col)= Mc;
        end %  Myosin in combo
    end
end




%% ---- Nucleation points and actin Initialization by calling 'InitiateActin' module
InitiateActin

%% Clear myosin from actin area

for i=1:length(cellCenters)
    cellCenTemp = cellCenters(i,:);
    cellCenter_x = cellCenTemp(1,2); %% col
    cellCenter_y = - cellCenTemp(1,1); %% row
    
    for x =cellCenter_x-CapR : cellCenter_x+ CapR
        y1= cellCenter_y - sqrt((CapR)^2- (x-cellCenter_x)^2);
        y2= cellCenter_y + sqrt((CapR)^2- (x-cellCenter_x)^2);
        for y = y1:y2
            row =round(-y);
            col = round(x);
       
            if  Combo(round(row),round(col))== Ac
                
                for r= row -5 : row +5
                    for c= col-5: col+5
                        Combo(r,c)= 0;    %  Myosin in combo
                        Combo(row,col)= Ac;    %  Actin in combo
                    end
                end
            end
        end
    end
end

%% ---Defining Myosin nodes properties ----- 

% % making Node IDs

% Mc_Node is an array that holds the information and flags for each node of
% myosin, like its active or passive, its position etc.

Mc_NodeCount=1;
NodeCountPerPixel=zeros(bound,bound);
disp('Making NodeIds')
for colNo=1:bound
    for rowNo=1:bound
        if Combo(rowNo,colNo)>0
            
            NodeCountPerPixel(rowNo, colNo) = NodeCountPerPixel(rowNo, colNo) + 1 ; % update the #nodes on the picked pixel
            Mc_NodeNo(rowNo,colNo)=Mc_NodeCount; % Mc_NodeNo will hold the node ID of node at that position. right now Count=ID
            
            Mc_Node(Mc_NodeCount,ROW)=rowNo; %Update Myosin array with rowNo
            Mc_Node(Mc_NodeCount,COL)=colNo;  %Update Myosin array with colNo
            Mc_Node(Mc_NodeCount,NO_OF_NEIGHBOR)=0; %Update Myosin array with NoOfNeighbour
            Mc_Node(Mc_NodeCount,COMBOVALUE)=Combo(rowNo,colNo);  %% Helpful when testing active passive nodes
            Mc_Node(Mc_NodeCount,MEAN_LENOFCONNEC)=0;  %% helpful for calculation of tn_rest dynamically for each simulation
            if Combo(rowNo, colNo)==Mc
                Mc_Node(Mc_NodeCount,MYO_CONC)=Myo_conc;%
            else
                Mc_Node(Mc_NodeCount,MYO_CONC)=0;
            end
            
            if rowNo> (bound-15) | rowNo<15 | colNo>(bound-15) | colNo<15
                Mc_Node(Mc_NodeCount,MOVEORNOT)=1;%Moveable NOde
            else
                Mc_Node(Mc_NodeCount,MOVEORNOT)=0;
            end
        
            Mc_NodeCount=Mc_NodeCount+1;
        end
    end
end
Mc_NodeCount=Mc_NodeCount-1;

%%%%%%---------Initialization Ends---------%%%%%

disp('Making Node Neighbors')
%% ------Making Node neighbors-------

neighborhood = zeros(Mc_NodeCount,MaxNeighAllowed);
%
AvailNeighRegistry=zeros(Mc_NodeCount,1);
RandNodes=randperm(Mc_NodeCount);

for i = 1:Mc_NodeCount
    NodeID= RandNodes(i) ; %randomly slect a myosin node
    AN=0;
    
    if(Mc_Node(NodeID,4)<MaxNeighAllowed) % check if the #neighbors is less than allowed
        
        rowNo = Mc_Node(NodeID,ROW);
        colNo = Mc_Node(NodeID,COL);
        
        
        % % search a circular area for neighbors
        for x= round(rowNo)-Dis_Thres_Neigh1:round(rowNo)+Dis_Thres_Neigh1
            y1= colNo - sqrt(Dis_Thres_Neigh1^2- (x-round(rowNo))^2);
            y2= colNo + sqrt(Dis_Thres_Neigh1^2- (x-round(rowNo))^2);
            %                  [NodeID x y1 y2]
            for y= y1:y2
                if Mc_NodeNo(round(x),round(y))>0 && Mc_NodeNo(round(x),round(y))~= NodeID %% if you found a neighbor
                    NeighborID = Mc_NodeNo(round(x),round(y));
                    
                                      
                    % % check if the neighbor you found is already
                    % listed as neighbor to avoid redundancy
                    
                    alreadyANeighbor = 0;
                    for k= 1:Mc_Node(NodeID,NO_OF_NEIGHBOR)
                        if NeighborID == neighborhood(NodeID,k)
                            alreadyANeighbor=1;
                            break;
                        end
                    end
                    
                    % % add the neighbor in the available neighor
                    % list to be used for selecting neighbors
                    if ((alreadyANeighbor==0) && Mc_Node(NeighborID,4)<(MaxNeighAllowed))
                        AN=AN+1;
                        AvailableNeighbors(NodeID, AN)=NeighborID;
                      
                    end

                end
            end
        end
        
        
        % % check here if Available Neighbors array is formed
        
        neighborAllowed = MaxNeighAllowed - Mc_Node(NodeID,4); %  get to know how many neibhors can node form 
        
        if AN>0
            
            if AN <= neighborAllowed  
                
                for i= 1:AN
                    
                    NeighborID=AvailableNeighbors(NodeID,i);
                    
                    %make connection
                    neighborhood(NodeID, Mc_Node(NodeID,NO_OF_NEIGHBOR)+1) = NeighborID; %add neighborId to neighborhood array of selected node
                    neighborhood(NeighborID, Mc_Node(NeighborID,NO_OF_NEIGHBOR)+1) = NodeID; % add selected nodeId to neighborhood array of neighbor
                    
                    %update neighborcount of both neighbor and node
                    Mc_Node(NeighborID,NO_OF_NEIGHBOR) = Mc_Node(NeighborID,NO_OF_NEIGHBOR)+1; 
                    Mc_Node(NodeID,NO_OF_NEIGHBOR) = Mc_Node(NodeID,NO_OF_NEIGHBOR)+1;
                    
                end
                
            else % this is the case when there are more available neighbors than required and so you have to select which neighbors you will make connection with
                R = randperm(AN,neighborAllowed);
                
                % repeat making connections 
                for i= 1:neighborAllowed
                    
                    NeighborID=AvailableNeighbors(NodeID,R(i));
                    
                    %make connection
                    neighborhood(NodeID, Mc_Node(NodeID,NO_OF_NEIGHBOR)+1) = NeighborID;
                    neighborhood(NeighborID, Mc_Node(NeighborID,NO_OF_NEIGHBOR)+1) = NodeID;
                    %update neighborcount
                    Mc_Node(NeighborID,NO_OF_NEIGHBOR) = Mc_Node(NeighborID,NO_OF_NEIGHBOR)+1;
                    Mc_Node(NodeID,NO_OF_NEIGHBOR) = Mc_Node(NodeID,NO_OF_NEIGHBOR)+1;
                    
                end
            end
        end
    end
    
     
end
%

%% Following lines calculate the mean resting length of connections
count  = 0
csvdata=zeros(Mc_NodeCount*6,7);

for NodeID = 1: Mc_NodeCount
    Connection_length=0;
    tn_Nb=0;
    for i = 1: Mc_Node(NodeID,NO_OF_NEIGHBOR)  % for all neighbors
        
        NbID = neighborhood(NodeID,i);  % get Neighbor iD
        
        tn_c_Nb= Mc_Node(NbID,COL)- Mc_Node(NodeID,COL);  % OUTWARD VECTOR since pulling from neighbor will be directed out
        
        tn_r_Nb= Mc_Node(NbID,ROW)- Mc_Node(NodeID,ROW);
        
        tn_Nb = 0.001+sqrt((tn_r_Nb)^2 + (tn_c_Nb)^2);
        
        Connection_length= Connection_length+tn_Nb;
        count  = count +1 
        csvdata(count,:)=[NodeID,NbID,Mc_Node(NodeID,ROW),Mc_Node(NodeID,COL),Mc_Node(NbID,ROW),Mc_Node(NbID,COL),(tn_Nb-3.5)];   
    
    end
    
    if Mc_Node(NodeID,NO_OF_NEIGHBOR)>0
        Mc_Node(NodeID,MEAN_LENOFCONNEC)=Connection_length/Mc_Node(NodeID,NO_OF_NEIGHBOR) ;
    else
        Mc_Node(NodeID,MEAN_LENOFCONNEC)=0;
    end
end
tn_rest= mean(Mc_Node(:,MEAN_LENOFCONNEC));


cd (absoluteFolderPath)
csvwrite('csvdata0.csv',csvdata(1:count,:));
cd ..

%
% %% deleting actin nodes if simulations are to be done without actin cap
% % 
% le= length(Ac_Node);
% for NodeID= le:-1:1
%         Ac_Node(NodeID,:)=[];
%         Ac_NodeCount=Ac_NodeCount-1;
% end
% Ac_NodeNo=zeros(bound,bound);

%% print an image of initial setup 

cd (absoluteFolderPath)
filename=strcat('initialImage');
fig = figure;
imagesc(Combo)
% imshow(Cap)
print(fig,filename,'-dpng');
% filename=strcat('initialMyo');
% fig = figure;
%
% imshow(Myo)
% print(fig,filename,'-dpng');
cd ..
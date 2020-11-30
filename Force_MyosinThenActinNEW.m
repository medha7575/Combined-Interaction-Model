%% This file calculates forces first on myosin nodes and then on actin nodes.
%% Myosin node experience forces from neighboring nodes and due to actin pushing.
%% Actin node experiences forces from actin polymerization directed outward from nucleation point to enable actin growth and from myosin pushing.



%% ---- This section calculates myosin forces -----
count  = 0
csvdata=zeros(Mc_NodeCount*6,7);

RandNodes=randperm(Mc_NodeCount);
for i = 1: Mc_NodeCount
    NodeID= RandNodes(i); %randomly pick up one myosin node
    %   %   Calculating forces
    if Mc_Node(NodeID,MOVEORNOT)==1 % check if it is a unmovable node
        %             disp('found a fixed node, so net force is zero')
        F_MyoNet(NodeID,ROW)= 0;  % set all forces on this node to be zero
        F_MyoNet(NodeID,COL)= 0;
        % disp('NodeID=')NodeID ;
    else
        % Force calculation of myosin  starts here sicne these are movable nodes
        
        for i = 1: Mc_Node(NodeID,NO_OF_NEIGHBOR)  % for all neighbors
            
            NbID = neighborhood(NodeID,i);  % pick one neighbor and get Neighbor iD
            if NbID ~=0    %Neigbor ID can not be zero,
                tn_c_Nb= Mc_Node(NbID,COL)- Mc_Node(NodeID,COL);  % OUTWARD VECTOR since pulling from neighbor will be directed out
                
                tn_r_Nb= Mc_Node(NbID,ROW)- Mc_Node(NodeID,ROW);
                
                tn_Nb = 0.001+sqrt((tn_r_Nb)^2 + (tn_c_Nb)^2);   % distance between node and neighbor (adding .001 to avoid numerical instaibility)
                
                 if mod(mcs,500)==0
                    count  = count +1 
                    csvdata(count,:)=[NodeID,NbID,Mc_Node(NodeID,ROW),Mc_Node(NodeID,COL),Mc_Node(NbID,ROW),Mc_Node(NbID,COL),(tn_Nb-3.5)];
                   
                end
                
                %% If there is tension in the spring
                %   ----Implementing Equation-2 ---------
                F_NbSp=0;
                if tn_Nb>tn_rest
                    
                    % % Spring force calculation
                    
                    
                    F_NbSp = KSp*(tn_Nb- tn_rest); % force generated due to spring = contstant times stretch
                    
                    F_r_NbSp(NodeID,i)= F_NbSp*(tn_r_Nb/tn_Nb); % outward force felt by the node due to tension
                    F_c_NbSp(NodeID,i)= F_NbSp*(tn_c_Nb/tn_Nb); % This is equal to magnitude*unit vector
                    
                end
                
                
                F_Spring(NodeID,ROW)=  F_Spring(NodeID,ROW) + F_r_NbSp(NodeID,i); %adding forces by all the neighbors
                F_Spring(NodeID,COL)=  F_Spring(NodeID,COL)+ F_c_NbSp(NodeID,i); % Net Spring force
                
                % % Inter-myosin force calculation
                %   ----Implementing Equation-1 ---------
                % % calculate pulling forces from  myosin neighbors
                F_Nb=0;
                AM_Neighbor= Mc_Node(NbID,MYO_CONC);% %fetch the myosin activity  on node and neighbor
                AM_BaseNode= Mc_Node(NodeID,MYO_CONC);%
                
                F_Nb =(K*AM_Neighbor*AM_BaseNode); ; % force magnitude
                
                F_r_Nb(NodeID,i)=F_Nb*(tn_r_Nb/tn_Nb);  % Actomyosin force on the node due to ith nieghbor
                F_c_Nb(NodeID,i)=F_Nb*(tn_c_Nb/tn_Nb);  % magnitude * unit vector
                
                
                F_am(NodeID,ROW)= F_am(NodeID,ROW)+ F_r_Nb(NodeID,i); % %adding forces by all the neighbors
                F_am(NodeID,COL)= F_am(NodeID,COL)+ F_c_Nb(NodeID,i); %  Net inter-myosin pulling force
                
            end
            
        end
        
        
        %%-----Following lines of code will calculate actin pushing forces on myosin due to the actin polymerization----%%%%%%%
         %   ----Implementing Equation-7 ---------
        totalPoints=0;
        points=0;
        
        colNo= Mc_Node(NodeID,COL); % the position of current node
        rowNo= Mc_Node(NodeID,ROW);
        for R= rowNo-SR : rowNo+SR   %%% searching around the Myosin node in search radius
            for C= colNo-SR:colNo+SR
                if Ac_NodeNo(round(R),round(C))> 0 % if there is an actin node in vicinity (search radius)
                    totalPoints= totalPoints+1;
                    points(totalPoints,1) = R; %% Save the position where actin node is present
                    points(totalPoints,2) = C;
                end
            end
        end
        
        %---Following lines will calculate forces from actin nodes found in the
        %---above determined positions
        
        Fcap_R=0;
        Fcap_C=0;
        
        for i= 1: totalPoints
            R= points(i,1) ;
            C= points(i,2);
            
            Act2Myo_r= rowNo - R;
            Act2Myo_c= colNo - C;  % vector directed towards the myosin node
            tn_Act2Myo= sqrt((Act2Myo_r)^2+(Act2Myo_c)^2);
            
            % accessing the features of actin - its ID, nucleation point and
            % direction
            ActID=  Ac_NodeNo(round(R),round(C));
            npfID=  Ac_Node(ActID,ASSCNPF);
            theta = Ac_Node(ActID,THETA);
            
            dir_C= cos(deg2rad(theta)); %% x is the col
            dir_R= -sin(deg2rad(theta));
            
            Fpush_R=Knpf*dir_R;
            Fpush_C=Knpf*dir_C ;  %% forces in x direction
            Fpush = sqrt((Fpush_R)^2+(Fpush_C)^2); % magnitude of actin pushing force
            
            cosPhi= (Act2Myo_r*Fpush_R+Act2Myo_c*Fpush_C)/(tn_Act2Myo*Fpush);  % calculate cosPhi
            
            if cosPhi>=0  % if cosPhi>0 --> apply forces else dont apply
                
                Fcap_R(i)=  Fpush_R*cosPhi;
                Fcap_C(i)=  Fpush_C*cosPhi;
            else
                
                Fcap_R(i)=  0;
                Fcap_C(i)=  0;
                
            end
            
        end
        
        Fcappush(NodeID,ROW)=sum(Fcap_R); % sum up the different actin node forces
        Fcappush(NodeID,COL)=sum(Fcap_C);
        
        
        %%%%%%----- Net FORCES APPLIED TO MYOSIN NODE -------%%%%%%%%%
        F_myoT(NodeID,ROW)= F_am(NodeID,ROW)+F_Spring(NodeID,ROW);   % these are only through myosin material
        F_myoT(NodeID,COL)= F_am(NodeID,COL)+F_Spring(NodeID,COL);
        
        F_MyoNet(NodeID,ROW)= F_am(NodeID,ROW)+F_Spring(NodeID,ROW)+ Fcappush(NodeID,ROW); % These are total forces from both materials
        F_MyoNet(NodeID,COL)= F_am(NodeID,COL)+F_Spring(NodeID,COL)+ Fcappush(NodeID,COL);
        
        
    end
    
end
if mod(mcs,500)==0
    cd (absoluteFolderPath)
    filename=strcat('csvdata_',num2str(mcs),'.csv');
    
    csvwrite(filename,csvdata(1:count,:));
    cd ..
end
%% --- This section calculates actin forces -----

%----Implementing Equation-5------
%%---Following lines calculate actin polymerization force directed outward
%%---from nucelation points.

RandNodes=randperm(Ac_NodeCount);
for i= 1:Ac_NodeCount
    NodeID= RandNodes(i); % randomly pick up an actin node
    
    npfID=  Ac_Node(NodeID,ASSCNPF);  %  corresponding nucleation point
    
    theta = Ac_Node(NodeID,THETA); % angle of the filament used to get its direction
    
    dir_C= cos(deg2rad(theta)); %% x is the col
    dir_R= -sin(deg2rad(theta)); %% taking -ve sign because the y axis of the image is pointed downwards
    
    if Ac_Node(NodeID,LEN) > Lth % cap the actin filament after a length threshold
        Fpoly(NodeID,ROW) =0 ;    % If reached certain length growth stops.
        Fpoly(NodeID,COL) =0 ;
    else
        
        Fpoly(NodeID,ROW)= Knpf*dir_R;  % polymerization force
        Fpoly(NodeID,COL)= Knpf*dir_C;
    end
    
    %%----Following lines calculate forces from the myosin node directed towards actin node----%%%%%%%
     %   ----Implementing Supplementary Equation-1 ---------
    totalPoints=0;
    points=0;
    noOfNodesPosition=0;
    
    colNo= Ac_Node(NodeID,COL); %get the position of actin node
    rowNo= Ac_Node(NodeID,ROW);
    for R= rowNo-SR : rowNo+SR   %%% searching around the actin node in search radius
        for C= colNo-SR:colNo+SR
            if Mc_NodeNo(round(R),round(C))> 0  %% if myosin node is found in vicinity
                totalPoints= totalPoints+1;
                points(totalPoints,1) = R; %% Save the position of myosin
                points(totalPoints,2) = C;
                noOfNodesPosition(totalPoints) = NodeCountPerPixel((round(R)),(round(C))); % check how many nodes are actually present at this location
            end
        end
    end
    
    %---Following lines will calculate forces from myosin nodes found in the
    %---above determined positions
    
    RhoMyo_R=0;
    RhoMyo_C=0;
    noOfNodesPosition;
    for i= 1: totalPoints
        R= points(i,1) ;
        C=points(i,2);
        
        Myo2Act_r= rowNo - R;
        Myo2Act_c= colNo - C;  % vector directed towards the actin from myo
        tn_Myo2Act= sqrt((Myo2Act_c)^2+(Myo2Act_r)^2);
        
        
        MyoID=  Mc_NodeNo(round(R),round(C));
        
        F_myoR=F_myoT(MyoID,ROW);  % % forces in -y direction
        F_myoC=F_myoT(MyoID,COL);  %% forces in x direction
        F_myo= sqrt((F_myoR)^2+(F_myoC)^2);  %% net magnitude of force from myosin based on actomyosin material forces
        
        cosPhi= (Myo2Act_r*F_myoR+Myo2Act_c*F_myoC)/(tn_Myo2Act*F_myo);
        
        if cosPhi>=0  % force will be applied only if cosPhi>0
            
            RhoMyo_R(i)= noOfNodesPosition(totalPoints)* F_myoR*cosPhi; %force magnitude times the # of nodes at the position
            RhoMyo_C(i)= noOfNodesPosition(totalPoints)* F_myoC*cosPhi; % gives an assesment of cumulative effect of all nodes in that position
        else
            RhoMyo_R(i)= 0;
            RhoMyo_C(i)=  0;
            
        end
        
    end
    Rho_dynR(NodeID,1)=sum(RhoMyo_R);  % summing up all myo node forces
    Rho_dynC(NodeID,1)=sum(RhoMyo_C);
    
    
    Rho_myo(NodeID,ROW)= Rho_dynR(NodeID,1);%
    Rho_myo(NodeID,COL)= Rho_dynC(NodeID,1);%
    
        
    %%%%%%----- Net FORCES APPLIED TO ACTIN NODE -------%%%%%%%%%
    F_AcNet(NodeID,ROW)= Fpoly(NodeID,ROW)+Rho_myo(NodeID,ROW);
    F_AcNet(NodeID,COL)= Fpoly(NodeID,COL)+Rho_myo(NodeID,COL);
    
end

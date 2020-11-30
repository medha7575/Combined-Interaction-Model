%% This file moves actin node to a new position depending of the force applied on the node.
% ---When actin filament length reaches a threshold, it stopes further
% growth of the respective actin node. When actin filament reaches a
% predefined length, it enables branching of that filament, by creating
% new actin nodes and nucleation points.

%% -----Implementing equation-6 -------

RandNodes=randperm(Ac_NodeCount);
for i = 1: Ac_NodeCount
    NodeID= RandNodes(i);
    
    npfID=  Ac_Node(NodeID,ASSCNPF);
    
    % % rounding to 2 digit to control spatial resolution
    
    NewrowNo=  round(Ac_Node(NodeID,ROW)+  F_AcNet(NodeID,ROW),2); % New actin barbed end position
    NewcolNo=  round(Ac_Node(NodeID,COL)+  F_AcNet(NodeID,COL),2);  % New actin barbed end position
    
    
    for R= NewrowNo-Aura : NewrowNo+Aura   % search all around the node for myosin
        for C= NewcolNo-Aura:NewcolNo+Aura
            
            % % % if there is a myosin node in the vicinity of the new
            % position, actin node is kept in the old position
            if Mc_NodeNo(round(R),round(C))>-1  %% check if node is a myosin node
                NewrowNo=Ac_Node(NodeID,ROW);
                NewcolNo=Ac_Node(NodeID,COL);
            end
            
        end
    end
    
    %%---Updating actin node features----
    Ac_NodeNo(round(Ac_Node(NodeID,ROW)),round(Ac_Node(NodeID,COL)))= -1;
  
    Ac_Node(NodeID,ROW)=NewrowNo;
    Ac_Node(NodeID,COL)=NewcolNo;
    Ac_Node(NodeID,LEN)= sqrt((NewrowNo- NPF(npfID,1))^2+(NewcolNo- NPF(npfID,2))^2); %length of filament
    Ac_Node(NodeID,AGE)= mcs- Ac_Node(NodeID,MCSATBIRTH);  %% age of actin updated
    
    Ac_NodeNo(round(Ac_Node(NodeID,ROW)),round(Ac_Node(NodeID,COL)))=NodeID;
    
    Combo(round(NewrowNo),round(NewcolNo))= Ac;  % Add the barbed end to the combo lattice
    
    %% %%%%%%%% Creating Actin Branches %%%%%%%%%%%%%%%%%
    if Ac_Node(NodeID,AGE)> AgeBr
        
        if Ac_Node(NodeID,BRANCHORNOT) ==0 %check if no branch has been generated out of this branch till now
            
            theta_Node= Ac_Node(NodeID,THETA); % fetching angle of the current node
            Rflip= rand(1);
            if Rflip>0.5
                flip = 1;
            else
                flip = -1;
            end
            
            theta_NPF= theta_Node+ flip* 70 ;% degree allignment of the new branch FIXED 70
            
            BrPoint_row = (NPF(npfID,ROW)+Ac_Node(NodeID,ROW))/2;
            BrPoint_col = (NPF(npfID,COL)+Ac_Node(NodeID,COL))/2;
            %
            dir_C= cos(deg2rad(theta_NPF)); % x is col
            dir_R= -sin(deg2rad(theta_NPF)); % y is negavtive for image
            
            Fnpf_c= Knpf*dir_C;
            Fnpf_r= Knpf*dir_R;
            
            % % rounding to 2 didgit to control spatial resolution:
            
            BRrowNo=round((BrPoint_row+Fnpf_r),2); % branched actin position
            BRcolNo=round((BrPoint_col+Fnpf_c),2);  % branched actin position
            
            
            if Mc_NodeNo(round(BRrowNo),round(BRcolNo))> -1  %% just check for presence of mysoin on that point
                disp('cant branch');
            else
                
                Ac_Node(NodeID,BRANCHORNOT) = 1;   %%%now branched
                %----creating new nucleation point for the branch----
                NPFCount=NPFCount+1;   %% new npf
                
                NPF(NPFCount,ROW)= BrPoint_row;
                NPF(NPFCount,COL)= BrPoint_col;
                
                Ac_NodeCount=Ac_NodeCount+1;
                
                % Updating Node array for the newborn actin which just
                % branched
                Ac_Node(Ac_NodeCount,ROW)=BRrowNo;
                Ac_Node(Ac_NodeCount,COL)=BRcolNo;
                Ac_Node(Ac_NodeCount,BRANCHORNOT)=0; %unbranched
                Ac_Node(Ac_NodeCount,THETA)=theta_NPF ;%direction of new actin filamnet
                Ac_Node(Ac_NodeCount,LEN)= sqrt((BRrowNo-NPF(NPFCount,1))^2+(BRcolNo-NPF(NPFCount,2))^2);
                Ac_Node(Ac_NodeCount,ASSCNPF)= NPFCount;  % associated NPF
                Ac_Node(Ac_NodeCount,MCSATBIRTH)= mcs;  % time of creation
                Ac_Node(Ac_NodeCount,AGE)= mcs-Ac_Node(Ac_NodeCount,MCSATBIRTH);  % age of actin
                
                Combo(round(Ac_Node(Ac_NodeCount,ROW)),round(Ac_Node(Ac_NodeCount,COL)))= Ac;
                Ac_NodeNo(round(Ac_Node(Ac_NodeCount,ROW)),round(Ac_Node(Ac_NodeCount,COL)))= Ac_NodeCount;
                
                NewActin=NewActin+1;
                
            end
        end
        
    end
end




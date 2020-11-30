      %% This file moves actin node to a new position depending of the force applied on the node.
      %----When moving a mysoin node, it checks for availability of the
      %----new position in two checks. 1. Node/pixel saturation 2. actin node
      %----occupancy.
      %% Lastly, this file applies the elevated myosin upon saturation. Thus testing the hypothesis.
      
    %% -----Implementing Equation-4-----  
        RandNodes=randperm(Mc_NodeCount);
    for i = 1: Mc_NodeCount
        NodeID= RandNodes(i);
    
        colNo= Mc_Node(NodeID,COL);
        rowNo= Mc_Node(NodeID,ROW);
       
        % % rounding to 2 didgit to control spatial resolution
        NewcolNo = round(colNo+ F_MyoNet(NodeID,COL),2);
        NewrowNo = round(rowNo+ F_MyoNet(NodeID,ROW),2);
        
        %%---check-1:  Node/pixel saturation
        
         noOfNodesAtNewPosition = NodeCountPerPixel((round(NewrowNo)),(round(NewcolNo)));
        
        if noOfNodesAtNewPosition == Max_myosin_nodes_per_pixels
                    NewcolNo= colNo;
                    NewrowNo= rowNo;
        end
        %%--- Check-2: Actin node occupancy: checking if there is actin in the way
        
        for R= NewrowNo-Aura:NewrowNo+Aura
            for C= NewcolNo-Aura:NewcolNo+Aura
                if Ac_NodeNo ( round(R),round(C))> -1
                    NewcolNo= colNo;
                    NewrowNo= rowNo;
                end
            end
        end
        
         %%---Updating myosin node features----
        Mc_NodeNo(round(Mc_Node(NodeID,ROW)),round(Mc_Node(NodeID,COL)))= -1;   %% reset the NodeNo in order to create it again afresh
        
        Mc_Node(NodeID,ROW)= NewrowNo; 
        Mc_Node(NodeID,COL)= NewcolNo;
        
        Mc_NodeNo(round(NewrowNo),round(NewcolNo))=NodeID;  % new place has you
        
        Combo(round(NewrowNo),round(NewcolNo))=Mc; % refresh combo
        
        NodeCountPerPixel((round(NewrowNo)),(round(NewcolNo))) = NodeCountPerPixel((round(NewrowNo)),(round(NewcolNo))) + 1;
        NodeCountPerPixel((round(rowNo)), (round(colNo))) = NodeCountPerPixel((round(rowNo)),(round(colNo))) - 1 ;
        
        
       %% if trying to reach a staurated point then Increase myosin
        
         if noOfNodesAtNewPosition>5 %%increasing myosin when density per pixel increases 5
        
         Mc_Node(NodeID,MYO_CONC)=0.6;
         
         end
        
        
    end
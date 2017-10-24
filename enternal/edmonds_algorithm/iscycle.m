function[dist,path]=iscycle(G,S,D)
   if size(G,1)>0
       MAXN=max([S,D,max( unique(G(:,1:2)) )]);
       G=[G;MAXN,MAXN,0];
    DG = sparse(G(:,1),G(:,2),G(:,3));
    
        
  
    [dist,path]=graphshortestpath(DG,S,D);
    %if there is no path from S to D try D to S
    if dist==Inf
        [dist,path]=graphshortestpath(DG,D,S);
    end 
   else
     dist=Inf;
     path=[];
   end


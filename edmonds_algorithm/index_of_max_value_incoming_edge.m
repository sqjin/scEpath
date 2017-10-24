function[EDGE_VALUE,EDGE_INDEX]=index_of_max_value_incoming_edge(G,VERTEX)
%first find incoming edges
if numel(G)==0
    EDGE_VALUE=-1;
    EDGE_INDEX=NaN;
    return;
end
%
INDICES=find(G(:,2)==VERTEX);
VALUES=(G(INDICES,3));
[EDGE_VALUE,LOC]=min(VALUES);
EDGE_INDEX=INDICES(LOC);

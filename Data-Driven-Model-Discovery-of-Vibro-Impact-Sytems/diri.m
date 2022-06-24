function z=diri(e,ben)
%za Nx2 cordinates of the nodes; e is the element numbers; ben: nodes on the bdd

z=zeros(max(size(e)),1);

for i=1:max(size(e))
    
    if ismember(e(i),ben)        
        z(i)=0;   % this does not have to be zero. it can be any constant. For the cantilever beam problem the left end dofs correspond to zero.    
    end        
    
end

end
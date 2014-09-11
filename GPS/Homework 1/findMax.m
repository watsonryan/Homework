function [Max]=findMax(ans)

k=1;
y=0;
for i=1:37
    y=y+1
       
        Max(y,1)=max(ans(y,:))
    
end
end
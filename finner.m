function super2= finner(xt,bin)
%super=[];
super2=[];
 n=length(xt);
   for i=1:1:n-1
    dxtt(i)=(xt(i+1,:)-xt(i,:))./bin;
    for j=1:1:bin
       xxt(j)= xt(i,:)+dxtt(i)*j;
    end
   % super(i,:)=[xt(i) xxt];
   super=[xt(i) xxt];
   super2=[super2 super(:,1:end-1)];
   end

end
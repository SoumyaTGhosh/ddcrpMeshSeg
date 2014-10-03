function [tables,t]=fast_cc(N,assig)
%computes table assignments from links.


% create an adjacency matrix.
An = (sparse(1:N,assig,ones(1,N),N,N));

%make it symmetric
An = max(An,An');
tables  = components(An);
t = max(tables);

end




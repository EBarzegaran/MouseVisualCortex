function M_OUT = Breakload2cells(M,dim, subdims)
% devides the first dimension of cell # dims into subdims, for unfolded
% dimension in PARAFAC

M_temp = reshape(M{dim},[squeeze(subdims) size(M{dim},2)]);
for d=1:numel(subdims)
    M_temp2 = M_temp;
    for d2 = setdiff(1:numel(subdims), d)
        M_temp2 = max(M_temp2,[], d2);
    end
    M_OUT{d} = squeeze( M_temp2 );
end
M_OUT = [M(1:dim-1) M_OUT M(dim+1:end)];
end
function idz = findz(vec)

idz = unique(find(vec(2:end).*vec(1:end-1) <= 0));

for j = 1:length(idz)
    if abs(vec(idz(j))) > abs(vec(idz(j)+1))
        idz(j) = idz(j)+1;
    end
end

end
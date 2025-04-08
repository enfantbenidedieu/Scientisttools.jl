
function intersect(a,b)
    ia = findall(in(b), a)
    return unique(view(a,ia))
end;
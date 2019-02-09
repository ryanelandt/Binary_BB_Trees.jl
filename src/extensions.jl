
function NumericalTricks.volume(eM::eMesh{T1,Tet}) where {T1}
    vol = 0.0
    for k = 1:n_tet(eM)
        vol += volume(eM.point[eM.tet[k]])
    end
    return vol
end

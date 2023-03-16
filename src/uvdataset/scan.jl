export compute_scans!


"""
    compute_scans!(dataset::UVDataSet, minΔt::Number=-1)

Identify Scan IDs from a given time series and minimum scan seperation.
Returns an integer array for the identified scan IDs.
- dataset::DimStack
- `minΔt`: minimum seperation between the scans in seconds.
"""
function compute_scans!(uvd::UVDataSet, minΔt::Number=-1)
    uvd[:baseline] = compute_scans(uvd[:baseline], minΔt)
    return
end


"""
    compute_scans(dataset::DimStack, minΔt::Number=-1)

Identify Scan IDs from a given time series and minimum scan seperation.
Returns an integer array for the identified scan IDs.
- dataset::DimStack
- `minΔt`: minimum seperation between the scans in seconds.
"""
function compute_scans(dataset::DimStack, minΔt::Number=-1)
    mjd = dataset[:mjd].data

    if minΔt > 0
        Δtthres = minΔt / 86400.0
    else
        Δtthres = -1
    end

    scanid = identify_scans(mjd, Δtthres)
    d = Dim{:data}(1:length(scanid))
    return append(dataset, DimArray(data=scanid, dims=(d), name=:scan))
end


"""
    identify_scans(t::AbstractArray, minΔt::Number=-1)

Identify Scan IDs from a given time series and minimum scan seperation.
Returns an integer array for the identified scan IDs.
- `t::AbstractArray`: time stamps
- `minΔt`: minimum seperation between the scans. should be in the same time unit with `t`.
    Default to 2 x minimum seperation
"""
function identify_scans(t::AbstractArray, minΔt::Number=-1)
    # get the unique set of the time stamps
    t̃, tidx = unique_ids(t)

    # measure the gaps of the time stamps
    Δt̃ = diff(t̃)

    # minimum time
    if minΔt < 0
        Δt̃thres = minimum(Δt̃) * 2
    else
        Δt̃thres = minΔt
    end

    # the array for the scan id
    nt̃ = length(t̃)
    s̃ = zeros(Int64, nt̃)
    s̃i = 1
    for i in 1:nt̃
        if i != 1
            if Δt̃[i-1] > Δt̃thres
                s̃i += 1
            end
        end
        s̃[i] = s̃i
    end

    return s̃[tidx]
end
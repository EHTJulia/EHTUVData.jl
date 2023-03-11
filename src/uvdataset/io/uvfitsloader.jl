export load_uvfits


"""
    load_uvfits(filename::AbstractString)::UVDataSet

load visibility data from the given uvfits file. 
Data will be loaded through astropy.io.fits module
loaded by PyCall.
"""
function load_uvfits(filename::AbstractString; ex=ThreadedEx())::UVDataSet
    # Load pyfits
    pf = pyimport("astropy.io.fits")

    # Open UVFITS file
    hdulist = pf.open(filename)

    # Load GroupHDU, HDUs for AIPS AN/FQ Tables
    ghdu, antab, fqtab = hdulist2hdus(hdulist)

    # Get the baseline-based UV data set (i.e. visibility) from GroupHDU
    blds = hdulist2bl(ghdu, ex=ex)

    # Get the frequency information
    blds = concat(blds, hdulist2freq(ghdu, antab, fqtab))

    # Get the antds
    antds = hdulist2ant(antab)

    # Get the metadata
    metadata = hdulist2metadata(ghdu, antab)

    # close HDUList
    hdulist.close()

    # group visds and antds in OrderedDict
    datasets = OrderedDict(:antenna => antds, :baseline => blds)

    # Form UVData
    uvdata = UVDataSet(datasets, metadata)

    # compute frequency and uvw
    compute_ν!(uvdata)
    compute_uvw!(uvdata)

    # output
    return uvdata
end


"""
    hdulist2hdus(hdulist)

Read the given hdulist (output of pf.io.open) and return
GroupHDU and HDUs of AIPS AN/FQ Tables.
"""
function hdulist2hdus(hdulist)
    # Number of HDUs
    nhdu = length(hdulist)

    # grouphdu
    ghdu = hdulist[0]

    # read other hdu
    antab = NaN
    fqtab = NaN
    for i in 1:(nhdu-1)
        hdu = hdulist[i]
        if string(hdu.header.get("EXTNAME")) == "AIPS AN"
            if typeof(antab) == Py
                println("Warning: there are multiple AN tables in the UVFITS file. The latest one will be loaded.")
            end
            antab = hdu
        elseif string(hdu.header.get("EXTNAME")) == "AIPS FQ"
            if typeof(fqtab) == Py
                println("Warning: there are multiple FQ tables in the UVFITS file. The latest one will be loaded.")
            end
            fqtab = hdu
        end
    end

    if typeof(antab) != Py
        @throwerror ValueError "The input UVFITS file does not have an AIPS FQ Table."
    end

    if typeof(fqtab) != Py
        @throwerror ValueError "The input UVFITS file does not have an AIPS FQ Table."
    end

    return ghdu, antab, fqtab
end


"""
    hdulist2bl(ghdu; ex=ThreadedEx())
"""
function hdulist2bl(ghdu; ex=ThreadedEx())
    # import numpy
    np = pyimport("numpy")

    # define useful shortcuts
    float64vector(x::Py) = pyconvert(Vector, np.asarray(x, dtype=np.float64))
    int32vector(x::Py) = pyconvert(Vector, np.asarray(x, dtype=np.int32))
    float32array(x::Py) = pyconvert(Array, np.asarray(x, dtype=np.float32))

    # group hdu
    ghdudata = float32array(ghdu.data.data)

    # size of visibility data
    ndata, ndec, nra, nspw, nch, npol, _ = size(ghdudata)

    if (ndec > 1) || (nra > 1)
        print("Warning: GroupHDU has more than a single coordinae (Nra, Ndec) = (", nra, ", ", ndec, "). We will pick up only the first one.")
    end

    # read visibilities
    Vre = ghdudata[:, 1, 1, :, :, :, 1] # Dim: data (time x baseline), spw, ch, pol 
    Vim = ghdudata[:, 1, 1, :, :, :, 2]
    Vwe = ghdudata[:, 1, 1, :, :, :, 3]

    # permutate dimensions
    Vre = permutedims(Vre, (3, 2, 1, 4)) # Dim: ch, spw, data (time x baseline), pol
    Vim = permutedims(Vim, (3, 2, 1, 4))
    Vwe = permutedims(Vwe, (3, 2, 1, 4))

    # compute visibility, sigma and flag
    Vcmp = ComplexF64.(Vre .+ 1im .* Vim)
    σV = (abs.(Float64.(Vwe))) .^ (-0.5)
    flag = Int8.(sign.(Vwe))

    # reset flags based on the value of sigma
    cond = σV .== NaN .|| isinf.(σV) .|| σV == 0
    idx = findall(cond)
    σV[idx] .= 0.0
    flag[idx] .= 0

    # Read random parameters
    paridxes = [-1 for i in 1:9]
    parnames = ghdu.data.parnames
    npar = length(parnames)
    pardata = DataFrame()
    for ipar in 1:npar
        parname = string(parnames[ipar-1])
        if occursin("UU", parname)
            paridxes[1] = ipar
            pardata[!, :usec] = float64vector(ghdu.data.par(ipar - 1))
        elseif occursin("VV", parname)
            paridxes[2] = ipar
            pardata[!, :vsec] = float64vector(ghdu.data.par(ipar - 1))
        elseif occursin("WW", parname)
            paridxes[3] = ipar
            pardata[!, :wsec] = float64vector(ghdu.data.par(ipar - 1))
        elseif occursin("DATE", parname)
            if paridxes[4] < 0
                paridxes[4] = ipar
                pardata[!, :mjd] = jd2mjd(float64vector(ghdu.data.par(ipar - 1)))
            elseif paridxes[5] < 0
                paridxes[5] = ipar
                pardata[!, :mjd] .+= float64vector(ghdu.data.par(ipar - 1))
            else
                println(ipar)
                @throwerror KeyError "Random Parameter have too many `DATE` columns."
            end
        elseif occursin("BASELINE", parname)
            paridxes[6] = ipar
            pardata[!, :baseline] = float64vector(ghdu.data.par(ipar - 1))
        elseif occursin("SOURCE", parname)
            paridxes[7] = ipar
            pardata[!, :source] = int32vector(ghdu.data.par(ipar - 1))
        elseif occursin("INTTIM", parname)
            paridxes[8] = ipar
            pardata[!, :inttime] = float64vector(ghdu.data.par(ipar - 1))
        elseif occursin("FREQSEL", parname)
            paridxes[9] = ipar
            pardata[!, :freqset] = int32vector(ghdu.data.par(ipar - 1))
        end
    end

    # check the loaded random parameters
    if paridxes[7] > 0
        if length(unique(pardata[:source])) > 1
            print("Warning: Group HDU contains data on more than a single source. ")
            print("It will likely cause a problem since this library assumes a single source UVFITS file.\n")
        end
    end

    if paridxes[8] < 0
        print("Warning: Group HDU do not have a random parameter for the integration time. ")
        print("It will be estimated with a minimal time interval of data.\n")
        dmjd = median(diff(sort(unique(pardata[!, :mjd]))))
        pardata[!, :dmjd] = fill(dmjd, ndata)
    else
        pardata[!, :dmjd] = pardata[!, :inttime] ./ 86400
    end

    if paridxes[9] > 0
        if length(unique(pardata[:freqset])) > 1
            print("Warning: Group HDU contains data on more than a single frequency setting. ")
            print("It will likely cause a problem since this library assumes a UVFITS file with a single frequency setup.\n")
        end
    end

    # antenna ID
    subarray = zeros(Int64, ndata)
    pardata[!, :antid1] = zeros(Int64, ndata)
    pardata[!, :antid2] = zeros(Int64, ndata)
    @floop ex for i in 1:ndata
        subid, blid = modf(pardata[i, :baseline])
        subarray[i] = Int64(100 * subid + 1)
        pardata[i, :antid1] = div(blid, 256)
        pardata[i, :antid2] = rem(blid, 256)
    end
    if length(unique(subarray)) > 1
        print("Warning: Group HDU contains data on more than a single subarray. ")
        print("It will likely cause a problem since this library assumes a UVFITS file with a single subarray.\n")
    end

    # polarization
    dp = pyconvert(Int64, ghdu.header.get("CDELT3"))
    ipref = pyconvert(Int64, ghdu.header.get("CRPIX3"))
    pref = pyconvert(Int64, ghdu.header.get("CRVAL3"))
    polids = (dp*(1-ipref)+pref):dp:(dp*(npol-ipref)+pref)
    pol = [uvfits_polid2name[string(polid)] for polid in polids]

    # form dimensional arrays
    c = Dim{:ch}(collect(1:nch))
    s = Dim{:spw}(collect(1:nspw))
    d = Dim{:data}(collect(1:ndata))
    p = Dim{:pol}(collect(1:npol))

    # ch, spw, data (time x baseline), pol
    blds = DimStack(
        DimArray(data=Vcmp, dims=(c, s, d, p), name=:visibility),
        DimArray(data=σV, dims=(c, s, d, p), name=:sigma),
        DimArray(data=flag, dims=(c, s, d, p), name=:flag),
        DimArray(data=pol, dims=(p), name=:polarization),
        DimArray(data=pardata[!, :usec], dims=(d), name=:usec),
        DimArray(data=pardata[!, :vsec], dims=(d), name=:vsec),
        DimArray(data=pardata[!, :wsec], dims=(d), name=:wsec),
        DimArray(data=pardata[!, :mjd], dims=(d), name=:mjd),
        DimArray(data=pardata[!, :dmjd], dims=(d), name=:Δmjd),
        DimArray(data=pardata[!, :antid1], dims=(d), name=:antid1),
        DimArray(data=pardata[!, :antid2], dims=(d), name=:antid2),
    )

    return blds
end

"""
    hdulist2freq(ghdu, antab, fqtab)
"""
function hdulist2freq(ghdu, antab, fqtab)
    # pyimport
    np = pyimport("numpy")

    # define useful shortcuts
    vector(x::Py) = pyconvert(Vector, x)

    # Reference Frequency
    reffreq = pyconvert(Float64, antab.header.get("FREQ"))

    # Get data dimension
    _, _, _, nspw, nch, _, _ = pyconvert(Tuple, ghdu.data.data.shape)

    # Check FREQSEL
    nfreqset = length(fqtab.data["FRQSEL"])
    if nfreqset > 1
        println("Input FQ Table has more than a single Frequency setting. This library only handles a UVFITS file with a single frequency setting. The first setting will be loaded.")
    end

    function arraylize(input)
        if isa(input, Number)
            return [input]
        else
            return input
        end
    end

    # Get frequency settings
    spwfreq = vector(np.asarray(fqtab.data["IF FREQ"][0], dtype=np.float64).reshape([nspw]))
    chbw = vector(np.asarray(fqtab.data["CH WIDTH"][0], dtype=np.float64).reshape([nspw]))
    sideband = vector(np.asarray(fqtab.data["SIDEBAND"][0], dtype=np.float64).reshape([nspw]))

    # Axis
    s = Dim{:spw}(collect(1:nspw))

    fqds = DimStack(
        DimArray(data=reffreq .+ spwfreq, dims=s, name=:νspw),
        DimArray(data=chbw, dims=s, name=:Δνch),
        DimArray(data=sideband, dims=s, name=:sideband),
    )
    return fqds
end

"""
    hdulist2ant(antab)
"""
function hdulist2ant(antab)
    # Load numpy
    np = pyimport("numpy")

    # define useful shortcuts
    float64vector(x::Py) = pyconvert(Vector, np.asarray(x, dtype=np.float64))
    int32vector(x::Py) = pyconvert(Vector, np.asarray(x, dtype=np.int32))
    float64array(x::Py) = pyconvert(Array, np.asarray(x, dtype=np.float64))

    # Get the antenna infromation
    nant = length(antab.data)

    # Get the antenna name
    antname = [string(name) for name in antab.data["ANNAME"]]
    xyz = float64array(antab.data["STABXYZ"])

    # Check polarization labels
    pola = np.unique(antab.data["POLTYA"])
    polb = np.unique(antab.data["POLTYB"])
    if length(pola) > 1
        @warn "Mixed polarization feed: POLTYB have more than a single polarization across the array."
    end
    if length(polb) > 1
        @warn "Mixed polarization feed: POLTYB have more than a single polarization across the array."
    end
    feed1 = [string(pol) for pol in antab.data["POLTYA"]]
    feed2 = [string(pol) for pol in antab.data["POLTYB"]]

    # Parse Field Rotation Information
    #   See AIPS MEMO 117
    #      0: ALT-AZ, 1: Eq, 2: Orbit, 3: X-Y, 4: Naismith-R, 5: Naismith-L
    #      6: Manual
    mntsta = int32vector(antab.data["MNTSTA"])
    fr_pa_coeff = ones(nant)
    fr_el_coeff = zeros(nant)
    fr_offset = zeros(nant)
    for i in 1:nant
        if mntsta[i] == 0  # azel
            fr_pa_coeff[i] = 1
            fr_el_coeff[i] = 0
        elseif mntsta[i] == 1  # Equatorial
            fr_pa_coeff[i] = 0
            fr_el_coeff[i] = 0
        elseif mntsta[i] == 4  # Nasmyth-R
            fr_pa_coeff[i] = 1
            fr_el_coeff[i] = 1
        elseif mntsta[i] == 5  # Nasmyth-L
            fr_pa_coeff[i] = 1
            fr_el_coeff[i] = -1
        end
    end

    # Antenna Type
    anttype = [:ground for i in 1:nant]

    # Set dimensions
    a = Dim{:ant}(collect(1:nant))

    # create data set
    antds = DimStack(
        DimArray(data=antname, dims=a, name=:name),
        DimArray(data=xyz[:, 1], dims=a, name=:x),
        DimArray(data=xyz[:, 2], dims=a, name=:y),
        DimArray(data=xyz[:, 3], dims=a, name=:z),
        DimArray(data=anttype, dims=a, name=:type),
        DimArray(data=fr_pa_coeff, dims=a, name=:fr_pa_coeff),
        DimArray(data=fr_el_coeff, dims=a, name=:fr_el_coeff),
        DimArray(data=fr_offset, dims=a, name=:fr_offset),
        DimArray(data=feed1, dims=a, name=:feed1),
        DimArray(data=feed2, dims=a, name=:feed2),
    )

    return antds
end

"""
    hdulist2metadata(ghdu, antab)
"""
function hdulist2metadata(ghdu, antab)
    # Output metadata
    metadata = OrderedDict()

    # Useful short cuts
    float64(x::Py) = pyconvert(Float64, x)
    any(x::Py) = pyconvert(Any, x)

    # Initialize
    for key in keys(uvdataset_metadata_default)
        metadata[key] = uvdataset_metadata_default[key]
    end

    metadata[:instrument] = string(antab.header.get("ARRNAM"))
    metadata[:observer] = string(ghdu.header.get("OBSERVER"))

    metadata[:source] = string(ghdu.header.get("OBJECT"))
    metadata[:ra] = deg2rad(float64(ghdu.header.get("CRVAL6")))
    metadata[:dec] = deg2rad(float64(ghdu.header.get("CRVAL7")))

    equinox = any(ghdu.header.get("EQUINOX"))
    epoch = any(ghdu.header.get("EPOCH"))
    if isa(equinox, Nothing) == false
        if isa(equinox, Number)
            metadata[:equinox] = float(equinox)
        elseif isa(equinox, String)
            metadata[:equinox] = parse(Float64, replace(replace(uppercase(equinox), "J" => ""), "B" => ""))
        else
            metadata[:equinox] = -1
        end
        metadata[:coordsys] = "fk5"
    elseif isa(epoch, Nothing) == false
        if isa(epoch, Number)
            metadata[:equinox] = float(epoch)
        elseif isa(equinox, String)
            metadata[:equinox] = parse(Float64, replace(replace(uppercase(epoch), "J" => ""), "B" => ""))
        else
            metadata[:equinox] = -1
        end
    else
        metadata[:equinox] = -1
    end

    if metadata[:equinox] < 0
        metadata[:coordsys] = "icrs"
    end

    return metadata
end
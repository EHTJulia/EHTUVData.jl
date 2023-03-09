export save_uvfits!


function save_uvfits!(uvd::UVDataSet, filename::AbstractString; ex=ThreadedEx())
    # open hdulist
    hdulist = [
        uvdataset2ghdu(uvd, ex=ex),
        uvdataset2antab(uvd, ex=ex),
        uvdataset2fqtab(uvd, ex=ex),
    ]

    py"""
    import astropy.io.fits as pf

    def writeto(filename, hdulist, overwrite):
        pf.HDUList(hdulist).writeto(filename, overwrite=overwrite)
    """
    py"writeto"(filename, hdulist, overwrite=true)

    return
end


function uvdataset2ghdu(uvd::UVDataSet; ex=ThreadedEx())
    # load pyfits
    copy!(pyfits, pyimport_conda("astropy.io.fits", "astropy"))

    # shortcut to each data set
    blds = uvd[:baseline]

    # get data size
    nch, nspw, ndata, npol = size(blds)

    # Time: jdref, mjdmin, mjd1, mjd2
    #    the first mjd
    mjdmin = floor(minimum(blds[:mjd].data))
    #    the reference date (for some reason Jan 1st of the observing year)
    Δmjd1 = zeros(ndata)
    Δmjd2 = zeros(ndata)
    @floop ex for i = 1:ndata
        Δmjd2[i], Δmjd1[i] = modf(blds[:mjd].data[i])
    end

    # Frequency
    Δνch = blds[:Δνch].data[1]
    νref = minimum(blds[:ν].data)

    # Source Parameters
    obsra = rad2deg(uvd.metadata[:ra])
    obsdec = rad2deg(uvd.metadata[:dec])

    # uvfits dimension ndata, ndec, nra, nspw, nch, npol, 3
    # uvdata dimension nch, nspw, ndata, npol

    # Visibility Data
    visarr = zeros(Float32, ndata, 1, 1, nspw, nch, npol, 3)
    @floop ex for ipol = 1:npol, idata = 1:ndata, ispw = 1:nspw, ich = 1:nch
        # get flag and uncertainty
        flag = blds[:flag].data[ich, ispw, idata, ipol]
        sigma = blds[:sigma].data[ich, ispw, idata, ipol]

        if sigma == 0 || isinf(sigma) == true || flag == 0
            continue
        end

        vcmp = blds[:visibility].data[ich, ispw, idata, ipol]
        visarr[idata, 1, 1, ispw, ich, ipol, 1] = real(vcmp)
        visarr[idata, 1, 1, ispw, ich, ipol, 2] = imag(vcmp)

        weight = 1 / sigma^2

        if flag > 0
            visarr[idata, 1, 1, ispw, ich, ipol, 3] = weight
        else
            visarr[idata, 1, 1, ispw, ich, ipol, 3] = -weight
        end
    end

    # create group HDU
    ghdu = pyfits.GroupsHDU(
        pyfits.GroupData(
            input=visarr,
            parnames=[
                "UU---SIN",
                "VV---SIN",
                "WW---SIN",
                "BASELINE",
                "DATE",
                "DATE",
                "INTTIM"
            ],
            pardata=[
                Float32.(blds[:usec].data),
                Float32.(blds[:vsec].data),
                Float32.(blds[:wsec].data),
                Float32.(blds[:antid1].data .* 256 .+ blds[:antid2].data),
                mjd2jd(Δmjd1),
                Δmjd2,
                Float32.(blds[:Δmjd].data .* 86400)
            ],
            bscale=1.0,
            bzero=0.0,
            bitpix=-32
        )
    )

    # PTYPE Header Cards
    for i in 1:7
        ghdu.header.set(format("PZERO{}", i), 0.0, "")
        ghdu.header.set(format("PSCAL{}", i), 1.0, "")
    end

    # CTYPE Header Cards
    #   Complex
    ghdu.header.set("CTYPE2", "COMPLEX", "real, imaginary, weight")
    ghdu.header.set("CRPIX2", 1.0, "")
    ghdu.header.set("CRVAL2", 1.0, "")
    ghdu.header.set("CDELT2", 1.0, "")
    ghdu.header.set("CROTA2", 0.0, "")
    #   Polarizaiton
    firstpol = uvfits_polname2id[blds[:polarization].data[1]]
    ghdu.header.set("CTYPE3", "STOKES", "stokes parameter")
    ghdu.header.set("CRPIX3", 1.0, "")
    ghdu.header.set("CRVAL3", float(firstpol), "")
    ghdu.header.set("CDELT3", float(sign(firstpol)), "")
    ghdu.header.set("CROTA3", 0.0, "")
    #   FREQ
    ghdu.header.set("CTYPE4", "FREQ", "frequency (spectral channel) in Hz")
    ghdu.header.set("CRPIX4", 1.0, "")
    ghdu.header.set("CRVAL4", νref, "")
    ghdu.header.set("CDELT4", blds[:Δνch].data[1], "")
    ghdu.header.set("CROTA4", 0.0, "")
    #   IF
    ghdu.header.set("CTYPE5", "IF", "spectral window number")
    ghdu.header.set("CRPIX5", 1.0, "")
    ghdu.header.set("CRVAL5", 1.0, "")
    ghdu.header.set("CDELT5", 1.0, "")
    ghdu.header.set("CROTA5", 0.0, "")
    #   RA
    ghdu.header.set("CTYPE6", "RA", "right ascention of the phase center in degree")
    ghdu.header.set("CRPIX6", 1.0, "")
    ghdu.header.set("CRVAL6", obsra, "")
    ghdu.header.set("CDELT6", 1.0, "")
    ghdu.header.set("CROTA6", 0.0, "")
    #   DEC
    ghdu.header.set("CTYPE7", "DEC", "declination of the phase center in degree")
    ghdu.header.set("CRPIX7", 1.0, "")
    ghdu.header.set("CRVAL7", obsdec, "")
    ghdu.header.set("CDELT7", 1.0, "")
    ghdu.header.set("CROTA7", 0.0, "")

    # Other entries
    #   Time Stamp
    ghdu.header.set("MJD", mjdmin, "")
    ghdu.header.set("DATE-OBS", Dates.format(mjd2datetime(mjdmin), DateFormat("yyyy-mm-dd")), "observational date")
    #   Observation Info
    ghdu.header.set("TELESCOP", uvd.metadata[:instrument], "telescope name")
    ghdu.header.set("INSTRUME", uvd.metadata[:instrument], "instrument name")
    ghdu.header.set("OBSERVER", uvd.metadata[:observer], "observer name")
    ghdu.header.set("OBJECT", uvd.metadata[:source], "source name")
    #   Target Coordinates
    ghdu.header.set("OBSRA", obsra, "phase center: right ascention in degree")
    ghdu.header.set("OBSDEC", obsdec, "phase center: declination in degree")
    ghdu.header.set("EQUINOX", uvd.metadata[:equinox], "equinox of source coordinates")
    #   Frequency related metadata
    ghdu.header.set("VELREF", 3, ">256 radio, 1 lsr, 2 hel, 3 obs")
    ghdu.header.set("ALTRVAL", 0.0, "altenate FREQ/VEL ref value")
    ghdu.header.set("ALTRPIX", 1.0, "altenate FREQ/VEL ref pixel")
    ghdu.header.set("RESTFREQ", νref, "rest frequency in Hz")
    #   Scaling of visbility products
    ghdu.header.set("BSCALE", 1.0, "real = tape * bscale + bzero")
    ghdu.header.set("BZERO", 0.0, "")
    ghdu.header.set("BUNIT", "JY", "unit of flux")

    # Add comments to the exsiting block
    #   Parameter Type
    ghdu.header.set("NAXIS1", nothing, "no stand image just group")
    ghdu.header.set("PTYPE1", nothing, "u baseline coordinates in seconds")
    ghdu.header.set("PTYPE2", nothing, "v baseline coordinates in seconds")
    ghdu.header.set("PTYPE3", nothing, "w baseline coordinates in seconds")
    ghdu.header.set("PTYPE4", nothing, "julian date 1")
    ghdu.header.set("PTYPE5", nothing, "julian date 2")
    ghdu.header.set("PTYPE6", nothing, "integration time in seconds")
    return ghdu
end


function uvdataset2antab(uvd::UVDataSet; ex=ThreadedEx())
    # load pyfits
    copy!(pyfits, pyimport_conda("astropy.io.fits", "astropy"))
    copy!(numpy, pyimport_conda("numpy", "numpy"))

    # quick shortcuts to data set
    antds = uvd[:antenna]
    blds = uvd[:baseline]

    # get data size
    nch, nspw, ndata, npol = size(blds)
    nant = length(antds[:x])
    nfeed = length(antds[:feed])

    # Columns
    #   ANNAME
    c1 = pyfits.Column(
        name="ANNAME", format="8A", unit=" ",
        array=numpy.asarray(uvd[:antenna][:name].data, dtype="|S8")
    )
    #   STABXYZ
    stabxyz = zeros(Float64, nant, 3)
    stabxyz[:, 1] = antds[:x]
    stabxyz[:, 2] = antds[:y]
    stabxyz[:, 3] = antds[:z]
    c2 = pyfits.Column(
        name="STABXYZ", format="3D", unit="METERS",
        array=stabxyz
    )
    #   ORBPARM
    c3 = pyfits.Column(
        name="ORBPARM", format="1D", unit=" ",
        array=zeros(Float64, nant)
    )
    #   NOSTA
    c4 = pyfits.Column(
        name="NOSTA", format="1J", unit=" ",
        array=Int32.(1:nant)
    )
    #   MNTSTA
    mntsta = zeros(Int32, nant)
    for i in 1:nant
        frcoeff = (antds[:fr_pa_coeff].data[i], antds[:fr_el_coeff].data[i])
        if frcoeff == (1, 0)
            mntsta[i] = 0
        elseif frcoeff == (0, 0)
            mntsta[i] = 1
        elseif frcoeff == (1, 1)
            mntsta[i] = 4
        elseif frcoeff == (1, -1)
            mntsta[i] = 5
        else
            @warn "Invalid Feed Rotation Type for Antenna ", i, ": ", frcoeff, ". We will use AZEL."
        end

        if antds[:fr_offset].data[i] != 0
            @warn "Feed Rotation Offset Angle is not zerofor Antenna ", i, ". This will be ignored."
        end
    end
    c5 = pyfits.Column(
        name="MNTSTA", format="1J", unit=" ",
        array=mntsta
    )
    #   STAXOF
    c6 = pyfits.Column(
        name="STAXOF", format="1E", unit="METERS",
        array=zeros(Float32, nant)
    )
    #   POLTYA
    c7 = pyfits.Column(
        name="POLTYA", format="1A", unit="",
        array=fill(antds[:feed].data[1], nant)
    )
    #   POLTYA
    c8 = pyfits.Column(
        name="POLAA", format="1E", unit="DEGREES",
        array=zeros(Float32, nant)
    )
    #   POLTYA
    c9 = pyfits.Column(
        name="POLCALA", format="1E", unit=" ",
        array=zeros(Float32, nant)
    )
    #   POLTYA
    c10 = pyfits.Column(
        name="POLTYB", format="1A", unit=" ",
        array=fill(antds[:feed].data[2], nant)
    )
    #   POLTYA
    c11 = pyfits.Column(
        name="POLAB", format="1E", unit="DEGREES",
        array=zeros(Float32, nant)
    )
    #   POLTYA
    c12 = pyfits.Column(
        name="POLCALB", format="1E", unit=" ",
        array=zeros(Float32, nant)
    )

    # Make HDU
    hdu = pyfits.BinTableHDU.from_columns(
        pyfits.ColDefs([c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12])
    )

    # Add Headers
    hdu.header.set("EXTNAME", "AIPS AN", "extension name")
    hdu.header.set("EXTVER", Int16(1), "subarray number")
    hdu.header.set("EXTLEVEL", Int32(1), "")
    hdu.header.set("ARRAYX", Float32(0.0), "x coordinates of array center (meters)")
    hdu.header.set("ARRAYY", Float32(0.0), "y coordinates of array center (meters)")
    hdu.header.set("ARRAYZ", Float32(0.0), "z coordinates of array center (meters)")
    hdu.header.set("GSTIA0", Float32(0.0), "GST at 0h on reference date (degrees)")
    hdu.header.set("DEGPDY", Float32(0.0), "earth's rotation rate (degrees/day)")
    hdu.header.set("FREQ", Float32(minimum(blds[:ν])), "reference frequency (Hz)")
    hdu.header.set("RDATE", Dates.format(mjd2datetime(minimum(blds[:mjd])), DateFormat("yyyy-mm-dd")), "reference date")
    hdu.header.set("POLARX", Float32(0.0), "x coodinates of North Pole (arcseconds)")
    hdu.header.set("POLARY", Float32(0.0), "y coodinates of North Pole (arcseconds)")
    hdu.header.set("UT1UTC", Float32(0.0), "UT1 - UTC (sec)")
    hdu.header.set("DATUTC", Float32(0.0), "time system - UTC (sec)")
    hdu.header.set("TIMESYS", "UTC", "time system")
    hdu.header.set("ARRNAM", uvd.metadata[:instrument], "array name")
    hdu.header.set("XYZHAND", "RIGHT", "handness of station coordinates")
    hdu.header.set("FRAME", "????", "coordinate frame")
    hdu.header.set("NUMORB", Int16(0), "number of orbital parameters in table")
    hdu.header.set("NO_IF", Int16(nspw), "number IFs")
    hdu.header.set("NOPCAL", Int16(0), "number of polarization calibration values")
    hdu.header.set("POLTYPE", "VLBI", "type of polarization calibration")
    hdu.header.set("FREQID", Int16(1), "frequency set up number")
    #
    hdu.header.set("TTYPE1", nothing, "antenna name")
    hdu.header.set("TTYPE2", nothing, "antenna station coordinates")
    hdu.header.set("TTYPE3", nothing, "orbital parameters")
    hdu.header.set("TTYPE4", nothing, "antenna number")
    hdu.header.set("TTYPE5", nothing, "mount type")
    hdu.header.set("TTYPE6", nothing, "axis offset")
    hdu.header.set("TTYPE7", nothing, "feed A: polarization letter")
    hdu.header.set("TTYPE8", nothing, "feed A: position angle")
    hdu.header.set("TTYPE9", nothing, "feed A: calibration parameters")
    hdu.header.set("TTYPE10", nothing, "feed B: polarization letter")
    hdu.header.set("TTYPE11", nothing, "feed B: position angle")
    hdu.header.set("TTYPE12", nothing, "feed B: calibration parameters")

    return hdu
end


function uvdataset2fqtab(uvd::UVDataSet; ex=ThreadedEx())
    # load pyfits
    copy!(pyfits, pyimport_conda("astropy.io.fits", "astropy"))
    copy!(numpy, pyimport_conda("numpy", "numpy"))

    # quick shortcuts to datasets
    blds = uvd[:baseline]

    # number of spectral
    nch, nspw, _, _ = size(blds)

    # columns
    c1 = pyfits.Column(
        name="FRQSEL", format="1J", unit=" ",
        array=Int32.([1])
    )
    c2 = pyfits.Column(
        name="IF FREQ", format=format("{}D", nspw), unit="HZ",
        array=Float64.(reshape(blds[:νspw].data, 1, nspw) .- minimum(blds[:ν]))
    )
    c3 = pyfits.Column(
        name="CH WIDTH", format=format("{}E", nspw), unit="HZ",
        array=Float32.(reshape(blds[:Δνch].data, 1, nspw))
    )
    c4 = pyfits.Column(
        name="TOTAL BANDWIDTH", format=format("{}E", nspw), unit="HZ",
        array=Float32.(reshape(blds[:Δνch].data, 1, nspw))
    )
    c5 = pyfits.Column(
        name="SIDEBAND", format=format("{}J", nspw), unit=" ",
        array=Int32.(reshape(blds[:sideband].data, 1, nspw))
    )
    hdu = pyfits.BinTableHDU.from_columns(
        pyfits.ColDefs([c1, c2, c3, c4, c5])
    )

    # header for columns
    hdu.header.set("TTYPE1", nothing, "frequency setup ID number")
    hdu.header.set("TTYPE2", nothing, "frequency offset")
    hdu.header.set("TTYPE3", nothing, "spectral channel separation")
    hdu.header.set("TTYPE4", nothing, "total width of spectral wndow")
    hdu.header.set("TTYPE5", nothing, "sideband")

    # Add Headers
    hdu.header.set("EXTNAME", "AIPS FQ", "extension name")
    hdu.header.set("EXTVER", Int32(1), "extension version")
    hdu.header.set("EXTLEVEL", Int32(1), "")
    hdu.header.set("NO_IF", Int16(nspw), "number IFs")

    return hdu
end
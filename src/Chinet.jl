module Chinet

struct PlateFormat
    size :: NTuple{2,Int} #size of the plate
end

Base.size(pf::PlateFormat) = pf.size

plate96 = PlateFormat((8,12))
plate384 = PlateFormat((16,24))

#wells
struct Well
    target
    sample
end

target(w::Well) = w.target
sample(w::Well) = w.sample

struct AssayPlate
    format :: PlateFormat
    wells :: Matrix{Well}
    #ensure that the size of wells will fit in the plate format
    function AssayPlate(format::PlateFormat,wells::Matrix{Well})
        @assert all(size(wells) .<= size(PlateFormat)) "too many wells to fit in the specified plate format"
        new(format,wells)
    end
end

#try to lay out a plate in the least annoying format possible
function buildassay(samples::Vector,controls::Vector,targets::Vector;replicates=3,formats=[plate96,plate384])
    push!(samples,nothing) #also need no template control
    numsamples = length(samples)
    #we want all the plates to be the same format, but we want the minimum number of plates. We prefer plateformats
    #in the order provided in `formats`
    nplates = 1
    while true
        #keep trying to squeeze into `nplates` plates until everything fits
    end
end

end # module Chinet

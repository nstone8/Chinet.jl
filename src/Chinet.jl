module Chinet

using Counters
import FileIO

export plate96, plate384, buildassay, target, sample, format, wells
export mastermix, sampledilutions

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

Base.show(io::IO,w::Well) = print(io,"$(w.sample) | $(w.target)")

target(w::Well) = w.target
sample(w::Well) = w.sample

struct AssayPlate
    format :: PlateFormat
    wells :: Matrix{Well}
    #ensure that the size of wells will fit in the plate format
    function AssayPlate(format::PlateFormat,wells::Matrix{Well})
        @assert all(size(wells) .<= size(format)) "too many wells to fit in the specified plate format"
        new(format,wells)
    end
end

format(ap::AssayPlate) = ap.format
wells(ap::AssayPlate) = ap.wells

"""
```julia
subdivide(nplates,entries)
```
Subdivide the array `entries` into `nplates` arrays of as close to equal length as possible
"""
function subdivide(nplates::Int,entries)
    divided=[]
    curpos=0
    for remainingplates in nplates:-1:1
        numtopop=round(Int,(length(entries) - curpos)/remainingplates)
        push!(divided,entries[curpos+1:curpos+numtopop])
        curpos+=numtopop
    end
    return divided
end

#try to lay out a plate in the least annoying format possible
"""
```julia
buildassay(samples,targets;replicates=3,formats=[plate96,plate384])
```
Lay out qpcr plates. Use the smallest possible number of plates, preferring plate formats in the order given in `formats`
"""
function buildassay(samples::Vector,targets::Vector;replicates=3,formats=[plate96,plate384])
    #we want all the plates to be the same format, but we want the minimum number of plates. We prefer plateformats
    #in the order provided in `formats`
    nplates = 1
    while true
        #keep trying to squeeze into `nplates` plates until everything fits
        #first subdivide our samples across plates
        samplesbyplate=subdivide(nplates,samples)
        #build the assayplates
        wellsvec=map(samplesbyplate) do sbp
            allsamples=vcat(sbp,nothing) #need no template control
            repsamples=vcat((repeat([s],replicates) for s in allsamples)...)
            return [Well(targets[ti],repsamples[si]) for si in 1:length(repsamples), ti in 1:length(targets)]
        end
        #check to see if any formats work, if so build the assayplate and return, else increase the number of plates
        for f in formats
            for orientation in [:unpermuted,:permuted] #not using the loop variable, just trying to illustrate what I'm doing
                if all(
                    map(wellsvec) do wv
                        all(size(wv) .<= size(f))
                    end
                    )

                    #everything fits, return assayplates
                    toreturn = [AssayPlate(f,wv) for wv in wellsvec]
                    #return a scalar if there's only one plate
                    if length(toreturn) == 1
                        return toreturn[1]
                    else
                        return toreturn
                    end
                end
                #try permuting all of the well matrices
                wellsvec=permutedims.(wellsvec)
            end
        end
    end
end

function Base.convert(::Type{Vector{NamedTuple}},ap::AssayPlate)
    wells=permutedims(ap.wells)
    map(1:size(wells)[1]) do row
        NamedTuple((Symbol(colname) => colval) for (colname,colval) in zip(reverse(collect('A':'Z')[1:size(wells)[2]]),wells[row,:]))
    end
end

"""
```julia
mastermix(vassay,assayplates...;overage=0.1,additionalcomponents...)
```
Calculate a recipe for the master mixes required for targets in `assayplate`.
`vassay` is the volume of primers for each reaction, overage is pipetting overage
(10% default) and additional components are the additional ingredients.

# Example
```julia
mastermix(0.5,ap,overage=0.1,mastermix=5,water=3.5)
```
"""
function mastermix(vassay,aps::AssayPlate...;overage=0.1,additionalcomponents...)
    #need to calculate the recipes for assay and sample mixes
    #need to provide the volumes of potentially arbitrary ingredients in the
    #master mix (assay volume required)
    adcomps=keys(additionalcomponents) |> collect
    targetcounts = sum(aps) do ap
        #count the targets in ap
        target.(wells(ap)) |> counter
    end

    #return a Vector{NamedTuple} with one entry per target
    return map(keys(targetcounts) |> collect) do t
        (target = string(t), assay_volume = vassay*targetcounts[t]*(1+overage),
         (Symbol(adc) => additionalcomponents[adc] * targetcounts[t]*(1+overage) for adc in adcomps)...)
    end
end

"""
```julia
sampledilutions(vsample,qsample,concsamples,aps...;overage=0.1)
```
Calculate the dilutions for all the samples in the `AssayPlate`s `aps`. `vsample`
is the total volume of cDNA per reaction, `qsample` is the desired cDNA concentration,
`concsamples` is the current concentration of each cDNA sample.
"""
function sampledilutions(vsample,qsample,concsamples,aps...;overage=0.1)
    #also need concentrations of sample, current concentrations and final volumes
    #per reaction

    #get a count of the number of reactions per sample
    samplecounts = sum(aps) do ap
        #count the samples in ap
        sample.(wells(ap)) |> counter
    end
    allsamples = keys(samplecounts) |> collect
    #make sure we have an entry in concsamples for each sample, assume it is dict-like
    for s in allsamples
        if isnothing(s) #we won't require a conc for NTC
            continue
        end
        @assert s in keys(concsamples) "no concentration provided for $s"        
    end
    #check that we have no extra entries
    @assert length(keys(concsamples) |> collect) == length(filter(allsamples) do as
                                                               !isnothing(as) #don't want to worry about NTC
                                                           end) "unused sample concentrations provided"

    #calculate the 'desired' sample concentration
    finalconc = qsample/vsample
    for (s,c) in concsamples
        @assert c > finalconc "the source concentration of $s is too low to achieve the provided quantity"
    end

    #return a Vector with a NamedTuple for each sample giving the amount of cDNA and water to use
    return map(allsamples) do s
        #scale for number of reactions
        rxnscale=samplecounts[s]*(overage+1)
        if isnothing(s)
            return (sample = string(s), vol_cDNA = zero(rxnscale)*vsample, vol_water = vsample*rxnscale)
        end
        #per reaction
        vcdna = qsample / concsamples[s]
        vwater = vsample - vcdna
        return (sample = string(s), vol_cDNA = vcdna*rxnscale, vol_water = vwater*rxnscale)
    end
end

#add a method to FileIO.save for AssayPlates
function FileIO.save(filename,ap::AssayPlate)
    FileIO.save(filename,convert(Vector{NamedTuple},ap))
end

end # module Chinet

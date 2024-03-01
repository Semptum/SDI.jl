module Utils
using Statistics
using ImageFiltering
using Peaks
using Unitful

import ..Images: open_image_series

function auto_bright_dark(image, b = 0.75, d = 0.5)
    return image .> quantile(image[:],b), image .< quantile(image[:],d)
end

function SDI_contrast(image, bright, dark)
    B,D = sum(image[bright]), sum(image[dark])
    return (B-D)/(B+D)
end

function contrast_sequence(image_sequence::Array{T, 3}, bright, dark) where {T}
    N = size(image_sequence)[3]
    contrasts = Vector{Float64}(undef, N)
    for i in 1:N
        contrasts[i] = SDI_contrast(image_sequence[:,:,i], bright, dark)
    end
    return contrasts
end

function smooth(values::Vector; size = 3)
    ker = ImageFiltering.Kernel.gaussian((3,))
    imfilter(values, ker)
end

function extrema(values::Vector, isolation = 15; index0::Int = -1)
    maxima = Peaks.argmaxima(values, isolation)
    minima = Peaks.argminima(values, 20)
    extrema = sort!(vcat(maxima,minima))
    if index0>0 # If we actually get an index
        insert!(extrema,1,index0)
    end
    return extrema
end


function angles(n)
    return collect(0:n-1).*pi
end

function index_to_lambda(lambda_n_point::Dict)
    i2l = zeros(Int64, sum(values(lambda_n_point)))
    lambda = sorted(collect(keys(lambda_n_point)))
    i = 1
    k = 1
    for l in lambda
        i2l[i:i+lambda_n_point[l]] = k
        k+=1
        i+=lambda_n_point[l]
    end
    return i2l
end

function xy_from_size(shape)
    x = collect(1:shape[2])
    y = collect(1:shape[1])
    xy = vcat(x,y)
    return xy
end

function lambdas_to_i(lambdas::Vector{T}) where {T}
    V = unique(lambdas)
    D = Dict(zip(sort(V),1:length(V)))
    result = Vector{Int64}(undef,0)
    for i in lambdas
        append!(result, D[i])
    end
    return result
end

function concat_values(values::Vector{T}, Ind, Len) where {T}
    cVal = Vector{T}(undef,0)
    cLen = cumsum(Len)
    for i in sort(unique(Ind))
        if length(cVal)==0
            append!(cVal,values[Ind.==i])
        else
            append!(cVal,cLen[i-1].+values[Ind.==i])
        end
    end
    return cVal
end

function load_data(data_prefix, data_range, param; n_cutoff = 10)
    inversions = Vector{Float64}(undef,0) # Times of inversions in picoseconds, concatenated
    angles = Vector{Float64}(undef,0) # Phase at moments of inversion (increments of pi), concatenated
    contrasts = Vector{Float64}(undef,0) # Actual contrast vectors, concatenated
    integer_inversions = Vector{Int64}(undef,0) # Index at which occurs the inversions, concatenated
    raw_contrasts = Vector{Float64}(undef,0) # Actual non-smoothed contrast vectors, concatenated
    
    lambdas = Vector{Float64}(undef,0) # Wavelengths, concatenated
    indices = Vector{Int64}(undef,0) # Number of acquisition, incrementing from 1, concatenated
    lengths = Vector{Int64}(undef,0) # Length of each contrast vector
    steps = Vector{Float64}(undef,0) # step sizes

    k = 1
    for scan in data_range
        shot_range,lambda,folder,step, delay0, Bright, Dark = scan # Loads the parameters that were manually indicated
        if delay0>0 # If we started the scan with negative delay
            index0 = Int(round(delay0/step))
        else
            index0 = -1
        end
        img = open_image_series(joinpath(data_prefix,folder), shot_range) # Loads all the images into a big array
        B_auto,D_auto = auto_bright_dark(img[:,:,1])
        if isnothing(Bright)
            Bright = B_auto
        end
        if isnothing(Dark)
            Dark = D_auto
        end
        raw_contrast = contrast_sequence(img,Bright,Dark) # Calculate contrast curve
        contrast = smooth(raw_contrast) # Smooth out the contrast curve
        ex_int = extrema(contrast; index0 = index0) # Find the indices of the extrema of the contrast
        ex_int = ex_int[1:min(n_cutoff, length(ex_int))] # If necessary, do a cutoff
        ex = ex_int * step # Moments of time of the extrema

        append!(inversions,ex)
        append!(angles, pi.* (0:length(ex)-1))
        append!(lambdas, lambda*ones(length(ex)))
        append!(indices,k*ones(Int64,length(ex)))
        append!(steps,step*ones(Int64,length(ex)))
        append!(contrasts, contrast)
        append!(raw_contrasts, raw_contrast)
        append!(lengths, length(contrast))
        append!(integer_inversions, ex_int)
        k+=1
    end
    cInv = concat_values(integer_inversions, indices, lengths)

    return Dict(
        "Inv"=>inversions, 
        "Ang"=>angles,
        "Lam"=>lambdas,
        "Ind"=>indices,
        "Con"=>contrasts,
        "Len"=>lengths,
        "iInv"=>integer_inversions,
        "Ste"=>steps,
        "rCon"=>raw_contrasts,
        "cInv"=>cInv, # Concatenated inversion moments (for plots)
    )
end


n_c(λ::Quantity) = 1.1e21u"cm^-3"*(1u"µm"/λ)^2 |> u"cm^-3"
n_c(λ::Number) = 1.1e21u"cm^-3"*(1000/λ)^2 # Assuming it's in nanometers

end
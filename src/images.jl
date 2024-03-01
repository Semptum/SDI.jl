module Images

export open_image

using TiffImages
using FixedPointNumbers

function open_image(filename)
    Float64.(reinterpret(FixedPointNumbers.N0f16,tiff_open(filename)))
end

function tiff_open(filename)
    TiffImages.load(filename)
end

function open_image_series(folder_path, range; base_name = "tir", extension = ".tiff")
    img_name = string(base_name, first(range), extension)
    filename = joinpath(folder_path, img_name)
    image_size = size(TiffImages.load(filename, mmap = true))
    img = zeros(Float64, (image_size[1], image_size[2], length(range)))
    k=1
    for i in range
        img_name = string(base_name, i, extension)
        filename = joinpath(folder_path, img_name)
        initial_gc_state = GC.enable(false)
        @views img[:,:, k] = open_image(filename)
        GC.enable(initial_gc_state)
        k+=1
    end
    return img
end


end
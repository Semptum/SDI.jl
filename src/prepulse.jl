module PrePulse

import ..Images: open_image
import ..Utils: xy_from_size

export fluence_airy, fluence_gauss, fluence_image, fit_airy, fit_gauss, fit_image

using TiffImages
using FixedPointNumbers
using LsqFit
using ImageFiltering
using Bessels
using Statistics


function load_prepulse(filename)
    open_image(filename)
end

function twoD_Airy(xy, p; size = (1024,1024))
    x = xy[1:size[1]]
    y = xy[size[1]+1:end]
    I0, xo, yo, sigma_x, sigma_y, offset = p
    R = sqrt.(((x'.-xo)/sigma_x).^2 .+ ((y.-yo)/sigma_y).^2)
    g = I0 .* (2 * besselj1.(R)./(R.+1e-5)).^2 .+ offset
    return g[:]
end

function twoD_Gaussian(xy, p; size = (1024,1024))
    amplitude, xo, yo, sigma_x, sigma_y, theta, offset = p
    a = (cos(theta)^2)/(2*sigma_x^2) + (sin(theta)^2)/(2*sigma_y^2)
    b = -(sin(2*theta))/(4*sigma_x^2) + (sin(2*theta))/(4*sigma_y^2)
    c = (sin(theta)^2)/(2*sigma_x^2) + (cos(theta)^2)/(2*sigma_y^2)

    x = xy[1:size[1]]
    y = xy[size[1]+1:end]
    #+ 2 .* b .* (x' .- xo) .* (y .- yo)
    exp_ins = - (a.*((x' .- xo).^2)  .+ c * ((y .- yo).^2))
    g = offset .+ amplitude .* exp.(exp_ins )
    return g[:]
end


function fit_gauss(image)
    H,W = size(image)
    xy = xy_from_size((H,W))
    p0 = Float64.([maximum(image), W/2, H/2, W/4, H/4, 0, 0])
    fit_gauss = LsqFit.curve_fit(twoD_Gaussian, xy, image[:], p0)
    return fit_gauss.param
end

function fit_airy(image)
    H,W = size(image)
    xy = xy_from_size((H,W))
    p0 = Float64.([maximum(image), W/2, H/2, W/4, H/4, 0])
    fit_gauss = LsqFit.curve_fit(twoD_Airy, xy, image[:], p0)
    return fit_gauss.param
end

function fit_image(image)
    center = find_center(image)
    shape = size(image)
    x = collect(1:shape[2])
    y = collect(1:shape[1])
    R = sqrt.(((x'.-center[2])).^2 .+ ((y.-center[1])).^2)
    return image.-mean(image[R.>0.9*maximum(R)]), R
end

function find_center(image)
    ker = ImageFiltering.Kernel.gaussian((50,50))
    argmax(imfilter(image,ker))
end

function fluence_airy(param, distance_px)
    I0, xo, yo, sigma_x, sigma_y, offset = param
    R = distance_px/sqrt(sigma_x*sigma_y)
    g = (2 * besselj1.(R)./(R.+1e-5)).^2
end

function fluence_gauss(param, distance_px)
    amplitude, xo, yo, sigma_x, sigma_y, theta, offset = param
    a = (cos(theta)^2)/(2*sigma_x^2) + (sin(theta)^2)/(2*sigma_y^2)
    c = (sin(theta)^2)/(2*sigma_x^2) + (cos(theta)^2)/(2*sigma_y^2)
    exp(- (sqrt(a*c)*distance_px^2))
end

function fluence_image(param, distance_px)
    image, R = param
    mean(image[distance_px-10 .< R .< distance_px+10])/mean(image[R.<20])
end

function scaling_factor(distance_sat; focal = 51.25, spacing = 4, wavelength = 800)
    real_distance_sat = focal*wavelength/spacing
    distance_sat/real_distance_sat # scaling factor from nanometers to pixels
end

function satellite_distance_factor(scale; focal = 51.25, spacing = 4)
    real_distance_sat = focal/spacing
    real_distance_sat*scale
end

end
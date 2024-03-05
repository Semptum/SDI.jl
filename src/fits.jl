module Fits
using LsqFit
using Unitful


import ..Utils: lambdas_to_i, n_c
import ..PrePulse: satellite_distance_factor, scaling_factor

"""
We suppose an isothermal model and fix a given sound speed. 
We let the ratio of expansion at the positions of the satellites be free.

The parameters are the expansion factors (between 0 and 1) for each wavelength
"""
struct FixedSpeed_FreeRatios end

"""
We suppose an isothermal model and a sound speed. 
The satellites expand with an expansion factor of sqrt(F)


The parameter is the sound speed.
"""
struct FreeSpeed_FixedFluenceRatios end


"""
We suppose an isothermal model and a sound speed. 
We let the ratio of expansion at the positions of the satellites be free.

The parameters are the expansion factors (between 0 and 1) for each wavelength as well as the sound speed.
"""
struct FreeSpeed_FreeRatios end


"""
We suppose an isothermal model and a sound speed. 
The satellites do not expand

The parameter is the sound speed.
"""
struct FreeSpeed_NullRatios end


"""
We suppose an isothermal model and a sound speed. 

The parameter is the sound speed as well as an activation fluence
"""
struct FreeSpeed_OffsetFluenceRatios end


"""
We suppose an isothermal model and a sound speed. 

The parameter is the activation fluence
"""
struct FixedSpeed_OffsetFluenceRatios end


"""
We suppose that the satellite expansion speed is proportional to the square root of the fluence

The parameters are the speeds at different wavelengths
"""
struct FreeSpeeds_NullFluenceRatios end

"""
We suppose that the satellite expansion speed is proportional to the square root of the fluence

The parameters are the speeds at different wavelengths
"""
struct FreeSpeeds_FixedFluenceRatios end


"""
We suppose that the satellite expansion speed is constant for all wavelengths.

The parameters are the speeds at different wavelengths as well as that single ratio.
"""
struct FreeSpeeds_ConstantRatios end


"""
We suppose the isothermal model. We suppose that the satellite expansion speed is constant for all wavelengths.

The parameters are the speeds at different wavelengths
"""
struct FreeSpeed_ConstantRatios end

"""
We suppose the isothermal model. We suppose that the satellite expansion speed is proportional to the fluence but starting from a certain wavelength it becomes 0.

The parameters are the speeds at different wavelengths
"""
struct FreeSpeed_LambdaCutoffFluenceRatios end

"""
We suppose that the satellite expansion speed is proportional to the square root of the fluence

The parameters are the speeds at different wavelengths but this time with a polynomial fit. First linear terms, then 2nd degree
"""
struct FreePolynomialSpeeds_FixedFluenceRatios end #Not yet re-implemented, not very useful


function fit_speed(inversions, angles, lambdas, param, ::FixedSpeed_FreeRatios; c_s = 22.0)
    θᵢ = param["incidence_angle"]
    #p_prepulse = param["prepulse_fit_param"]
    #fit_prepulse = param["prepulse_fit_method"]
    #Fluences = fit_prepulse.(Ref(p_prepulse),lambdas)
    n_e_0 = param["n_e_0"]*n_c(800)
    lambdas2i = lambdas_to_i(lambdas)

    Geom = 4*π*cos(θᵢ)./(lambdas*1u"nm")
    F_iso = log.(n_e_0 ./(exp(1)*n_c.(lambdas)*cos(θᵢ)^2))

    model(t, exp_ratio) = Geom .* F_iso .* (1 .- exp_ratio[lambdas2i]) .* c_s*1u"nm/ps" .*t*1u"ps"

    n = maximum(lambdas2i)
    upper = ones(n)
    lower = zeros(n)

    fit = curve_fit(model, inversions, angles, rand(n), lower=lower, upper=upper)
    
    return fit.param, model(inversions,fit.param), model, fit
end

function fit_speed(inversions, angles, lambdas, param, ::FreeSpeed_FreeRatios)
    θᵢ = param["incidence_angle"]
    #p_prepulse = param["prepulse_fit_param"]
    #fit_prepulse = param["prepulse_fit_method"]
    #Fluences = fit_prepulse.(Ref(p_prepulse),lambdas)
    n_e_0 = param["n_e_0"]*n_c(800)
    lambdas2i = lambdas_to_i(lambdas)
    n = maximum(lambdas2i)

    Geom = 4*π*cos(θᵢ)./(lambdas*1u"nm")
    F_iso = log.(n_e_0 ./(exp(1)*n_c.(lambdas)*cos(θᵢ)^2))

    function model(t, p)
        exp_ratio = p[1:n]
        c_s = p[n+1]
        Geom .* F_iso .* (1 .- exp_ratio[lambdas2i]) .*c_s *1u"nm/ps" .*t*1u"ps"
    end

    upper = vcat(ones(n),300)
    lower = zeros(n+1)

    fit = curve_fit(model, inversions, angles, rand(n+1), lower=lower, upper=upper)
    
    return fit.param, model(inversions,fit.param), model, fit
end

function fit_speed(inversions, angles, lambdas, param, ::FreeSpeed_NullRatios)
    θᵢ = param["incidence_angle"]
    #p_prepulse = param["prepulse_fit_param"]
    #fit_prepulse = param["prepulse_fit_method"]
    #Fluences = fit_prepulse.(Ref(p_prepulse),lambdas)
    n_e_0 = param["n_e_0"]*n_c(800)
    lambdas2i = lambdas_to_i(lambdas)
    n = maximum(lambdas2i)

    Geom = 4*π*cos(θᵢ)./(lambdas*1u"nm")
    F_iso = log.(n_e_0 ./(exp(1)*n_c.(lambdas)*cos(θᵢ)^2))

    function model(t, c_s)
        Geom .* F_iso .*c_s *1u"nm/ps" .* t * 1u"ps"
    end

    upper = [300.0]
    lower = [0.0]

    fit = curve_fit(model, inversions, angles, rand(1), lower=lower, upper=upper)
    
    return fit.param, model(inversions,fit.param), model, fit
end


function fit_speed(inversions, angles, lambdas, param, ::FreeSpeed_ConstantRatios)
    θᵢ = param["incidence_angle"]
    #p_prepulse = param["prepulse_fit_param"]
    #fit_prepulse = param["prepulse_fit_method"]
    #Fluences = fit_prepulse.(Ref(p_prepulse),lambdas)
    n_e_0 = param["n_e_0"]*n_c(800)
    lambdas2i = lambdas_to_i(lambdas)
    n = maximum(lambdas2i)

    Geom = 4*π*cos(θᵢ)./(lambdas*1u"nm")
    F_iso = log.(n_e_0 ./(exp(1)*n_c.(lambdas)*cos(θᵢ)^2))

    function model(t, p)
        c_s, sat_ratio = p
        Geom .* F_iso .* (1 - sat_ratio) .*c_s *1u"nm/ps" .* t * 1u"ps"
    end

    upper = [300.0, 1.0]
    lower = [0.0, 0.0]

    fit = curve_fit(model, inversions, angles, rand(2), lower=lower, upper=upper)
    
    return fit.param, model(inversions,fit.param), model, fit
end



function fit_speed(inversions, angles, lambdas, param, ::FreeSpeed_OffsetFluenceRatios)
    θᵢ = param["incidence_angle"]
    p_prepulse = param["prepulse_fit_param"]
    fluence_function = param["prepulse_fit_method"]
    d_sat_800 = param["distance_sat"]
    d_sat_px_per_nm = satellite_distance_factor(scaling_factor(d_sat_800))
    Fluences = fluence_function.(Ref(p_prepulse),lambdas*d_sat_px_per_nm)

    n_e_0 = param["n_e_0"]*n_c(800)
    lambdas2i = lambdas_to_i(lambdas)
    n = maximum(lambdas2i)

    Geom = 4*π*cos(θᵢ)./(lambdas*1u"nm")
    F_iso = log.(n_e_0 ./(exp(1)*n_c.(lambdas)*cos(θᵢ)^2))

    function model(t, p)
        F0 = p[1]
        c_s = p[2]
        Geom .* F_iso .* (1 .- sqrt.(Fluences.+F0)) .*c_s *1u"nm/ps" .* t * 1u"ps"
    end

    upper = [1,300.0]
    lower = [0.0, 0.0]

    fit = curve_fit(model, inversions, angles, rand(2), lower=lower, upper=upper)
    
    return fit.param, model(inversions,fit.param), model, fit
end

function fit_speed(inversions, angles, lambdas, param, ::FreeSpeed_FixedFluenceRatios)
    θᵢ = param["incidence_angle"]
    p_prepulse = param["prepulse_fit_param"]
    fluence_function = param["prepulse_fit_method"]
    d_sat_800 = param["distance_sat"]
    d_sat_px_per_nm = satellite_distance_factor(scaling_factor(d_sat_800))
    Fluences = fluence_function.(Ref(p_prepulse),lambdas*d_sat_px_per_nm)

    n_e_0 = param["n_e_0"]*n_c(800)
    lambdas2i = lambdas_to_i(lambdas)
    n = maximum(lambdas2i)

    Geom = 4*π*cos(θᵢ)./(lambdas*1u"nm")
    F_iso = log.(n_e_0 ./(exp(1)*n_c.(lambdas)*cos(θᵢ)^2))

    function model(t, c_s)
        Geom .* F_iso .* (1 .- sqrt.(Fluences)) .*c_s *1u"nm/ps" .* t * 1u"ps"
    end

    upper = [300.0]
    lower = [0.0]

    fit = curve_fit(model, inversions, angles, rand(1), lower=lower, upper=upper)
    
    return fit.param, model(inversions,fit.param), model, fit
end

function fit_speed(inversions, angles, lambdas, param, ::FixedSpeed_OffsetFluenceRatios; c_s = 22)
    θᵢ = param["incidence_angle"]
    p_prepulse = param["prepulse_fit_param"]
    fluence_function = param["prepulse_fit_method"]
    d_sat_800 = param["distance_sat"]
    d_sat_px_per_nm = satellite_distance_factor(scaling_factor(d_sat_800))
    Fluences = fluence_function.(Ref(p_prepulse),lambdas*d_sat_px_per_nm)

    n_e_0 = param["n_e_0"]*n_c(800)
    lambdas2i = lambdas_to_i(lambdas)
    n = maximum(lambdas2i)

    Geom = 4*π*cos(θᵢ)./(lambdas*1u"nm")
    F_iso = log.(n_e_0 ./(exp(1)*n_c.(lambdas)*cos(θᵢ)^2))

    function model(t, F0)
        Geom .* F_iso .* (1 .- sqrt.(Fluences.+F0)) .*c_s *1u"nm/ps" .* t * 1u"ps"
    end

    upper = [1.0]
    lower = [0.0]

    fit = curve_fit(model, inversions, angles, rand(1), lower=lower, upper=upper)
    
    return fit.param, model(inversions,fit.param), model, fit
end

function fit_speed(inversions, angles, lambdas, param, ::FreeSpeeds_FixedFluenceRatios)
    θᵢ = param["incidence_angle"]
    p_prepulse = param["prepulse_fit_param"]
    fluence_function = param["prepulse_fit_method"]
    d_sat_800 = param["distance_sat"]
    d_sat_px_per_nm = satellite_distance_factor(scaling_factor(d_sat_800))
    Fluences = fluence_function.(Ref(p_prepulse),lambdas*d_sat_px_per_nm)

    n_e_0 = param["n_e_0"]*n_c(800)
    lambdas2i = lambdas_to_i(lambdas)
    n = maximum(lambdas2i)

    Geom = 4*π*cos(θᵢ)./(lambdas*1u"nm")

    function model(t, speeds)
        Geom .* (1 .- sqrt.(Fluences)) .*speeds[lambdas2i] *1u"nm/ps" .* t * 1u"ps"
    end

    upper = 300*ones(n)
    lower = zeros(n)

    fit = curve_fit(model, inversions, angles, rand(n), lower=lower, upper=upper)
    
    return fit.param, model(inversions,fit.param), model, fit
end

function fit_speed(inversions, angles, lambdas, param, ::FreeSpeeds_NullFluenceRatios)
    θᵢ = param["incidence_angle"]

    n_e_0 = param["n_e_0"]*n_c(800)
    lambdas2i = lambdas_to_i(lambdas)
    n = maximum(lambdas2i)

    Geom = 4*π*cos(θᵢ)./(lambdas*1u"nm")

    function model(t, speeds)
        Geom .*speeds[lambdas2i] *1u"nm/ps" .* t * 1u"ps"
    end

    upper = 300*ones(n)
    lower = zeros(n)

    fit = curve_fit(model, inversions, angles, rand(n), lower=lower, upper=upper)
    
    return fit.param, model(inversions,fit.param), model, fit
end


function fit_speed(inversions, angles, lambdas, param, ::FreeSpeeds_ConstantRatios)
    θᵢ = param["incidence_angle"]

    n_e_0 = param["n_e_0"]*n_c(800)
    lambdas2i = lambdas_to_i(lambdas)
    n = maximum(lambdas2i)

    Geom = 4*π*cos(θᵢ)./(lambdas*1u"nm")

    function model(t, p)
        speeds = p[1:end-1]
        sat_ratio = p[end]
        Geom .* (1-sat_ratio) .*speeds[lambdas2i] *1u"nm/ps" .* t * 1u"ps"
    end

    upper = vcat(300*ones(n), 1.0)
    lower = zeros(n+1)

    fit = curve_fit(model, inversions, angles, rand(n+1), lower=lower, upper=upper)
    
    return fit.param, model(inversions,fit.param), model, fit
end



function fit_speed(inversions, angles, lambdas, param, ::FreeSpeed_LambdaCutoffFluenceRatios; cutoff = 740)
    θᵢ = param["incidence_angle"]
    p_prepulse = param["prepulse_fit_param"]
    fluence_function = param["prepulse_fit_method"]
    d_sat_800 = param["distance_sat"]
    d_sat_px_per_nm = satellite_distance_factor(scaling_factor(d_sat_800))
    Fluences = fluence_function.(Ref(p_prepulse),lambdas*d_sat_px_per_nm)

    n_e_0 = param["n_e_0"]*n_c(800)
    lambdas2i = lambdas_to_i(lambdas)

    Geom = 4*π*cos(θᵢ)./(lambdas*1u"nm")
    F_iso = log.(n_e_0 ./(exp(1)*n_c.(lambdas)*cos(θᵢ)^2))

    function model(t, c_s)
        Geom .* F_iso .* (1 .- sqrt.(Fluences).*(lambdas.<cutoff)) .*c_s *1u"nm/ps" .* t * 1u"ps"
    end

    upper = [300.0]
    lower = [0.0]

    fit = curve_fit(model, inversions, angles, rand(1), lower=lower, upper=upper)
    
    return fit.param, model(inversions,fit.param), model, fit
end


end
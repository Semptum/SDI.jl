module Plotting

import ..PrePulse: fit_image, 
                    fit_gauss, 
                    fit_airy, 
                    fluence_image,
                    fluence_gauss,
                    fluence_airy,
                    scaling_factor,
                    satellite_distance_factor


import ..Fits:  FixedSpeed_FreeRatios,
                FreeSpeed_FixedFluenceRatios,
                FreeSpeed_FreeRatios,
                FreeSpeed_NullRatios,
                FreeSpeed_OffsetFluenceRatios,
                FixedSpeed_OffsetFluenceRatios,
                FreeSpeeds_NullFluenceRatios,
                FreeSpeeds_FixedFluenceRatios,
                FreeSpeeds_ConstantRatios,
                FreeSpeed_ConstantRatios,
                FreePolynomialSpeeds_FixedFluenceRatios,
                fit_speed

using Plots
using LsqFit
                
function plot_prepulse_fit(image, param, figure_folder; filename = "prepulse_fit.pdf")
    p_im = fit_image(image)
    p_g = fit_gauss(image)
    p_a = fit_airy(image)
    lambdas = param["lambdas"]
    scaling_sat_px_per_nm = satellite_distance_factor(scaling_factor(param["distance_sat"]))
    plot(20:5:1000,
        sqrt.(
                fluence_image.(
                    Ref(p_im),
                    (20:5:1000)*scaling_sat_px_per_nm
                )
            ),
    label = "From image")

    plot!(20:5:1000,
    sqrt.(
            fluence_gauss.(
                Ref(p_g),
                (20:5:1000)*scaling_sat_px_per_nm
            )
        ),
    label = "Gaussian fit")

    plot!(20:5:1000,
    sqrt.(
            fluence_airy.(
                Ref(p_a),
                (20:5:1000)*scaling_sat_px_per_nm
            )
        ),
    label = "Airy fit")
    vline!(lambdas, label = "Measurement wavelengths")
    xlabel!("Wavelength [nm]")
    ylabel!("Fluence [relative]")
    title!("Fluence seen by satellites of wavelength λ")
    savefig(joinpath(figure_folder, filename))
end

function plot_Extrema(Loaded_Data, figure_folder; filename = "Extrema.pdf")
    cInv, rCon, Con, Lam = Loaded_Data["cInv"],Loaded_Data["rCon"], Loaded_Data["Con"], Loaded_Data["Lam"]
    plot(rCon, label="Unsmoothed contrast")
    #scatter!(cInv,Con[cInv], label="Extrema (from smoothed contrast)")
    k = 1
    for l in sort(unique(Lam))
        rng = k:k+sum(Lam.==l)-1
        scatter!(cInv[rng],Con[cInv][Lam.==l], color = Int.(l), label = "$l nm")
        k += sum(Lam.==l)
    end
    xlabel!("Time/Acquisitions (concatendated)")
    ylabel!("Contrast")
    title!("Contrast plot and extracted extrema")
    savefig(joinpath(figure_folder,filename))
end

function plot_FixedSpeed_FreeRatios(Loaded_Data, figure_folder, param; filename = "FixedSpeed_FreeRatios.pdf", c_s = 22.0)
    Inv,Ang,Lam, Steps, cInv = Loaded_Data["Inv"],Loaded_Data["Ang"],Loaded_Data["Lam"], Loaded_Data["Ste"], Loaded_Data["cInv"]
    p,A,model,fit = fit_speed(Inv, Ang, Lam, param, 
        FixedSpeed_FreeRatios(); c_s = c_s)
    s = satellite_distance_factor(scaling_factor(param["distance_sat"]))
    conf95 = [i[2]-i[1] for i in confidence_interval(fit, 0.05)]
    plot(300:5:1000,
        sqrt.(
                fluence_image.(
                    Ref(param["prepulse_fit_param"]),
                    (300:5:1000)*s
                )
            ),
    label = "√F measured by camera")
    vline!(sort(unique(Lam)), label = "Measurement wavelengths")
    scatter!(sort(unique(Lam)), p, yerror=conf95, label = "Parameters obtained by fit, cₛ = $(round(c_s, sigdigits=2))")
    xlabel!("Wavelength [nm]")
    ylabel!("Satellite expansion speed [relative]")
    title!("Expansion ratio of satellites compared to theory")
    savefig(joinpath(figure_folder,filename))
    #
    plot(Steps.*cInv, A, label="Fit")
    k = 1
    for l in sort(unique(Lam))
        rng = k:k+sum(Lam.==l)-1
        scatter!(cInv[rng].*Steps[Lam.==l],Ang[Lam.==l], color = Int.(l), label = "$l nm")
        k += sum(Lam.==l)
    end
    xlabel!("Concatenated measurements")
    ylabel!("Angle [rad]")
    title!("Visualisation of fit")
    savefig(joinpath(figure_folder,"vis_fit_"*filename))
end

function plot_FreeSpeed_FixedFluenceRatios(Loaded_Data, figure_folder, param; filename = "FreeSpeed_FixedFluenceRatios.pdf")
    Inv,Ang,Lam, Steps, cInv = Loaded_Data["Inv"],Loaded_Data["Ang"],Loaded_Data["Lam"], Loaded_Data["Ste"], Loaded_Data["cInv"]
    p,A,model,fit = fit_speed(Inv, Ang, Lam, param, 
        FreeSpeed_FixedFluenceRatios())
    conf95 = [i[2]-i[1] for i in confidence_interval(fit, 0.05)]
    print("$(p[1]) ± $(conf95[1]/2) nm/ps")
    plot(Steps.*cInv, A, label="Fit")
    k = 1
    for l in sort(unique(Lam))
        rng = k:k+sum(Lam.==l)-1
        scatter!(cInv[rng].*Steps[Lam.==l],Ang[Lam.==l], color = Int.(l), label = "$l nm")
        k += sum(Lam.==l)
    end
    xlabel!("Concatenated measurements")
    ylabel!("Angle [rad]")
    title!("Visualisation of fit, cₛ = $(round(p[1],digits = 2)) ± $(round(conf95[1]/2, sigdigits=2)) nm/ps")
    savefig(joinpath(figure_folder,"vis_fit_"*filename))
end

function plot_FreeSpeed_FreeRatios(Loaded_Data, figure_folder, param; filename = "FreeSpeed_FreeRatios.pdf")
    Inv,Ang,Lam, Steps, cInv = Loaded_Data["Inv"],Loaded_Data["Ang"],Loaded_Data["Lam"], Loaded_Data["Ste"], Loaded_Data["cInv"]
    p,A,model,fit = fit_speed(Inv, Ang, Lam, param, 
        FreeSpeed_FreeRatios())
    s = satellite_distance_factor(scaling_factor(param["distance_sat"]))
    conf95 = [i[2]-i[1] for i in confidence_interval(fit, 0.05)]
    plot(300:5:1000,
        sqrt.(
                fluence_image.(
                    Ref(param["prepulse_fit_param"]),
                    (300:5:1000)*s
                )
            ),
    label = "√F measured by camera")
    vline!(sort(unique(Lam)), label = "Measurement wavelengths")
    scatter!(sort(unique(Lam)), p[1:end-1], yerror=conf95[1:end-1], label = "Parameters obtained by fit, cₛ = $(round(p[end], sigdigits=3))")
    xlabel!("Wavelength [nm]")
    ylabel!("Satellite expansion speed [relative]")
    title!("Expansion ratio of satellites compared to theory")
    savefig(joinpath(figure_folder,filename))
    #
    plot(Steps.*cInv, A, label="Fit")
    k = 1
    for l in sort(unique(Lam))
        rng = k:k+sum(Lam.==l)-1
        scatter!(cInv[rng].*Steps[Lam.==l],Ang[Lam.==l], color = Int.(l), label = "$l nm")
        k += sum(Lam.==l)
    end
    xlabel!("Concatenated measurements")
    ylabel!("Angle [rad]")
    title!("Visualisation of fit, cₛ = $(round(p[end],digits = 2)) ± $(round(conf95[end]/2, sigdigits=2)) nm/ps")
    savefig(joinpath(figure_folder,"vis_fit_"*filename))
end

function plot_FreeSpeed_NullRatios(Loaded_Data, figure_folder, param; filename = "FreeSpeed_NullRatios.pdf")
    Inv,Ang,Lam, Steps, cInv = Loaded_Data["Inv"],Loaded_Data["Ang"],Loaded_Data["Lam"], Loaded_Data["Ste"], Loaded_Data["cInv"]
    p,A,model,fit = fit_speed(Inv, Ang, Lam, param, 
        FreeSpeed_NullRatios())
    conf95 = [i[2]-i[1] for i in confidence_interval(fit, 0.05)]
    #
    plot(Steps.*cInv, A, label="Fit")
    k = 1
    for l in sort(unique(Lam))
        rng = k:k+sum(Lam.==l)-1
        scatter!(cInv[rng].*Steps[Lam.==l],Ang[Lam.==l], color = Int.(l), label = "$l nm")
        k += sum(Lam.==l)
    end
    xlabel!("Concatenated measurements")
    ylabel!("Angle [rad]")
    title!("Visualisation of fit, cₛ = $(round(p[end],digits = 2)) ± $(round(conf95[end]/2, sigdigits=2)) nm/ps")
    savefig(joinpath(figure_folder,"vis_fit_"*filename))
end

function plot_FreeSpeed_OffsetFluenceRatios(Loaded_Data, figure_folder, param; filename = "FreeSpeed_OffsetFluenceRatios.pdf", c_s = 22.0)
    Inv,Ang,Lam, Steps, cInv = Loaded_Data["Inv"],Loaded_Data["Ang"],Loaded_Data["Lam"], Loaded_Data["Ste"], Loaded_Data["cInv"]
    p,A,model,fit = fit_speed(Inv, Ang, Lam, param, 
        FreeSpeed_OffsetFluenceRatios())
    conf95 = [i[2]-i[1] for i in confidence_interval(fit, 0.05)]
    #
    plot(Steps.*cInv, A, label="Fit")
    k = 1
    for l in sort(unique(Lam))
        rng = k:k+sum(Lam.==l)-1
        scatter!(cInv[rng].*Steps[Lam.==l],Ang[Lam.==l], color = Int.(l), label = "$l nm")
        k += sum(Lam.==l)
    end
    xlabel!("Concatenated measurements")
    ylabel!("Angle [rad]")
    title!("Vis. of fit, cₛ = $(round(p[end],digits = 2)) ± $(round(conf95[end]/2, sigdigits=2)) nm/ps, F₀ = $(round(p[1],digits = 2)) ± $(round(conf95[1]/2, sigdigits=2))")
    savefig(joinpath(figure_folder,"vis_fit_"*filename))
end

function plot_FixedSpeed_OffsetFluenceRatios(Loaded_Data, figure_folder, param; filename = "FixedSpeed_OffsetFluenceRatios.pdf", c_s = 22.0)
    Inv,Ang,Lam, Steps, cInv = Loaded_Data["Inv"],Loaded_Data["Ang"],Loaded_Data["Lam"], Loaded_Data["Ste"], Loaded_Data["cInv"]
    p,A,model,fit = fit_speed(Inv, Ang, Lam, param, 
        FixedSpeed_OffsetFluenceRatios();c_s=c_s)
    conf95 = [i[2]-i[1] for i in confidence_interval(fit, 0.05)]
    #
    plot(Steps.*cInv, A, label="Fit")
    k = 1
    for l in sort(unique(Lam))
        rng = k:k+sum(Lam.==l)-1
        scatter!(cInv[rng].*Steps[Lam.==l],Ang[Lam.==l], color = Int.(l), label = "$l nm")
        k += sum(Lam.==l)
    end
    xlabel!("Concatenated measurements")
    ylabel!("Angle [rad]")
    title!("Vis. of fit, cₛ=$c_s nm/ps, F₀ = $(round(p[1],digits = 2)) ± $(round(conf95[1]/2, sigdigits=2))")
    savefig(joinpath(figure_folder,"vis_fit_"*filename))
end

function plot_FreeSpeeds_NullFluenceRatios(Loaded_Data, figure_folder, param; filename = "FreeSpeeds_NullFluenceRatios.pdf")
    Inv,Ang,Lam, Steps, cInv = Loaded_Data["Inv"],Loaded_Data["Ang"],Loaded_Data["Lam"], Loaded_Data["Ste"], Loaded_Data["cInv"]
    p,A,model,fit = fit_speed(Inv, Ang, Lam, param, 
        FreeSpeeds_NullFluenceRatios())
    conf95 = [i[2]-i[1] for i in confidence_interval(fit, 0.05)]
    scatter(sort(unique(Lam)), p, yerror=conf95 ,label = "Reflection point speed from fit")
    xlabel!("Wavelength [nm]")
    ylabel!("Velocity [nm/ps]")
    title!("Reflection surface speed without satellite movement")
    savefig(joinpath(figure_folder,filename))
    #
    plot(Steps.*cInv, A, label="Fit")
    k = 1
    for l in sort(unique(Lam))
        rng = k:k+sum(Lam.==l)-1
        scatter!(cInv[rng].*Steps[Lam.==l],Ang[Lam.==l], color = Int.(l), label = "$l nm")
        k += sum(Lam.==l)
    end
    xlabel!("Concatenated measurements")
    ylabel!("Angle [rad]")
    title!("Visualisation of fit")
    savefig(joinpath(figure_folder,"vis_fit_"*filename))
end

function plot_FreeSpeeds_FixedFluenceRatios(Loaded_Data, figure_folder, param; filename = "FreeSpeeds_FixedFluenceRatios.pdf")
    Inv,Ang,Lam, Steps, cInv = Loaded_Data["Inv"],Loaded_Data["Ang"],Loaded_Data["Lam"], Loaded_Data["Ste"], Loaded_Data["cInv"]
    p,A,model,fit = fit_speed(Inv, Ang, Lam, param, 
        FreeSpeeds_FixedFluenceRatios())
    conf95 = [i[2]-i[1] for i in confidence_interval(fit, 0.05)]
    scatter(sort(unique(Lam)), p, yerror=conf95 ,label = "Reflection point speed from fit")
    xlabel!("Wavelength [nm]")
    ylabel!("Velocity [nm/ps]")
    title!("Reflection surface speed with √F satellite speed")
    savefig(joinpath(figure_folder,filename))
    #
    plot(Steps.*cInv, A, label="Fit")
    k = 1
    for l in sort(unique(Lam))
        rng = k:k+sum(Lam.==l)-1
        scatter!(cInv[rng].*Steps[Lam.==l],Ang[Lam.==l], color = Int.(l), label = "$l nm")
        k += sum(Lam.==l)
    end
    xlabel!("Concatenated measurements")
    ylabel!("Angle [rad]")
    title!("Visualisation of fit")
    savefig(joinpath(figure_folder,"vis_fit_"*filename))
end

function plot_FreeSpeeds_ConstantRatios(Loaded_Data, figure_folder, param; filename = "FreeSpeeds_ConstantRatios.pdf")
    Inv,Ang,Lam, Steps, cInv = Loaded_Data["Inv"],Loaded_Data["Ang"],Loaded_Data["Lam"], Loaded_Data["Ste"], Loaded_Data["cInv"]
    p,A,model,fit = fit_speed(Inv, Ang, Lam, param, 
        FreeSpeeds_ConstantRatios())
    conf95 = [i[2]-i[1] for i in confidence_interval(fit, 0.05)]
    scatter(sort(unique(Lam)), p[1:end-1], yerror=conf95[1:end-1] ,label = "Reflection point speed from fit")
    xlabel!("Wavelength [nm]")
    ylabel!("Velocity [nm/ps]")
    title!("Surface speed, expansion ratio of $(round(p[end], sigdigits=2))±$(round(conf95[end]; sigdigits=2))")
    savefig(joinpath(figure_folder,filename))
    #
    plot(Steps.*cInv, A, label="Fit")
    k = 1
    for l in sort(unique(Lam))
        rng = k:k+sum(Lam.==l)-1
        scatter!(cInv[rng].*Steps[Lam.==l],Ang[Lam.==l], color = Int.(l), label = "$l nm")
        k += sum(Lam.==l)
    end
    xlabel!("Concatenated measurements")
    ylabel!("Angle [rad]")
    title!("Visualisation of fit")
    savefig(joinpath(figure_folder,"vis_fit_"*filename))
end

function plot_FreeSpeed_ConstantRatios(Loaded_Data, figure_folder, param; filename = "FreeSpeed_ConstantRatios.pdf")
    Inv,Ang,Lam, Steps, cInv = Loaded_Data["Inv"],Loaded_Data["Ang"],Loaded_Data["Lam"], Loaded_Data["Ste"], Loaded_Data["cInv"]
    p,A,model,fit = fit_speed(Inv, Ang, Lam, param, 
        FreeSpeed_ConstantRatios())
    conf95 = [i[2]-i[1] for i in confidence_interval(fit, 0.05)]
    #
    plot(Steps.*cInv, A, label="Fit")
    k = 1
    for l in sort(unique(Lam))
        rng = k:k+sum(Lam.==l)-1
        scatter!(cInv[rng].*Steps[Lam.==l],Ang[Lam.==l], color = Int.(l), label = "$l nm")
        k += sum(Lam.==l)
    end
    xlabel!("Concatenated measurements")
    ylabel!("Angle [rad]")
    title!("Visualisation of fit, cₛ = $(round(p[1],digits = 2)) ± $(round(conf95[1]/2, sigdigits=2)) nm/ps")
    savefig(joinpath(figure_folder,"vis_fit_"*filename))
end
end
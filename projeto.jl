#imports
using SpecialFunctions #to use sinc function
using FFTW #to use fft
using Plots
using Statistics
using Statistics  # to use corrcoef

function graph1(t0,tf,f_samp,fc)
    #time values
    t = t0:1/f_samp:tf
    
    c_t = c(fc,t)

    plot(t, c_t, xlabel="Time (s)", ylabel="Amplitude",
        title="Carrier c(t)", ylim=(-1.1, 1.1),
        legend=false)
end

function graph2(t0,tf,f_samp,fc)
    #calculate sampling
    N = round(Int, (tf - t0) * f_samp)

    #time values
    t = range(t0, stop=tf, length=N)

    #carrier fourier transform
    C_f = fftshift(fft(c(fc,t))) / N
    freq = range(-f_samp/2, stop=f_samp/2, length=N)

    #specifying the subplots frequency ranges
    freq_range_neg = (freq .>= -2e6) .& (freq .<= -1.8e6)
    freq_range_pos = (freq .>= 1.8e6) .& (freq .<= 2e6)

    #spectrum magnitude in the range of interest
    spectrum_neg = abs.(C_f[freq_range_neg])
    spectrum_pos = abs.(C_f[freq_range_pos])

    plot1 = plot(freq[freq_range_neg] ./ 1e6, spectrum_neg,
                xlabel="Frequency(MHz)", ylabel="Magnitude", 
                title="Negative Part of the Spectrum", legend=false)
    
    plot2 = plot(freq[freq_range_pos] ./ 1e6, spectrum_pos,
                xlabel="Frequency (MHz)", ylabel="Magnitude",
                title="Positive Part of the Spectrum", legend=false)

    plot(plot1, plot2, layout=(1, 2), size=(850,450))

end

function graph3(t0,tf,f_samp)
    #time values
    t = range(0, stop=tf, length=round(Int, tf * f_samp))

    # sinc argument
    x = (t .- 100e-6) .* 1e6

    #message signal
    m_t = msg(x)

    #selecting the portion of the message signal within the 90-110µs interval
    t_interval = t[(t .>= t0) .& (t .<= tf)]
    m_t_interval = m_t[(t .>= t0) .& (t .<= tf)]

    plot(t_interval .* 1e6, m_t_interval, xlabel="Time (µs)", ylabel="Amplitude", title="Message Signal in the 90-110µs Interval")

end

function graph4(t0,tf,f_samp)
    #time values
    t = range(t0, stop=tf, length=round(Int, tf * f_samp))
    
    #message argument
    x = (t .- 100e-6) .* 1e6

    #getting message spectrum
    M_f = fftshift(fft(msg(x))) / length(msg(x))
    frequencies = range(-f_samp/2, stop=f_samp/2, length=length(M_f))
    freq_range = (frequencies .>= -2e6) .& (frequencies .<= 2e6)
    spectrum = abs.(M_f[freq_range])
    
    plot(frequencies[freq_range] ./ 1e6, spectrum, xlabel="Frequency (MHz)", ylabel="Magnitude", title="Baseband Message Signal Spectrum")
    
    #calculating half-power bandwidth (-3dB)
    peak_magnitude = maximum(spectrum)
    half_power_point = peak_magnitude / √2
    indices_above_half_power = findall(x -> x >= half_power_point, spectrum)
    bandwidth_half_power = frequencies[maximum(indices_above_half_power)] - frequencies[minimum(indices_above_half_power)]
    
    #converting into MHz
    bandwidth_half_power_mhz = bandwidth_half_power / 1e6
    
    println("Largura de meia potência (-3dB): $bandwidth_half_power_mhz MHz")
    
end

function graph5(t0,tf,f_samp,fc)
    #time values
    t = range(0, stop=tf, length=round(Int, tf * f_samp))
    
    #sinc argment
    x = (range(0, stop= tf * 1e6, length=round(Int, tf * f_samp)) .- 100)
    
    #modulation of the message m(t) with the carrier c(t)
    s = msg(x) .* c(fc,t)
    
    #selecting the portion of the modulated signal within the interval
    t_interval = t[(t .>= t0) .& (t .<= tf)]
    s_t_interval = s[(t .>= t0) .& (t .<= tf)]
    
    plot(t_interval .* 1e6, s_t_interval, xlabel="Time (µs)", ylabel="Amplitude", title="Modulated Signal in the 90-110µs Interval")

end

function graph6(t0,tf,f_samp,fc)
    #time values
    t = range(t0, stop=tf, length=round(Int, tf * f_samp))

    #sinc argment
    x = (range(t0, stop= tf * 1e6, length=round(Int, tf * f_samp)) .- 100)

    #modulation of the message m(t) with the carrier c(t)
    s = msg(x) .* c(fc,t)



    #getting message spectrum
    S_f = fftshift(fft(s) / length(s))
    freq = range(-f_samp/2, stop=f_samp/2, length=length(S_f))

    # define frequenci range
    freq_range = (freq .>= -5e6) .& (freq .<= 5e6)

    spectrum = abs.(S_f[freq_range])

    plot(freq[freq_range] ./ 1e6, spectrum, xlabel="Frequency (MHz)", ylabel="Magnitude",title="Spectrum of the Modulated Message Signal")

end


function graph7(t0, tf, f_samp, fc)

    #time values
    t = range(t0, stop = tf, length=round(Int, tf * f_samp))

    #sinc argment
    x = (range(t0, stop= tf * 1e6, length=round(Int, tf * f_samp)) .- 100)

    #modulated signal
    s = msg(x) .* c(fc,t)

    # demodulated signal
    ds = s .* c(fc,t)

    # getting spectrum
    D_f = fftshift(fft(ds)/length(ds))
    freq = range(-f_samp/2, stop=f_samp/2, length=length(D_f))
    freq_range = (freq .>= -6e6 ).& (freq .<=6e6)

    spectrum = abs.(D_f[freq_range])

    plot(freq[freq_range] ./ 1e6, spectrum, xlabel="Frequency (MHz)", ylabel="Magnitude",title="Spectrum of the Demodulated Message Signal")

end

function graph8(t0, tf, f_samp, fc)
    #time values
    t = range(t0, stop = tf, length=round(Int, tf * f_samp))

    #sinc argment
    x = (range(t0, stop= tf * 1e6, length=round(Int, tf * f_samp)) .- 100)

    #modulated signal
    s = msg(x) .* c(fc,t)

    # demodulated signal
    ds = s .* c(fc,t)

    # getting spectrum
    D_f = fftshift(fft(ds)/length(ds))
    freq = range(-f_samp/2, stop=f_samp/2, length=length(D_f))
    freq_range = (freq .>= -6e6 ).& (freq .<=6e6)

    #Defining filter
    filter = (freq .>= -2e6) .& (freq .<=2e6)

    #Appling filter
    D_f_filtered = D_f .* filter

    # Spectrum filtered
    spectrum_filtered = abs.(D_f_filtered[freq_range])

    # Plotting filtered spectrum
    plot(freq[freq_range] ./ 1e6, spectrum_filtered, xlabel="Frequency (MHz)", ylabel="Magnitude",title="Filtred spectrum")

end

function graph9(t0, tf, f_samp, fc)
    #time values
    t = range(t0, stop = tf, length=round(Int, tf * f_samp))

    #sinc argment
    x = (range(t0, stop= tf * 1e6, length=round(Int, tf * f_samp)) .- 100)

    #modulated signal
    s = msg(x) .* c(fc,t)

    # demodulated signal
    ds = s .* c(fc,t)

    # getting spectrum
    D_f = fftshift(fft(ds)/length(ds))
    freq = range(-f_samp/2, stop=f_samp/2, length=length(D_f))
    freq_range = (freq .>= -6e6 ).& (freq .<=6e6)

    #Defining filter
    filter = (freq .>= -2e6) .& (freq .<=2e6)

    #Appling filter
    D_f_filtered = D_f .* filter

    # Aplying IFFT to recover signal
    recovered = real(ifft(ifftshift(D_f_filtered)) * length(D_f_filtered))

    # Comparing
    original = msg(x)

    # calculating correlation
    similarity = cor(original, recovered)
    println("Correlation coefficient: ", similarity)

    # plotting both original and recovered
    plot(t .*1e6, recovered, label="Recovered Message", xlabel="Time (µs)", ylabel="Amplitude", title="Recovered message")
    plot!(t .* 1e6, original, label="Original Message")

end

#main block
begin
    # samp freq of 50MHz
    f_samp = 50e6

    t0 = 0
    t0_3_5 = 90e-6

    tf1 = 5e-6
    tf2 = 200e-6
    tf3_5 = 110e-6
    
    t_end = 200e-6

    msg(t) = sinc.(t);
    c(fc,t) = cos.(2*π*fc.*t);

    graph3(t0_3_5,tf3_5,f_samp)
    savefig("graph3.png")

    graph4(t0,t_end,f_samp)
    savefig("graph4.png")


    for fc in [2e6,0.5e6]
        graph1(t0,tf1,f_samp,fc)
        savefig("graph1_fc_$(fc/1e6)MHz.png")
        
        graph2(t0,tf2,f_samp,fc)
        savefig("graph2_fc_$(fc/1e6)MHz.png")
        
        graph5(t0_3_5,tf3_5,f_samp,fc)
        savefig("graph5_fc_$(fc/1e6)MHz.png")

        graph6(t0, t_end, f_samp, fc)
        savefig("graph6_fc_$(fc / 1e6)MHz.png")

        graph7(t0, t_end, f_samp, fc)
        savefig("graph7_fc_$(fc / 1e6)MHz.png")

        
        graph8(t0, t_end, f_samp, fc)
        savefig("graph8_fc_$(fc / 1e6)MHz.png")

        graph9(t0, t_end, f_samp, fc)
        savefig("graph9_fc_$(fc / 1e6)MHz.png")

    end

end

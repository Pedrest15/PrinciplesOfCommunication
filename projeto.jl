#imports
using SpecialFunctions #to use sinc function
using FFTW #to use fft
using Plots

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
    t = range(0, stop=t_end, length=round(Int, t_end * f_samp))

    # sinc argument
    x = (t .- 100e-6) .* 1e6

    #message signal
    m_t = msg(x)

    #selecting the portion of the message signal within the 90-110µs interval
    t_interval = t[(t .>= t0) .& (t .<= tf)]
    m_t_interval = m_t[(t .>= t0) .& (t .<= tf)]

    plot(t_interval .* 1e6, m_t_interval, xlabel="Time (µs)", ylabel="Amplitude", title="Message Signal in the 90-110µs Interval")

end

function graph5(t0,tf,f_samp,fc)
    #time values
    t = range(0, stop=t_end, length=round(Int, t_end * f_samp))
    
    #sinc argment
    x = (range(0, stop=t_end * 1e6, length=round(Int, t_end * f_samp)) .- 100)
    
    #modulation of the message m(t) with the carrier c(t)
    s = msg(x) .* c(fc,t)
    
    #selecting the portion of the modulated signal within the interval
    t_interval = t[(t .>= t0) .& (t .<= tf)]
    s_t_interval = s[(t .>= t0) .& (t .<= tf)]
    
    plot(t_interval .* 1e6, s_t_interval, xlabel="Time (µs)", ylabel="Amplitude", title="Modulated Signal in the 90-110µs Interval")

end

#main block
begin
    # samp freq of 50MHz
    f_samp = 50e6
    t_end = 200e-6

    msg(t) = sinc.(t);
    c(fc,t) = cos.(2*π*fc.*t);

    for fc in [2e6,0.5e6]
        t0 = 0
        t0_3_5 = 90e-6

        tf1 = 5e-6
        tf2 = 200e-6
        tf3_5 = 110e-6 

        graph1(t0,tf1,f_samp,fc)
        savefig("graph1_fc_$(fc/1e6)MHz.png")
        
        graph2(t0,tf2,f_samp,fc)
        savefig("graph2_fc_$(fc/1e6)MHz.png")
        
        graph3(t0_3_5,tf3_5,f_samp)
        savefig("graph3_fc_$(fc/1e6)MHz.png")
        
        graph5(t0_3_5,tf3_5,f_samp,fc)
        savefig("graph5_fc_$(fc/1e6)MHz.png")
    end
end
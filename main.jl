# Imports
using SpecialFunctions  # to use sinc function
using FFTW  # to use fft
using Plots
using Statistics  # to use corrcoef

#  Funcoes auxiliares

# define o intervalo de tempo
function time_values(t0, tf, f_samp)
    N = round(Int, (tf - t0) * f_samp)
    t = range(t0, stop=tf, length=N)
    return t
end

# define o intrevalo de freq
function frequency_values(N, f_samp)
    freq = range(-f_samp/2, stop=f_samp/2, length=N)
    return freq
end

# fft para obter o espectro
function compute_spectrum(signal, f_samp)
    N = length(signal)
    spectrum = fftshift(fft(signal)) / N
    frequencies = frequency_values(N, f_samp)
    return spectrum, frequencies
end

# Plots

# portadora
function plot_carrier(t0, tf, f_samp, fc)
    t = time_values(t0, tf, f_samp)
    c_t = c(fc, t)
    plot(t, c_t, xlabel="Time (s)", ylabel="Amplitude",
         title="Carrier c(t)", ylim=(-1.1, 1.1), legend=false)
end

# Plota o espectro da portadora
function plot_spectrum(t0, tf, f_samp, signal_fn, fc, title_str, freq_limits=(-2e6, 2e6))
    t = time_values(t0, tf, f_samp)
    signal = signal_fn(fc, t)
    spectrum, frequencies = compute_spectrum(signal, f_samp)
    freq_range = (frequencies .>= freq_limits[1]) .& (frequencies .<= freq_limits[2])
    spectrum_mag = abs.(spectrum[freq_range])
    plot(frequencies[freq_range] ./ 1e6, spectrum_mag,
         xlabel="Frequency (MHz)", ylabel="Magnitude",
         title=title_str, legend=false)
end

# Plota a mensagem
function plot_message_signal(t0, tf, f_samp)
    t = time_values(0, tf, f_samp)
    x = (t .- 100e-6) .* 1e6
    m_t = msg(x)
    t_interval = t[(t .>= t0) .& (t .<= tf)]
    m_t_interval = m_t[(t .>= t0) .& (t .<= tf)]
    plot(t_interval .* 1e6, m_t_interval, xlabel="Time (µs)",
         ylabel="Amplitude", title="Message Signal in the 90-110µs Interval")
end

# Plota o espectro da mensagem
function plot_message_spectrum(t0, tf, f_samp)
    t = time_values(t0, tf, f_samp)
    x = (t .- 100e-6) .* 1e6
    spectrum, frequencies = compute_spectrum(msg(x), f_samp)
    freq_range = (frequencies .>= -2e6) .& (frequencies .<= 2e6)
    spectrum_mag = abs.(spectrum[freq_range])
    plot(frequencies[freq_range] ./ 1e6, spectrum_mag,
         xlabel="Frequency (MHz)", ylabel="Magnitude",
         title="Baseband Message Signal Spectrum")

    peak_magnitude = maximum(spectrum_mag)
    half_power_point = peak_magnitude / √2
    indices_above_half_power = findall(x -> x >= half_power_point, spectrum_mag)
    bandwidth_half_power = frequencies[maximum(indices_above_half_power)] - frequencies[minimum(indices_above_half_power)]
    bandwidth_half_power_mhz = bandwidth_half_power / 1e6
    println("Largura de meia potência (-3dB): $bandwidth_half_power_mhz MHz")
end

# PLota o sinal modulado
function plot_modulated_signal(t0, tf, f_samp, fc, title_str)
    t = time_values(0, tf, f_samp)
    x = (range(0, stop=tf * 1e6, length=length(t)) .- 100)
    s = msg(x) .* c(fc, t)
    t_interval = t[(t .>= t0) .& (t .<= tf)]
    s_t_interval = s[(t .>= t0) .& (t .<= tf)]
    plot(t_interval .* 1e6, s_t_interval, xlabel="Time (µs)",
         ylabel="Amplitude", title=title_str)
end

# Plota o sinal filtrado
function plot_filtered_spectrum(t0, tf, f_samp, fc)
    t = time_values(t0, tf, f_samp)
    x = (range(t0, stop=tf * 1e6, length=length(t)) .- 100)
    s = msg(x) .* c(fc, t)
    ds = s .* c(fc, t)
    spectrum, frequencies = compute_spectrum(ds, f_samp)
    freq_range = (frequencies .>= -6e6) .& (frequencies .<= 6e6)
    filter = abs.(frequencies) .<= 2e6
    #filter = (frequencies .>= -2e6) .& (frequencies .<= 2e6)
    D_f_filtered = spectrum .* filter
    spectrum_filtered = abs.(D_f_filtered[freq_range])
    plot(frequencies[freq_range] ./ 1e6, spectrum_filtered,
         xlabel="Frequency (MHz)", ylabel="Magnitude",
         title="Filtered Spectrum")
end

# Plota o sinal recuperado
function plot_recovered_message(t0, tf, f_samp, fc)
    t = time_values(t0, tf, f_samp)
    x = (range(t0, stop=tf * 1e6, length=length(t)) .- 100)
    s = msg(x) .* c(fc, t)
    ds = s .* c(fc, t)
    spectrum, frequencies = compute_spectrum(ds, f_samp)
    filter = (frequencies .>= -2e6) .& (frequencies .<= 2e6)
    D_f_filtered = spectrum .* filter
    recovered = real(ifft(ifftshift(D_f_filtered)) * length(D_f_filtered))
    original = msg(x)
    similarity = cor(original, recovered)
    println("Correlation coefficient for $fc Hz: ", similarity)
    plot(t .* 1e6, recovered, label="Recovered Message", xlabel="Time (µs)",
         ylabel="Amplitude", title="Recovered Message")
    plot!(t .* 1e6 ,original, label="Original Message")
end

# Main block
begin
    # frequencia de amostragem
    f_samp = 50e6

    # intervalos de tempo
    t0 = 0
    
    # tempo incial para as questoes 3 e 5
    t0_3_5 = 90e-6

    # tempo final questao 1
    tf1 = 5e-6

    # tempo final para questao 2
    tf2 = 200e-6

    #tempo finala para as questoes 3 e 5
    tf3_5 = 110e-6

    # tempo final para as demais
    t_end = 200e-6

    # seta mensagem como sinc
    msg(t) = sinc.(t)
    #define a portadora
    c(fc, t) = cos.(2*π*fc.*t)

    # Questao 3 - plota a mensagem
    plot_message_signal(t0_3_5, tf3_5, f_samp)
    savefig("graph3.png")

    # Questao 4 - plota o espectro da mensagem
    plot_message_spectrum(t0, t_end, f_samp)
    savefig("graph4.png")

    # variando a frequencia de operacao entre 0.5MHz e 2MHz
    for fc in [2e6, 0.5e6]
        # Questao 1 - Plota a portadora
        plot_carrier(t0, tf1, f_samp, fc)
        savefig("graph1_fc_$(fc / 1e6)MHz.png")

        # Questao 2 - Plota o espectro da portadora
        plot_spectrum(t0, tf2, f_samp, c, fc, "Carrier Spectrum", (-2e6, 2e6))
        savefig("graph2_fc_$(fc / 1e6)MHz.png")

        # Questao 5 - Plota a mensagem modulada
        plot_modulated_signal(t0_3_5, tf3_5, f_samp, fc, "Modulated Signal in the 90-110µs Interval")
        savefig("graph5_fc_$(fc / 1e6)MHz.png")

        # Questao 6 - Plota o espectro da mensagem modulada
        plot_spectrum(t0, t_end, f_samp, (fc, t) -> msg((t .- 100e-6) .* 1e6) .* c(fc, t), fc, "Spectrum of the Modulated Message Signal", (-5e6, 5e6))
        savefig("graph6_fc_$(fc / 1e6)MHz.png")

        # Questao 7 - Plota o sinal demodulado
        plot_spectrum(t0, t_end, f_samp, (fc, t) -> msg((t .- 100e-6) .* 1e6) .* c(fc, t) .* c(fc, t), fc, "Spectrum of the Demodulated Message Signal", (-6e6, 6e6))
        savefig("graph7_fc_$(fc / 1e6)MHz.png")
        
        # Questao 8 - Plota o sinal demodulado filtrado
        plot_filtered_spectrum(t0, t_end, f_samp, fc)
        savefig("graph8_fc_$(fc / 1e6)MHz.png")

        # Questão  9 - Plota o sinal recuperado
        plot_recovered_message(t0, t_end, f_samp, fc)
        savefig("graph9_fc_$(fc / 1e6)MHz.png")
    end
end

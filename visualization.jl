module visual                                    #注意现在程序的可视化结构体定义在visual模块中了

export plot_real, plot_momentum, plot_HHG, plot_Ground, visualization, generate_Animi, plot_probability    #HHG之后还可以用来读取完数据之后画图.

using ..TDSE
using PyPlot, FFTW, CSV, DataFrames

import ..TDSE.Potential_Matrix.ω_0

const T = 2pi/ω_0

                                                                         #可以在其它部分导入模块后定义两个上面的结构,一个画实空间,一个画动量空间
function plot_real(P::Parameter, Wave::wave_function)
    plot(P.sampling, abs2.(Wave.real_space),  
    color="teal")
    grid()
    ylabel("amplitude") 
    xlabel("position/a.u")
    title("Real space")
    xlim(-20,20)
    fig.shoot.snap()
end

function plot_probability(P::Parameter)
    local df = CSV.read("Time_evolution_wave_function.csv", DataFrame)

    PyPlot.axes(yscale = "log")
    plot(df.x, df.Probability, color = "black")
    grid()
    title("E_0 = 0.1, ω=0.148, t=16T")
    ylabel("Probability")
    xticks(collect(LinRange(-P.scope, P.scope, 11)))
    xlabel("X(a.u)")

end

function plot_momentum(P::Parameter, Wave::wave_function)
    plot(P.frequency_space, 2 * abs.(Wave.momentum_space)./P.N ,         #这里有改进空间, 可以把第一个变量放到P中,不用每次画图都计算一遍
    color="indigo")
    grid()
    title("Momentum space")
    xlabel("Momentum")
    ylabel("amplitude") 
end

function plot_HHG()
    local df = CSV.read("High_Harmonic_Generation.csv", DataFrame)
    local L = length(df.t)
    local xlocal = collect(0:T:16T)
    local Harmonic_order = 1/(df.t[end]-df.t[1]) .* (LinRange(0, L-1, L))
    
    subplot(2,1,1)
    plot(df.t, df.a_t, color = "maroon")
    xticks(xlocal, 1:length(xlocal))
    ylabel("a(t)")
    xlabel("cycles")
    title("The a(t) of the Wave packet")
    grid()
    subplot(2,1,2)
    semilogy(Harmonic_order, 2/L * abs2.(fft(df.a_t))
    , color = "navy")
    ylabel(L"$log_{10}(a(\omega)^2)$")
    grid()
    xticks(ω_0/(2pi) * collect(0:2:20),collect(0:2:20))
    xlabel("Hamonic order")
    xlim(0,20 * ω_0/(2pi))
    show()
end


function plot_Ground()
    local df = CSV.read("Ground_Wave_Func.csv", DataFrame)

    df.wave =@. real(parse(Complex{Float64}, df.wave))

    plot(df.x, abs2.(df.wave), color="black")
    grid()
    title("Ground state of wave function")
    xlabel("x")
    ylabel("amplitude")

end

end
module visual                                    #注意现在程序的可视化结构体定义在visual模块中了

export plot_real, plot_momentum, plot_HHG, plot_Ground, visualization, generate_Animi, plot_probability    #HHG之后还可以用来读取完数据之后画图.

using ..Mystruct
using PyCall, PyPlot, FFTW, CSV, DataFrames
Snap = pyimport("celluloid")
import ..Mystruct.Potential_Matrix.ω_0

const T = 2pi/ω_0

Base.@kwdef mutable struct visualization         #用来可视化的参数,这里因为只有在这里用到,干脆就直接定义在这里了,其它计算部分导入
    Canvas:: Figure = figure()                   #该部分会比较浪费时间了,而且对于包引用也比较方便一些,如果定义在Varible_initial里还需要在其中引入
    shoot:: PyObject = Snap.Camera(Canvas)       #上面的三个包的导入也很麻烦,这里Camera相当于把对应的函数的对象设定为Canvas
    Animi:: PyObject = shoot.animate()
end

                                                                         #可以在其它部分导入模块后定义两个上面的结构,一个画实空间,一个画动量空间

function plot_real(P::Parameter, Wave::wave_function, fig::visualization)
    plot(P.sampling, real.(Wave.real_space),  
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

function plot_momentum(P::Parameter, Wave::wave_function, fig::visualization)
    plot(P.frequency_space, 2 * abs.(Wave.momentum_space)./P.N ,         #这里有改进空间, 可以把第一个变量放到P中,不用每次画图都计算一遍
    color="indigo")
    grid()
    title("Momentum space")
    xlabel("Momentum")
    ylabel("amplitude") 
    #title("Momentum space wave function")
    #xlim(-5,5)
    fig.shoot.snap()
    
end

function generate_Animi(fig::visualization)
    fig.Animi = fig.shoot.animate()
    fig.Animi.save("Time_Wave_resolution.gif",writer="pillow")
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

    plot(df.x, df.wave, color="black")
    grid()
    title("Ground state of wave function")
    xlabel("x")
    ylabel("amplitude")

end

end
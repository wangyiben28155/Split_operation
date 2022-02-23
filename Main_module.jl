module TDSE
export wave_function, Parameter, plot_real, plot_momentum, plot_HHG, plot_Ground, plot_probability, 
    V_Matrix, T_Matrix, positive_fft!, inverse_fft!, Calculation!, visualization, ITR!            #大模块输出的部分,可以给小模块,也可以给外部调用


using Distributions
import Base.@kwdef

const L, x_num, step_t = (200.0, 20001, 13001)
const Δt = -0.05im

@kwdef mutable struct wave_function{T<:AbstractFloat}      #用来计算迭代的参数
    real_space::Vector{Complex{T}} = complex(collect(pdf(Truncated(Normal(0, 30), -Inf64, Inf64), LinRange(-L, L, x_num))))
    momentum_space::Vector{Complex{T}} = zeros(eltype(real_space), length(real_space))
    Time_dipole::Vector{T} = zeros(Float64, step_t)
    Time::Union{T,Complex{T}} = 0.0
end


@kwdef struct Parameter{T<:AbstractFloat}               #用来控制计算参数的
    N::Int64 = x_num                                    #划分的格点的总数,后面做离散傅里叶变换的时候会用得到
    scope::T = L                                        #确定波函数的计算范围为-scope到+scope
    Δx::T = 2 * scope / (N - 1)                         #波函数的离散的空间间隔
    Δt::Union{T,Complex{T}} = Δt                        #划分的时间间隔, 尝试时间迭代区间, 因为考虑到虚时演化, 所以类型设定为复数
    Step_t::Int64 = step_t
    sampling::Vector{T} = collect(LinRange(-scope, scope, N))
    frequency_space::Vector{T} = 1 / (2 * scope) * collect(LinRange(-N / 2, N / 2 - 1, N))
end

include("Numerical_Diff_DiscreteFunc.jl")
include("Potential.jl")                                        #注意这里文件插入的顺序也是有关的,后面的文件对应模块能调用前面的,而前面的
include("Resolution.jl")                                #不能调用后面的模块,所以尽量把比较基本的模块写在前面比较好
include("visualization.jl")
include("physical_quantity.jl")
include("Real_Evolution.jl")
include("virtual_evolution.jl")                                       #计算迭代的函数放到最后,因为它也是最核心的步骤


using .Potential_Matrix, .Resolution, .visual, .Calculation,            #这里using也是为了在模块的顶端可以export,对其它无影响
    .physical_quantity, .Ground_state, .Discrete_Func_diff
end

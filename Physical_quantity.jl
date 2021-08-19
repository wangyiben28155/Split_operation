module physical_quantity

export quantity_a, eigen_energy, Normalization

using ..Mystruct
import ..Mystruct.Potential_Matrix.V_0
import ..Mystruct.Potential_Matrix.E

using PyCall, FiniteDifferences
Ig = pyimport("scipy.integrate")

function Normalization(P::Parameter, Wave::wave_function)
    return sqrt(Ig.trapz(abs2.(Wave.real_space),P.sampling)) 
end

function quantity_a(P::Parameter, Wave::wave_function)
    local Operator =@. -central_fdm(5, 1)(V_0, P.sampling) + E(Wave.Time)   #这里得到的是-∂V/∂r+E算符的函数(因为这个算符不包含对波函数的微分,
                                                                    #所以直接写成函数形式)
    return Ig.trapz((@. abs2(Wave.real_space) * Operator), P.sampling)
end


function Cubic_spline(x::Vector{T},y::Vector{T}) where T<: AbstractFloat        #三次样条插值,等拿回数值分析之后再写这部分,主要是用来求本征能量的,下面函数的1/2的部分
    

end

function eigen_energy(P::Parameter, Wave::wave_function)
    local V_operator = @. V_0(P.sampling)
    local Operator_on_Wave = @. (V_operator * Wave.real_space -  1/2)

    return Ig.trapz(@. (conj(Wave.real_space) * Operator_on_Wave),P.sampling)
end

end
module physical_quantity

export quantity_a, eigen_energy, Normalization

using ..Mystruct
using ..Mystruct.Discrete_Func_diff
import ..Mystruct.Potential_Matrix.V_0
import ..Mystruct.Potential_Matrix.E


using NumericalIntegration


@inline function Normalization(P::Parameter, Wave::wave_function)
    return sqrt.(integrate(P.sampling, abs2.(Wave.real_space), SimpsonEven()))
end

function quantity_a(P::Parameter, Wave::wave_function)
    local Operator = -Five_Point(V_0.(P.sampling), dL = P.Δx) .+ E(Wave.Time)   #这里得到的是-∂V/∂r+E算符的函数(因为这个算符不包含对波函数的微分,
    #所以直接写成函数形式)
    return integrate(P.sampling, @. abs2(Wave.real_space) * Operator)
end


function eigen_energy(P::Parameter, Wave::wave_function)
    local V_operator = V_0.(P.sampling)
    local Operator_on_Wave = (V_operator .* abs.(Wave.real_space) .- (1 / 2) * Derivative_2(abs.(Wave.real_space), dL = P.Δx))

    return integrate(P.sampling, abs.(Wave.real_space) .* Operator_on_Wave, SimpsonEven())
end

end
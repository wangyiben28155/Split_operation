module Ground_state

export ITR!                              #虚时演化

using ..TDSE
using ..TDSE.Resolution
using ..TDSE.Potential_Matrix
using ..TDSE.physical_quantity

using CSV, DataFrames

function record(P::Parameter, Wave::wave_function)                                      #这个函数虽然和function_1里的函数同名但是作用域是隔离的
    local df = DataFrame()

    df.x = P.sampling
    df.wave = Wave.real_space 
    CSV.write("Ground_Wave_Func.csv", df)
end

function ITR!(P::Parameter, Wave::wave_function)
    local count::Int64 = 0
    local T_0 = T_Matrix(P)                     #而这里由于势能不含时,所以把两个矩阵都当参数传入,避免每次循环计算一遍
    local V_0 = V_0_Matrix(P)
    local Normalizer::Float64 = 0.0
    local Energy::Vector{Float64} = [0.0, 1.0]

    positive_fft!(Wave)                         #对动量空间波函数进行初始化

    if imag(P.Δt) != 0.0                        #这里确保是虚时演化

        while true

            if abs(Energy[2] - Energy[1]) >= 1e-7                      #至少先让程序运行一次,进而判断是否满足
                Energy[1] = Energy[2]
            
                resolution!(P, Wave, Momentum_operator = T_0, V_operator = V_0)       #这里对波函数进行一次更新,因为动能算符和势能算符都是不含时的
            
                Normalizer = Normalization(P, Wave)                                  #因为要重复使用,这里应该需要储存到一个变量里比较好.
                Wave.real_space /= Normalizer                                        #把实空间的波函数归一化
                Wave.momentum_space /= Normalizer                                    #因为波函数归一化后,动量的波函数和实空间不是一一对应的了

                Energy[2] = eigen_energy(P, Wave)
                #这里为了节省计算资源,对动量波函数通过矢量除法进行归一化
                count += 1
            else
                Wave.Time = -count * P.Δt
                println("The Energy is ", eigen_energy(P, Wave))
                break
            end

        end

        record(P, Wave)                         #这里结束后说明迭代得到了稳定的波函数
        println("迭代次数:$(count),虚时间$(Wave.Time)")

    else
        return @error "the Time shold be a imaginary number"
    end
end

end
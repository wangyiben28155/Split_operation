module Calculation  #这个模块可以对波函数进行迭代, 比喻的话, 这个模块相当于餐厅(根据引用流)
# 对计算的波函数做迭代
export calculation


using ..Mystruct
using ..Mystruct.Resolution
using ..Mystruct.visual
using ..Mystruct.physical_quantity
import ..Mystruct.Potential_Matrix.T_Matrix

using CSV, DataFrames

function record(P::Parameter, Wave::wave_function)                    #将计算的数据写入csv文件中
    local df_1 = DataFrame()
    local df_2 = DataFrame()

    df_1.x = P.sampling
    df_1.Probability = abs2.(Wave.real_space)

    df_2.a_t = Wave.Time_dipole
    df_2.t = LinRange(P.Δt, Wave.Time, P.Step_t)

    CSV.write("High_Harmonic_Generation.csv", df_2)
    CSV.write("Time_evolution_wave_function.csv", df_1)
end


function calculation(P::Parameter, Wave::wave_function, fig::visualization)
    local T_0_Matrix = T_Matrix(P)                         #因为已经将演化过程函数化了,有些参数如果不想每次循环重新计算,就在初始的顶端
                                                         #定义好常矩阵方便传参,没传参的默认每次循环计算一遍,这里的V矩阵因为含时所以每次都需要重新计算

    positive_fft!(Wave)                                  #初始化动量空间波函数
    if imag(P.Δt) == 0.0

        for i in 1:P.Step_t                                   #共演化Step_t次
            resolution(P, Wave, Momentum_operator = T_0_Matrix)                                                        

            Wave.Time_dipole[i] = quantity_a(P, Wave)                                  #充分利用这个计算过程进行参数采集, 
                                                                                   #这里之所以要i+1因为我们第一次初始的情况也计算在内
        #    plot_real(P, Wave, fig)
        #    plot_momentum(P, Wave, fig)
            Wave.Time = i*P.Δt                                                         #经过上面的运算相当于经过了Δt的时间,运行时间+Δt
        end

        record(P, Wave)                                                                   #将计算得到的数据保存到csv文件中去

        #println("Calculation is finished! Now is generating the gif.")
        #generate_Animi(fig)
        #println("The animation is stored in you program directory.")
    else
        return throw(DomainError("The time is not a real number"))
    end

end

end
module RunHHg
include("Main_module.jl"); using .Mystruct      #这个文件是用来直接计算高次谐波的, 主要求解演化过后的波函数的概率密度,以及随时间的加速度,以及我们关系的高次谐波
using CSV, DataFrames                           #用来读取数据

df = CSV.read("Ground_Wave_Func.csv", DataFrame)    #加载基态波函数,之后会自动计算动量空间的波函数

P = Parameter()        #初始化三个数据结构
gif = visualization()
Wave = wave_function()
Wave.real_space = @. real(parse(Complex{Float64}, df.wave))    #读取过程的类型转化,其实当时存储的时候完全可以只取其模来当作基态的波函数,但实际上这个相位因子并不大重要 


calculation(P, Wave, gif)
end

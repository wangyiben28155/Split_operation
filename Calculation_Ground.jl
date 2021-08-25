module RunGround
include("Main_module.jl"); using .Mystruct      #我们直接把所有的函数模块都放在Mystruct里面, 
                                                #这个文件是用来计算基态波函数的,并储存到csv文件当中,用于后续的演化过程的计算
P = Parameter() 
gif = visualization()
Wave = wave_function()

ITR(P, Wave)
plot_Ground()
end
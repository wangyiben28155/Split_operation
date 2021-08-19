module Resolution

export positive_fft!, inverse_fft!, resolution

using ..Mystruct
using ..Mystruct.Potential_Matrix

using FFTW, SparseArrays                         #在模块内使用快速傅里叶变换的包

function positive_fft!(Wave::wave_function)
    Wave.momentum_space = fftshift(fft(Wave.real_space))           #按照书中的值的范围为-N/2+1到N/2
end

function inverse_fft!(Wave::wave_function)
    Wave.real_space = ifft(ifftshift(Wave.momentum_space))
end

function resolution(P::Parameter, Wave::wave_function;                        #这里还是把两个算符都当作默认参量使用,方便进行固定或者改变
     Momentum_operator::SparseMatrixCSC = T_Matrix(P), V_operator::SparseMatrixCSC = V_Matrix(P, Wave)) #这里因为这个函数其实返回的是固定的矩阵,不像势能因为电场有个时间变化的项干扰

        Wave.momentum_space = Momentum_operator * Wave.momentum_space                          
        inverse_fft!(Wave)
        Wave.real_space = V_operator * Wave.real_space
        positive_fft!(Wave)
        Wave.momentum_space = Momentum_operator * Wave.momentum_space
        inverse_fft!(Wave)                                                                      #这里已经形成一个闭环,实空间的波函数和动量空间的是一一对应的
end

end
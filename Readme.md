# 几点说明

- 想要运行项目的文件, 在$REPL$中输入include("文件名.jl")或者在cmd,powershell中运行julia 文件名.jl即可,注意julia的REPL的路径要调整到当前文件的路径(用cd("当前路径")注意双斜杠), 在cmd中同理

- 使用的包可以直接用以下的命令全部安装即可

  ]activate 当前项目路径

  instantiate

- 在REPL中要设置python的路径ENV["PYTHON"]="python的路径"(注意要双斜杠)为想要使用的python的路径,然后用

  ]build PyCall即可
  
- 然后在python中安装matplotlib(这是julia的PyPlot调用的backend), celluloid, 和scipy即可, python包的安装网上应该有很多教程, 这里不再赘述.


using CSV, DataFrames, Plots, Statistics

var = 16
dfA = CSV.read("A_$var.csv", DataFrame)
dfB = CSV.read("B_$var.csv", DataFrame)

# расчет коэффициентов регрессионной прямой по МНК
function ols_coeff(x::Vector,y::Vector)
    x_mean = mean(x)
    y_mean = mean(y)
    b1 = sum((x.-x_mean).*(y.-y_mean))/sum((x.-x_mean).^2)
    b0 = y_mean-b1*x_mean
    return b0, b1
end

# расчет RSS
function RSS_calc(y_init::Vector,y_model::Vector)
    RSS=sum((y_init.-y_model).^2)
    return RSS
end

# расчет R2
function R2_calc(y_init::Vector,y_model::Vector)
    RSS=sum((y_init.-y_model).^2)
    SST=sum((y_init.-mean(y_init)).^2)
    R2=1-(RSS/SST)
    return R2
end

bo,b1=ols_coeff(dfA.x, dfA.y)
yA=bo.+dfA.x*b1
RSS_A=round(RSS_calc(dfA.y,yA), digits =2)
R2_A=round(R2_calc(dfA.y,yA), digits = 2)
scatter(dfA.x, dfA.y, label = "исходные данные")
plot!(dfA.x, yA, label = "модель",lw = 3, xlabel = "x", ylabel = "y",
title = ("Выборка A, вариант $var"), annotate = ([30,30], [10,5],
[text("RSS=$RSS_A", 10, :left), text("R2=$R2_A", 10, :left)]))
savefig("Модель А.png")

bo,b1=ols_coeff(dfB.x, dfB.y)
yB=bo.+dfB.x*b1
RSS_B=round(RSS_calc(dfB.y,yB))
R2_B=round(R2_calc(dfB.y,yB), digits = 2)
scatter(dfB.x, dfB.y, label = "исходные данные")
plot!(dfB.x, yB, label = "модель",lw = 3, xlabel = "x", ylabel = "y",
title = ("Выборка B, вариант $var"), annotate = ([20,20], [70,5],
[text("RSS=$RSS_B", 10, :left), text("R2=$R2_B", 10, :left)]))
savefig("Модель B.png")

residual_A=dfA.y-yA
#гистограмма
histogram(residual_A, bins=21, normalize=:pdf, title="Гистограмма
остатков", label = "Выборка А", xlabel = "x", ylabel = "P(x)" )
savefig("Гистограмма А.png")
#скаттерограмма
pop!(residual_A) #убираем последний элемент из массива
residual_x=residual_A #данные по оси x
residual_A=dfA.y-yA #восстанавливаем исходный массив
popfirst!(residual_A) #убираем первый элемент из массива
residual_y=residual_A #данные по оси y
scatter(residual_x, residual_y, title="Скаттерограмма остатков", label
= "Выборка А", xlabel = "residual", ylabel = "residual-1")
hline!([0,0], linestyle=:dash, lw=2, label="")
savefig("Скатеррограмма А.png")

residual_B=dfB.y-yB
#гистограмма
histogram(residual_B, bins=21, normalize=:pdf, title="Гистограмма
остатков", label = "Выборка B", xlabel = "x", ylabel = "P(x)" )
savefig("Гистограмма В.png")
#скаттерограмма
pop!(residual_B) #убираем последний элемент из массива
residual_x=residual_B #данные по оси x
residual_B=dfB.y-yB #восстанавливаем исходный массив
popfirst!(residual_B) #убираем первый элемент из массива
residual_y=residual_B #данные по оси y
scatter(residual_x, residual_y, title="Скаттерограмма остатков", label
= "Выборка B", xlabel = "residual", ylabel = "residual-1")
hline!([0,0], linestyle=:dash, lw=2, label="")
savefig("Скаттерограмма B.png")
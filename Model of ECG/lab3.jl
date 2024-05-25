using Plots
struct WaveParams
name::String
A::Float64
b1::Float64
b2::Float64
mu::Float64
end
all_waves = Vector{WaveParams}()
push!(all_waves, WaveParams("P",0.11,0.014,0.014,0.399))
push!(all_waves, WaveParams("Q",-0.004,0.008,0.008,0.45))
push!(all_waves, WaveParams("R",1.453,0.008,0.008,0.474))
push!(all_waves, WaveParams("S",-1.053,0.007,0.007,0.495))
push!(all_waves, WaveParams("ST",0.063,0.04,0.04,0.574))
push!(all_waves, WaveParams("T",0.52,0.056,0.024,0.7))
depression_waves = Vector{WaveParams}()
push!(depression_waves, WaveParams("P",0.04,0.03,0.03,0.203))
push!(depression_waves, WaveParams("Q",0,0.066,0.066,0.266))
push!(depression_waves, WaveParams("R",0.64,0.016,0.026,0.296))
push!(depression_waves, WaveParams("S",-0.1,0.03,0.03,0.4))
push!(depression_waves, WaveParams("ST",-0.23,0.15,0.2,0.45))
push!(depression_waves, WaveParams("T",0.06,0.1,0.08,0.7))
function model_cmpx(waves::Vector{WaveParams}, t0m::Float64,
Fs::Number)
time_grid = collect(0:1/Fs:t0m) # временная сетка
out = fill(0.0,length(time_grid))
for point in waves
# числитель формулы для Z(t)
num = (time_grid .- point.mu).^2
# знаменатель считается в зависимости от mu
# так задается ассиметричность волны
den = Float64[]
for t in time_grid
if t<point.mu
push!(den, point.b1^2)
else
push!(den, point.b2^2)
end
end
delta_A = -0.07:0.01:0.07 # допуски на вариабельность амплитуды
Amp = point.A*(1+rand(delta_A))

# в соответствии с формулой (2) результаты суммируем
out .+= Amp*exp.(-num./den)
end
return out
end
normal_ecg = Float64[]
depression_ecg = Float64[]
Fs=500; T=1/Fs;
HR=60 # ЧСС
t0 = 60/HR
V = -0.2:0.01:0.2 # допуски на вариабельность ритма
N = 20 # число комплексов
# ЭКГ в норме
for m = 1:N
t0m = t0*(1+rand(V))
out = model_cmpx(all_waves, t0m, Fs)
append!(normal_ecg,out)
end
# ЭКГ с депрессией ST сегмента
for m = 1:N
t0m = t0*(1+rand(V))
out = model_cmpx(depression_waves, t0m, Fs)
append!(depression_ecg,out)
end
# сигнал с помехами
time_normal=collect(0:1/Fs:(length(normal_ecg)-1)/Fs) # вектор времени для синусоид
drift=2*sin.(2*pi*0.05*time_normal) # дрейф нуля
net=0.3*sin.(2*pi*50*time_normal) # сетевая наводка
signal=normal_ecg.+drift.+net # получившийся сигнал
time_depression=collect(0:1/Fs:(length(depression_ecg)-1)/Fs)
# графики
plotly()
plot(time_normal, normal_ecg, label="ecg", title="ЭКГ в норме",
xlabel="t, с", ylabel="U, мВ")
plot(time_depression, depression_ecg, label="ecg", title="ЭКГ с
депрессией ST", xlabel="t, с", ylabel="U, мВ")
plot(time_normal, signal, title="ЭКГ с помехами",label="ecg",
xlabel="t, с", ylabel="U, мВ")
# в удобном масштабе
plot(time_normal, normal_ecg, label="ecg", title="ЭКГ в норме",
xlabel="t, с", ylabel="U, мВ", xlim=(5.5,9), ylim=(-2,2))
plot(time_depression, depression_ecg, label="ecg", title="ЭКГ с
депрессией ST", xlabel="t, с", ylabel="U, мВ", xlim=(5.7,8.3), ylim=(-
1,1))
plot(time_normal, signal, title="ЭКГ с помехами",label="ecg",
xlabel="t, с", ylabel="U, мВ", xlim=(5.5,9), ylim=(-0.5,3.5))
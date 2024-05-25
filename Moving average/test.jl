using CSV
using DataFrames
using Plots

# фильтр скользящего среднего
mutable struct SlideMeanFilter{T}
    buf::Vector{T} # кольцевой буфер
    k::Int # состояние фильтра
    need_restart::Bool # маркер инциализации фильтра
    
    function SlideMeanFilter{T}(window::Int) where T
        new(fill(T(0), window-1), 1, true) # буфер заполнится нулями
    end
end

# фильтр взвеш скользящего среднего
mutable struct VzveshSlideMeanFilter{T}
    buf_k::Vector{T} # кольцевой буфер
    k::Int # состояние фильтра
    need_restart::Bool # маркер инциализации фильтра
    buf_w::Vector{Int64} #буфер весов фильтра
    
    function VzveshSlideMeanFilter{T}(window::Int) where T
        new(fill(T(0), window-1), 1, true, collect(Int64, 1:window-1)) # буферы заполняются
    end
end

function SMF(obj::SlideMeanFilter{T}, inp::AbstractVector{T}, out::AbstractVector) where T
    buf, k = obj.buf, obj.k
    if obj.need_restart # инициализация на первой точке
        # заполняем буфер первой точкой
        fill!(buf, inp[1])
        # отмечаем, что инициализация уже не требуется
        obj.need_restart = false
    end

    for i in 1:length(inp)
        sum_x = inp[i]
        # сумма всех элементов в буфере + 1 новая точка
        for xi in buf
            sum_x += xi
        end
        #расчет среднего
        window = length(buf) + 1
        y = sum_x / window
        # в буфер записываем новую точку
        buf[k] = inp[i]
        k += 1
        # проверка, не кончился ли буфер
        if k > length(buf)
            k = 1
        end
        # фиксируем состояние фильтра
        obj.k = k
        push!(out, y)
    end
    return out
end

function VSMF(obj::VzveshSlideMeanFilter{T}, inp::AbstractVector{T}, out::AbstractVector) where T
    buf_k, k, buf_w = obj.buf_k, obj.k, obj.buf_w
    if obj.need_restart # инициализация на первой точке
        # заполняем буфер первой точкой
        fill!(buf_k, inp[1])
        # отмечаем, что инициализация уже не требуется
        obj.need_restart = false
    end

    for i in 1:length(inp)
        #расчет знаменателя
        window = length(buf_k) + 1
        sum_w = 0
        for i in 1:window
            sum_w += i
        end 
        #расчет числителя   
        sum = inp[i]*window
        for j in 1:length(buf_k)
            sum += buf_k[j]*buf_w[j]
        end
        #расчет среднего
        y = sum/sum_w
        #перезапись буфера весов
        a = fill(0, length(buf_w))
        for i in 1:length(buf_w)
            if i==length(buf_w)
                a[1] = buf_w[i]
            else
                a[i+1] = buf_w[i]
            end
        end
        buf_w = a
        # в буфер записываем новую точку
        buf_k[k] = inp[i]
        k += 1
        # проверка, не кончился ли буфер
        if k > length(buf_k)
            k = 1
        end
        # фиксируем состояние фильтра
        obj.k = k
        push!(out, y)
    end
    return out
end


#чтение файла
data_start = CSV.read("hr10.csv", DataFrame)
#частота дискритизации
Fs = 0.1
#получение отсчетов
samples = collect(Int64, 1:length(data_start[!,1]))
#перевод в секунды
time_samples = samples.*Int64(1/Fs)
#перевод в вариант
time_samples_v16_s = filter(x -> x>=16*6*60, time_samples)
#перевод в часы
time_samples_v16_h = time_samples_v16_s./3600 
#получение значений для варианта
inp = data_start[Int(6*96/(Fs*10)):length(data_start[!,1]),1]
#построение графика
min_data = minimum(inp)
max_data = maximum(inp)
time_data_v16 = time_samples_v16_h[end] - time_samples_v16_h[1]
plot_v16 = plot(time_samples_v16_h, inp,title="ЧСС_в16", label=["Тренд ЧСС"], linewidth=0.8, legend = false)
xlims!(time_samples_v16_h[1], time_samples_v16_h[end])
ylims!(min_data, max_data) 

#фильтрация
window = 3
slide_flt = SlideMeanFilter{eltype(inp)}(window)
v_slide_flt = VzveshSlideMeanFilter{eltype(inp)}(window)
data_out_smf_1 = []
data_out_vsmf_1 = []
data_out_smf_1 = SMF(slide_flt, inp, data_out_smf_1)
data_out_vsmf_1 = VSMF(v_slide_flt, inp, data_out_vsmf_1)

window = 6
slide_flt = SlideMeanFilter{eltype(inp)}(window)
v_slide_flt = VzveshSlideMeanFilter{eltype(inp)}(window)
data_out_smf_2 = []
data_out_vsmf_2 = []
data_out_smf_2 = SMF(slide_flt, inp, data_out_smf_2)
data_out_vsmf_2 = VSMF(v_slide_flt, inp, data_out_vsmf_2)

window = 12
slide_flt = SlideMeanFilter{eltype(inp)}(window)
v_slide_flt = VzveshSlideMeanFilter{eltype(inp)}(window)
data_out_smf_3 = []
data_out_vsmf_3 = []
data_out_smf_3 = SMF(slide_flt, inp, data_out_smf_3)
data_out_vsmf_3 = VSMF(v_slide_flt, inp, data_out_vsmf_3)



#выбираем небольшой по времени участок для визуализации
istart = 1500 #начало участка
len = 30 #длина участка
data_visual = inp[istart:istart+len-1]
out_smf_visual_1 = data_out_smf_1[istart:istart+len-1]
out_smf_visual_2 = data_out_smf_2[istart:istart+len-1]
out_smf_visual_3 = data_out_smf_3[istart:istart+len-1]
out_vsmf_visual_1 = data_out_vsmf_1[istart:istart+len-1]
out_vsmf_visual_2 = data_out_vsmf_2[istart:istart+len-1]
out_vsmf_visual_3 = data_out_vsmf_3[istart:istart+len-1]

#вывод графиков
plot(data_visual, marker = :circle, title = "Скользящее среднее", label = "Исходный", xlabel = "Отсчеты", ylabel = "ЧСС, уд./мин", xlim = (0, 30), ylim = (73, 90))
plot!(out_smf_visual_1, marker = :circle, label = "Window=3", xlim = (0, 30), ylim = (73, 90))
plot!(out_smf_visual_2, label = "Window=6", marker = :circle)
plot!(out_smf_visual_3, label = "Window=12", marker = :circle)
savefig("Скользящее среднее")

plot(data_visual, marker = :circle, title = "Взвешанное скользящее среднее", label = "Исходный", xlabel = "Отсчеты", ylabel = "ЧСС, уд./мин", xlim = (0, 30), ylim = (73, 90))
plot!(out_vsmf_visual_1, marker = :circle, label = "Window=3", xlim = (0, 30), ylim = (73, 90))
plot!(out_vsmf_visual_2, label = "Window=6", marker = :circle)
plot!(out_vsmf_visual_3, label = "Window=12", marker = :circle)
savefig("Взвешенное скользящее среднее")
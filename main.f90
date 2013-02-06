program main
use signal_generator
use omp_lib
use ifport
!------------------------------------------------------------------------------
!--------Входные данные--------------------------------------------------------
    real*4 :: min_s(3) = [-3., -3., 1.5]                     !Минимум координат в сетке    
    real*4 :: max_s(3) = [3., 3., 4.5]                       !Макс координат в сетке
    real*4,    parameter :: speed  = 4                       !Скорость в среде
    integer*4, parameter :: ticks_per_sec = 500              !Количество отсчетов в секунду
    integer*4, parameter :: num_d_x = 5                      !Размерность датчиков по Ox
    integer*4, parameter :: num_d_y = 5                      !Размерность датчиков по Oy
    integer*4, parameter :: height = 1000                    !Высота окна
    integer*4, parameter :: split  = 40                      !Количество усредняемых значений
    integer*4, parameter :: num_s  = 32                      !Размерность сетки    
    integer*4, parameter :: num = num_d_x * num_d_y          !Общее количество датчиков
!------------------------------------------------------------------------------
!--------Вспомогательные переменные--------------------------------------------
    integer*4, dimension(:),        allocatable :: intgodog    !Целая часть годографа
    real*4,    dimension(:),        allocatable :: avret       !Усредненное значение по окну и точке        
    real*4,    dimension(:),        allocatable :: ret         !Возвращаемое значение                               
    real*4,    dimension(:),        allocatable :: godograf    !Годограф
    real*4,    dimension(:),        allocatable :: delta_godog !Дробная часть годографа    
    real*4,    dimension(:),        allocatable :: f_i         !F(i)        
    real*4,    dimension(:,:),      allocatable :: w1          !1 окно    
    real*4,    dimension(:,:),      allocatable :: w2          !2 окно    
    real*4,    dimension(:,:,:),    allocatable :: coord_d     !Координаты датчиков  
    real*4,    dimension(:,:,:,:),  allocatable :: res         !Результат

    real*4 point(3)                                            !Координаты точек
    real*4 r(3)                                                !Радиус вектор между точками
    real*4 delta(3)                                            !Размеры сетки

    integer*4 num_part                                         !Число частей    
    integer*4 i,j,k,l,m,t,s                                    !Вспомогательные переменные
    integer time1, time2                                       !Время работы алгоритма
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!----Генерация данных----------------------------------------------------------
    type (struct) inpData                                  !Источник сигнала
    real*4, dimension(:,:), allocatable :: out_data        !Массив выходных данных
    !---Input data----
    inpData%rxCount         = (/ num_d_x, num_d_y/)        !Размеры сетки приемников
    inpData%txCount         = (/ num_s, num_s, num_s/)     !Размеры сетки источников
    inpData%windowHeight    = height * 2                   !Высота окна
    inpData%txCoord         = (/ 0, 0, 3/)                 !Координаты источника
    inpData%X               = (/-3, 3/)                    !(min X, max X) сетки источников
    inpData%Y               = (/-3, 3/)                    !(min Y, max Y) сетки источников
    inpData%Z               = (/ 1.5, 4.5/)                !(min Z, max Z) сетки источников
    inpData%dT              = 0.002                        !Время между отсчетами
    inpData%alpha           = 70.0                         !Кооофиценты для
    inpData%omega           = 20.0                         !импульса Пузырева
    inpData%noiseCoef       = 1                            !Коэф. шума
    inpData%startTime       = 0.2                          !Время возн. сигнала
    inpData%speed           = speed                        !Скорость в среде
    !-----------------        
    allocate(out_data(0:inpData%windowHeight-1, inpData%rxCount(1) * inpData%rxCount(2)))
    !Генерация    
    call SG_generateData(inpData, out_data)                
!----Конец генерации данных----------------------------------------------------
!------------------------------------------------------------------------------

!----Выделение памяти под массивы----------------------------------------------
!------------------------------------------------------------------------------    
    ALLOCATE (res(height/split, 0:num_s, 0:num_s, 0:num_s))
    ALLOCATE (coord_d(0:num_d_x-1, 0:num_d_y-1, 3))
    ALLOCATE (w1(0:height-1, num_d_x*num_d_y))                
    ALLOCATE (w2(0:height-1, num_d_x*num_d_y))                
    ALLOCATE (avret(height/split))                
    ALLOCATE (ret(1:height))                
    ALLOCATE (godograf(num))                         
    ALLOCATE (intgodog(num))                         
    ALLOCATE (delta_godog(num))                    
    ALLOCATE (f_i(0:height))

    delta = max_s-min_s
    num_part = height/split
    !Вывод частоты сетки
    print * ,"Delta: ",delta/(num_s-1)
    coord_d=0

    !Заполнение координат сетки
    do i = 0,num_d_x-1
        do j = 0,num_d_y-1
            coord_d(i,j,1) = min_s(1) + i*delta(1)/(num_d_x-1)
            coord_d(i,j,2) = min_s(2) + j*delta(2)/(num_d_y-1)
            coord_d(i,j,3) = 0
        end do
    end do

    !Получение сгенерированных данных
    do i = 1,num
        w1(:,i) = out_data(:height-1,i)
        w2(:,i) = out_data(height:,i)
    end do

    !Начало отсчета    
    time1 = time()
    
    !Распараллеленый цикл
    !!$OMP PARALLEL DO schedule(dynamic,1) collapse (3) private ( point, godograf, f_i, ret, intgodog, delta_godog, avret, l,m, t) NUM_THREADS(8)
    do i = 0,num_s-1
        do j = 0,num_s-1
            do k = 0,num_s-1
                !Получение координат обрабатываемой точки
                point = [i,j,k] * delta / (num_s - 1) + min_s    
                
                !Вычисление годографа        
                godograf = 0
                do l = 0,num_d_x-1
                    do m = 0,num_d_y-1
                        !Высчитываем радиус вектор
                        r(:) = coord_d(l,m,:) - point(:)
                        !Высчитываем время достижения датчика                        
                        godograf(l*num_d_y+m+1) = sqrt(r(1) ** 2 + r(2) ** 2 + r(3) ** 2)  / speed * ticks_per_sec            
                        r = 0
                    end do
                end do

                !Нормируем годограф до нуля
                godograf = godograf - minval(godograf)                
                !Вычисление целой части годографа
                intgodog = floor(godograf)
                !Вычисление дробной части годографа
                delta_godog = godograf - intgodog
                ret = 0
                !Выполняем суммирование
                do t = 1,num
                    !Получение данных для когерентного суммирования                    
                    f_i(0:height-1-intgodog(t)) = w1(intgodog(t):(height-1),t)
                    f_i(height-intgodog(t):height) = w2(0:intgodog(t),t)        
                       !Выполняем суммирование
                    ret = ret + f_i(0:height-1)* (1 - delta_godog(t)) + delta_godog(t)*f_i(1:height) 
                end do

                !Производим усреднение
                avret = 0
                do m = 0,num_part-1
                    avret(m+1) = sqrt(sum(ret(m*split+1:(m+1)*split) ** 2)) / split
                end do
            
                !Запись результата усреднения
                res(:,i,j,k) = avret
                avret = 0    
            end do
        end do
    end do
    !!$OMP END PARALLEL DO
    !Конец параллельного цикла

    !Высчитываем время работы 
    time2 = time()     
    time2 = time2 - time1
    
    !Вывод результата
    do s= 1,num_part
        print * , "Delta loc: ", maxloc(res(s,:,:,:)) * delta/(num_s-1) + min_s - inpData%txCoord, &
                  "Value: ",     maxval(res(s,:,:,:))
      end do
    print * , "Elapsed time: ", time2, "seconds"


!------------------------------------------------------------------------------
!----Освобождение памяти-------------------------------------------------------
    DEALLOCATE(out_data)     
    DEALLOCATE (res)
    DEALLOCATE (coord_d)
    DEALLOCATE (w1)                
    DEALLOCATE (w2)            
    DEALLOCATE (godograf)                 
    DEALLOCATE (intgodog)                     
    DEALLOCATE (delta_godog)                  
    DEALLOCATE (f_i)
end program

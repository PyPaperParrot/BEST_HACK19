# BEST_HACK19
Anunaki's team repository for BEST HACK'19 (data science track)
# Задача
Создать математическую модель летящего груза, как материальной точки для определения начальных условий сброса: координат (X, Z) и направления вектора скорости (угол alpha - изображен положительным на рисунке). При этом надо обеспечить ряд следующих условий:
  
  1. Возможность задания начальных услови:
  
    H_0 - высоты сброса не более 1400 м;
  
    V_0 - начальная скорость груза не более 250 м/с;
  
    m - масса;
  
    F_a - аэродинамическая сила на грузе, направление действия противоположно вектору скорости;
  
    База данных (БД) ветров - зависимость направления ветра от высоты.
  
  2. Обеспечить попадание груза в зону приземления, приветствуется большая точность.

# Работа с приложением
Запуск приложения необходимо производить в ОС Ubuntu.
Все файлы проекта находятся в zip-архиве BEST HACK19.tar.giz. Скачайте данный архив и извлеките в одну директорию все содержащиеся в нем файлы:
    
    1) F.csv -- файл с таблицей 1 (Зависимость аэродинамической силы от скорости движения груза);
    2) Wind.csv -- файл с таблицей 2 (Распределение поля ветров по высоте, где: Y - высота над уровнем земли; Wx, Wz - проекции ветра на ось X и Z соответственно);
    3) main.cpp -- файл с кодом на языке C++, в котором описана математическая модель летящего груза и выполнена основная задача;
    4) trajectory.txt -- файл содержит время, скорость, три координаты по осям в каждый момент времени;
    5) trajectory.py -- файл содержит программу на python для визуализации траектории; 
    6) BEST.sh - файл последовательно запускает скрипиты .cpp и .py;
    
# Запуск
Существует несколько способов запустить приложение:
  1. Правой кнопкой мыши щелкните в директории с файлами и выберите "Открыть в терминале". Далее в командной строке введите "bash BEST.sh". После чего запуститься файл main.cpp, который загрузит необходимые для моделирования данные из файлов F.csv и Wind.csv

    

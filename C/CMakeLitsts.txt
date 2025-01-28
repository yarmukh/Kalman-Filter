cmake_minimum_required(VERSION 3.10)

# Название проекта
project(KalmanFilter)

# Установка стандарта C
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED True)

# Добавление исполняемого файла
add_executable(KalmanFilter Source.c )

target_link_libraries(KalmanFilter m)
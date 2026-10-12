cmake_minimum_required(VERSION 3.16)
project(factorizer LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_AUTOUIC ON)

find_package(Qt6 REQUIRED COMPONENTS Quick QuickControls2)

qt_standard_project_setup()
# Устраняем предупреждение QTP0001
qt_policy(SET QTP0001 NEW)

qt_add_executable(appfactorizer
    WIN32
    main.cpp
    factorizerbackend.cpp

    # Заголовочные файлы:
    factorizerbackend.h
    u128.hpp
    ubig.hpp
    u128_utils.h
    u128_intrinsics.hpp
    siqs.hpp
    generic_factorizer.hpp
)

qt_add_qml_module(appfactorizer
    URI factorizer
    VERSION 1.0
    QML_FILES Main.qml
    SOURCES
        factorizerbackend.h
        factorizerbackend.cpp
        generic_factorizer.hpp
        siqs.hpp
)

target_link_libraries(appfactorizer PRIVATE
    Qt6::Quick
    Qt6::QuickControls2
)
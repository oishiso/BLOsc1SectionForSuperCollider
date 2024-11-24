cmake_minimum_required(VERSION 3.15)

project(BLOsc1section)

# SuperColliderのインクルードディレクトリを設定
set(SUPERCOLLIDER_PLUGIN_INTERFACE_PATH "/Users/so/supercollider/include/plugin_interface")
set(SUPERCOLLIDER_COMMON_PATH "/Users/so/supercollider/include/common")

include_directories(${SUPERCOLLIDER_PLUGIN_INTERFACE_PATH} ${SUPERCOLLIDER_COMMON_PATH})

# BLOsc1sectionライブラリを作成
add_library(BLOsc1section MODULE BLOsc1section.cpp)


# コンパイルオプションを設定
target_compile_options(BLOsc1section PRIVATE -arch arm64 -std=c++17)
target_compile_definitions(BLOsc1section PRIVATE SC_USE_FASTCALL=0)

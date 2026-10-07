cmake_minimum_required(VERSION 3.10)
project(frida1 Fortran)

enable_language(Fortran)

# ========== Compiler settings ==========
set(CMAKE_Fortran_COMPILER gfortran)
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O -A -W -s -w -ffixed-line-length-132 -std=legacy")

# Optional: static linking of Fortran and system libs
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")

# ========== Source Groups ==========

# L group
set(LL
    l1.f l2.f l3.f l4.f l5.f
    l6.f l0x11.f
)

# G group
set(GG
    g1.f g2.f
)

# I group
set(II
    i00.f i01.f i10.f i20.f i23.f i25.f
    i30.f i32.f i40.f i41.f i42.f i43.f
    i50.f i60.f i66.f i67.f
    i70.f i71.f i72.f i73.f i74.f i75.f
    i76.f i80.f i87.f i95.f i99.f
    i77.f i81.f
)

# All objects
set(SOURCES ${II} ${GG} ${LL})

# ========== Executable ==========
add_executable(frida1 ${SOURCES})

# ========== Include Directories (if needed) ==========
# include_directories(/usr/local/include)

# ========== Libraries ==========
find_library(MATH_LIB m)
target_link_libraries(frida1 ${MATH_LIB})

target_link_libraries(frida1
#    -L/usr/local/lib -lget
#    -L/usr/local/lib -lz -ljpeg
#    -L/usr/local/lib/hdf/lib -lhdf5 -lmfhdf -ldf
#    -L/opt/absoft10.2/lib -lV77 -lU77 -lNeXus
#    -L/usr/lib/gcc/x86_64-linux-gnu/ -lgfortran -lm -lnag
)

# Optional: m32 for 32-bit build (ensure 32-bit libs installed)
# set_target_properties(frida1 PROPERTIES COMPILE_OPTIONS "-m32" LINK_FLAGS "-m32")

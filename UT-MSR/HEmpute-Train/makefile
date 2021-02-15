all: HEmpute_Train

CC = g++
comp_flags = -c -O3 -Wall
gsl_flags = -lgsl -lgslcblas -lz
exec_name = bin/HEmpute_Train
LIB_DIR = src

# Define pattern rule for building object files.
%.o: %.cpp
	${CC} ${comp_flags} $< -o $@

objs = \
${LIB_DIR}/hmpt_main.o\
${LIB_DIR}/hmpt_nucleotide.o \
${LIB_DIR}/hmpt_utils.o\
${LIB_DIR}/hmpt_xlog_math.o \
${LIB_DIR}/hmpt_ansi_thread.o \
${LIB_DIR}/hmpt_nomenclature.o \
${LIB_DIR}/hmpt_ansi_string.o \
${LIB_DIR}/hmpt_annot_region_tools.o \
${LIB_DIR}/hmpt_variation_tools.o \
${LIB_DIR}/hmpt_imputation_utils.o 

HEmpute_Train: ${objs}
	${CC} -O3 -pthread ${gsl_flags} -o ${exec_name} ${objs}

clean:
	rm -f ${objs} ${exec_name}

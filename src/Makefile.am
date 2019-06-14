## Process this file with automake to produce Makefile.in

lib_LTLIBRARIES = libGI.la
libGI_la_SOURCES = gi_blas.h gi_container.h gi_context.h gi_cutter.h gi_gl.h \
                   gi_image.h gi_math.h gi_memory.h gi_mesh.h \
                   gi_multiresolution.h gi_numerics.h gi_parameterizer.h \
                   gi_sampler.h gi_thread.h \
                   gi_blas.c gi_container.c gi_context.c gi_cutter.c gi_gl.c \
                   gi_image.c gi_math.c gi_memory.c gi_mesh.c \
                   gi_multiresolution.c gi_numerics.c gi_parameterizer.c \
                   gi_sampler.c gi_thread.c
libGI_la_CPPFLAGS = -I$(srcdir)/../include
libGI_la_CFLAGS = $(SSECFLAGS)
libGI_la_LDFLAGS = -version-info $(LIBRARY_VERSION)

dist_noinst_DATA = Mainpage.dox


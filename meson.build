project('tspbwt', ['c', 'cpp'], version : '0.1', default_options : ['c_std=c99', 'cpp_std=c++20'])

boost_dep = dependency('boost')

kastore_inc = include_directories('tskit/c/subprojects/kastore')
kastore = static_library('kastore', 'tskit/c/subprojects/kastore/kastore.c')
kastore_dep = declare_dependency(link_with : kastore, include_directories: kastore_inc)

cc = meson.get_compiler('c')
m_dep = cc.find_library('m', required : false)

tskit_inc = [kastore_inc, include_directories('tskit/c')]
lib_deps = [m_dep, kastore_dep, boost_dep, declare_dependency(include_directories:tskit_inc)]

lib_sources = [
    'tskit/c/tskit/core.c', 'tskit/c/tskit/tables.c', 'tskit/c/tskit/trees.c',
    'tskit/c/tskit/genotypes.c', 'tskit/c/tskit/stats.c', 'tskit/c/tskit/convert.c', 'tskit/c/tskit/haplotype_matching.c']
lib_headers = [
    'tskit/c/tskit/core.h', 'tskit/c/tskit/tables.h', 'tskit/c/tskit/trees.h',
    'tskit/c/tskit/genotypes.h', 'tskit/c/tskit/stats.h', 'tskit/c/tskit/convert.h', 'tskit/c/tskit/haplotype_matching.h']
tskit_lib = static_library('tskit', sources: lib_sources, dependencies: lib_deps)
tskit_dep = declare_dependency(include_directories:tskit_inc, link_with: tskit_lib)

executable('TSPBWT_v0.1', 'src/tspbwt.cpp', dependencies : tskit_dep, install : true)
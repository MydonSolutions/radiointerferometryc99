project('radiointerferometryc99', 'c',
  version: '0.7.5',
  default_options: [
    'c_std=c99',
  ]
)
src_lst =[]
inc_lst =[]
inc_lst = [
	include_directories('include')
]
dep_lst = [
  dependency('erfa'),
]
cc = meson.get_compiler('c')
m_dep = cc.find_library('m', required : true)
dep_lst += [
  m_dep
]

subdir('src')
subdir('include')

lib_radiointerferometry = library(
    'radiointerferometryc99',
    src_lst,
    include_directories : inc_lst,
    dependencies: dep_lst,
    install: true
)
lib_radiointerferometry_dep = declare_dependency(
  include_directories: inc_lst,
  dependencies: dep_lst,
  link_with: lib_radiointerferometry,
)

build_dir = meson.current_build_dir()
py = import('python').find_installation('python3', required: false)
subdir('tests')

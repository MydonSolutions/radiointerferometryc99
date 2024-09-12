iers_filepath = '@0@/iers_finals2000A.snippet.all'.format(meson.current_source_dir())

test('iers_read', executable(
  'iers_read', ['iers_read.c'],
	dependencies: lib_radiointerferometry_dep,
	),
	args : [iers_filepath],
	is_parallel: false
)

test('posangle', executable(
  'posangle', ['posangle.c'],
	dependencies: lib_radiointerferometry_dep,
	),
	args : [iers_filepath],
	is_parallel: false
)

if py.found()
	# ATA-like accumulation.
	test(
		'posangle_reference',
		py,
		args: [files('posangle_reference.py'), '@0@/libradiointerferometryc99.so'.format(build_dir), iers_filepath],
		is_parallel: false,
		timeout: 0,
	)
endif

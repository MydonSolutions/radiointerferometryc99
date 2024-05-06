test('iers_read', executable(
  'iers_read', ['iers_read.c'],
	dependencies: lib_radiointerferometry_dep,
	),
	args : files(['iers_finals2000A.snippet.all']),
	is_parallel: false
)

test('posangle', executable(
  'posangle', ['posangle.c'],
	dependencies: lib_radiointerferometry_dep,
	),
	args : files(['iers_finals2000A.snippet.all']),
	is_parallel: false
)
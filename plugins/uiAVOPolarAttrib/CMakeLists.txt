SET(OD_MODULE_DEPS AVOPolarAttrib uiODMain)
SET(OD_IS_PLUGIN yes)
SET(OD_MODULE_SOURCES
	uiavopolarattribpi.cc
	uiavopolarattrib.cc)
SET( OD_PLUGIN_ALO_EXEC ${OD_MAIN_EXEC} )
OD_INIT_MODULE()

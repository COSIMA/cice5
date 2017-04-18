
access-om01:
	bld/build.sh nci auscom 3600x2700
access-om1:
	bld/build.sh nci auscom 360x300
access-cm1:
	bld/build.sh nci access 360x300
access-om025:
	bld/build.sh nci auscom 1440x1080
access-cm025:
	bld/build.sh nci access 1440x1080

clean:
	rm -rf build_*

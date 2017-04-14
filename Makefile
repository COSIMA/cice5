
access-om01:
	bld/build.sh nci access-om 3600x2700
access-om1:
	bld/build.sh nci access-om 360x300
access-cm1:
	bld/build.sh nci access-cm 360x300
access-om025:
	bld/build.sh nci access-om 1440x1080
access-cm025:
	bld/build.sh nci access-cm 1440x1080

clean:
	rm -rf build_*

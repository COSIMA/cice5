platform ?= nci

1deg:
	bld/build.sh $(platform) auscom 360x300 debug
025deg:
	bld/build.sh $(platform) auscom 1440x1080
01deg:
	bld/build.sh $(platform) auscom 3600x2700

clean:
	rm -rf build_*

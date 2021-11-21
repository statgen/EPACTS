## Helloworld C++: an example project using Automake

## Place generated object files (.o) into the same directory as their source
## files, in order to avoid collisions when non-recursive make is used.
ACLOCAL_AMFLAGS = -I m4 ${ACLOCAL_FLAGS}
AUTOMAKE_OPTIONS = foreign
SUBDIRS = src data scripts

docker-centos:
	docker build -t epacts/master/centos:latest \
		-f docker/centos/Dockerfile .

docker-ubuntu:
	docker build -t epacts/master/ubuntu:latest \
		-f docker/ubuntu/Dockerfile .

install-exec-hook:
	mkdir --p $(libdir)
	R CMD INSTALL -l $(libdir) $(datadir)/${PACKAGE}/epactsR_3.3.2.tar.gz
	R CMD INSTALL -l $(libdir) $(datadir)/${PACKAGE}/mmSKAT_0.95.tar.gz

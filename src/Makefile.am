## Process this file with automake to generate Makefile.in
##

## Programs

progdir = @bindir@/@PACKAGE@

prog_PROGRAMS     = darith dcirc dclip dcont dcor dgdef dgdist dinit dnadd dnanal dplot drank dspot \
dsymm dtinfo dtmem dtscl dvar fdplot tback tboot tfilt tgen tnadd tplot ddisc dline tarith tgauss tmolly

darith_SOURCES = darith.cc
dcirc_SOURCES  = dcirc.cc
dclip_SOURCES  = dclip.cc
dcont_SOURCES  = dcont.cc
dcor_SOURCES   = dcor.cc
ddisc_SOURCES  = ddisc.cc
dgdef_SOURCES  = dgdef.cc
dgdist_SOURCES = dgdist.cc
dinit_SOURCES  = dinit.cc
dline_SOURCES  = dline.cc
dnadd_SOURCES  = dnadd.cc
dnanal_SOURCES = dnanal.cc
dplot_SOURCES  = dplot.cc
drank_SOURCES  = drank.cc
dspot_SOURCES  = dspot.cc
dsymm_SOURCES  = dsymm.cc
dtinfo_SOURCES = dtinfo.cc
dtmem_SOURCES  = dtmem.cc
dtscl_SOURCES  = dtscl.cc
dvar_SOURCES   = dvar.cc
fdplot_SOURCES = fdplot.cc
tarith_SOURCES = tarith.cc
tback_SOURCES  = tback.cc
tboot_SOURCES  = tboot.cc
tfilt_SOURCES  = tfilt.cc
tgauss_SOURCES = tgauss.cc
tgen_SOURCES   = tgen.cc 
tmolly_SOURCES = tmolly.cc
tnadd_SOURCES  = tnadd.cc
tplot_SOURCES  = tplot.cc

DATE = $(shell date)

install-data-hook: $(progdir)/dadd $(progdir)/dmul $(progdir)/ddiv $(progdir)/dsub $(progdir)/tadd $(progdir)/tmul $(progdir)/tdiv $(progdir)/tsub $(progdir)/dset
	echo "This is $(PACKAGE)-$(VERSION), built on $(DATE)" > $(progdir)/VERSION
	@echo "Installed link targets"

# link targets

$(progdir)/dadd: $(progdir)/darith
	cd $(progdir); ln -sf darith $@

$(progdir)/ddiv: $(progdir)/darith
	cd $(progdir); ln -sf darith $@

$(progdir)/dmul: $(progdir)/darith
	cd $(progdir); ln -sf darith $@

$(progdir)/dsub: $(progdir)/darith
	cd $(progdir); ln -sf darith $@

$(progdir)/dset: $(progdir)/darith
	cd $(progdir); ln -sf darith $@

$(progdir)/tadd: $(progdir)/tarith
	cd $(progdir); ln -sf tarith $@

$(progdir)/tdiv: $(progdir)/tarith
	cd $(progdir); ln -sf tarith $@

$(progdir)/tmul: $(progdir)/tarith
	cd $(progdir); ln -sf tarith $@

$(progdir)/tsub: $(progdir)/tarith
	cd $(progdir); ln -sf tarith $@

# back to normal

INCLUDES = -I../include -I../.

LDADD = libtomog.la

## Library

lib_LTLIBRARIES = libtomog.la 

libtomog_la_SOURCES = trm_trail.cc trm_dmap.cc optr.cc


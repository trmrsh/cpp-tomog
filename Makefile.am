## Process this file with automake to generate Makefile.in

SUBDIRS = src doc include

EXTRA_DIST = pgplot_libs

export: dist
	scp $(PACKAGE)-$(VERSION).tar.gz $(WEB_SERVE):/$(WEB_PATH)/software/.
	ssh $(WEB_SERVE) "cd $(WEB_PATH)/software ; \rm -f $(PACKAGE).tar.gz ; ln -s $(PACKAGE)-$(VERSION).tar.gz $(PACKAGE).tar.gz; \rm -fr $(PACKAGE); mkdir $(PACKAGE)"
	scp -r doc/html $(WEB_SERVE):$(WEB_PATH)/software/$(PACKAGE)/.

make-html:
	cd doc; make html

ALIASES   = Tomog

.PHONY : $(ALIASES)

DATE      = $(shell date)

progdir = @bindir@/@PACKAGE@

prog_DATA = $(ALIASES)

$(ALIASES):
	@echo '#' > $(ALIASES)
	@echo '# This file defines aliases for the $(PACKAGE) v$(VERSION) software.' >> $(ALIASES)
	@echo '# Source it to define them.' >> $(ALIASES)
	@echo '#' >> $(ALIASES)
	@echo 'alias dadd      $(progdir)/dadd'       >> $(ALIASES)
	@echo 'alias dcirc     $(progdir)/dcirc'      >> $(ALIASES)
	@echo 'alias dclip     $(progdir)/dclip'      >> $(ALIASES)
	@echo 'alias dcont     $(progdir)/dcont'      >> $(ALIASES)
	@echo 'alias dcor      $(progdir)/dcor'       >> $(ALIASES)
	@echo 'alias ddisc     $(progdir)/ddisc'      >> $(ALIASES)
	@echo 'alias ddiv      $(progdir)/ddiv'       >> $(ALIASES)
	@echo 'alias dgdef     $(progdir)/dgdef'      >> $(ALIASES)
	@echo 'alias dgdist    $(progdir)/dgdist'     >> $(ALIASES)
	@echo 'alias dinit     $(progdir)/dinit'      >> $(ALIASES)
	@echo 'alias dline     $(progdir)/dline'      >> $(ALIASES)
	@echo 'alias dmul      $(progdir)/dmul'       >> $(ALIASES)
	@echo 'alias dnadd     $(progdir)/dnadd'      >> $(ALIASES)
	@echo 'alias dnanal    $(progdir)/dnanal'     >> $(ALIASES)
	@echo 'alias dplot     $(progdir)/dplot'      >> $(ALIASES)
	@echo 'alias drank     $(progdir)/drank'      >> $(ALIASES)
	@echo 'alias dset      $(progdir)/dset'       >> $(ALIASES)
	@echo 'alias dspot     $(progdir)/dspot'      >> $(ALIASES)
	@echo 'alias dsub      $(progdir)/dsub'       >> $(ALIASES)
	@echo 'alias dsymm     $(progdir)/dsymm'      >> $(ALIASES)
	@echo 'alias dtinfo    $(progdir)/dtinfo'     >> $(ALIASES)
	@echo 'alias dtmem     $(progdir)/dtmem'      >> $(ALIASES)
	@echo 'alias dtscl     $(progdir)/dtscl'      >> $(ALIASES)
	@echo 'alias dvar      $(progdir)/dvar'       >> $(ALIASES)
	@echo 'alias fdplot    $(progdir)/fdplot'     >> $(ALIASES)
	@echo 'alias tadd      $(progdir)/tadd'       >> $(ALIASES)
	@echo 'alias tback     $(progdir)/tback'      >> $(ALIASES)
	@echo 'alias tboot     $(progdir)/tboot'      >> $(ALIASES)
	@echo 'alias tdiv      $(progdir)/tdiv'       >> $(ALIASES)
	@echo 'alias tfilt     $(progdir)/tfilt'      >> $(ALIASES)
	@echo 'alias tgen      $(progdir)/tgen'       >> $(ALIASES)
	@echo 'alias tmolly    $(progdir)/tmolly'     >> $(ALIASES)
	@echo 'alias tmul      $(progdir)/tmul'       >> $(ALIASES)
	@echo 'alias tnadd     $(progdir)/tnadd'      >> $(ALIASES)
	@echo 'alias tplot     $(progdir)/tplot'      >> $(ALIASES)
	@echo 'alias tsub      $(progdir)/tsub'       >> $(ALIASES)
	@echo 'alias tgauss    $(progdir)/tgauss'     >> $(ALIASES)
	@echo ' ' >> $(ALIASES)
	@echo 'echo " "' >> $(ALIASES)
	@echo 'echo "Welcome to $(PACKAGE) version $(VERSION), built $(DATE)"' >> $(ALIASES)
	@echo 'echo " "' >> $(ALIASES)
	@echo 'echo "Point a browser at ${prefix}/html/$(PACKAGE)/index.html for help."' >> $(ALIASES)
	@echo 'echo " "' >> $(ALIASES)

# to create packages
PACKAGE_FILES = Manifest.text Manifest.bin
CHECK_SUMS    = SHA1SUMS MD5SUMS
DIR           := $(shell pwd)
CUR_DIR       := $(notdir $(DIR))

# to upload packages
FTP_USER      = anonymous
FTP_PASSWD    = redstar_@users.sourceforge.net
FTP_DIR       = /incoming
FTP_HOST      = upload.sourceforge.net

# to upload homepage
WEB_SERVER    = shell.sourceforge.net
WEB_PATH      = /home/groups/m/ma/matrixss/htdocs
SCP_USER      = redstar_
WEB_FILES     = doc/www/*.html README PackageInfo.g ChangeLog 
WEBMAN_PATH   = /home/groups/m/ma/matrixss/htdocs/manual
WEBMAN_FILES  = htm/*.htm

# to perform test/benchmark
TEST_FILE     = tst/test.g
BENCH_FILE    = tst/bench.g
PROFILE_FILE  = tst/profile.g

DATE = $(shell date +%Y%m%d)
# package default name
ifeq ($(PACKAGE),)
	PACKAGE = matrixss-0.9.$(DATE)
	DEBPKG  = gap-matrix-schreiersims_0.9.$(DATE)
endif

# output package files
PACKAGES = tmp/$(PACKAGE).tar.gz tmp/$(PACKAGE)-win.zip \
           tmp/$(PACKAGE).zoo tmp/$(PACKAGE).tar.bz2

# checksum files regexp
SHA1 := $(word 1, $(CHECK_SUMS))
MD5  := $(word 2, $(CHECK_SUMS))
CHECKS := $(SHA1)\|$(MD5)

# for making html manual
HTML_MAN = ../../etc/convert.pl

all:

test:
	gap -q < $(TEST_FILE)

bench:
	gap -o 1g -K 2g -q < $(BENCH_FILE)

profile:
	gap -q < $(PROFILE_FILE)

release: upload www

checksums: $(CHECK_SUMS)

$(CHECK_SUMS): ChangeLog manual
	md5sum `grep --invert-match --extended-regexp --regexp=$(CHECKS) Manifest.text` > $(MD5)
	md5sum -b `grep --invert-match --extended-regexp --regexp=$(CHECKS) Manifest.bin` >> $(MD5)
	sha1sum --text `grep --invert-match --extended-regexp --regexp=$(CHECKS) Manifest.text` > $(SHA1)
	sha1sum --binary `grep --invert-match --extended-regexp --regexp=$(CHECKS) Manifest.bin` >> $(SHA1)

NEWS: doc/www/news.html
	html2text < $< > $@

ChangeLog: NEWS
	echo -e "Matrix Schreier-Sims ChangeLog\n" | cvs2cl --prune --separate-header --gmt --usermap AUTHORS --header -

deb: $(PACKAGES)
	cvs tag -F $(DEB_TAG) `grep --invert-match --extended-regexp --regexp=$(CHECKS) Manifest.text` Makefile doc/code.msk doc/config.matrixss doc/report/schreiersims.bib debian/*
	cp tmp/$(PACKAGE).tar.gz tmp/gap-matrix-schreiersims/$(DEBPKG).orig.tar.gz
	cvs-buildpackage -R $(DIR)/tmp -rfakeroot -kredstar -Mmatrixss

manual: 
	$(MAKE) -C doc $@
	$(HTML_MAN) -c -t -i -n matrixss doc htm

package: manual $(PACKAGES)

$(PACKAGES): ChangeLog checksums
	mkdir --parents tmp/gap-matrix-schreiersims
#	mkdir --parents tmp/doc
#	cp ../../doc/gapmacro.tex tmp/doc/
	cat $(PACKAGE_FILES) | sed --expression='s/^/$(CUR_DIR)\//' | tar --verify --verbose --create --directory .. --file tmp/$(PACKAGE).tar --files-from - 
	bzip2 --keep --verbose --best --force tmp/$(PACKAGE).tar
	bzip2 --test --verbose tmp/$(PACKAGE).tar.bz2
	gzip --best --force --verbose tmp/$(PACKAGE).tar
	gzip --verbose --test tmp/$(PACKAGE).tar.gz
	cd .. && cat $(CUR_DIR)/Manifest.text | sed --expression='s/^/$(CUR_DIR)\//' | zip -v -l -9 $(CUR_DIR)/tmp/$(PACKAGE)-win.zip -@ && cd $(CUR_DIR)
	cd .. && cat $(CUR_DIR)/Manifest.bin | sed --expression='s/^/$(CUR_DIR)\//' | zip -v -9 $(CUR_DIR)/tmp/$(PACKAGE)-win.zip -@ && cd $(CUR_DIR)
	unzip -t tmp/$(PACKAGE)-win.zip
	cd .. && cat $(CUR_DIR)/Manifest.text | sed --expression='s/^\([^!/]\)/$(CUR_DIR)\/\1/' | perl -e '{ while(<>) { print; print "!TEXT!\n/END\n"; } }' | zoo achPI $(CUR_DIR)/tmp/$(PACKAGE).zoo && cd $(CUR_DIR)
	cd .. && cat $(CUR_DIR)/Manifest.bin | sed --expression='s/^\([^!/]\)/$(CUR_DIR)\/\1/' | zoo ahPI $(CUR_DIR)/tmp/$(PACKAGE).zoo && cd $(CUR_DIR)
	zoo xN tmp/$(PACKAGE).zoo

upload:
	ncftpput -u $(FTP_USER) -p $(FTP_PASSWD) $(FTP_HOST) $(FTP_DIR) $(PACKAGES) tmp/gap-matrix-schreiersims/$(DEBPKG)*

www: ChangeLog
	scp $(WEB_FILES) $(SCP_USER)@$(WEB_SERVER):$(WEB_PATH)
	scp $(WEBMAN_FILES) $(SCP_USER)@$(WEB_SERVER):$(WEBMAN_PATH)

distclean:
	rm --force $(PACKAGES) $(CHECK_SUMS) ChangeLog NEWS

clean:
	rm --force $(PACKAGES) $(CHECK_SUMS) ChangeLog NEWS htm/*.htm
	rm --recursive --force tmp *~ *.bak
	$(MAKE) -C doc clean

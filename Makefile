# to create packages
PACKAGE_FILES = Manifest
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

# to perform test/benchmark
TEST_FILE     = tst/test.g
BENCH_FILE    = tst/bench.g
PROFILE_FILE  = tst/profile.g

# package default name
ifeq ($(PACKAGE),)
	PACKAGE = matrixss
endif

# output package files
PACKAGES = $(PACKAGE).tar.gz $(PACKAGE).zip $(PACKAGE).zoo $(PACKAGE).tar.bz2

# checksum files regexp
SHA1 := $(word 1, $(CHECK_SUMS))
MD5  := $(word 2, $(CHECK_SUMS))
CHECKS := $(SHA1)\|$(MD5)

all:

test:
	gap -q < $(TEST_FILE)

bench:
	gap -q < $(BENCH_FILE)

profile:
	gap -q < $(PROFILE_FILE)

release: upload www

checksums: $(CHECK_SUMS)

$(CHECK_SUMS): ChangeLog
	md5sum `grep --invert-match --extended-regexp --regexp=$(CHECKS) Manifest` > $(MD5)
	sha1sum --text `grep --invert-match --extended-regexp --regexp=$(CHECKS) Manifest` > $(SHA1)

ChangeLog: 
	echo -e "Matrix Schreier-Sims ChangeLog\n" | cvs2cl --prune --separate-header --gmt --usermap AUTHORS --header -

package: $(PACKAGES)

$(PACKAGES): ChangeLog checksums
	cat $(PACKAGE_FILES) | sed --expression='s/^/$(CUR_DIR)\//' | tar --verify --verbose --create --directory .. --file $(PACKAGE).tar --files-from - 
	bzip2 --keep --verbose --best --force $(PACKAGE).tar
	bzip2 --test --verbose $(PACKAGE).tar.bz2
	gzip --best --force --verbose $(PACKAGE).tar
	gzip --verbose --test $(PACKAGE).tar.gz
	cd .. && cat $(CUR_DIR)/$(PACKAGE_FILES) | sed --expression='s/^/$(CUR_DIR)\//' | zip -v -9 $(CUR_DIR)/$(PACKAGE).zip -@ && cd $(CUR_DIR)
	unzip -t $(PACKAGE).zip
	cd .. && cat $(CUR_DIR)/$(PACKAGE_FILES) | sed --expression='s/^\([^!/]\)/$(CUR_DIR)\/\1/' | perl -e '{ while(<>) { print; print "!TEXT!\n/END\n"; } }' | zoo achPI $(CUR_DIR)/$(PACKAGE).zoo && cd $(CUR_DIR)
	zoo xN $(PACKAGE).zoo

upload: package
	ncftpput -u $(FTP_USER) -p $(FTP_PASSWD) $(FTP_HOST) $(FTP_DIR) $(PACKAGES)

www: ChangeLog
	scp $(WEB_FILES) $(SCP_USER)@$(WEB_SERVER):$(WEB_PATH)

distclean:
	rm --force $(PACKAGES) $(CHECK_SUMS) ChangeLog

clean:
	rm --force $(PACKAGES) $(CHECK_SUMS) ChangeLog
	rm --recursive --force *~ *.bak
	$(MAKE) -C doc/report clean

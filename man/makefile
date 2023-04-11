all: aflock.txt mflock.txt chromflock.txt

VERSION_MAJOR != awk -F' ' '{if($$2=="CF_VERSION_MAJOR") print $$3}' ../src/cf_version.h | sed 's/"//g'
VERSION_MINOR != awk -F' ' '{if($$2=="CF_VERSION_MINOR") print $$3}' ../src/cf_version.h | sed 's/"//g'
VERSION_PATCH != awk -F' ' '{if($$2=="CF_VERSION_PATCH") print $$3}' ../src/cf_version.h | sed 's/"//g'
VERSION="$(VERSION_MAJOR).$(VERSION_MINOR).$(VERSION_PATCH)"

.SECONDARY: chromflock.1

%.txt: %.1
	man ./$< | col -b > $@

%.1: %.md
	pandoc $< -s \
                 --standalone \
	         --to man \
                 -o $@
	sed -i "s/ABCDE/$(VERSION)/g" $@

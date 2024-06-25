all: joss

joss: paper.pdf

paper.pdf: paper.md paper.bib
	@docker run --rm --volume $$PWD:/data --user $$(id -u):$$(id -g) --env JOURNAL=joss openjournals/inara

%:
	cd src && $(MAKE) $*

commit:
	@git commit -am "."
	@git push

clean:
	rm -rf paper.jats paper.pdf media

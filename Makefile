.PHONY: default clean cleanup publish render

default: render

cleanup:
	rm -rf site_libs/

clean: cleanup
	rm -rf .quarto/_freeze
	rm -rf _freeze/

publish:
	quarto publish gh-pages

render:
	quarto render

preview:
	quarto preview

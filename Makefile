.PHONY: install test lint validate stub clean

install:
	pip install -e ".[dev]"

test:
	pytest -q

lint:
	ruff check src tests

validate:
	whale validate-study
	whale validate-status

stub:
	nextflow run workflows/methylation -profile test -stub

clean:
	rm -rf work results .nextflow .nextflow.log*

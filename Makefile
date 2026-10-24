.PHONY: install test lint validate stub real-run clean

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

# Every real command, in Docker, on the simulated fixtures; then check the calls.
real-run:
	nextflow run workflows/methylation -profile test,docker
	python3 tests/check_real_run.py results/test/TEST2026_FIXTURE_WGBS

clean:
	rm -rf work results .nextflow .nextflow.log*

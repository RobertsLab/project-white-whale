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
	nextflow run workflows/snp_mask -profile test -stub

# Every real command, in Docker, on the simulated fixtures; then check the calls,
# and that the SNP mask finds C/T SNPs planted in a copy of the reference.
real-run:
	nextflow run workflows/methylation -profile test,docker
	python3 tests/check_real_run.py results/test/TEST2026_FIXTURE_WGBS
	python3 tests/plant_snps.py results/test/snp_reference
	nextflow run workflows/snp_mask -profile test,docker --reference_dir results/test/snp_reference
	python3 tests/plant_snps.py --check results/test/snp_reference results/test/snp_mask

clean:
	rm -rf work results .nextflow .nextflow.log*

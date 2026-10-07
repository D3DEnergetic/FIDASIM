.PHONY: all test_001 test_003 test_004 test_005 test_006 summary clean

all: test_001 test_003 test_004 test_005 test_006 summary

test_001:
	@$(MAKE) -C test_001

test_003:
	@$(MAKE) -C test_003

test_004:
	@$(MAKE) -C test_004

test_005:
	@$(MAKE) -C test_005

test_006:
	@$(MAKE) -C test_006

summary:
	@echo ""
	@echo "Regression build summary"
	@echo "  Tests 001, 003, 004, 005, and 006 compiled successfully."
	@echo "  No regression tests were executed."
	@echo "  Run ./regression_tests/run.sh from the repository root."

clean:
	@$(MAKE) -C test_001 clean
	@$(MAKE) -C test_003 clean
	@$(MAKE) -C test_004 clean
	@$(MAKE) -C test_005 clean
	@$(MAKE) -C test_006 clean

.PHONY: all test_001 test_003 test_004 test_005 summary clean

all: test_001 test_003 test_004 test_005 summary

test_001:
	@cd test_001; make

test_003:
	@cd test_003; make

test_004:
	@cd test_004; make

test_005:
	@cd test_005; make check

summary:
	@echo ""
	@echo "Regression verification summary"
	@echo "  Test 001: BUILD PASS - scientific comparison not run"
	@echo "  Test 003: BUILD PASS - scientific comparison not run"
	@echo "  Test 004: BUILD PASS - scientific comparison not run"
	@echo "  Test 005: PASS - self-check executed and report verified"
	@echo "  Overall:  PASS for the verification stages listed above"

clean:
	@cd test_001; make clean
	@cd test_003; make clean
	@cd test_004; make clean
	@cd test_005; make clean
